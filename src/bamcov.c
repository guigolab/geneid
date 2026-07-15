/*************************************************************************
*   Module: bamcov -- per-split RNA-seq coverage from an indexed BAM.     *
*                                                                        *
*   htslib wrapper: open + header + index, then a region query that walks *
*   only the reads overlapping [start,end), accumulates per-base depth    *
*   (CIGAR M/=/X; split at N and D), and reports it as run-length         *
*   intervals through the shared bwIntervalCB callback -- the same shape  *
*   the bigWig reader emits, so the downstream sr[] fill is identical.    *
*                                                                        *
*   Built only under WITH_HTSLIB (see the Makefile); nothing else in      *
*   geneid references htslib.                                            *
*                                                                        *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "bamcov.h"

struct BamCov {
  htsFile*    fp;
  sam_hdr_t*  hdr;
  hts_idx_t*  idx;
};

/* Flags that exclude a read from coverage. */
#define BAMCOV_SKIP (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | \
                     BAM_FQCFAIL | BAM_FDUP)

/* Transcription strand of a read under the given library mode. RF (dUTP):
   read1/single-end antisense, read2 sense. FR: the reverse. */
static char readTxnStrand(const bam1_t* b, int libMode)
{
  int rev   = (b->core.flag & BAM_FREVERSE) != 0;
  int read2 = (b->core.flag & BAM_FPAIRED) && (b->core.flag & BAM_FREAD2);
  int fwd;                                  /* is the transcription strand '+' ? */
  if (libMode == BAMLIB_FR)
    fwd = read2 ? rev : !rev;               /* read1/SE sense, read2 antisense */
  else
    fwd = read2 ? !rev : rev;               /* RF: read1/SE antisense, read2 sense */
  return fwd ? '+' : '-';
}

BamCov* bamOpen(const char* path)
{
  BamCov* bc = (BamCov*) calloc(1, sizeof(BamCov));
  if (bc == NULL) return NULL;

  bc->fp = hts_open(path, "r");
  if (bc->fp == NULL) { free(bc); return NULL; }

  /* Sniff: accept only BAM (a text HSP GFF / bigWig is not, -> fall back). */
  const htsFormat* fmt = hts_get_format(bc->fp);
  if (fmt == NULL || fmt->format != bam) { bamClose(bc); return NULL; }

  bc->hdr = sam_hdr_read(bc->fp);
  if (bc->hdr == NULL) { bamClose(bc); return NULL; }

  bc->idx = sam_index_load(bc->fp, path);   /* needs a .bai/.csi next to the BAM */
  if (bc->idx == NULL) { bamClose(bc); return NULL; }

  return bc;
}

void bamClose(BamCov* bc)
{
  if (bc == NULL) return;
  if (bc->idx != NULL) hts_idx_destroy(bc->idx);
  if (bc->hdr != NULL) sam_hdr_destroy(bc->hdr);
  if (bc->fp  != NULL) hts_close(bc->fp);
  free(bc);
}

/* Sum mapped reads over every reference from the index meta -- the same numbers
   `samtools idxstats` prints, obtained without touching alignment records.
   hts_idx_get_stat returns <0 for a reference with no index bin, which covers
   both a reference that simply has no reads (common: a chr21-only subset BAM
   still carries the full header) and a statless index. We skip such references
   (count them as 0) and fail (return -1) only when NO reference carries stats,
   so a subset BAM still yields the true library size. */
long bamMappedReads(BamCov* bc)
{
  int nref, tid, anystat = 0;
  long total = 0;

  if (bc == NULL || bc->idx == NULL) return -1;
  nref = hts_idx_nseq(bc->idx);
  for (tid = 0; tid < nref; tid++) {
    uint64_t mapped = 0, unmapped = 0;
    if (hts_idx_get_stat(bc->idx, tid, &mapped, &unmapped) < 0)
      continue;                 /* reference with no reads / no per-ref stats */
    anystat = 1;
    total += (long) mapped;
  }
  return anystat ? total : -1;
}

long bamCoverageQuery(BamCov* bc, const char* chrom, long start, long end,
                      char wantStrand, int libMode, bwIntervalCB cb, void* userData)
{
  int tid = sam_hdr_name2tid(bc->hdr, chrom);
  if (tid < 0) return 0;                       /* chrom not in this BAM */
  if (end <= start) return 0;

  long span = end - start;
  int* depth = (int*) calloc((size_t) span, sizeof(int));
  if (depth == NULL) return -1;

  hts_itr_t* iter = sam_itr_queryi(bc->idx, tid, start, end);
  if (iter == NULL) { free(depth); return -1; }

  bam1_t* b = bam_init1();
  if (b == NULL) { hts_itr_destroy(iter); free(depth); return -1; }

  long ret;
  while ((ret = sam_itr_next(bc->fp, iter, b)) >= 0) {
    if (b->core.flag & BAMCOV_SKIP) continue;
    if (wantStrand && readTxnStrand(b, libMode) != wantStrand) continue;

    long refpos = b->core.pos;                 /* 0-based leftmost ref coord */
    uint32_t* cig = bam_get_cigar(b);
    int n = b->core.n_cigar, i;
    for (i = 0; i < n; i++) {
      int op  = bam_cigar_op(cig[i]);
      int len = bam_cigar_oplen(cig[i]);
      int type = bam_cigar_type(op);           /* bit1 = consumes query, bit2 = ref */
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        /* read-supported reference bases: add depth over [refpos, refpos+len) */
        long s = refpos > start ? refpos : start;
        long e = refpos + len < end ? refpos + len : end;
        long p;
        for (p = s; p < e; p++) depth[p - start]++;
      }
      if (type & 2) refpos += len;             /* N, D, M/=/X advance the reference */
    }
  }
  bam_destroy1(b);
  hts_itr_destroy(iter);

  if (ret < -1) { free(depth); return -1; }    /* truncated/corrupt iteration */

  /* Run-length encode the non-zero depth into intervals and report them. */
  long count = 0, i = 0;
  while (i < span) {
    int v = depth[i];
    if (v == 0) { i++; continue; }
    long j = i;
    while (j < span && depth[j] == v) j++;
    cb(start + i, start + j, (float) v, userData);
    count++;
    i = j;
  }

  free(depth);
  return count;
}

long bamJunctionQuery(BamCov* bc, const char* chrom, long start, long end,
                      bamJunctionCB cb, void* userData)
{
  int tid = sam_hdr_name2tid(bc->hdr, chrom);
  if (tid < 0) return 0;

  hts_itr_t* iter = sam_itr_queryi(bc->idx, tid, start, end);
  if (iter == NULL) return -1;

  bam1_t* b = bam_init1();
  if (b == NULL) { hts_itr_destroy(iter); return -1; }

  long count = 0, ret;
  while ((ret = sam_itr_next(bc->fp, iter, b)) >= 0) {
    if (b->core.flag & BAMCOV_SKIP) continue;

    /* Junction strand: XS (genomic transcript strand; STAR/HISAT2) if present,
       else minimap2 ts (transcript strand RELATIVE TO THE READ) mapped to the
       genome via the read orientation, else '.' -- the caller then infers the
       strand from the splice motif. */
    char strand = '.';
    uint8_t* xs = bam_aux_get(b, "XS");
    if (xs != NULL) {
      char v = bam_aux2A(xs);
      if (v == '+' || v == '-') strand = v;
    }
    if (strand == '.') {
      uint8_t* ts = bam_aux_get(b, "ts");
      if (ts != NULL) {
        char v = bam_aux2A(ts);
        if (v == '+' || v == '-')
          strand = (b->core.flag & BAM_FREVERSE) ? (v == '+' ? '-' : '+') : v;
      }
    }

    long refpos = b->core.pos;
    uint32_t* cig = bam_get_cigar(b);
    int n = b->core.n_cigar, i;
    for (i = 0; i < n; i++) {
      int op  = bam_cigar_op(cig[i]);
      int len = bam_cigar_oplen(cig[i]);
      if (op == BAM_CREF_SKIP) {              /* N = intron gap */
        cb(refpos, refpos + len, strand, userData);
        count++;
      }
      if (bam_cigar_type(op) & 2) refpos += len;
    }
  }
  bam_destroy1(b);
  hts_itr_destroy(iter);

  if (ret < -1) return -1;                    /* truncated/corrupt iteration */
  return count;
}
