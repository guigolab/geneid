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

long bamCoverageQuery(BamCov* bc, const char* chrom, long start, long end,
                      bwIntervalCB cb, void* userData)
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

    /* Junction strand comes from the XS tag (spliced-transcript strand, set by
       STAR/HISAT2). Without it a junction cannot be placed on a strand. */
    uint8_t* xs = bam_aux_get(b, "XS");
    char strand = (xs != NULL) ? bam_aux2A(xs) : '.';
    if (strand != '+' && strand != '-') continue;

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
