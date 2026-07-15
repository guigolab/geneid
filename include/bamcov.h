/*************************************************************************
*   Module: bamcov -- per-split RNA-seq COVERAGE from an indexed BAM,     *
*   via htslib. Region queries over the .bai/.csi index compute depth for *
*   only the current fragment -- never the whole file -- so geneid can    *
*   take a STAR / minimap2 BAM directly instead of a precomputed bigWig.  *
*                                                                        *
*   Depth rule: reference positions consumed by CIGAR M/=/X are counted   *
*   (+1 per overlapping read); N (intron) and D (deletion) advance the    *
*   reference without adding depth, so reads are "split" at introns like  *
*   `bedtools genomecov -split`. Reads flagged unmapped/secondary/        *
*   supplementary/qcfail/duplicate are skipped. Unstranded (every read    *
*   counts); stranded library handling is a later addition.              *
*                                                                        *
*   Compiled only when geneid is built WITH_HTSLIB=1; the default build   *
*   never references htslib.                                             *
*                                                                        *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#ifndef BAMCOV_H
#define BAMCOV_H

#include "bigwig.h"   /* bwIntervalCB: the coverage callback, shared with bigWig */

typedef struct BamCov BamCov;

/* Library strandedness for deriving a read's transcription strand from its FLAG.
   NONE = unstranded (every read counts on both strands). RF = "reverse"/dUTP
   (read1 & single-end reads are antisense, read2 sense) -- the common Illumina
   protocol. FR = "forward" (read1 & single-end sense, read2 antisense). */
#define BAMLIB_NONE 0
#define BAMLIB_RF   1
#define BAMLIB_FR   2

/* Open an indexed, coordinate-sorted BAM (reads header + .bai/.csi index).
   Returns NULL on error, or when the file is not a BAM (so this doubles as the
   -S format sniff) or has no index. */
BamCov* bamOpen(const char* path);

void bamClose(BamCov* bc);

/* Invoke cb for each maximal run of constant non-zero depth overlapping
   [start,end) on `chrom` (0-based half-open genomic coords, value = depth).
   wantStrand selects which reads contribute by transcription strand (derived
   from the FLAG via libMode): '+' or '-' for a single strand, or 0 for
   unstranded (all reads, libMode ignored). Returns the number of runs
   reported, or -1 on error; unknown chrom -> 0. */
long bamCoverageQuery(BamCov* bc, const char* chrom, long start, long end,
                      char wantStrand, int libMode, bwIntervalCB cb, void* userData);

/* Called once per spliced-read junction (a CIGAR N gap) overlapping the query.
   start/end are the intron's 0-based half-open reference coordinates; strand is
   the transcription strand from the read's XS tag ('+' or '-'). */
typedef void (*bamJunctionCB)(long start, long end, char strand, void* userData);

/* Invoke cb for every CIGAR-N junction of reads overlapping [start,end) on
   `chrom`. Only reads carrying an XS strand tag contribute (a junction needs a
   strand to become an intron); reads flagged unmapped/secondary/supplementary/
   qcfail/dup are skipped. Returns the number of junctions reported, or -1 on
   error; unknown chrom -> 0. */
long bamJunctionQuery(BamCov* bc, const char* chrom, long start, long end,
                      bamJunctionCB cb, void* userData);

/* Total mapped reads across all references, read from the index meta (like
   `samtools idxstats`) -- no read scan. Used to set MRM (millions of reads
   mapped) for the rpkm report when -N is not given. Returns the count, or -1
   if the index carries no per-reference stats. */
long bamMappedReads(BamCov* bc);

#endif
