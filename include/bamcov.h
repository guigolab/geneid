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

/* Open an indexed, coordinate-sorted BAM (reads header + .bai/.csi index).
   Returns NULL on error, or when the file is not a BAM (so this doubles as the
   -S format sniff) or has no index. */
BamCov* bamOpen(const char* path);

void bamClose(BamCov* bc);

/* Invoke cb for each maximal run of constant non-zero depth overlapping
   [start,end) on `chrom` (0-based half-open genomic coords, value = depth).
   Returns the number of runs reported, or -1 on error; unknown chrom -> 0. */
long bamCoverageQuery(BamCov* bc, const char* chrom, long start, long end,
                      bwIntervalCB cb, void* userData);

#endif
