/*************************************************************************
*   Module: bigwig -- minimal random-access reader for the UCSC bigWig    *
*   (BBI) binary format. Range queries over the R-tree index decode only  *
*   the data sections overlapping [start,end) -- never the whole file --   *
*   so geneid can pull just the RNA-seq coverage for the current split.   *
*                                                                        *
*   bigWig shares the BBI container with bigBed (same header layout,      *
*   chrom B+-tree and R-tree index); only the leaf data blocks differ --   *
*   bigWig stores wiggle sections (bedGraph/varStep/fixedStep) of typed   *
*   per-base values instead of BED records.                              *
*                                                                        *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#ifndef BIGWIG_H
#define BIGWIG_H

typedef struct BigWig BigWig;

/* Called once per stored interval overlapping the query. start/end are the
   interval's 0-based half-open genomic coordinates (as stored, NOT clipped to
   the query); value is the signal for every base in [start,end). */
typedef void (*bwIntervalCB)(long start, long end, float value, void* userData);

/* Open a bigWig file (reads + validates the header and the chrom B+-tree).
   Returns NULL on error (not a bigWig, unreadable, unsupported byte order). */
BigWig* bwOpen(const char* path);

void bwClose(BigWig* bw);

/* Invoke cb for every stored interval overlapping [start,end) on `chrom`
   (0-based half-open: an interval touching exactly the query edge does not
   overlap). Returns the number of intervals reported, or -1 on I/O/decode
   error. An unknown chrom name is not an error -> returns 0. */
long bwQuery(BigWig* bw, const char* chrom, long start, long end,
             bwIntervalCB cb, void* userData);

#endif
