/*************************************************************************
*   Module: bigbed -- minimal random-access reader for the UCSC bigBed   *
*   (BBI) binary format. Range queries over the R-tree index return only *
*   the data blocks overlapping [start,end) -- never the whole file --    *
*   so geneid can pull just the evidence for the current split.          *
*                                                                        *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#ifndef BIGBED_H
#define BIGBED_H

typedef struct BigBed BigBed;

/* Called once per record overlapping the query. `rest` is the tab-joined
   trailing BED fields after chrom/start/end (i.e. name, score, strand, ...);
   start/end are the record's 0-based half-open coordinates. */
typedef void (*bbRecordCB)(long start, long end, const char* rest, void* userData);

/* Open a bigBed file (reads + validates the header and the chrom B+-tree).
   Returns NULL on error (not a bigBed, unreadable, unsupported byte order). */
BigBed* bbOpen(const char* path);

void bbClose(BigBed* bb);

/* Invoke cb for every record overlapping [start,end) on `chrom` (0-based).
   Returns the number of records reported, or -1 on I/O/decode error. An
   unknown chrom name is not an error -> returns 0. */
long bbQuery(BigBed* bb, const char* chrom, long start, long end,
             bbRecordCB cb, void* userData);

#endif
