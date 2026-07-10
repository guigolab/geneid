/* Test driver for the bamcov junction reader: print CIGAR-N junctions overlapping
   a range. Usage: juncdump <file.bam> <chrom> <start> <end>
   Output: one line per junction occurrence -> "start\tend\tstrand". */
#include <stdio.h>
#include <stdlib.h>
#include "bamcov.h"

static void rec(long s, long e, char strand, void* ud){
  (void)ud;
  printf("%ld\t%ld\t%c\n", s, e, strand);
}

int main(int argc, char** argv){
  if (argc != 5){ fprintf(stderr, "usage: %s file.bam chrom start end\n", argv[0]); return 2; }
  BamCov* bc = bamOpen(argv[1]);
  if (!bc){ fprintf(stderr, "cannot open BAM (or no index): %s\n", argv[1]); return 1; }
  long n = bamJunctionQuery(bc, argv[2], atol(argv[3]), atol(argv[4]), rec, NULL);
  bamClose(bc);
  if (n < 0){ fprintf(stderr, "query error\n"); return 1; }
  return 0;
}
