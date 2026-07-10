/* Test driver for the bamcov reader: print coverage intervals overlapping a range.
   Usage: bamdump <file.bam> <chrom> <start> <end>
   Output: one line per run -> "start\tend\tdepth". */
#include <stdio.h>
#include <stdlib.h>
#include "bamcov.h"

static void rec(long s, long e, float v, void* ud){
  (void)ud;
  printf("%ld\t%ld\t%g\n", s, e, v);
}

int main(int argc, char** argv){
  if (argc != 5){ fprintf(stderr, "usage: %s file.bam chrom start end\n", argv[0]); return 2; }
  BamCov* bc = bamOpen(argv[1]);
  if (!bc){ fprintf(stderr, "cannot open BAM (or no index): %s\n", argv[1]); return 1; }
  long n = bamCoverageQuery(bc, argv[2], atol(argv[3]), atol(argv[4]), rec, NULL);
  bamClose(bc);
  if (n < 0){ fprintf(stderr, "query error\n"); return 1; }
  return 0;
}
