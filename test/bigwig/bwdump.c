/* Test driver for the bigwig reader: print intervals overlapping a range.
   Usage: bwdump <file.bw> <chrom> <start> <end>
   Output: one line per interval -> "start\tend\tvalue". */
#include <stdio.h>
#include <stdlib.h>
#include "bigwig.h"

static void rec(long s, long e, float v, void* ud){
  (void)ud;
  printf("%ld\t%ld\t%g\n", s, e, v);
}

int main(int argc, char** argv){
  if (argc != 5){ fprintf(stderr, "usage: %s file.bw chrom start end\n", argv[0]); return 2; }
  BigWig* bw = bwOpen(argv[1]);
  if (!bw){ fprintf(stderr, "cannot open bigWig: %s\n", argv[1]); return 1; }
  long n = bwQuery(bw, argv[2], atol(argv[3]), atol(argv[4]), rec, NULL);
  bwClose(bw);
  if (n < 0){ fprintf(stderr, "query error\n"); return 1; }
  return 0;
}
