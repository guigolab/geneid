/* Test driver for the bigbed reader: print records overlapping a range.
   Usage: bbdump <file.bb> <chrom> <start> <end>
   Output: one line per record -> "start\tend\trest" (rest = tabbed BED tail). */
#include <stdio.h>
#include <stdlib.h>
#include "bigbed.h"

static void rec(long s, long e, const char* rest, void* ud){
  (void)ud;
  printf("%ld\t%ld\t%s\n", s, e, rest);
}

int main(int argc, char** argv){
  if (argc != 5){ fprintf(stderr, "usage: %s file.bb chrom start end\n", argv[0]); return 2; }
  BigBed* bb = bbOpen(argv[1]);
  if (!bb){ fprintf(stderr, "cannot open bigBed: %s\n", argv[1]); return 1; }
  long n = bbQuery(bb, argv[2], atol(argv[3]), atol(argv[4]), rec, NULL);
  bbClose(bb);
  if (n < 0){ fprintf(stderr, "query error\n"); return 1; }
  return 0;
}
