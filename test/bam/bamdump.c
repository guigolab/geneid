/* Test driver for the bamcov reader: print coverage intervals overlapping a range.
   Usage: bamdump <file.bam> <chrom> <start> <end> [wantStrand libmode]
     wantStrand: + | - | . (unstranded, default)   libmode: rf | fr | none
   Output: one line per run -> "start\tend\tdepth". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bamcov.h"

static void rec(long s, long e, float v, void* ud){
  (void)ud;
  printf("%ld\t%ld\t%g\n", s, e, v);
}

int main(int argc, char** argv){
  if (argc != 5 && argc != 7){
    fprintf(stderr, "usage: %s file.bam chrom start end [wantStrand libmode]\n", argv[0]);
    return 2;
  }
  char want = 0;
  int libMode = BAMLIB_NONE;
  if (argc == 7){
    if (argv[5][0] == '+' || argv[5][0] == '-') want = argv[5][0];
    if (!strcmp(argv[6], "rf")) libMode = BAMLIB_RF;
    else if (!strcmp(argv[6], "fr")) libMode = BAMLIB_FR;
  }
  BamCov* bc = bamOpen(argv[1]);
  if (!bc){ fprintf(stderr, "cannot open BAM (or no index): %s\n", argv[1]); return 1; }
  long n = bamCoverageQuery(bc, argv[2], atol(argv[3]), atol(argv[4]), want, libMode, rec, NULL);
  bamClose(bc);
  if (n < 0){ fprintf(stderr, "query error\n"); return 1; }
  return 0;
}
