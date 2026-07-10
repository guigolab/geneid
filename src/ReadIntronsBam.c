/*************************************************************************
*   Module: ReadIntronsBam                                               *
*                                                                        *
*   Per-split intron ingest from an indexed BAM: for the current         *
*   fragment [l1,l2] of a locus, collect the spliced-read junctions      *
*   (CIGAR N gaps, via bamcov's bamJunctionQuery), tally identical       *
*   junctions into a read count, and commit each as an Intron evidence   *
*   feature through the shared AddEvidenceExon -- the same routing the    *
*   text/-R and bigBed Intron features get (which waives the intron-      *
*   length penalty for junction-supported introns).                      *
*                                                                        *
*   A BAM junction spanning reference [p,p+len) maps exactly onto a       *
*   bigBed Intron record (start=p, end=p+len -> geneid 1-based p+1..p+len)*
*   so this mirrors ReadExonsBigBed. Junction strand comes from the read  *
*   XS tag (bamcov); the read count becomes the evidence score.          *
*                                                                        *
*   htslib-free: all BAM decoding lives in bamcov.c. Built only under     *
*   WITH_HTSLIB (see the Makefile).                                      *
*                                                                        *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#include "geneid.h"
#include "bamcov.h"

extern float EvidenceFactor;
extern float EvidenceEW;
extern int FWD, RVS;

/* One spliced-read junction occurrence (0-based half-open reference coords). */
typedef struct { long start, end; char strand; } juncRec;

typedef struct {
  juncRec* v;
  long     n, cap;
} juncList;

static void collectCB(long start, long end, char strand, void* ud){
  juncList* L = (juncList*)ud;
  if (L->n >= L->cap){
    L->cap = L->cap ? L->cap * 2 : 256;
    L->v = (juncRec*)realloc(L->v, L->cap * sizeof(juncRec));
    if (!L->v) printError("Not enough memory: BAM junction buffer");
  }
  L->v[L->n].start  = start;
  L->v[L->n].end    = end;
  L->v[L->n].strand = strand;
  L->n++;
}

/* Sort junctions by (start, end, strand) so identical ones are adjacent. */
static int cmpJunc(const void* a, const void* b){
  const juncRec* x = (const juncRec*)a;
  const juncRec* y = (const juncRec*)b;
  if (x->start != y->start) return (x->start > y->start) - (x->start < y->start);
  if (x->end   != y->end)   return (x->end   > y->end)   - (x->end   < y->end);
  return (x->strand > y->strand) - (x->strand < y->strand);
}

/* Populate external->evidence[0] with the junction-derived Intron features whose
   acceptor (start+1) falls in this fragment's OWNED range (ownedLo, ownedHi] --
   the same non-overlapping per-fragment assignment ReadExonsBigBed makes -- so
   each distinct junction is committed in exactly one fragment. Returns the
   number of evidence exon slots produced (including frame/intron replicas). */
long ReadIntronsBam(BamCov* bc, packExternalInformation* external, dict* d,
                    char* Locus, long l1, long l2, long ownedLo, long ownedHi){
  packEvidence* ev = external->evidence[0];
  juncList L = { NULL, 0, 0 };
  long lastAcceptor = -INFI;
  long i;

  ev->nvExons = 0;
  ev->nvSites = 0;

  /* Reads overlapping [l1,l2+1) (0-based); ownership (below) is in 1-based
     acceptor units, matching the GFF/bigBed evidence assignment. */
  if (bamJunctionQuery(bc, Locus, l1, l2 + 1, collectCB, &L) < 0)
    printError("BAM junction query failed");

  qsort(L.v, L.n, sizeof(juncRec), cmpJunc);

  /* Walk the sorted junctions, tallying each distinct (start,end,strand) run. */
  i = 0;
  while (i < L.n){
    long j = i + 1;
    while (j < L.n && cmpJunc(&L.v[i], &L.v[j]) == 0) j++;
    long readCount = j - i;                    /* reads supporting this junction */

    long begin  = L.v[i].start + 1;            /* 0-based -> 1-based intron start */
    long end    = L.v[i].end;                  /* 0-based half-open end == 1-based inclusive */
    char strand = L.v[i].strand;
    char lineCopy[MAXLINE];

    i = j;

    if (begin <= ownedLo || begin > ownedHi) continue;   /* owned by another fragment */
    if ((strand == '+' && !FWD) || (strand == '-' && !RVS)) continue;

    float score = (float) readCount * EvidenceFactor + EvidenceEW;
    int three = 1;                             /* introns carry no frame -> 3 copies */

    sprintf(lineCopy, "%s\t.\tIntron\t%ld\t%ld\t%ld\t%c\t.\t.",
            Locus, begin, end, readCount, strand);

    strcpy((ev->vExons + ev->nvExons)->Type, sINTRON);
    (ev->vExons + ev->nvExons)->Score  = score;
    (ev->vExons + ev->nvExons)->Strand = strand;
    strcpy((ev->vExons + ev->nvExons)->Group, NOGROUP);
    (ev->vSites + ev->nvSites)->Position     = begin;
    (ev->vSites + ev->nvSites + 1)->Position = end;

    AddEvidenceExon(external, 0, d, &three, U2, U2, &lastAcceptor, lineCopy);
  }

  free(L.v);
  external->i1vExons = 0;
  external->i2vExons = ev->nvExons;
  external->ivExons  = ev->nvExons;
  return ev->nvExons;
}
