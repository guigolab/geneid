/*************************************************************************
*   Module: ReadExonsBigBed                                              *
*                                                                        *
*   Per-split evidence ingest from a bigBed file: for the current        *
*   fragment [l1,l2] of a locus, range-query the bigBed (src/bigbed.c),  *
*   turn each record into an evidence feature via the shared             *
*   AddEvidenceExon commit, and hand the fragment's slice to SortExons   *
*   -- so only the evidence overlapping this split is ever read.         *
*                                                                        *
*   BED record layout (Tyler's convention): chrom/start/end in fields    *
*   1-3, then the tab-joined "rest": name(4)=group, score(5), strand(6), *
*   and an optional type(7) = a geneid feature (Intron / First /         *
*   Internal / Terminal / Single). No BED12 blocks. Coordinates are      *
*   0-based half-open -> geneid 1-based (start+1 .. end).                *
*                                                                        *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#include "geneid.h"
#include "bigbed.h"

extern float EvidenceFactor;
extern float EvidenceEW;
extern int FWD, RVS;

/* One record collected from a range query, before grouping/commit. */
typedef struct { long start, end; char rest[MAXLINE]; } bbRec;

typedef struct {
  bbRec* v;
  long   n, cap;
} bbRecList;

static void collectCB(long start, long end, const char* rest, void* ud){
  bbRecList* L = (bbRecList*)ud;
  if (L->n >= L->cap){
    L->cap = L->cap ? L->cap * 2 : 256;
    L->v = (bbRec*)realloc(L->v, L->cap * sizeof(bbRec));
    if (!L->v) printError("Not enough memory: bigBed record buffer");
  }
  L->v[L->n].start = start;
  L->v[L->n].end   = end;
  strncpy(L->v[L->n].rest, rest, MAXLINE - 1);
  L->v[L->n].rest[MAXLINE - 1] = '\0';
  L->n++;
}

static int cmpByStart(const void* a, const void* b){
  long x = ((const bbRec*)a)->start, y = ((const bbRec*)b)->start;
  return (x > y) - (x < y);
}

/* Populate external->evidence[0] with the evidence features whose acceptor
   (start+1) falls in this fragment's OWNED range (ownedLo, ownedHi] -- the same
   non-overlapping per-fragment assignment SearchEvidenceExons makes for GFF, so
   an evidence record is committed in exactly one fragment. Returns the number of
   evidence exon slots produced (including frame/intron replicas). */
long ReadExonsBigBed(BigBed* bb, packExternalInformation* external, dict* d,
                     char* Locus, long l1, long l2, long ownedLo, long ownedHi){
  packEvidence* ev = external->evidence[0];
  bbRecList L = { NULL, 0, 0 };
  long lastAcceptor = -INFI;
  long i;

  ev->nvExons = 0;
  ev->nvSites = 0;

  /* l1/l2 are 0-based sequence indices (see geneid.c); bigBed is 0-based too,
     so the fragment maps straight through. Ownership (below) is in 1-based
     acceptor units to match GFF's SearchEvidenceExons. */
  if (bbQuery(bb, Locus, l1, l2, collectCB, &L) < 0)
    printError("bigBed range query failed");

  qsort(L.v, L.n, sizeof(bbRec), cmpByStart);

  for (i = 0; i < L.n; i++){
    long begin = L.v[i].start + 1;          /* 0-based -> 1-based acceptor */
    long end   = L.v[i].end;                 /* 0-based half-open end == 1-based inclusive */
    char* name; char* sScore; char* sStrand; char* sType;
    char strand; char* type; float score; int three;

    if (begin <= ownedLo || begin > ownedHi) continue;  /* owned by another fragment */

    /* Split the tab-joined BED tail: name score strand [type] */
    name    = strtok(L.v[i].rest, "\t");
    sScore  = strtok(NULL, "\t");
    sStrand = strtok(NULL, "\t");
    sType   = strtok(NULL, "\t\n");
    if (!sStrand) continue;                  /* need at least name/score/strand */
    strand = sStrand[0];
    if (strand != '+' && strand != '-') continue;
    /* Only the requested strand pass(es) of this fragment care about it */
    if ((strand == '+' && !FWD) || (strand == '-' && !RVS)) continue;

    type = (sType && sType[0]) ? sType : sEXON;   /* generic Exon expansion is TODO */
    /* score: '.' or missing -> MAXSCORE; else scaled like the GFF path */
    if (!sScore || !strcmp(sScore, ".") || sscanf(sScore, "%f", &score) != 1)
      score = MAXSCORE;
    else
      score = score * EvidenceFactor + EvidenceEW;

    /* Fill the staging slot, then let the shared commit validate + expand it. */
    strcpy((ev->vExons + ev->nvExons)->Type, type);
    (ev->vExons + ev->nvExons)->Score  = score;
    (ev->vExons + ev->nvExons)->Strand = strand;
    if (name && strcmp(name, ".") && strcmp(name, ""))
      strcpy((ev->vExons + ev->nvExons)->Group, name);
    else
      strcpy((ev->vExons + ev->nvExons)->Group, NOGROUP);
    (ev->vSites + ev->nvSites)->Position     = begin;
    (ev->vSites + ev->nvSites + 1)->Position = end;

    three = 1;                                /* bigBed carries no frame -> 3 copies */
    AddEvidenceExon(external, 0, d, &three, U2, U2, &lastAcceptor, L.v[i].rest);
  }

  free(L.v);
  external->i1vExons = 0;
  external->i2vExons = ev->nvExons;
  external->ivExons  = ev->nvExons;
  return ev->nvExons;
}
