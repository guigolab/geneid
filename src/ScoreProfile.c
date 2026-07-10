/*************************************************************************
*                                                                        *
*   Module: ScoreProfile                                                 *
*                                                                        *
*   Score a splice/signal profile at ONE given position, replicating the *
*   order-0/1/2 PWM/Markov math of BuildDonors.c / BuildAcceptors.c but  *
*   for a single site instead of a range scan. Used by the annotation-   *
*   scoring mode (-J) to score a FORCED evidence splice site under each  *
*   enabled profile (so its donor/acceptor score and intron subtype can  *
*   be reported), WITHOUT re-running or altering the ab-initio builders.  *
*                                                                        *
*   This file is part of the geneid distribution                         *
*************************************************************************/

#include "geneid.h"

/* Function TRANS: char -> integer such that A=0, C=1, G=2 and T/U=3 (N=4) */
extern int TRANS[];

/* Which optional splice profiles the param actually populated. The gparam struct
   is malloc'd so every profile pointer is non-NULL, but only the ones flagged
   here hold real values -- gate on these, exactly like manager.c's builders. */
extern int U12GTAG, U12ATAC, U2GCAG, U2GTA, U2GTG, U2GTY;

/* U12 branch-point scorer (defined in BuildU12Acceptors.c). Scans upstream of the
   acceptor at `positionAcc` (in the frame of `s`), fills splicesite->PositionBP
   (offset relative to the acceptor), and returns the best branch-point score. */
extern float ComputeU12BranchProfile(char* s, long positionAcc, long limitRight,
                                     profile* p, site* splicesite);

/* Raw profile score at signal position `pos` (0-based index into `s`, the same
   frame BuildDonors/BuildAcceptors store in site->Position). The scored window
   starts `p->offset` bases before the signal, exactly as in the builders:
   there a site's Position is left+is+p->offset while the window base is
   sequence+left+is, i.e. Position-p->offset. Returns afactor+bfactor*loglik.
   The UTR PeakEdgeScore term used by the builders is deliberately NOT applied
   here -- this is a pure signal score for reporting. Callers must guarantee the
   window [pos-p->offset, pos-p->offset+p->dimension) lies inside the loaded
   sequence (annotation-scoring guards fragment edges before calling). */
/* Raw log-likelihood sum of profile p at signal position `pos`, WITHOUT the
   afactor/bfactor rescaling (i.e. the builders' inner `score` before line
   `score = afactor + bfactor*score`). This is what the acceptor path reports as
   acceptor_profile_score, and what the branch-point contribution is added to
   before rescaling. ScoreProfileAt() applies the rescaling on top. */
static float ScoreProfileRawAt(char* s, long pos, profile* p)
{
  int i, index;
  float score = 0.0;
  char* w = s + (pos - p->offset);       /* window base (0-based) */

  if (p->order == 0)
    {
      for (i = 0; i < p->dimension; i++)
	{
	  index = TRANS[(int)(*(w + i))];
	  score += (index >= p->dimensionTrans) ? -(float)INFI
	                                        : p->transitionValues[i][index];
	}
    }
  else if (p->order == 1)
    {
      for (i = 0; i < p->dimension; i++)
	{
	  index = 5 * TRANS[(int)(*(w + i - 1))] + TRANS[(int)(*(w + i))];
	  score += (index >= p->dimensionTrans) ? -(float)INFI
	                                        : p->transitionValues[i][index];
	}
    }
  else
    {
      for (i = 0; i < p->dimension; i++)
	{
	  index = OligoToInt(w + i - p->order, p->order + 1, 5);
	  score += (index >= p->dimensionTrans) ? -(float)INFI
	                                        : p->transitionValues[i][index];
	}
    }

  return score;
}

float ScoreProfileAt(char* s, long pos, profile* p)
{
  return p->afactor + (p->bfactor * ScoreProfileRawAt(s, pos, p));
}

/* True if profile p's scoring window for a signal at 0-based index `pos` lies
   fully inside a sequence of length L (window = [pos-offset, pos-offset+dim)). */
static int windowInBounds(long pos, profile* p, long L)
{
  long base = pos - p->offset;
  return (base >= 0) && (base + p->dimension <= L);
}

/* Score profile p at pos, or -INFI if p is absent or its window is off-sequence.
   A profile whose invariant-dinucleotide mask rejects the site returns a hugely
   negative score, so "best across the enabled profiles" recovers the site's
   dinucleotide class automatically (a GC donor only scores under U2gcag, etc.). */
static float scoreGuarded(char* s, long pos, profile* p, long L, int enabled)
{
  if (!enabled || !p || !windowInBounds(pos, p, L)) return -(float)INFI;
  return ScoreProfileAt(s, pos, p);
}

/* Classification thresholds (log-odds units). U12 is chosen over U2 only when it
   both clears FLOOR and beats the best U2 profile by MARGIN -- this curbs U12
   over-calling on GT-AG donors (GC-AG stays U2; AT-AC has no U2 competitor so it
   is U12 whenever the U12 profile matches). PROFILE_MATCH is the score above
   which a profile's dinucleotide mask has accepted the site. Defaults are
   neutral (U12 must beat U2 and score >= 0); a param may tune them later. */
#define U12_CLASSIFY_FLOOR   0.0f
#define U12_CLASSIFY_MARGIN  0.0f
#define PROFILE_MATCH        (-500.0f)

/* Classify a donor at 0-based index `pos` in `s` under every enabled+present
   donor profile, and set the site's Score/subtype/class. The intron's U2-vs-U12
   type is read from the donor subtype downstream (PrintGIntron). */
static void ClassifyDonor(char* s, long pos, gparam* gp, long L, site* d)
{
  float bU2 = -(float)INFI; const char* bU2sub = sU2; float sc;
  sc = scoreGuarded(s,pos,gp->DonorProfile,L,1);           if (sc > bU2){ bU2=sc; bU2sub=sU2;    }
  sc = scoreGuarded(s,pos,gp->U2gtaDonorProfile,L,U2GTA);  if (sc > bU2){ bU2=sc; bU2sub=sU2gta; }
  sc = scoreGuarded(s,pos,gp->U2gtgDonorProfile,L,U2GTG);  if (sc > bU2){ bU2=sc; bU2sub=sU2gtg; }
  sc = scoreGuarded(s,pos,gp->U2gtyDonorProfile,L,U2GTY);  if (sc > bU2){ bU2=sc; bU2sub=sU2gty; }
  sc = scoreGuarded(s,pos,gp->U2gcagDonorProfile,L,U2GCAG);if (sc > bU2){ bU2=sc; bU2sub=sU2gcag;}

  float bU12 = -(float)INFI; const char* bU12sub = sU2; short bU12cls = U2;
  sc = scoreGuarded(s,pos,gp->U12gtagDonorProfile,L,U12GTAG); if (sc > bU12){ bU12=sc; bU12sub=sU12gtag; bU12cls=U12gtag; }
  sc = scoreGuarded(s,pos,gp->U12atacDonorProfile,L,U12ATAC); if (sc > bU12){ bU12=sc; bU12sub=sU12atac; bU12cls=U12atac; }

  int u2ok = (bU2 > PROFILE_MATCH), u12ok = (bU12 > PROFILE_MATCH);
  if (u12ok && (!u2ok || (bU12 >= U12_CLASSIFY_FLOOR && bU12 >= bU2 + U12_CLASSIFY_MARGIN)))
    { d->Score = bU12; strcpy(d->subtype, bU12sub); d->class = bU12cls; }
  else if (u2ok)
    { d->Score = bU2;  strcpy(d->subtype, bU2sub);  d->class = U2; }
  else /* non-canonical: report the generic donor score, keep U2 label */
    { d->Score = scoreGuarded(s,pos,gp->DonorProfile,L,1); strcpy(d->subtype, sU2); d->class = U2; }
}

/* Classify an acceptor at `pos` under the U2 and (if present) U12 acceptor
   profiles. Sets ScoreAccProfile + Score (BP/PPT added separately) + subtype +
   class. The U12 acceptor PWM alone is used here; its branch-point contribution
   is layered on later. */
/* Classify an acceptor at `pos`. `forceU12cls` lets the caller impose the U12
   subclass of the intron this acceptor closes: ab initio types an intron by its
   DONOR (PrintGIntron reads the donor subtype) and pairs it with a U12 acceptor
   even when the acceptor's own U12 signal is weak, so when the paired donor is
   U12 we score the acceptor under that same U12 profile rather than re-deciding
   U2-vs-U12 from the acceptor alone. forceU12cls == U2 means "decide from the
   acceptor" (used when the paired donor is U2 or unknown). */
static void ClassifyAcceptor(char* s, long pos, gparam* gp, long L, site* a,
                             short forceU12cls)
{
  float bU2 = scoreGuarded(s,pos,gp->AcceptorProfile,L,1); float sc;
  float bU12 = -(float)INFI; const char* bU12sub = sU2; short bU12cls = U2;
  profile* bU12prof = NULL;
  sc = scoreGuarded(s,pos,gp->U12gtagAcceptorProfile,L,U12GTAG); if (sc > bU12){ bU12=sc; bU12sub=sU12gtag; bU12cls=U12gtag; bU12prof=gp->U12gtagAcceptorProfile; }
  sc = scoreGuarded(s,pos,gp->U12atacAcceptorProfile,L,U12ATAC); if (sc > bU12){ bU12=sc; bU12sub=sU12atac; bU12cls=U12atac; bU12prof=gp->U12atacAcceptorProfile; }

  /* Donor-driven override: score under the U12 profile matching the paired donor. */
  if (forceU12cls == U12gtag) { bU12sub = sU12gtag; bU12cls = U12gtag; bU12prof = gp->U12gtagAcceptorProfile; }
  else if (forceU12cls == U12atac) { bU12sub = sU12atac; bU12cls = U12atac; bU12prof = gp->U12atacAcceptorProfile; }

  int u2ok = (bU2 > PROFILE_MATCH), u12ok = (bU12 > PROFILE_MATCH);
  int forced = (forceU12cls == U12gtag || forceU12cls == U12atac) && bU12prof;
  a->ScoreBP = 0.0; a->PositionBP = 0; a->ScorePPT = 0.0; a->PositionPPT = 0;
  /* Auto classification (U12 vs U2) uses the profile-only scores (bU2/bU12),
     matching the validated per-site call; a forced U12 class wins outright. */
  if (forced || (u12ok && (!u2ok || (bU12 >= U12_CLASSIFY_FLOOR && bU12 >= bU2 + U12_CLASSIFY_MARGIN))))
    {
      /* U12 acceptor: report the same score composition BuildU12Acceptors uses --
         acceptor_profile_score = raw profile sum (pre-afactor), plus a U12 branch
         point, combined through the profile's afactor/bfactor:
           acceptor_score = afactor + bfactor*(raw + bp_score).                 */
      float raw = ScoreProfileRawAt(s, pos, bU12prof);
      float bp  = ComputeU12BranchProfile(s, pos, L, gp->U12BranchPointProfile, a);
      a->ScoreAccProfile = raw;
      a->ScoreBP = bp;
      a->Score = bU12prof->afactor + bU12prof->bfactor * (raw + bp);
      strcpy(a->subtype, bU12sub);
      a->class = bU12cls;
    }
  else
    { /* U2 acceptor: profile score only (branch point off for U2 in these params). */
      a->ScoreAccProfile = bU2;
      a->Score = bU2;
      strcpy(a->subtype, sU2);
      a->class = U2;
    }
}

/* Score+classify ONE forced evidence exon's real splice boundary(ies), filling
   the REAL donor/acceptor profile scores on its (otherwise dummy) evidence sites
   so they report meaningful donor_score/acceptor_score and an intron subtype.
   Report-only: never feeds the DP (the annotation is forced regardless).

   An evidence exon labels its two sites GENOMICALLY: Acceptor = the left (lower-
   coordinate) boundary, Donor = the right. Which boundary is a real splice
   DONOR vs ACCEPTOR depends on strand and exon type (transcription orientation):
     First    -> transcription donor only
     Internal -> both
     Terminal -> transcription acceptor only
     Single   -> neither
   On '+' the transcription donor is the genomic-right site (Donor field) scored
   on the forward Sequence at Position-COFFSET; on '-' it is the genomic-left
   site (Acceptor field) scored on RSequence at L-Position (roles swap). The
   filled site objects are exactly the ones PrintExons/PrintGIntron read back
   (they apply the same strand swap when printing). gp is the (single) isochore. */
/* Pass 1 of the join: score the transcription DONOR and START of one exon (the
   sites a downstream exon's acceptor is paired against). On '+' the donor is the
   genomic-right site (Donor field) and the start the genomic-left (Acceptor); on
   '-' the roles swap onto RSequence. */
static void ScoreEvidenceDonorStart(exonGFF* e, char* Sequence, char* RSequence,
                                    long L, gparam* gp)
{
  int isFirst    = !strcmp(e->Type, sFIRST);
  int isInternal = !strcmp(e->Type, sINTERNAL);
  int isSingle   = !strcmp(e->Type, sSINGLE);
  int hasDonor    = isFirst || isInternal;      /* transcription donor    */
  int hasStart    = isFirst || isSingle;        /* transcription start (ATG) */
  site* leftSite  = e->Acceptor;                /* genomic-left  boundary */
  site* rightSite = e->Donor;                   /* genomic-right boundary */

  if (e->Strand == '+')
    {
      if (hasDonor) ClassifyDonor(Sequence, rightSite->Position - COFFSET, gp, L, rightSite);
      if (hasStart) leftSite->Score = scoreGuarded(Sequence, leftSite->Position - COFFSET, gp->StartProfile, L, 1);
    }
  else
    {
      if (hasDonor) ClassifyDonor(RSequence, L - leftSite->Position, gp, L, leftSite);
      if (hasStart) rightSite->Score = scoreGuarded(RSequence, L - rightSite->Position, gp->StartProfile, L, 1);
    }
}

/* The class of the donor that pairs with e's acceptor (the intron e closes): it
   lives on e's transcription-UPSTREAM neighbour e->PreviousExon, on the field
   that ScoreEvidenceDonorStart wrote -- Donor on '+', Acceptor on '-'. Returns U2
   when there is no real evidence predecessor. */
static short PairedDonorClass(exonGFF* e)
{
  exonGFF* prev = e->PreviousExon;
  if (!prev || prev->Strand == '*' || !prev->evidence) return U2;
  return (e->Strand == '+') ? prev->Donor->class : prev->Acceptor->class;
}

/* Pass 2 of the join: score the transcription ACCEPTOR of one exon, told the
   paired donor's class so a U12 intron gets a U12 acceptor (donor-driven). */
static void ScoreEvidenceAcceptor(exonGFF* e, short pairedDonorCls,
                                  char* Sequence, char* RSequence, long L, gparam* gp)
{
  int isInternal = !strcmp(e->Type, sINTERNAL);
  int isTerminal = !strcmp(e->Type, sTERMINAL);
  int hasAcceptor = isInternal || isTerminal;   /* transcription acceptor */
  if (!hasAcceptor) return;

  if (e->Strand == '+')
    ClassifyAcceptor(Sequence, e->Acceptor->Position - COFFSET, gp, L, e->Acceptor, pairedDonorCls);
  else /* '-' : transcription acceptor is the genomic-right site, on RSequence */
    ClassifyAcceptor(RSequence, L - e->Donor->Position, gp, L, e->Donor, pairedDonorCls);
}

/* Annotation-scoring mode (-J): score+classify the evidence splice sites of the
   OPTIMAL gene chain that is about to be printed. `chain` is genes->GOptim; it is
   walked backward via PreviousExon to the Ghost sentinel (Strand '*').

   Walking the PRINTED chain (not the loaded evidence->vExons) is what makes this
   correct on MULTI-SPLIT sequences: there BackupGenes/backupExon deep-copies each
   selected exon AND its sites into the dumpster BEFORE this runs, and the printed
   genes reference those copies -- so classifying the originals would never reach
   the output. The copies retain the evidence flag and genomic site Positions, so
   scoring them here fills the exact objects PrintExons/PrintGIntron read back. On
   a single-split sequence GOptim points straight at the evidence exons, so the
   same walk works. Only forced-evidence exons are (re)scored; ab-initio exons in
   a mixed -R solution keep their builder-assigned scores. gp = single isochore. */
void ScoreEvidenceSites(exonGFF* chain, char* Sequence, char* RSequence,
                        long L, gparam* gp)
{
  exonGFF* e;
  /* Two passes so an acceptor can be scored with its paired donor's class known:
     pass 1 fills every donor (+ start); pass 2 scores acceptors join-aware. */
  for (e = chain; e != NULL && e->Strand != '*'; e = e->PreviousExon)
    if (e->evidence)
      ScoreEvidenceDonorStart(e, Sequence, RSequence, L, gp);
  for (e = chain; e != NULL && e->Strand != '*'; e = e->PreviousExon)
    if (e->evidence)
      ScoreEvidenceAcceptor(e, PairedDonorClass(e), Sequence, RSequence, L, gp);
}
