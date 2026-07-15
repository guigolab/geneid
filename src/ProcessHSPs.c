/*************************************************************************
*                                                                        *
*   Module: ProcessHSPs                                                  *
*                                                                        *
*                                                                        *
*   This file is part of the geneid distribution                         *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                *
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           * 
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/

#include "geneid.h"
#include "bigwig.h"
#ifdef WITH_HTSLIB
#include "bamcov.h"
extern int BAMSTRAND;   /* -y library type: BAMLIB_NONE (unstranded) / RF / FR */
#endif

extern float MRM;
extern int UTR;
extern int SRP;
extern float NO_SCORE;
extern int VRB;         /* -v verbose: gates the Stage-1 coverage-background diagnostic */


/* Projection of HSPs: save the maximum for each nucleotide */
/* Requirement: HSPs must sorted by Position1 */
void HSPScan(packExternalInformation* external,
			 packHSP* hsp, 
			 int Strand, 
			 long l1, long l2)
{
  short x;
  short frameStart, frameEnd;
  long i,j;
  float scoreHSP;
  
  if (Strand == FORWARD)
    {
	  frameStart = 0; 
	  frameEnd = FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = l1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] || 
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] || 
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <= l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			}
		}
	}
  else
	{
	  /* HSPs in REVERSE strand */
	  frameStart = FRAMES; 
	  frameEnd = 2*FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  /* HSPs in reverse strand are reverse-sorted by Position2 */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j < l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos2; 
					  j >= hsp->sPairs[x][i]->Pos1 && j >= l1;
					  j--)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			}
		}
	}
}
/* Additive coverage accumulation, shared by the text RNA-seq path (ReadScan)
   and, later, the frameless bigWig/BAM fill: summed depth is capped at COV,
   while raw read support is tracked separately for the rpkm report. Extracted
   verbatim from ReadScan so the -u -S path stays byte-identical. */
static void CoverAdd(packExternalInformation* external, short x, long idx,
                     float scoreHSP, float rawScore)
{
  if (external->sr[x][idx] == NO_SCORE)
    external->sr[x][idx] = scoreHSP;
  else
    external->sr[x][idx] = MIN(external->sr[x][idx] + scoreHSP, COV);

  /* readcount[] (raw read support, for the rpkm report) is only allocated under
     -u; without it, coverage still fills sr[] to score exons but keeps no rpkm. */
  if (UTR)
    external->readcount[x][idx] = external->readcount[x][idx] + rawScore;
}

/* Projection of RNA-seq reads: accumulate summed depth for each nucleotide */
/* Requirement: HSPs must sorted by Position1 */
void ReadScan(packExternalInformation* external,
			 packHSP* hsp, 
			 int Strand, 
			 long l1, long l2)
{
  short x;
  short frameStart, frameEnd;
  long i,j;
  float scoreHSP;
  
  if (Strand == FORWARD)
    {
	  frameStart = 0; 
	  frameEnd = FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->readcount[x][i] = 0.0;
		  
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
				  /* / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = l1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  CoverAdd(external, x, j-l1, scoreHSP, hsp->sPairs[x][i]->Score);
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; 
/* 				    /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  CoverAdd(external, x, j-l1, scoreHSP, hsp->sPairs[x][i]->Score);
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; 
/* 				    /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <= l2;
					  j++)
					{
					  CoverAdd(external, x, j-l1, scoreHSP, hsp->sPairs[x][i]->Score);
					}
				}
			}
		}
	}
  else
	{
	  /* HSPs in REVERSE strand */
	  frameStart = FRAMES; 
	  frameEnd = 2*FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  /* HSPs in reverse strand are reverse-sorted by Position2 */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->readcount[x][i] = 0.0;
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
/* 				  / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j < l2;
					  j++)
					{
					  CoverAdd(external, x, j-l1, scoreHSP, hsp->sPairs[x][i]->Score);
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
/* 				  /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  CoverAdd(external, x, j-l1, scoreHSP, hsp->sPairs[x][i]->Score);
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
/* 				  / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos2; 
					  j >= hsp->sPairs[x][i]->Pos1 && j >= l1;
					  j--)
					{
					  CoverAdd(external, x, j-l1, scoreHSP, hsp->sPairs[x][i]->Score);
					}
				}
			}
		}
	}
}

/* Preprocessing of HSPs projections */
void HSPScan2(packExternalInformation* external,
			  packHSP* hsp, 
			  int Strand, 
			  long l1, long l2)
{
  short x;
  long i;
  float previousScore;
  float previousReadCount;
  short frameStart, frameEnd;
  
  if (Strand == FORWARD)
    {
	  frameStart = 0; 
	  frameEnd = FRAMES;
	}
  else
	{
	  frameStart = FRAMES; 
	  frameEnd = 2*FRAMES;
	}
  
  for(x=frameStart; x < frameEnd; x++)
    {
      previousScore = 0.0;
      previousReadCount = 0.0;
      /* Screening the whole sequence to accumulate the sum in every base */
      for (i=l1; i<=l2; i++)
		{
		  /* Accumulating step */
		  if (UTR){
		    external->sr[x][i-l1] = previousScore + log(external->sr[x][i-l1] + 1);
		    previousScore = external->sr[x][i-l1];
		    external->readcount[x][i-l1] = previousReadCount + external->readcount[x][i-l1];
		    previousReadCount = external->readcount[x][i-l1];
		  }else{
		    external->sr[x][i-l1] = previousScore + external->sr[x][i-l1];
		    previousScore = external->sr[x][i-l1];
		  }
		}
    }  
}


/* --- RNA-seq expression-scoring redesign, Stage 1 -------------------------- *
 * Estimate a robust background coverage level (lambda_bg) for the current
 * fragment/strand from the raw per-base coverage in sr[], read BEFORE HSPScan2
 * turns sr[] into a prefix sum. We use the MEDIAN of covered positions, not the
 * mean: RNA-seq coverage has a heavy hyper-expressed tail (rRNA, pileups) that
 * inflates the mean ~16x on real data, so a mean-based null would sit far above
 * normal genes and penalise them. The median tracks the typical background/low-
 * expression level and scales with sequencing depth, so deeper data raises the
 * null and the signal together (the null stays calibrated).
 *
 * Stage 1 only REPORTS this value (verbose diagnostic, no scoring change); in
 * Stage 2 it becomes the null of a per-base log-likelihood-ratio coverage term.
 * sr[] holds depth/MRM capped at COV (see CoverAdd), so read depth ~= sr*MRM;
 * we report in depth units. Covered positions are sr[] != NO_SCORE. A strand's
 * three frame planes are identical copies, so we scan only the first. */
static void ReportCoverageBackground(packExternalInformation* external,
                                     int Strand, long l1, long l2)
{
  short frame = (Strand == FORWARD) ? 0 : FRAMES;
  long  len   = l2 - l1 + 1;
  long  i, covered = 0, run, half, medianDepth = 0;
  double lambdaBg;
  /* Depth histogram for the median. Depth is small for the vast majority of
     bases; a fixed cap with an overflow bin keeps this O(len) and allocation
     free. The median is a low quantile, so depths >= HISTCAP (all above it)
     only need counting, not exact binning. */
  enum { HISTCAP = 1024 };
  long hist[HISTCAP + 1];
  char mess[MAXSTRING];

  if (!VRB) return;              /* pure diagnostic; nothing else consumes it yet */

  for (i = 0; i <= HISTCAP; i++) hist[i] = 0;
  for (i = 0; i < len; i++) {
    long depth;
    if (external->sr[frame][i] == NO_SCORE) continue;   /* uncovered */
    depth = (long)(external->sr[frame][i] * MRM + 0.5);
    if (depth < 0) depth = 0;
    if (depth > HISTCAP) depth = HISTCAP;
    hist[depth]++;
    covered++;
  }

  if (covered == 0) {
    sprintf(mess, "Coverage background [%ld-%ld] %s: no covered bases",
            l1, l2, (Strand == FORWARD) ? "fwd" : "rvs");
    printMess(mess);
    return;
  }

  half = (covered + 1) / 2;
  run = 0;
  for (i = 0; i <= HISTCAP; i++) {
    run += hist[i];
    if (run >= half) { medianDepth = i; break; }
  }
  lambdaBg = medianDepth + 1.0;   /* + pseudocount */
  sprintf(mess,
          "Coverage background [%ld-%ld] %s: covered %ld/%ld (%.1f%%), "
          "median depth %ld, lambda_bg %.1f",
          l1, l2, (Strand == FORWARD) ? "fwd" : "rvs",
          covered, len, 100.0 * (double) covered / (double) len,
          medianDepth, lambdaBg);
  printMess(mess);
}

/* Management function to score and filter exons */
void ProcessHSPs(long l1,
                long l2,
                int Strand,
		packExternalInformation* external,
                packHSP* hsp
                )
{

  /* Fill in the temporary HSP arrays (pre-processing) */
  /* GENIS hack */
  if (SRP)
	{
	  if (UTR){
	    printMess("Preprocessing read information: step 1");
	    ReadScan(external,hsp,Strand,l1,l2);
	    ReportCoverageBackground(external, Strand, l1, l2);
	  }else{
	    printMess("Preprocessing homology information: step 1");
	    HSPScan(external,hsp,Strand,l1,l2);
	  }

	  printMess("Preprocessing homology information: step 2");
	  HSPScan2(external,hsp,Strand,l1,l2);
	}
}

/* ------------------------------------------------------------------------- *
 *  bigWig RNA-seq coverage: the same sr[]/readcount[] fill as the text       *
 *  ReadScan path, but sourced from a per-fragment bigWig range query instead *
 *  of the preloaded HSP list. Fills sr[] to score exons whether or not -u is  *
 *  set; readcount[] (rpkm) is only touched under -u (see CoverAdd).          *
 * ------------------------------------------------------------------------- */

/* One coverage interval in the CURRENT strand's coordinate frame (genomic for
   FORWARD, RSequence for REVERSE). */
typedef struct { long s, e; float v; } covIv;

/* Growable buffer + mapping context for the bwQuery callback. */
typedef struct {
  covIv* a;
  long   n, cap;
  int    strand;
  long   L;          /* LengthSequence, for the REVERSE genomic<->RSequence flip */
} covBuf;

/* bwQuery reports genomic intervals [s,e) 0-based half-open; store each in the
   1-based position frame the sr[] fill uses (the text HSP path stores GFF
   1-based coords). FORWARD: 0-based [s,e) -> 1-based half-open [s+1, e+1).
   REVERSE: the manager runs on RSequence, so genomic 1-based p maps to
   RSequence coord L-p+1; genomic 0-based [s,e) = 1-based [s+1,e], which reverses
   to RSequence 1-based [L-e+1, L-s], i.e. half-open [L-e+1, L-s+1). */
static void covCollect(long s, long e, float v, void* ud)
{
  covBuf* b = (covBuf*) ud;
  long rs, re;

  if (b->strand == FORWARD) { rs = s + 1;        re = e + 1; }
  else                      { rs = b->L - e + 1; re = b->L - s + 1; }

  if (b->n == b->cap) {
    b->cap = b->cap ? b->cap * 2 : 64;
    b->a = (covIv*) realloc(b->a, b->cap * sizeof(covIv));
    if (b->a == NULL) printError("Not enough memory: bigWig coverage buffer");
  }
  b->a[b->n].s = rs;
  b->a[b->n].e = re;
  b->a[b->n].v = v;
  b->n++;
}

/* Fill sr[]/readcount[] over fragment [l1,l2] from a frameless coverage stream.
   The signal has no reading frame, so it is replicated identically into all
   three frame planes of the strand (0..2 FWD, 3..5 RVS) -- the same 3-copy
   replication the text path gets on disk (frame '.' in ReadHSP), done here in
   memory. sr[]/readcount[] are per-fragment scratch (reset every call, like
   ReadScan), so fragments never double-count across the OVERLAP band. */
static void FillCoverageFrameless(packExternalInformation* external, int Strand,
                                  long l1, long l2, covIv* iv, long niv)
{
  short frameStart = (Strand == FORWARD) ? 0 : FRAMES;
  short frameEnd   = frameStart + FRAMES;
  short x;
  long  i, k, j;

  for (x = frameStart; x < frameEnd; x++) {
    for (i = 0; i < l2 - l1 + 1; i++) {
      external->sr[x][i] = NO_SCORE;
      if (UTR) external->readcount[x][i] = 0.0;   /* readcount[] is -u-only (see CoverAdd) */
    }
    for (k = 0; k < niv; k++) {
      float scoreHSP = (RREADS / MRM) * iv[k].v;
      long a = (iv[k].s < l1)     ? l1     : iv[k].s;   /* clip to [l1, l2]  */
      long b = (iv[k].e > l2 + 1) ? l2 + 1 : iv[k].e;   /* half-open upper   */
      for (j = a; j < b; j++)
        CoverAdd(external, x, j - l1, scoreHSP, iv[k].v);
    }
  }
}

/* bigWig counterpart of ProcessHSPs (see geneid.h). */
void ProcessCoverageBigWig(long l1, long l2, int Strand,
                           packExternalInformation* external,
                           long LengthSequence)
{
  BigWig* bw = (Strand == FORWARD) ? external->bwPlus : external->bwMinus;
  covBuf buf;
  long gS, gE;   /* genomic half-open query range for this fragment */

  buf.a = NULL; buf.n = 0; buf.cap = 0;
  buf.strand = Strand; buf.L = LengthSequence;

  /* 0-based genomic query range for this fragment, widened by one base each
     side so an interval abutting the fragment edge is still returned; the exact
     1-based mapping + clip to [l1,l2] happens in covCollect/FillCoverageFrameless. */
  if (Strand == FORWARD) {
    gS = l1 - 1;                   /* FWD fill frame (1-based) ~ genomic+1 */
    gE = l2 + 2;
  } else {
    gS = LengthSequence - 2 - l2;  /* RSeq [l1,l2] -> genomic ~[L-1-l2, L-1-l1] */
    gE = LengthSequence - l1 + 1;
  }
  if (gS < 0) gS = 0;

  printMess("Preprocessing bigWig coverage: step 1");
  if (bw != NULL && external->curLocus != NULL)
    bwQuery(bw, external->curLocus, gS, gE, covCollect, &buf);

  FillCoverageFrameless(external, Strand, l1, l2, buf.a, buf.n);
  free(buf.a);
  ReportCoverageBackground(external, Strand, l1, l2);

  printMess("Preprocessing bigWig coverage: step 2");
  HSPScan2(external, NULL, Strand, l1, l2);
}

/* BAM counterpart of ProcessCoverageBigWig (see geneid.h). Identical shape --
   query the current fragment, map/clip via covCollect + FillCoverageFrameless,
   run step 2 -- with bamCoverageQuery in place of bwQuery. Unstranded: both
   strands see the same per-base depth. */
void ProcessCoverageBam(long l1, long l2, int Strand,
                        packExternalInformation* external,
                        long LengthSequence)
{
#ifdef WITH_HTSLIB
  covBuf buf;
  long gS, gE;   /* genomic half-open query range for this fragment */

  buf.a = NULL; buf.n = 0; buf.cap = 0;
  buf.strand = Strand; buf.L = LengthSequence;

  if (Strand == FORWARD) {
    gS = l1 - 1;
    gE = l2 + 2;
  } else {
    gS = LengthSequence - 2 - l2;
    gE = LengthSequence - l1 + 1;
  }
  if (gS < 0) gS = 0;

  /* Unstranded: count every read on both strands. Stranded (-y rf|fr): this
     pass wants only the reads whose transcription strand matches the frame
     being filled (FORWARD -> '+', REVERSE -> '-'). */
  char wantStrand = 0;
  if (BAMSTRAND != BAMLIB_NONE)
    wantStrand = (Strand == FORWARD) ? '+' : '-';

  printMess("Preprocessing BAM coverage: step 1");
  if (external->bam != NULL && external->curLocus != NULL)
    bamCoverageQuery(external->bam, external->curLocus, gS, gE,
                     wantStrand, BAMSTRAND, covCollect, &buf);

  FillCoverageFrameless(external, Strand, l1, l2, buf.a, buf.n);
  free(buf.a);
  ReportCoverageBackground(external, Strand, l1, l2);

  printMess("Preprocessing BAM coverage: step 2");
  HSPScan2(external, NULL, Strand, l1, l2);
#else
  (void) l1; (void) l2; (void) Strand; (void) external; (void) LengthSequence;
  printError("BAM input requires building geneid with WITH_HTSLIB=1");
#endif
}







