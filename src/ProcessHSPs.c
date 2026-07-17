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
#include <float.h>
#include "bigwig.h"
#ifdef WITH_HTSLIB
#include "bamcov.h"
extern int BAMSTRAND;   /* -y library type: BAMLIB_NONE (unstranded) / RF / FR */
#endif

extern int UTR;
extern int SRP;
extern float NO_SCORE;
extern int VRB;         /* -v verbose: gates the Stage-1 coverage-background diagnostic */
extern int EXPRLLR;     /* -L: Poisson per-base expression LLR coverage scoring */
extern float LLRK, LLRW;/* -L fold-change k (>1); -Q weight/scale of the LLR term */
extern int LLRMINLIBS;  /* -K: libraries that must support a base (order statistic) */


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
				  scoreHSP = (RREADS/COVNORM) * hsp->sPairs[x][i]->Score;
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
				  scoreHSP = (RREADS/COVNORM) * hsp->sPairs[x][i]->Score; 
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
				  scoreHSP = (RREADS/COVNORM) * hsp->sPairs[x][i]->Score; 
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
				  scoreHSP = (RREADS/COVNORM) * hsp->sPairs[x][i]->Score;
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
				  scoreHSP = (RREADS/COVNORM) * hsp->sPairs[x][i]->Score;
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
				  scoreHSP = (RREADS/COVNORM) * hsp->sPairs[x][i]->Score;
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
  /* -L expression LLR: replace the legacy per-base coverage term with a Poisson
     two-state log-likelihood ratio, covered depth vs the background rate
     covLambdaBg, expressed in units of that background:

         term = LLRW * ( (depth/lambda_bg)*log(k) - (k-1) )

     i.e. the Poisson LLR (depth*log k - (k-1)*lambda_bg) DIVIDED BY lambda_bg.
     The division is what makes the term depth-invariant: the raw LLR is linear
     in depth, so a 2x deeper library would silently double the coverage term's
     weight against the (depth-independent) coding and site scores. Dividing by
     lambda_bg leaves the term a function of fold-enrichment over background
     only, which is a property of the transcript, not of how deep we sequenced.
     (This only holds because lambda_bg is itself depth-linear -- see the
     global trimmed-mean estimator in SetCoverageBackground.)

     It crosses zero at depth/lambda_bg = (k-1)/log k, a pure enrichment
     threshold: 1.44x background for k=2, 1.82x for k=3. Uncovered bases score
     -LLRW*(k-1), a bounded penalty rather than one that grows with depth.

     Gated on a valid covLambdaBg (>0), so the protein-homology path and
     coverage-free fragments keep the legacy term -> default off (-L absent) is
     byte-identical. */
  /* covLambdaBg > 0 means the coverage fill already wrote finished per-base LLR
     values into sr[] (FillCoverageLLR): the transform has to happen per library,
     before the MAX, so it cannot be done here. All that is left is the running
     sum. Legacy paths (protein homology, text HSP, -L absent) keep their
     per-base term below, so the default stays byte-identical. */
  int preLLR = (EXPRLLR && external->covLambdaBg > 0.0);

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
		  if (preLLR){
		    /* sr[] already holds this base's LLR: accumulate only. */
		    external->sr[x][i-l1] = previousScore + external->sr[x][i-l1];
		    previousScore = external->sr[x][i-l1];
		    if (UTR){
		      external->readcount[x][i-l1] = previousReadCount + external->readcount[x][i-l1];
		      previousReadCount = external->readcount[x][i-l1];
		    }
		  }else if (UTR){
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

/* --- RNA-seq expression-scoring redesign: global background (lambda_bg) ---- *
 * lambda_bg is the null of the per-base LLR term (-L): the coverage level an
 * unexpressed base is expected to show. It is estimated ONCE PER SEQUENCE AND
 * STRAND by sampling evenly spaced windows and taking a trimmed mean of the
 * per-base depth, then cached (covLambdaGlobal/covLambdaLocus).
 *
 * WHY GLOBAL, NOT PER FRAGMENT. Two biases pull in opposite directions and
 * cannot both be avoided within one fragment:
 *   - Anything conditioned on "covered" (e.g. the median of covered positions)
 *     is DEPTH-biased: as depth falls, low-coverage bases drop below 1 and leave
 *     the covered set, so the statistic falls only ~sqrt(depth). Measured on
 *     pancreas chr21: 6.0 -> 4.0 -> 3.0 across a 4x depth cut, which hardened the
 *     effective enrichment threshold ~2x on shallow data.
 *   - Anything over ALL positions is DENSITY-biased per fragment: most of a
 *     chromosome is gene desert, so the mean collapses to ~0 there and every
 *     covered base then looks enormously enriched. Measured per fragment on
 *     chr21: median 0.010, max 20.1 -- a ~400x swing in the null.
 * Estimating over the whole sequence removes the density term by construction
 * (density is then a fixed property of the sequence) and leaves a mean that is
 * exactly linear in depth -- the property the LLR needs, since it compares
 * coverage to lambda_bg as a RATIO. Absolute accuracy barely matters: a constant
 * multiplicative bias is absorbed by the fold-change k (the term crosses zero at
 * c/lambda_bg = (k-1)/log k), which is why ~1% sampling error is irrelevant.
 * Trimming the top LLR_TRIM fraction drops the hyper-expressed tail (rRNA,
 * pileups) and the most-expressed genes, none of which are background.
 *
 * Depth values come from the same range-query machinery the per-fragment fill
 * uses, so this works for BAM and bigWig alike. The text HSP path has no range
 * query and so cannot supply a global null: there the LLR stays off. */

#define BG_HISTCAP 1024      /* depth histogram cap; the trim discards the tail anyway */

/* Depth histogram accumulator for the background sampler. Intervals arrive as
   [s,e) with a constant depth v, so each contributes (e-s) positions. */
typedef struct { long* hist; long covered; } bgBuf;

static void bgCollect(long s, long e, float v, void* ud)
{
  bgBuf* b = (bgBuf*) ud;
  long d = (long)(v + 0.5);
  long n = e - s;

  if (n <= 0) return;
  if (d < 0) d = 0;
  if (d > BG_HISTCAP) d = BG_HISTCAP;
  b->hist[d] += n;
  b->covered += n;
}

/* Trimmed-mean depth over evenly spaced sample windows of this sequence/strand. */
static float SampleCoverageBackground(packExternalInformation* external,
                                      BamCov* bam, int Strand, long seqLen)
{
  enum { NWIN = 100, WINLEN = 10000 };
  long hist[BG_HISTCAP + 1];
  bgBuf b;
  long i, sampled = 0, keep, run, cap, kept, sum;
  double lam;

  if (external->curLocus == NULL || seqLen <= 0) return -1.0;

  for (i = 0; i <= BG_HISTCAP; i++) hist[i] = 0;
  b.hist = hist; b.covered = 0;

  for (i = 0; i < NWIN; i++) {
    long gS = (long)((double) seqLen * (double) i / (double) NWIN);
    long gE = gS + WINLEN;
    long before = b.covered;

    if (gE > seqLen) gE = seqLen;
    if (gE <= gS) continue;

    if (bam != NULL) {
#ifdef WITH_HTSLIB
      /* Stranded libraries get a per-strand null; unstranded sees every read. */
      char wantStrand = 0;
      if (BAMSTRAND != BAMLIB_NONE)
        wantStrand = (Strand == FORWARD) ? '+' : '-';
      bamCoverageQuery(bam, external->curLocus, gS, gE,
                       wantStrand, BAMSTRAND, bgCollect, &b);
#endif
    } else {
      BigWig* bw = (Strand == FORWARD) ? external->bwPlus : external->bwMinus;
      if (bw == NULL) return -1.0;
      bwQuery(bw, external->curLocus, gS, gE, bgCollect, &b);
    }

    sampled += (gE - gS);
    hist[0] += (gE - gS) - (b.covered - before);   /* the window's uncovered bases */
  }

  if (sampled <= 0) return -1.0;

  /* Depth cap = the (1 - LLR_TRIM) quantile over the sampled positions. */
  keep = (long)((1.0 - LLR_TRIM) * (double) sampled);
  if (keep < 1) keep = 1;
  run = 0; cap = BG_HISTCAP;
  for (i = 0; i <= BG_HISTCAP; i++) {
    run += hist[i];
    if (run >= keep) { cap = i; break; }
  }

  sum = 0; kept = 0;
  for (i = 0; i <= cap; i++) { sum += i * hist[i]; kept += hist[i]; }
  lam = kept ? (double) sum / (double) kept : 0.0;
  return (float) lam;
}

/* Fill this sequence/strand's per-LIBRARY background cache (covLambdaGlobal[si][b]),
   computing it on first use for the sequence. Each library keeps its OWN null:
   lambda_bg is tissue biology, not depth -- across 5 human total-RNA libraries it
   spans 14x per read -- so one pooled null would be wrong for every library in the
   pool (see MAXBAMS in geneid.h). */
static void SetCoverageBackground(packExternalInformation* external,
                                  int Strand, long seqLen)
{
  int si = (Strand == FORWARD) ? 0 : 1;
  int b, nlib;
  char mess[MAXSTRING];

  external->covLambdaBg = -1.0;
  if (!EXPRLLR && !VRB) return;
  if (external->curLocus == NULL) return;

  /* New sequence -> drop the cached values. */
  if (strcmp(external->covLambdaLocus, external->curLocus) != 0) {
    int s2;
    for (s2 = 0; s2 < 2; s2++)
      for (b = 0; b < MAXBAMS; b++)
        external->covLambdaGlobal[s2][b] = -1.0;
    strncpy(external->covLambdaLocus, external->curLocus, MAXSTRING - 1);
    external->covLambdaLocus[MAXSTRING - 1] = '\0';
  }

  /* bigWig (or the single-bigWig path) has no library array: slot 0, bam NULL. */
  nlib = (external->nBams > 0) ? external->nBams : 1;
  for (b = 0; b < nlib; b++) {
    if (external->covLambdaGlobal[si][b] >= 0.0) continue;
    external->covLambdaGlobal[si][b] =
      SampleCoverageBackground(external,
                               (external->nBams > 0) ? external->bams[b] : NULL,
                               Strand, seqLen);
    if (VRB) {
      sprintf(mess, "Coverage background %s %s lib %d/%d: lambda_bg %.4f (global, sampled)",
              external->curLocus, (Strand == FORWARD) ? "fwd" : "rvs",
              b + 1, nlib, external->covLambdaGlobal[si][b]);
      printMess(mess);
    }
  }

  /* Legacy single-source path keeps using covLambdaBg directly. */
  if (nlib == 1 && external->covLambdaGlobal[si][0] >= LLR_MINLAMBDA)
    external->covLambdaBg = external->covLambdaGlobal[si][0];
}

/* Query ONE library's coverage for this fragment as per-base DEPTH into dst[]
   (0 where uncovered), saturating exactly like CoverAdd does (COV*COVNORM). */
static void QueryLibDepth(packExternalInformation* external, BamCov* bam, int Strand,
                          long l1, long l2, long LengthSequence, float* dst)
{
  covBuf buf;
  long gS, gE, i, k, j, len = l2 - l1 + 1;

  buf.a = NULL; buf.n = 0; buf.cap = 0;
  buf.strand = Strand; buf.L = LengthSequence;

  if (Strand == FORWARD) { gS = l1 - 1; gE = l2 + 2; }
  else                   { gS = LengthSequence - 2 - l2; gE = LengthSequence - l1 + 1; }
  if (gS < 0) gS = 0;

  for (i = 0; i < len; i++) dst[i] = 0.0;

  if (bam != NULL) {
#ifdef WITH_HTSLIB
    char wantStrand = 0;
    if (BAMSTRAND != BAMLIB_NONE)
      wantStrand = (Strand == FORWARD) ? '+' : '-';
    if (external->curLocus != NULL)
      bamCoverageQuery(bam, external->curLocus, gS, gE,
                       wantStrand, BAMSTRAND, covCollect, &buf);
#endif
  } else {
    BigWig* bw = (Strand == FORWARD) ? external->bwPlus : external->bwMinus;
    if (bw != NULL && external->curLocus != NULL)
      bwQuery(bw, external->curLocus, gS, gE, covCollect, &buf);
  }

  for (k = 0; k < buf.n; k++) {
    long a = (buf.a[k].s < l1)     ? l1     : buf.a[k].s;
    long z = (buf.a[k].e > l2 + 1) ? l2 + 1 : buf.a[k].e;
    for (j = a; j < z; j++) {
      float d = dst[j - l1] + buf.a[k].v;
      dst[j - l1] = (d > (float)(COV * COVNORM)) ? (float)(COV * COVNORM) : d;
    }
  }
  free(buf.a);
}

/* -L fill: write the per-base LLR straight into sr[], combining libraries by MAX.
 *
 * The LLR must be computed PER LIBRARY and only then combined, because each
 * library has its own lambda_bg: you cannot max (or sum) raw depths across
 * libraries whose backgrounds differ 14x. MAX gives union semantics -- "expressed
 * in ANY tissue", which is what annotation wants -- and makes a redundant library
 * NEUTRAL (max of a duplicate changes nothing) where merging made it actively
 * harmful (measured: a 2nd brain region adding 0 new genes still cost gSN
 * .178->.168 purely by contributing background). Summing LLRs instead would be
 * AND-ish and would kill exactly the tissue-specific genes extra tissues are for.
 *
 * sr[] then holds the finished per-base term, so HSPScan2 only prefix-sums it.
 * Returns 0 when no library has usable signal here (caller keeps the legacy term). */
static int FillCoverageLLR(packExternalInformation* external, int Strand,
                           long l1, long l2, long LengthSequence)
{
  short frameStart = (Strand == FORWARD) ? 0 : FRAMES;
  short frameEnd   = frameStart + FRAMES;
  int   si   = (Strand == FORWARD) ? 0 : 1;
  long  len  = l2 - l1 + 1;
  double logk = log((double) LLRK);
  double km1  = (double) LLRK - 1.0;
  int    nlib = (external->nBams > 0) ? external->nBams : 1;
  int    b, nvalid = 0, useSecond;
  long   i;
  short  x;

  if (UTR)
    for (i = 0; i < len; i++) external->readcount[frameStart][i] = 0.0;

  /* Track the best (covComb) and 2nd best (covComb2) per-base LLR across the
     libraries, so -K can pick which order statistic to score with. */
  for (i = 0; i < len; i++) {
    external->covComb[i]  = -FLT_MAX;
    external->covComb2[i] = -FLT_MAX;
  }

  for (b = 0; b < nlib; b++) {
    float lam = external->covLambdaGlobal[si][b];
    if (lam < LLR_MINLAMBDA) continue;          /* this library has nothing usable here */

    QueryLibDepth(external, (external->nBams > 0) ? external->bams[b] : NULL,
                  Strand, l1, l2, LengthSequence, external->covTmp);

    for (i = 0; i < len; i++) {
      float llr = (float) ((double) LLRW *
                  (((double) external->covTmp[i] / (double) lam) * logk - km1));
      if (llr > external->covComb[i]) {
        external->covComb2[i] = external->covComb[i];
        external->covComb[i]  = llr;
      } else if (llr > external->covComb2[i]) {
        external->covComb2[i] = llr;
      }
      if (UTR) external->readcount[frameStart][i] += external->covTmp[i];
    }
    nvalid++;
  }

  if (nvalid == 0) return 0;                    /* no usable library */

  /* -K n: score with the n-th highest library. With fewer libraries than -K, fall
     back to the max, so one library behaves the same however -K is set. */
  useSecond = (LLRMINLIBS >= 2 && nvalid >= 2);

  /* The signal has no reading frame: replicate into the strand's 3 frame planes. */
  for (x = frameStart; x < frameEnd; x++)
    for (i = 0; i < len; i++) {
      external->sr[x][i] = useSecond ? external->covComb2[i] : external->covComb[i];
      if (UTR && x != frameStart)
        external->readcount[x][i] = external->readcount[frameStart][i];
    }

  external->covLambdaBg = 1.0;   /* >0 = sr[] already holds the per-base LLR */
  return 1;
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
	    /* The text HSP path has no range query, so it cannot supply the global
	       background the LLR needs; leave the legacy per-base term in place. */
	    external->covLambdaBg = -1.0;
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
      float scoreHSP = (RREADS / COVNORM) * iv[k].v;
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
  SetCoverageBackground(external, Strand, LengthSequence);

  if (EXPRLLR && FillCoverageLLR(external, Strand, l1, l2, LengthSequence))
    {
      free(buf.a);
    }
  else
    {
      external->covLambdaBg = -1.0;
      if (bw != NULL && external->curLocus != NULL)
        bwQuery(bw, external->curLocus, gS, gE, covCollect, &buf);
      FillCoverageFrameless(external, Strand, l1, l2, buf.a, buf.n);
      free(buf.a);
    }

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
  SetCoverageBackground(external, Strand, LengthSequence);

  /* -L: score each library against its OWN background and combine by MAX (this
     is the only path that supports several libraries; they are never merged).
     Otherwise fall back to the legacy single-source depth fill. */
  if (EXPRLLR && FillCoverageLLR(external, Strand, l1, l2, LengthSequence))
    {
      free(buf.a);
    }
  else
    {
      external->covLambdaBg = -1.0;            /* legacy per-base term in HSPScan2 */
      if (external->bam != NULL && external->curLocus != NULL)
        bamCoverageQuery(external->bam, external->curLocus, gS, gE,
                         wantStrand, BAMSTRAND, covCollect, &buf);
      FillCoverageFrameless(external, Strand, l1, l2, buf.a, buf.n);
      free(buf.a);
    }

  printMess("Preprocessing BAM coverage: step 2");
  HSPScan2(external, NULL, Strand, l1, l2);
#else
  (void) l1; (void) l2; (void) Strand; (void) external; (void) LengthSequence;
  printError("BAM input requires building geneid with WITH_HTSLIB=1");
#endif
}







