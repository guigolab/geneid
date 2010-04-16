/*************************************************************************
*                                                                        *
*   Module: ScoreExons                                                   *
*                                                                        *
*   Score(exon) = reliability measure about coding potential regions     *
*                                                                        *
*   This file is part of the geneid 1.3 distribution                     *
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

extern float MRM;
extern int UTR;
extern int scanORF;
extern float EW;
extern float U12EW;
extern int SRP;
extern float NO_SCORE;
extern int U12GTAG;
extern int U12ATAC;
extern float RSSMARKOVSCORE;
extern int RSS;

/* Matrix to translate characters to numbers. borrowed from jwf */
int TRANS[] = {
  /* Control characters */    
  4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,
  /* Punctuation and digits */
  4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,
  /* Capitals */     /*  A=0; C=1; G=2; T=3; other = 4  */
  4,0,4,1,4,4,4,2,   /* @,A-G: A,C and G found */
  4,4,4,4,4,4,4,4,   /* H-O   */
  4,4,4,4,3,3,4,4,   /* P-W: T and U found */
  4,4,4,4,4,4,4,4,   /* X-Z,etc */
  /* Lower case */
  4,0,4,1,4,4,4,2,   /*  @,A-G  */
  4,4,4,4,4,4,4,4,   /*   H-O   */
  4,4,4,4,3,3,4,4,   /*   P-W   */
  4,4,4,4,4,4,4,4    /* X-Z,etc */
};

/* Translation from string into integer: for oligonucleotides with this length */
/* s: sequence; ls: length of sequence; cardinal: length of the alphabet */
long OligoToInt(char* s, int ls, int cardinal)
{ 
  long index;
  int weight;
  short i;
  
  index = 0;
  weight = 1;
  
  for ( i=ls-1; i>=0; --i )
    {
      index += weight*TRANS[(int)s[i]];
      weight *= cardinal;
    }
  return(index);
}


/* Select the isochore trained to work with this G+C content */
int SelectIsochore(float percent, gparam** isochores)
{
  int i;
  int stop;

  /* Translation from 0.xy to XY */
  percent = PERCENT * percent;

  /* Isochore key */
  i = 0;

  /* Value is between the current isochore range? */
  stop = (percent >= isochores[i]->leftValue &&
	      percent <= isochores[i]->rightValue);

  while(!stop)
    {
      i++;

      /* Value is between the current isochore range? */
      stop = (percent >= isochores[i]->leftValue &&
			  percent <= isochores[i]->rightValue);
    }

  return(i);
}

/* Compute the percentage of G+C nucleotides on a DNA sequence */
float ComputeGC(packGC* GCInfo, long inigc, long endgc)
{
  float percentGC;

  /* %GC = number of C|G divided by the number of "useful" nucleotides */
  /* The idea is to skip the N's in the computing */
  /* Accumulated sum technique: rest between both positions */
  percentGC = ((float)(GCInfo->GC[endgc] - GCInfo->GC[inigc]))
	/ ((float)(endgc-inigc+1 - (GCInfo->N[endgc] - GCInfo->N[inigc])));

  return (percentGC);
}

/* Counting the frequency of C/Gs or Ns found until reaching very position */
void GCScan(char* s, packGC* GCInfo, long l1, long l2)
{
  long i;

  /* Initializing array values: setting first nucleotide */
  switch(*(s+l1))
	{
	case 'C': 
	  GCInfo->GC[0] = 1;
	  GCInfo->N[0] = 0;
	  break;
	case 'G': 
	  GCInfo->GC[0] = 1;
	  GCInfo->N[0] = 0;
	  break;
	case 'A': 
	  GCInfo->GC[0] = 0;
	  GCInfo->N[0] = 0;
	  break;
	case 'T': 
	  GCInfo->GC[0] = 0;
	  GCInfo->N[0] = 0;
	  break;
	default: 
	  GCInfo->GC[0] = 0;
	  GCInfo->N[0] = 1;
	  break;
	}          

  /* 2. Pre-processing the fragment to get the accumulated sum of values */
  for (i = l1+1; i <= l2; i++)
    {
      switch(*(s+i))
		{
		case 'C': 
		  /* Increasing GC counter, preserve N counter */
		  GCInfo->GC[i-l1] = GCInfo->GC[i-l1-1] + 1;
		  GCInfo->N[i-l1] = GCInfo->N[i-l1-1];
		  break;
		case 'G': 
		  /* Increase GC counter, preserve N counter */
		  GCInfo->GC[i-l1] = GCInfo->GC[i-l1-1] + 1;
		  GCInfo->N[i-l1] = GCInfo->N[i-l1-1];
		  break;
		case 'A': 
		  /* Preserve GC counter, preserve N counter */
		  GCInfo->GC[i-l1] = GCInfo->GC[i-l1-1];
		  GCInfo->N[i-l1] = GCInfo->N[i-l1-1];
		  break;
		case 'T': 
		  /* Preserve GC counter, preserve N counter */
		  GCInfo->GC[i-l1] = GCInfo->GC[i-l1-1];
		  GCInfo->N[i-l1] = GCInfo->N[i-l1-1];
		  break;
		default: 
		  /* Preserve GC counter, increase N counter */
		  GCInfo->GC[i-l1] = GCInfo->GC[i-l1-1];
		  GCInfo->N[i-l1] = GCInfo->N[i-l1-1] + 1;
		  break;
		}          
	}
}

/* Compute the coding potential statistic for a sequence: pre-processing */
/* There are Initial (penta) and Transition (hexa) score matrices: */
/* Transition scores for every position are accumulated sums */
void MarkovScan(char* sequence,
                gparam* gp,
                float* OligoDistIni[3], 
                float* OligoDistTran[3],
                long l1, long l2) 
{
  int OligoLength_1;
  long i;
  int intword;
  short x,cp;
  float previousScore;

  /* Pentanucleotides score: initial values for Markov chains */
  for (i=l1; (i<=l2 && *(sequence+i+gp->OligoLength - 1)) ; i++)
    {
      /* Indexing the initialMarkov with the oligonucleotide */
      intword=OligoToInt(sequence+i,gp->OligoLength,4);

      /* Assign the pentanucleotide score depending on the codon position */
      if (intword>=gp->OligoDim)
		for (x=0;x<3;x++)
		  OligoDistIni[x][i-l1] = NULL_OLIGO_SCORE;
      else 
		for (x=0;x<3;x++) 
		  OligoDistIni[x][i-l1] = gp->OligoLogsIni[x][intword];
    }    
  
  /* Hexanucleotides score: transition values for Markov chains */
  /* Accumulated sum is stored for every position and codon position */
  OligoLength_1=gp->OligoLength+1;
  
  /* For every codon position computing the accumulated sum of scores */
  for (x=0; x<FRAMES; x++)
    {
      previousScore = 0.0;
      /* Codon position and frame are different properties of exons */
      cp = (3-x) % 3;
      
      /* Screening the whole sequence to accumulate the sum in every base */
      for (i=l1; (i<=l2 && *(sequence+i+OligoLength_1 - 1)) ; i++)
		{
		  /* Indexing the array with the oligonucleotide identifier */
		  intword=OligoToInt(sequence+i,OligoLength_1,4);

		  if (intword>=gp->OligoDim_1)
			OligoDistTran[x][i-l1]= previousScore + NULL_OLIGO_SCORE;
		  else
			{
			  /* Accumulating step */
			  OligoDistTran[x][i-l1]=
				previousScore + gp->OligoLogsTran[cp][intword];
			}
		  previousScore = OligoDistTran[x][i-l1];

		  /* Shifting one nucleotide: increase codon position */
		  cp = (cp + 1) % 3;
		}
    }
}

/* /\* Projection of HSPs: save the maximum for each nucleotide *\/ */
/* /\* Requirement: HSPs must sorted by Position1 *\/ */
/* void HSPScan(packExternalInformation* external, */
/* 			 packHSP* hsp,  */
/* 			 int Strand,  */
/* 			 long l1, long l2) */
/* { */
/*   short x; */
/*   short frameStart, frameEnd; */
/*   long i,j; */
/*   float scoreHSP; */
  
/*   if (Strand == FORWARD) */
/*     { */
/* 	  frameStart = 0;  */
/* 	  frameEnd = FRAMES; */

/* 	  /\* For each frame and strand, preprocess homology information *\/ */
/* 	  for(x=frameStart; x < frameEnd; x++) */
/* 		{ */
/* 		  /\* Reset arrays: NO_SCORE values *\/ */
/* 		  for(i=0; i < l2-l1+1; i++) */
/* 			external->sr[x][i] = NO_SCORE; */
		  
		  
/* 		  if (hsp != NULL) */
/* 			{ */
/* 			  /\* A. Skip HSPs out of this range: [l1,l2] *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1;  */
/* 				   i++) */
/* 				; */
			  
/* 			  /\* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = hsp->sPairs[x][i]->Score /  */
/* 					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  /\* Save the projection of HSPs into the array: including negative HSPs *\/ */
/* 				  for(j = l1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (scoreHSP > external->sr[x][j-l1] ||  */
/* 						  external->sr[x][j-l1] == NO_SCORE) */
/* 						external->sr[x][j-l1] = scoreHSP; */
/* 					} */
/* 				} */
			  
/* 			  /\* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = hsp->sPairs[x][i]->Score /  */
/* 					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  /\* Save the projection of HSPs into the array: including negative HSPs *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (scoreHSP > external->sr[x][j-l1] ||  */
/* 						  external->sr[x][j-l1] == NO_SCORE) */
/* 						external->sr[x][j-l1] = scoreHSP; */
/* 					} */
/* 				} */
			  
/* 			  /\* Update partial counter: previous HSPs are useless for next split *\/ */
/* 			  external->iSegments[x] = i; 	   */
			  
/* 			  /\* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = hsp->sPairs[x][i]->Score /  */
/* 					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
/* 				  /\* Update the array sr with some positions of current HSPs *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <= l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (scoreHSP > external->sr[x][j-l1] || */
/* 						  external->sr[x][j-l1] == NO_SCORE) */
/* 						external->sr[x][j-l1] = scoreHSP; */
/* 					} */
/* 				} */
/* 			} */
/* 		} */
/* 	} */
/*   else */
/* 	{ */
/* 	  /\* HSPs in REVERSE strand *\/ */
/* 	  frameStart = FRAMES;  */
/* 	  frameEnd = 2*FRAMES; */

/* 	  /\* For each frame and strand, preprocess homology information *\/ */
/* 	  /\* HSPs in reverse strand are reverse-sorted by Position2 *\/ */
/* 	  for(x=frameStart; x < frameEnd; x++) */
/* 		{ */
/* 		  /\* Reset arrays: NO_SCORE values *\/ */
/* 		  for(i=0; i < l2-l1+1; i++) */
/* 			external->sr[x][i] = NO_SCORE; */
		  
/* 		  if (hsp != NULL) */
/* 			{ */
/* 			  /\* A. Skip HSPs out of this range: [l1,l2] *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2;  */
/* 				   i++) */
/* 				; */
			  
/* 			  /\* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = hsp->sPairs[x][i]->Score /  */
/* 					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j < l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (scoreHSP > external->sr[x][j-l1] || */
/* 						  external->sr[x][j-l1] == NO_SCORE) */
/* 						external->sr[x][j-l1] = scoreHSP; */
/* 					} */
/* 				} */
			  
/* 			  /\* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = hsp->sPairs[x][i]->Score /  */
/* 					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (scoreHSP > external->sr[x][j-l1] || */
/* 						  external->sr[x][j-l1] == NO_SCORE) */
/* 						external->sr[x][j-l1] = scoreHSP; */
/* 					} */
/* 				} */
			  
/* 			  /\* Update partial counter: previous HSPs are useless for next split *\/ */
/* 			  external->iSegments[x] = i; 	   */
			  
/* 			  /\* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = hsp->sPairs[x][i]->Score /  */
/* 					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
/* 				  /\* Update the array sr with some positions of current HSPs *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos2;  */
/* 					  j >= hsp->sPairs[x][i]->Pos1 && j >= l1; */
/* 					  j--) */
/* 					{ */
/* 					  if (scoreHSP > external->sr[x][j-l1] || */
/* 						  external->sr[x][j-l1] == NO_SCORE) */
/* 						external->sr[x][j-l1] = scoreHSP; */
/* 					} */
/* 				} */
/* 			} */
/* 		} */
/* 	} */
/* } */
/* /\* Projection of HSPs: save the maximum for each nucleotide *\/ */
/* /\* Requirement: HSPs must sorted by Position1 *\/ */
/* void ReadScan(packExternalInformation* external, */
/* 			 packHSP* hsp,  */
/* 			 int Strand,  */
/* 			 long l1, long l2) */
/* { */
/*   short x; */
/*   short frameStart, frameEnd; */
/*   long i,j; */
/*   float scoreHSP; */
  
/*   if (Strand == FORWARD) */
/*     { */
/* 	  frameStart = 0;  */
/* 	  frameEnd = FRAMES; */

/* 	  /\* For each frame and strand, preprocess homology information *\/ */
/* 	  for(x=frameStart; x < frameEnd; x++) */
/* 		{ */
/* 		  /\* Reset arrays: NO_SCORE values *\/ */
/* 		  for(i=0; i < l2-l1+1; i++) */
/* 			external->sr[x][i] = NO_SCORE; */
/* 		  /\* Reset arrays: NO_SCORE values *\/ */
/* 		  for(i=0; i < l2-l1+1; i++) */
/* 			external->readcount[x][i] = 0.0; */
		  
		  
/* 		  if (hsp != NULL) */
/* 			{ */
/* 			  /\* A. Skip HSPs out of this range: [l1,l2] *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1;  */
/* 				   i++) */
/* 				; */
			  
/* 			  /\* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; */
/* 				  /\* / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); *\/ */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  /\* Save the projection of HSPs into the array: including negative HSPs *\/ */
/* 				  for(j = l1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (external->sr[x][j-l1] == NO_SCORE){ */
/* 					    external->sr[x][j-l1] = scoreHSP; */
/* 					  }else{ */
/* 					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV); */
/* 					  } */
/* 					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score; */
/* 					} */
/* 				} */
			  
/* 			  /\* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;  */
/* /\* 				    /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); *\/ */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  /\* Save the projection of HSPs into the array: including negative HSPs *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (external->sr[x][j-l1] == NO_SCORE){ */
/* 					    external->sr[x][j-l1] = scoreHSP; */
/* 					  }else{ */
/* 					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV); */
/* 					  } */
/* 					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score; */
/* 					} */
/* 				} */
			  
/* 			  /\* Update partial counter: previous HSPs are useless for next split *\/ */
/* 			  external->iSegments[x] = i; 	   */
			  
/* 			  /\* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;  */
/* /\* 				    /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); *\/ */
				  
/* 				  /\* Update the array sr with some positions of current HSPs *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <= l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (external->sr[x][j-l1] == NO_SCORE){ */
/* 					    external->sr[x][j-l1] = scoreHSP; */
/* 					  }else{ */
/* 					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV); */
/* 					  } */
/* 					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score; */
/* 					} */
/* 				} */
/* 			} */
/* 		} */
/* 	} */
/*   else */
/* 	{ */
/* 	  /\* HSPs in REVERSE strand *\/ */
/* 	  frameStart = FRAMES;  */
/* 	  frameEnd = 2*FRAMES; */

/* 	  /\* For each frame and strand, preprocess homology information *\/ */
/* 	  /\* HSPs in reverse strand are reverse-sorted by Position2 *\/ */
/* 	  for(x=frameStart; x < frameEnd; x++) */
/* 		{ */
/* 		  /\* Reset arrays: NO_SCORE values *\/ */
/* 		  for(i=0; i < l2-l1+1; i++) */
/* 			external->sr[x][i] = NO_SCORE; */
/* 		  /\* Reset arrays: NO_SCORE values *\/ */
/* 		  for(i=0; i < l2-l1+1; i++) */
/* 			external->readcount[x][i] = 0.0; */
		  
/* 		  if (hsp != NULL) */
/* 			{ */
/* 			  /\* A. Skip HSPs out of this range: [l1,l2] *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2;  */
/* 				   i++) */
/* 				; */
			  
/* 			  /\* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) *\/ */
/* 			  for (i = external->iSegments[x];  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; */
/* /\* 				  / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); *\/ */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j < l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (external->sr[x][j-l1] == NO_SCORE){ */
/* 					    external->sr[x][j-l1] = scoreHSP; */
/* 					  }else{ */
/* 					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV); */
/* 					  } */
/* 					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score; */
/* 					} */
/* 				} */
			  
/* 			  /\* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; */
/* /\* 				  /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); *\/ */
				  
/* 				  /\* For each position in the HSP update the array sr *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos1;  */
/* 					  j <= hsp->sPairs[x][i]->Pos2 && j <l2; */
/* 					  j++) */
/* 					{ */
/* 					  if (external->sr[x][j-l1] == NO_SCORE){ */
/* 					    external->sr[x][j-l1] = scoreHSP; */
/* 					  }else{ */
/* 					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV); */
/* 					  } */
/* 					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score; */
/* 					} */
/* 				} */
			  
/* 			  /\* Update partial counter: previous HSPs are useless for next split *\/ */
/* 			  external->iSegments[x] = i; 	   */
			  
/* 			  /\* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) *\/ */
/* 			  for (;  */
/* 				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1;  */
/* 				   i++) */
/* 				{ */
/* 				  /\* Score value *\/ */
/* 				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; */
/* /\* 				  / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); *\/ */
				  
/* 				  /\* Update the array sr with some positions of current HSPs *\/ */
/* 				  for(j = hsp->sPairs[x][i]->Pos2;  */
/* 					  j >= hsp->sPairs[x][i]->Pos1 && j >= l1; */
/* 					  j--) */
/* 					{ */
/* 					  if (external->sr[x][j-l1] == NO_SCORE){ */
/* 					    external->sr[x][j-l1] = scoreHSP; */
/* 					  }else{ */
/* 					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV); */
/* 					  } */
/* 					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score; */
/* 					} */
/* 				} */
/* 			} */
/* 		} */
/* 	} */
/* } */

/* /\* Preprocessing of HSPs projections *\/ */
/* void HSPScan2(packExternalInformation* external, */
/* 			  packHSP* hsp,  */
/* 			  int Strand,  */
/* 			  long l1, long l2) */
/* { */
/*   short x; */
/*   long i; */
/*   float previousScore; */
/*   float previousReadCount; */
/*   short frameStart, frameEnd; */
  
/*   if (Strand == FORWARD) */
/*     { */
/* 	  frameStart = 0;  */
/* 	  frameEnd = FRAMES; */
/* 	} */
/*   else */
/* 	{ */
/* 	  frameStart = FRAMES;  */
/* 	  frameEnd = 2*FRAMES; */
/* 	} */
  
/*   for(x=frameStart; x < frameEnd; x++) */
/*     { */
/*       previousScore = 0.0; */
/*       previousReadCount = 0.0; */
/*       /\* Screening the whole sequence to accumulate the sum in every base *\/ */
/*       for (i=l1; i<=l2; i++) */
/* 		{ */
/* 		  /\* Accumulating step *\/ */
/* 		  if (UTR){ */
/* 		    external->sr[x][i-l1] = previousScore + log(external->sr[x][i-l1] + 1); */
/* 		    previousScore = external->sr[x][i-l1]; */
/* 		    external->readcount[x][i-l1] = previousReadCount + external->readcount[x][i-l1]; */
/* 		    previousReadCount = external->readcount[x][i-l1]; */
/* 		  }else{ */
/* 		    external->sr[x][i-l1] = previousScore + external->sr[x][i-l1]; */
/* 		    previousScore = external->sr[x][i-l1]; */
/* 		  } */
/* 		} */
/*     }   */
/* } */

/* Homology to protein score: using homology information (blast HSPs) */
float ScoreHSPexon(exonGFF* exon, 
				   int Strand, 
				   packExternalInformation* external, 
				   long l1, long l2)
{
  int index;
  short trueFrame;
  long iniExon, endExon;
  float Score;

  iniExon = exon->Acceptor->Position - l1 + COFFSET;
  endExon = exon->Donor->Position - l1 + COFFSET ;

  if (Strand == FORWARD)
    index = 0; 
  else
    index = FRAMES;
  
  /* Frame blast definition: according to the sequence start */
  trueFrame = (exon->Acceptor->Position + COFFSET + exon->Frame) % 3;
  trueFrame += index;
  
  /* Access the sr array to obtain the homology score for current score */
  Score = external->sr[trueFrame][endExon] - 
	external->sr[trueFrame][(iniExon>0)?iniExon-1:iniExon];

  return(Score);
}
/* Homology to protein score: using homology information (blast HSPs) */
float GetReadCount(exonGFF* exon, 
				   int Strand, 
				   packExternalInformation* external, 
				   long l1, long l2)
{
  int index;
  short trueFrame;
  long iniExon, endExon;
  float Score;
/*   float kb = 1000.000; */
  iniExon = exon->Acceptor->Position - l1 + COFFSET;
  endExon = exon->Donor->Position - l1 + COFFSET ;

  if (Strand == FORWARD)
    index = 0; 
  else
    index = FRAMES;
  
  /* Frame blast definition: according to the sequence start */
  trueFrame = (exon->Acceptor->Position + COFFSET + exon->Frame) % 3;
  trueFrame += index;
  
  /* Access the sr array to obtain the homology score for current score */
  if ((endExon-iniExon + 1)>0){
      Score = (external->readcount[trueFrame][endExon] - 
	    external->readcount[trueFrame][(iniExon>0)?iniExon-1:iniExon]);
    }else{
      Score = 0.0;
    }
  return(Score);
}
/* Computing the score of a list of exons from (Sites, Markov, homology) */
long Score(exonGFF *Exons,
           long nExons,
           long l1,
           long l2,
           int Strand,
	   packExternalInformation* external,
           packHSP* hsp,
           gparam** isochores,
           packGC* GCInfo)
{
  long iniExon, endExon, i, j;
  int exonLen;
  long inigc, endgc;
  short frame;
  short codonPosition;
  long n;
  float scoreMarkov;
  float scoreHSP;
  float exonR;
  float scoreTotal;
  int OligoLength_1;
  float ExonWeight;
  float percentGC;
  int currentIsochore;
  paramexons* p;
  int OligoLength;
  gparam* gp;
  int rss = 0;
/*   char mess[MAXSTRING]; */
/*   float million = 1000000; */
/*   int u12correction; */
/*   u12correction = 0; */
  /* Number of survivor exons after scoring and filtering */
  n=0;
  /* For every exon computing scores: protein coding and homology info */
  for (i=0;i<nExons;i++) 
    {
      /* Measure G+C around the exon to select the best isochore */
      /* NOTE: Stop is not contained in Terminals, Singles and ORFs */
      iniExon=(Exons+i)->Acceptor->Position - l1;
      endExon=(Exons+i)->Donor->Position - l1;
      exonLen=endExon-iniExon+1;
      rss = 0;
      /* 0. Get G+C content of the region around the exon (local) */
      /* selecting the proper isochore to score the exon */
      if (iniExon <= ISOCONTEXT)
		inigc=0;
      else
		inigc=iniExon-ISOCONTEXT;

      if (endExon+ISOCONTEXT >= l2-l1)
		endgc=l2-l1;
      else
		endgc=endExon+ISOCONTEXT;

      percentGC = ComputeGC(GCInfo,inigc,endgc); 
     
      currentIsochore = SelectIsochore(percentGC, isochores);
      gp = isochores[currentIsochore];

      OligoLength = gp->OligoLength;
      OligoLength_1 = OligoLength+1;

	  /* Default selection of params: UTR exons */
	  p = gp->utr;

      /* Selecting parameters according to exon type */
      if (!(strcmp((Exons+i)->Type,sFIRST)))
		p = gp->Initial;
   
      if (!(strcmp((Exons+i)->Type,sINTERNAL)))
		p = gp->Internal;

      if (!(strcmp((Exons+i)->Type,sZEROLENGTH)))
		p = gp->Internal;

      if (!(strcmp((Exons+i)->Type,sTERMINAL)))
		p = gp->Terminal;

      if (!(strcmp((Exons+i)->Type,sSINGLE)))
		p = gp->Single;
      
      if (!(strcmp((Exons+i)->Type,sORF)))
		p = gp->Single;

/*       if (!(strcmp((Exons+i)->Type,sUTRFIRST))) */
/* 		p = gp->utr; */
      
      /* 1. Coding potential score: initial plus accumulated sums */
      /* Checkpoint for exons shorter than a minimum value */
      if ((exonLen < MINEXONLENGTH)
	  || (strcmp((Exons+i)->Type,sFIRST)&&strcmp((Exons+i)->Type,sINTERNAL)&&strcmp((Exons+i)->Type,sTERMINAL)&&strcmp((Exons+i)->Type,sSINGLE)&&strcmp((Exons+i)->Type,sORF))

	  ){
	scoreMarkov = MINSCORELENGTH;
	if (exonLen == 0){scoreMarkov = RSSMARKOVSCORE;}
      }
      else
	{
	  scoreMarkov = 0.0;
	  frame = (Exons+i)->Frame;
         
	  /* Translate frame to position into codon */
	  codonPosition = (3 - frame) % 3;
   
	  /* Assign initial probability: pentanucleotide */
	  scoreMarkov += gp->OligoDistIni[codonPosition][iniExon];
         
	  /* Which one of the three combinations? */
	  j = (iniExon + (3-codonPosition)) % 3;
   
	  /* Accumulating transition probabilities: hexanucleotides */    
	  scoreMarkov +=
	    gp->OligoDistTran[j][(endExon>OligoLength_1)?endExon-OligoLength_1+1 : endExon]
	    - gp->OligoDistTran[j][(iniExon)? iniExon - 1 : 0];
	}

      /* First cutoff: coding potential score */
      if (scoreMarkov >= p->OligoCutoff) 
	{
		  
	  /* 2. Homology to protein score */
	  scoreHSP = 0.0;
	  exonR = 0.0;
	  if (SRP)
            scoreHSP = ScoreHSPexon((Exons+i),Strand,external,l1,l2);
	  if (UTR){
	    exonR = GetReadCount((Exons+i),Strand,external,l1,l2);
/*             scoreHSP = scoreHSP - log(NTMAPPED/1000000); */
/* 	    sprintf(mess,"rpkm: %f len:%d type:%s",exonRPKM,exonLen,(Exons+i)->Type); */
/* 	    printMess(mess); */
	    /* scoreHSP = scoreHSP *(RREADS/MRM); */

	    
	  }
	  /* 3. Total (combined) score */
	  /* Don't want to count shared signals twice - in the case of half exons */
	  if (((Strand == FORWARD)&&(!strcmp((Exons+i)->Type,sUTRFIRSTHALF)||!strcmp((Exons+i)->Type,sUTR5INTERNALHALF)))
	      ||
	      ((Strand == REVERSE)&&(!strcmp((Exons+i)->Type,sUTRTERMINALHALF)||!strcmp((Exons+i)->Type,sUTR3INTERNALHALF))))
	    {
	      scoreTotal = 
		(p->siteFactor * ((Exons+i)->Acceptor->Score))
		+ (p->exonFactor * scoreMarkov) 
		+ (p->HSPFactor * scoreHSP);
	    }else{
	    if(((Strand == REVERSE)&&(!strcmp((Exons+i)->Type,sUTRFIRSTHALF)||!strcmp((Exons+i)->Type,sUTR5INTERNALHALF)))
	       ||
	       ((Strand == FORWARD)&&(!strcmp((Exons+i)->Type,sUTRTERMINALHALF)||!strcmp((Exons+i)->Type,sUTR3INTERNALHALF))))
	      {
		scoreTotal = 
		  (p->siteFactor * ((Exons+i)->Donor->Score))
		  + (p->exonFactor * scoreMarkov) 
		  + (p->HSPFactor * scoreHSP);
	      }else{
	      scoreTotal = 
		(p->siteFactor * ((Exons+i)->Acceptor->Score + (Exons+i)->Donor->Score))
		+ (p->exonFactor * scoreMarkov) 
		+ (p->HSPFactor * scoreHSP); 
	    }
	  }
	  /* Second cutoff- final score */
	  if (scoreTotal >= p->ExonCutoff) 
	    {
	      Exons[n]=Exons[i];

	      /* -E: increase/decrease current ExonWeight parameter */
	      if (EW == NOVALUE)
		ExonWeight = p->ExonWeight;
	      else
		ExonWeight = p->ExonWeight + EW;

/* 	      ExonWeight = p->ExonWeight; */
/* 	      if (EW != NOVALUE){ */
/* 		ExonWeight = ExonWeight + EW; */
/* 	      } */
	      if ((U12EW != NOVALUE)&&(((Exons+i)->Acceptor->class != U2) && ((Exons+i)->Donor->class != U2))){	    
		ExonWeight = ExonWeight + U12EW;			       
	      }
	      (Exons+n)->Score = scoreTotal + ExonWeight;
	      (Exons+n)->PartialScore = scoreMarkov;
	      (Exons+n)->HSPScore = scoreHSP;
	      (Exons+n)->R = exonR;
/* 	      sprintf(mess,"rpkm: %f",(Exons+n)->R); */
/* 	      printMess(mess); */
	      n++;
	    }
	}
    }
  return(n);
}

/* Management function to score and filter exons */
void ScoreExons(char *Sequence, 
                packExons* allExons, 
                long l1,
                long l2,
                int Strand,
		packExternalInformation* external,
                packHSP* hsp,
                gparam** isochores,
                int nIsochores,
                packGC* GCInfo)
{
  long i;
  char mess[MAXSTRING];

  /* Fill in the temporary Markov arrays (pre-processing) */
  for(i=0; i<nIsochores; i++)
    {
      /* Filling OligoDistIni and OligoDistTran */
      printMess("Preprocessing coding potential scores");
      MarkovScan(Sequence, isochores[i], 
				 isochores[i]->OligoDistIni, 
				 isochores[i]->OligoDistTran, 
				 l1, l2);
    }

  /* Fill in the temporary HSP arrays (pre-processing) */
  /* GENIS hack */
/*   if (SRP) */
/* 	{ */
/* 	  if (UTR){ */
/* 	    printMess("Preprocessing read information: step 1"); */
/* 	    ReadScan(external,hsp,Strand,l1,l2); */
/* 	  }else{ */
/* 	    printMess("Preprocessing homology information: step 1"); */
/* 	    HSPScan(external,hsp,Strand,l1,l2); */
/* 	  } */

/* 	  printMess("Preprocessing homology information: step 2"); */
/* 	  HSPScan2(external,hsp,Strand,l1,l2); */
/* 	} */
  
  /* 2. Ready for Score and Filter Exons? */
  printMess("Scoring and filtering exons");
  
  allExons->nInitialExons = Score(allExons->InitialExons,
								  allExons->nInitialExons,
								  l1, l2, Strand, 
								  external, hsp, 
								  isochores, GCInfo);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 

  allExons->nInternalExons=Score(allExons->InternalExons,
								 allExons->nInternalExons, 
								 l1, l2, Strand, 
								 external, hsp, 
								 isochores, GCInfo);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 
  if (RSS){
    allExons->nZeroLengthExons=Score(allExons->ZeroLengthExons,
				     allExons->nZeroLengthExons, 
				     l1, l2, Strand, 
				     external, hsp, 
				     isochores, GCInfo);
    sprintf(mess,"Zero-length Exons \t\t%8ld", allExons->nZeroLengthExons);
    printRes(mess); 
  }
  allExons->nTerminalExons=Score(allExons->TerminalExons,
								 allExons->nTerminalExons,
								 l1, l2, Strand, 
								 external, hsp, 
								 isochores, GCInfo);
  sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
  printRes(mess);  
    
  allExons->nSingles=Score(allExons->Singles,
						   allExons->nSingles,
						   l1, l2, Strand, 
						   external, hsp, 
						   isochores, GCInfo);
  sprintf(mess,"Singles \t\t%8ld", allExons->nSingles);
  printRes(mess);

  if (scanORF)
    {
      allExons->nORFs=Score(allExons->ORFs,
							allExons->nORFs,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"ORFs \t\t\t%8ld", allExons->nORFs);
      printRes(mess);
    }
  if (UTR)
    {
      allExons->nUtrInitialExons=Score(allExons->UtrInitialExons,
							allExons->nUtrInitialExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR Initial Exons \t%8ld", allExons->nUtrInitialExons);
      printRes(mess);

      allExons->nUtrInitialHalfExons=Score(allExons->UtrInitialHalfExons,
							allExons->nUtrInitialHalfExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR Initial Half Exons \t%8ld", allExons->nUtrInitialHalfExons);
      printRes(mess);

      allExons->nUtrInternalExons=Score(allExons->UtrInternalExons,
							allExons->nUtrInternalExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR Internal Exons \t%8ld", allExons->nUtrInternalExons);
      printRes(mess);

      allExons->nUtr5InternalHalfExons=Score(allExons->Utr5InternalHalfExons,
							allExons->nUtr5InternalHalfExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR 5' Int. Half Exons \t%8ld", allExons->nUtr5InternalHalfExons);
      printRes(mess);

      allExons->nUtr3InternalHalfExons=Score(allExons->Utr3InternalHalfExons,
							allExons->nUtr3InternalHalfExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR 3' Int. Half Exons \t%8ld", allExons->nUtr3InternalHalfExons);
      printRes(mess);

      allExons->nUtrTerminalHalfExons=Score(allExons->UtrTerminalHalfExons,
							allExons->nUtrTerminalHalfExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR Term. Half Exons \t%8ld", allExons->nUtrTerminalHalfExons);
      printRes(mess);

      allExons->nUtrTerminalExons=Score(allExons->UtrTerminalExons,
							allExons->nUtrTerminalExons,
							l1, l2, Strand, 
							external, hsp, 
							isochores, GCInfo);
      sprintf(mess,"UTR Terminal Exons \t%8ld", allExons->nUtrTerminalExons);
      printRes(mess);
    }
      
    
}







