/*************************************************************************
*                                                                        *
*   Module: ScoreExons                                                   *
*                                                                        *
*   Score(exon) = reliability measure about coding potential regions     *
*                                                                        *
*   This file is part of the geneid 1.1 distribution                     *
*                                                                        *
*     Copyright (C) 2001 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
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

extern int scanORF;
extern float EW;
extern int SRP;

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

/* Homology to protein score: using homology information (SR or HSPs) */
double ScoreSRexon(exonGFF* exon, int Strand, packSR* sr)
{
  int index;
  short trueFrame;
  long iniExon, endExon;
  int i;
  double Score;
  long lIntersection;
  long lSR;
  long a,b;

  iniExon = exon->Acceptor->Position + COFFSET;
  endExon = exon->Donor->Position + COFFSET;

  if (Strand == FORWARD)
    index = 0; 
  else
    index = FRAMES;
  
  /* Frame blast definition: according to the sequence start */
  trueFrame = (iniExon + exon->Frame) % 3;
  trueFrame += index;
  
  /* Scan all of the SR's in this frame/strand (HSPs are accepted) */
  i = 0;
  Score = 0;
  while(i < sr->nRegions[trueFrame])
    {
      /* Looking for match between the SR and the exon */
      a = MIN(endExon,sr->sRegions[trueFrame][i].Pos2);
      b = MAX(iniExon,sr->sRegions[trueFrame][i].Pos1);
      lIntersection = a - b + 1;
      
      if (lIntersection > 0)
		{
		  lSR = sr->sRegions[trueFrame][i].Pos2 - sr->sRegions[trueFrame][i].Pos1 + 1;
		  Score += sr->sRegions[trueFrame][i].Score * ((float)lIntersection / (float)lSR);
		}
      i++;
    }

  return(Score);
}

/* Computing the score of a list of exons (Markov, SR and total) */
long Score(exonGFF *Exons,
           long nExons,
           long l1,
           long l2,
           int Strand,
           packSR* sr,
           gparam** isochores,
           packGC* GCInfo)
{
  long iniexon, endexon, i, j;
  int exonlen;
  long inigc, endgc;
  short frame;
  short codonPosition;
  long n;
  double scoreMarkov;
  double scoreSR;
  double scoreTotal;
  int OligoLength_1;
  double ExonWeight;
  float percentGC;
  int currentIsochore;
  paramexons* p;
  int OligoLength;
  gparam* gp;

  /* Number of survivor exons after scoring and filtering */
  n=0;
  /* For every exon computing scores: protein coding and homology info */
  for (i=0;i<nExons;i++) 
    {
      /* Measure G+C around the exon to select the best isochore */
      /* NOTE: Stop is not contained in Terminals, Singles and ORFs */
      iniexon=(Exons+i)->Acceptor->Position - l1;
      endexon=(Exons+i)->Donor->Position - l1;
      exonlen=endexon-iniexon+1;

      /* 0. Get G+C content of the region around the exon (local) */
      /* selecting the proper isochore to score the exon */
      if (iniexon <= ISOCONTEXT)
		inigc=0;
      else
		inigc=iniexon-ISOCONTEXT;

      if (endexon+ISOCONTEXT >= l2-l1)
		endgc=l2-l1;
      else
		endgc=endexon+ISOCONTEXT;

      percentGC = ComputeGC(GCInfo,inigc,endgc); 
     
      currentIsochore = SelectIsochore(percentGC, isochores);
      gp = isochores[currentIsochore];

      OligoLength = gp->OligoLength;
      OligoLength_1 = OligoLength+1;

	  /* Default selection of params: First exons */
	  p = gp->Initial;

      /* Selecting parameters according to exon type */
      if (!(strcmp((Exons+i)->Type,sFIRST)))
		p = gp->Initial;
      if (!(strcmp((Exons+i)->Type,sINTERNAL)))
		p = gp->Internal;
      if (!(strcmp((Exons+i)->Type,sTERMINAL)))
		p = gp->Terminal;
      if (!(strcmp((Exons+i)->Type,sSINGLE)))
		p = gp->Single;
      if (!(strcmp((Exons+i)->Type,sORF)))
		p = gp->Single;
      
      /* 1. Coding potential score: initial plus accumulated sums */
      /* Checkpoint for exons shorter than a minimum value */
      if (exonlen < MINEXONLENGTH)
		scoreMarkov = MINSCORELENGTH;
      else
		{
		  scoreMarkov = 0.0;
		  frame = (Exons+i)->Frame;
         
		  /* Translate frame to position into codon */
		  codonPosition = (3 - frame) % 3;
   
		  /* Assign initial probability: pentanucleotide */
		  scoreMarkov += gp->OligoDistIni[codonPosition][iniexon];
         
		  /* Which one of the three combinations? */
		  j = (iniexon + (3-codonPosition)) % 3;
   
		  /* Accumulating transition probabilities: hexanucleotides */    
		  scoreMarkov +=
			gp->OligoDistTran[j][(endexon>OligoLength_1)?endexon-OligoLength_1+1 : endexon]
			- gp->OligoDistTran[j][(iniexon)? iniexon - 1 : 0];
		}

      /* First cutoff: coding potential score */
      if (scoreMarkov >= p->OligoCutoff) 
		{
		  /* 2. Homology to protein score */
		  scoreSR = 0;
		  if (SRP)
            scoreSR = ScoreSRexon((Exons+i),Strand,sr);

		  /* 3. Total (combined) score */
		  scoreTotal = 
			((1 - p->OligoWeight) * ((Exons+i)->Acceptor->Score 
									 + (Exons+i)->Donor->Score))
			+ (p->OligoWeight * scoreMarkov) + scoreSR; 
	  
		  /* Second cutoff- final score */
		  if (scoreTotal >= p->ExonCutoff) 
			{
			  Exons[n]=Exons[i];

			  /* -E: increase/decrease current ExonWeight parameter */
			  if (EW == NOVALUE)
				ExonWeight = p->ExonWeight;
			  else
				ExonWeight = p->ExonWeight + EW;

			  (Exons+n)->Score = scoreTotal + ExonWeight;
			  (Exons+n)->PartialScore = scoreMarkov;
			  (Exons+n)->SRScore = scoreSR;
   
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
                packSR* sr,
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
      printMess("Computing accumulated sum for coding potential score");
      MarkovScan(Sequence, isochores[i], 
				 isochores[i]->OligoDistIni, 
				 isochores[i]->OligoDistTran, 
				 l1, l2);
    }
  
  /* 2. Ready for Score and Filter Exons? */
  printMess("Scoring and filtering exons");
  
  allExons->nInitialExons = Score(allExons->InitialExons,
								  allExons->nInitialExons,
								  l1, l2, Strand, sr, 
								  isochores, GCInfo);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 
  
  allExons->nInternalExons=Score(allExons->InternalExons,
								 allExons->nInternalExons, 
								 l1, l2, Strand, sr, 
								 isochores, GCInfo);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 
  
  allExons->nTerminalExons=Score(allExons->TerminalExons,
								 allExons->nTerminalExons,
								 l1, l2, Strand, sr, 
								 isochores, GCInfo);
  sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
  printRes(mess);  
  
  allExons->nSingles=Score(allExons->Singles,
						   allExons->nSingles,
						   l1, l2, Strand, sr, 
						   isochores, GCInfo);
  sprintf(mess,"Singles \t\t%8ld", allExons->nSingles);
  printRes(mess);

  if (scanORF)
    {
      allExons->nORFs=Score(allExons->ORFs,
							allExons->nORFs,
							l1, l2, Strand, sr, 
							isochores, GCInfo);
      sprintf(mess,"ORFs \t\t\t%8ld", allExons->nORFs);
      printRes(mess);
    }
}







