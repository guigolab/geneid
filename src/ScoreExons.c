/*************************************************************************
*                                                                        *
*   Module: ScoreExons                                                   *
*                                                                        *
*   This module scores each exon using Markov models.                    *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
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

extern float EW;
extern int SRP;

/* variable to translate characters to numbers. borrowed from jwf */
int TRANS[] = {
    /* Control characters */     /*  A=0; C=1; G=2; T=3; other = 4  */
    4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,
    /* Punctuation and digits */
    4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,
    /* Capitals */
    4,0,4,1,4,4,4,2,   /*  @,A-G  */
    4,4,4,4,4,4,4,4,   /*   H-O   */
    4,4,4,4,3,3,4,4,   /*   P-W   */
    4,4,4,4,4,4,4,4,   /* X-Z,etc */
    /* Lower case */
    4,0,4,1,4,4,4,2,   /*  @,A-G  */
    4,4,4,4,4,4,4,4,   /*   H-O   */
    4,4,4,4,3,3,4,4,   /*   P-W   */
    4,4,4,4,4,4,4,4    /* X-Z,etc */
  };

int SearchIsochore(float percent, gparam** isochores)
{
  int i;
  int stop;

  i = 0;

  percent = percent * PERCENT;
  stop = (percent >= isochores[i]->leftValue && 
	  percent <= isochores[i]->rightValue);

  while(!stop)
    {
      i++;

      stop = (percent >= isochores[i]->leftValue && 
	      percent <= isochores[i]->rightValue);
    }

  return(i);
}

float MeasureSequence(long l1,long l2,char* s)
{
  long i;
  long Measure[5];
  float percent;

  for (i = 0; i < 5; i++)
    Measure[i] = 0;

  for (i = l1; i < l2; i++)
      Measure[TRANS[(int)(*(s+i))]]++;

  percent = (float)(Measure[TRANS['C']] + Measure[TRANS['G']])/((float)l2-l1+1);

  return(percent);
}

void CGScan(char* s, long* CG, long l1, long l2)
{
  long i;

  /* Initializing ... */
  if((*(s+l1)) == 'C' || (*(s+l1)) == 'G')
    CG[0]=1;
  else
    CG[0]=0;

  /* Scanning DNA split sequence */
  for (i = l1+1; i <= l2; i++)
    {
      if((*(s+i)) == 'C' || (*(s+i)) == 'G')
	CG[i-l1] = CG[i-l1-1] + 1;
      else
	CG[i-l1] = CG[i-l1-1];
    }
}

long OligoToInt(char *s, int ls) 
{ 
  long index;
  int weight;
  short i;
  
  index = 0;
  weight = 1;
  
  for ( i=ls-1; i>=0; --i ) 
    {
      index += weight*TRANS[(int)s[i]];
      weight *= 4;
    }
  return(index);
}

void MarkovScan(char *sequence, gparam* gp,
		float *OligoDistIni[3], 
		float *OligoDistTran[3],
		long l1, long l2) 
{
  int OligoLength_1;
  long i;
  int intword;
  short x,cp;
  float previousScore;

  /* Pentanucleotides score */
  for (i=l1; (i<=l2 && *(sequence+i+gp->OligoLength - 1)) ; i++)
    {
      intword=OligoToInt(sequence+i,gp->OligoLength);
      if (intword>=gp->OligoDim)
	for (x=0;x<3;x++)
	  OligoDistIni[x][i-l1] = NULL_OLIGO_SCORE;
      else 
	for (x=0;x<3;x++) 
	  OligoDistIni[x][i-l1] = gp->OligoLogsIni[x][intword];
    }    
  
  /* Hexanucleotides score */
  OligoLength_1=gp->OligoLength+1;
  
  for (x=0; x<FRAMES; x++)
    {
      previousScore = 0.0;
      cp = (3-x) % 3;
      
      for (i=l1; (i<=l2 && *(sequence+i+OligoLength_1 - 1)) ; i++)
	{
	  intword=OligoToInt(sequence+i,OligoLength_1);
	  if (intword>=gp->OligoDim_1)
	    OligoDistTran[x][i-l1]= previousScore + NULL_OLIGO_SCORE;
	  else
	    {
	      OligoDistTran[x][i-l1]=
		previousScore + gp->OligoLogsTran[cp][intword];
	    }
	  previousScore = OligoDistTran[x][i-l1];
	  cp = (cp + 1) % 3;
	}
    }
}

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
  
  if (!strcmp(exon->Type,"Terminal"))
    endExon += LENGTHCODON;

  if (Strand == FORWARD)
    index = 0; 
  else
    index = FRAMES;
  
  /* Frame about start of sequence */
  trueFrame = (iniExon + exon->Frame) % 3;
  trueFrame += index;

  /* Forgetting sr-regions before this exon */
  i = sr->iRegions[trueFrame];

  while(i < sr->nRegions[trueFrame] && 
	iniExon > sr->sRegions[trueFrame][i].Pos2)
    i ++;
  sr->iRegions[trueFrame] = i;

  /* Getting sr-score for this exon */
  Score = 0;
  while(i < sr->nRegions[trueFrame] && 
	endExon > sr->sRegions[trueFrame][i].Pos1)
    {
      /* Scoring with current sr */
      a = MIN(endExon,sr->sRegions[trueFrame][i].Pos2);
      b = MAX(iniExon,sr->sRegions[trueFrame][i].Pos1);
      lIntersection = a - b + 1;
      lSR = sr->sRegions[trueFrame][i].Pos2 - sr->sRegions[trueFrame][i].Pos1 + 1;
      Score += sr->sRegions[trueFrame][i].Score * ((float)lIntersection / (float)lSR);

      i++;
    }

  return(Score);
}

long Score(exonGFF *Exons, long nExons,
	   long l1, long l2, int Strand, packSR* sr, 
	   long* CG, gparam** isochores)
{
  long iniexon, endexon, i, j;
  int exonlen;
  long inigc, endgc;
  short frame;
  short codonPosition;
  long n=0;
  double scoreMarkov;
  double scoreSR;
  double scoreTotal;
  int OligoLength_1;
  double ExonWeight;
  float valueCG;
  int currentIsochore;
  paramexons* p;
  int OligoLength;

  if (SRP)
    {
      /* Reset similarity regions counters */
      for (i=0; i<STRANDS * FRAMES; i++)
	sr->iRegions[i] = 0;
    }

  /* For each exon computing scores */
  for (i=0;i<nExons;i++) 
    {
      /* Measure C+G for current exon: selecting isochore */
      iniexon=(Exons+i)->Acceptor->Position - l1;
      endexon=(Exons+i)->Donor->Position - l1;
      exonlen=endexon-iniexon+1;

      /* get GC content of the region around the exon */
      if (iniexon <= ISOCONTEXT)
	inigc=0;
      else
	inigc=iniexon-ISOCONTEXT;

      if (endexon+ISOCONTEXT >= l2-l1)
	endgc=l2-l1;
      else
	endgc=endexon+ISOCONTEXT;

      valueCG = ((float)(CG[endgc] - CG[inigc])) / ((float)(endgc - inigc + 1));

      currentIsochore = SearchIsochore(valueCG, isochores);

      OligoLength = isochores[currentIsochore]->OligoLength;
      OligoLength_1 = OligoLength+1;

      /* Selecting statistics according to exon type */
      if (!(strcmp((Exons+i)->Type,sFIRST)))
	p = isochores[currentIsochore]->Initial;
      if (!(strcmp((Exons+i)->Type,sINTERNAL)))
	p = isochores[currentIsochore]->Internal;
      if (!(strcmp((Exons+i)->Type,sTERMINAL)))
	p = isochores[currentIsochore]->Terminal;
      if (!(strcmp((Exons+i)->Type,sSINGLE)))
	p = isochores[currentIsochore]->Single;

      if (EW == NOVALUE)
	ExonWeight = p->ExonWeight;
      else
	ExonWeight = p->ExonWeight + EW;

      /* 1.markov score exon */
      scoreMarkov=0.0;
      frame=(Exons+i)->Frame;
      
      /* Translate frame to position into codon */
      codonPosition = (3 - frame) % 3;

      /* initial probability: Pentanucleotide */
      scoreMarkov+=
	isochores[currentIsochore]->OligoDistIni[codonPosition][iniexon]; 
      
      /* What kind of adding? */
      j = (iniexon + (3-codonPosition)) % 3;

      /* transition probability: Hexanucleotide */    
      scoreMarkov +=
	isochores[currentIsochore]->OligoDistTran[j][endexon-OligoLength_1+1] - 
	isochores[currentIsochore]->OligoDistTran[j][(iniexon)?iniexon-1:0];   

      if (SRP)
	{
	  /* SR score */
	  scoreSR = ScoreSRexon((Exons+i),Strand,sr);
	  scoreMarkov += scoreSR; 
	}

      /* First cutoff- Markov */
      if (scoreMarkov >= p->OligoCutoff) 
	{
	  /* 2.score about two sites of exon */
	  scoreTotal = 
	    ((1 - p->OligoWeight) * ((Exons+i)->Acceptor->Score 
				     + (Exons+i)->Donor->Score))
	    + (p->OligoWeight * scoreMarkov); 
	  
	  /* Second cutoff- Global score */
	  if (scoreTotal >= p->ExonCutoff) 
	    {
	      Exons[n]=Exons[i];
	      (Exons+n)->Score = scoreTotal  + ExonWeight;
	      (Exons+n)->PartialScore = scoreMarkov;

	      if (SRP)
		(Exons+n)->SRScore = scoreSR;

	      n++;
	    }
	}
    }
  return(n);
}


void ScoreExons(char *Sequence, 
		packExons* allExons, 
		long l1, long l2, int Strand, packSR* sr,
		gparam** isochores, int nIsochores)
{
  long i;
  char mess[MAXSTRING];
  long l;
  long* CG;
 
  /* Exons are between l1 and l2 */
  l = l2 - l1 + 1;

  /* Creating and filling CG array */
  if ((CG = (long *) calloc(l, sizeof(long))) == NULL)
    printError("Not enough space to hold CG array"); 

  printMess("Computing CG content");
  CGScan(Sequence, CG, l1, l2);

  /* Fill in the Markov arrays */
  for(i=0; i<nIsochores; i++)
    {
      /* Filling OligoDistIni and OligoDistTran */
      printMess("Computing Markov logscore profile");
      MarkovScan(Sequence, isochores[i], 
		 isochores[i]->OligoDistIni, 
		 isochores[i]->OligoDistTran, 
		 l1, l2);
    }
  
  /* 2. Ready for Score and Filter Exons? */
  printMess("Filtering exons using oligo scores");
  
  allExons->nInitialExons = Score(allExons->InitialExons,
				  allExons->nInitialExons,
				  l1, l2, Strand, sr, 
				  CG, isochores);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 
  
  allExons->nInternalExons=Score(allExons->InternalExons,
				 allExons->nInternalExons, 
				 l1, l2, Strand, sr, 
				 CG, isochores);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 
  
  allExons->nTerminalExons=Score(allExons->TerminalExons,
				 allExons->nTerminalExons,
				 l1, l2, Strand, sr, 
				 CG, isochores);
  sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
  printRes(mess);  
  
  allExons->nSingles=Score(allExons->Singles,
			   allExons->nSingles,
			   l1, l2, Strand, sr, 
			   CG, isochores);
  sprintf(mess,"Singles \t\t%8ld", allExons->nSingles);
  printRes(mess);

  
  /* Freeing CG memory */  
  free(CG);
}


