/*************************************************************************
*                                                                        *
*   Module: BuildAcceptors.c                                             *
*                                                                        *
*   Signal prediction by using a Position Weighted Array                 *
*   using branch point or Poly Pyrimidine Tract profiles, if provided    *
*                                                                        *
*   This file is part of the geneid 1.2 distribution                     *
*                                                                        *
*     Copyright (C) 2003 - Enrique BLANCO GARCIA                         *
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

/* Function TRANS: char -> integer such that A=0, C=1, G=2 and T/U=2 */
extern int TRANS[];

/* Maximum allowed number of generic sites */
extern long NUMU12SITES;

/* Additional profiles */
extern int BP;
extern int PPT;

float ComputeU12BranchProfile(char* s,
						  long positionAcc,
						  long limitRight,
						  profile* p,
						  site* splicesite)
{
  float maxScore;
  float score;
  int index;
  int Opt;
  long end;
  long i,j;
  maxScore = -INF;
/*      char mess[MAXSTRING];  */
  i = MAX(p->order,positionAcc - ACCEPTOR_CONTEXT);
  end = MIN(positionAcc - MIN_BPACC_DIST + p->dimension - p->offset,limitRight);
  Opt = MAX(0,positionAcc - OPT_U12BP_DIST);
  for (;
	   i + p->dimension <= end;
	   i++)
	{
	  /* Applying the additional profile */
	  score=0.0;
	  for (j=0;j < p->dimension;j++)
		{
		  /* i is the position inside the region */
		  /* 5 is used because there are A,C,G,T and N */
		  index = OligoToInt(s + i + j - p->order, p->order+1,5);
		  
		  if (index >= p->dimensionTrans)
			score = score + -INFI;
		  else
			score = score + p->transitionValues[j][index];
		}
/*  	  score = score + 1.5 - score * (abs(i + p->offset -Opt)/(1.5*ACCEPTOR_CONTEXT - OPT_U12BP_DIST));  */
/*   	  sprintf(mess, "penalty \t%2.2f", ((float)(abs(i + p->offset -Opt))/((float)(ACCEPTOR_CONTEXT - p->offset - OPT_U12BP_DIST)))*((float)(abs(i + p->offset -1
	  -Opt))/((float)(ACCEPTOR_CONTEXT - p->offset - OPT_U12BP_DIST))));
	  printRes(mess); */
	   
	  score = score - U12BP_PENALTY_SCALING_FACTOR * (((float)(abs(i + p->offset -Opt))/((float)(ACCEPTOR_CONTEXT - p->offset - OPT_U12BP_DIST)))*((float)(abs(i + p->offset
	  -Opt))/((float)(ACCEPTOR_CONTEXT - p->offset - OPT_U12BP_DIST)))); 
	  if (score >= maxScore){
		maxScore = score;
		splicesite->PositionBP = i + p->offset - positionAcc;
		}
	}
  
  /* Cutoff for BranchPoint and PPtracts are useless */
  /* if (maxScore < p->cutoff) */
  /* 	maxScore = 0.0; */

  return maxScore;
}

float ComputePPTProfile(char* s,
						  long positionAcc,
						  long limitRight,
						  profile* p,
						  site* splicesite)
{
  float maxScore;
  float score;
  int index;
  long end;
  long i,j;

    
  maxScore = -INF;

  i = MAX(p->order,positionAcc + splicesite->PositionBP + 1);
  i = MAX(i,positionAcc - PPT_ACC_MAXDIST - p->dimension);
  end = MIN(positionAcc,limitRight);
  for (;
	   i + p->dimension <= end;
	   i++)
	{
	  /* Applying the additional profile */
	  score=0.0;
	  for (j=0;j < p->dimension;j++)
		{
		  /* i is the position inside the region */
		  /* 5 is used because there are A,C,G,T and N */
		  index = OligoToInt(s + i + j - p->order, p->order+1,5);
		  
		  if (index >= p->dimensionTrans)
			score = score + -INFI;
		  else
			score = score + p->transitionValues[j][index];
		}
	  
	  if (score >= maxScore){
		maxScore = score;
		splicesite->PositionPPT = i + p->offset - positionAcc;
		}
	}
  
  /* Cutoff for BranchPoint and PPtracts are useless */
  /* if (maxScore < p->cutoff) */
  /* 	maxScore = 0.0; */

  return maxScore;
}

/* Search for acceptor splice sites, using additional profiles */
long  BuildU12Acceptors(char* s,
					 char* type,
					 profile* u12_p,
					 profile* u12bp,
					 profile* ppt,
					 site* st, 
					 long l1, 
					 long l2) 
{ 
  int i,j;
  char* sOriginal;
  float score;
  float scoreBP;
  float scorePPT;
  long ns,is;
  long left,right;
  int index;

  /* Final number of predicted signals (that type) */
  ns = 0;

  /* Back-up the origin of the sequence */
  sOriginal = s;

  /* 1. Searching sites between beginning of the sequence and p->offset */
  if (!l1)
    {
      for (is = 0; is < u12_p->offset && (ns<NUMU12SITES); is++)
		{ 	  
			if (ns<NUMU12SITES){
		  	  
			  scorePPT = 0.0;
			  scoreBP = 0.0;
			  score=0.0;
			  /* Applying part of the profile */
			  for (i=u12_p->offset-is, j=0; i < u12_p->dimension; i++,j++) 
				{
				  /* i is the position inside the region */
				  index = OligoToInt(s+j, u12_p->order+1,5);

				  if (index >= u12_p->dimensionTrans)
					score = score + -INFI;
				  else
					score = score + u12_p->transitionValues[i][index];
				}
			  /* Using additional profiles */
			  scoreBP = ComputeU12BranchProfile(sOriginal,u12_p->offset-is,l2,u12bp,&st[ns]);
			  
		  		if (PPT)
				scorePPT = ComputePPTProfile(sOriginal,u12_p->offset-is,l2,ppt,&st[ns]);
				
			  	score = score + scoreBP;
		  
			  if (score >= u12_p->cutoff) 
				{
				  st[ns].Position = is + u12_p->order;
				  st[ns].ScoreBP = scoreBP;
			  	  st[ns].ScorePPT = scorePPT;
				  st[ns].Score=score;
				  /* st[ns].PositionBP= 5; */
			      strcpy(st[ns].subtype,type);
				  ns++;
				}
			}
		}
    }
  
  
  
  /* 2. Normal processing: predicting using the whole profile */
  /* left and right are the true boundaries of prediction */
  left  = MAX(0+u12_p->order, l1 - u12_p->offset);
  right = l2 - u12_p->offset;
  s += left;
  is = 0;     
  /* Case A: Using Markov chain with order 0: PWM */
  if (u12_p->order == 0)
    {
      /* discovering splice sites with current profile */
      while (*(s+u12_p->dimension) && (is < right- left + 1) && (ns<NUMU12SITES))
		{ 	
			if (ns<NUMU12SITES){

		  	  scorePPT = 0.0;
			  scoreBP = 0.0;
			  /* is = 0..right */
		  	  score=0.0;
		  	  for (i=0;i<u12_p->dimension;i++)
				{
			  	/* i is the position inside the region */
			  	index = TRANS[(int)(*(s + i))];
			  	if (index >= u12_p->dimensionTrans)
					score = score + -INFI;
			  	else
					score = score + u12_p->transitionValues[i][index];
				}
			  
			  /* Using additional profiles */
 			  scoreBP = ComputeU12BranchProfile(sOriginal,left + is + u12_p->offset,l2,u12bp,&st[ns]);  
			  
		  		if (PPT)
				scorePPT = ComputePPTProfile(sOriginal,left + is + u12_p->offset,l2,ppt,&st[ns]);
				
			  score = score + scoreBP;
			
			  if (score >= u12_p->cutoff) 
				{
				  st[ns].Position = left + is + u12_p->offset;
				  st[ns].ScoreBP = scoreBP;
			  	  st[ns].ScorePPT = scorePPT;
				  st[ns].Score=score;
				  /* st[ns].PositionBP= 5; */
				  strcpy(st[ns].subtype,type);
				  ns++;
				}
			}			
		  is++;
		  s++;
		}
    }
	
  /* case B: Using Markov chain with order 1: dinucleotides */
  else if (u12_p->order == 1)
    {

      /* discovering splice sites with current profile */
      while (*(s+u12_p->dimension) && (is < right- left + 1) && (ns<NUMU12SITES))
		{ 		
			if (ns<NUMU12SITES){
			/*Do for U12GTAG*/
			  
			  scorePPT = 0.0;
		  	  scoreBP = 0.0;
			  /* is = 0..right */
		  	  score=0.0;
		  	  for (i=0;i<u12_p->dimension;i++)
				{
			  	/* i is the position inside the region */
			  	index = 5*TRANS[(int)(*(s + i -1))] + TRANS[(int)(*(s + i))];
			  	if (index >= u12_p->dimensionTrans)
					score = score + -INFI;
			  	else
					score = score + u12_p->transitionValues[i][index];
				}
			  
			  /* Using additional profiles */
			  scoreBP = ComputeU12BranchProfile(sOriginal,left + is + u12_p->offset,l2,u12bp,&st[ns]);
			  
		  		if (PPT)
				scorePPT = ComputePPTProfile(sOriginal,left + is + u12_p->offset,l2,ppt,&st[ns]);
				
			  score = score + scoreBP;

			  if (score >= u12_p->cutoff) 
				{
				  st[ns].Position = left + is + u12_p->offset;
				  st[ns].ScoreBP = scoreBP;
			  	  st[ns].ScorePPT = scorePPT;
				  st[ns].Score=score;
				  /* st[ns].PositionBP= 5; */
				  strcpy(st[ns].subtype,type);
				  ns++;
				}
			}
		  is++;
		  s++;
		}
    }
  /* case C: Using Markov chain with order > 1 */
  else
    {
      /* discovering splice sites with current profile */
      while (*(s+u12_p->dimension) && (is < right- left + 1) && (ns<NUMU12SITES))
		{ 
			if (ns<NUMU12SITES){
			  
			  scorePPT = 0.0;
		  	  scoreBP = 0.0;
			  /* is = 0..right */
		  	  score=0.0;
		  	  for (i=0;i<u12_p->dimension;i++)
				{
			  	/* i is the position inside the region */
			  	/* 5 is used because there are A,C,G,T and N */
			    index = OligoToInt(s + i - u12_p->order, u12_p->order+1,5);
			  	if (index >= u12_p->dimensionTrans)
					score = score + -INFI;
			  	else
					score = score + u12_p->transitionValues[i][index];
				}
			  
			  /* Using additional profiles */
			  scoreBP = ComputeU12BranchProfile(sOriginal,left + is + u12_p->offset,l2,u12bp,&st[ns]);

		  		if (PPT)
				scorePPT = ComputePPTProfile(sOriginal,left + is + u12_p->offset,l2,ppt,&st[ns]);
				
			  score = score + scoreBP;
		  
			  if (score >= u12_p->cutoff) 
				{
				  st[ns].Position = left + is + u12_p->offset;
				  st[ns].ScoreBP = scoreBP;
			  	  st[ns].ScorePPT = scorePPT;
				  st[ns].Score=score;
				  /* st[ns].PositionBP= 5; */
				  strcpy(st[ns].subtype,type);
				  ns++;
				}
			}	
		  is++;
		  s++;
		}
    }
  
  if (ns >= NUMU12SITES)
    printError("Too many predicted sites: decrease RU12SITES parameter");
  
  return(ns);
}

 
  
