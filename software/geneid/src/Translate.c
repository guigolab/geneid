/*************************************************************************
*                                                                        *
*   Module: Translate                                                    *
*                                                                        *
*   Translate nucleotids to aminoacids.                                  *
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

/*  $Id: Translate.c,v 1.3 2001-04-23 13:11:49 eblanco Exp $  */

#include "geneid.h"

/* Translate a exon to aminoacids */
int Translate(long p1, long p2, short fra, short rmd,
              char *s, dict* dAA, char sAux[]) 
{
   int nAA;
   long i;
   char codon[LENGTHCODON+1];

   nAA = 0;
   /* Exon Frame */
   for(i=0; i<fra; i++)
     {
       sAux[i] = s[p1+i] + 32;
       nAA++;
     }

   /* Exon */
   for(i=p1+ fra; i<= p2 - rmd; i=i+3)
     {
     codon[0] = s[i];
     codon[1] = s[i+1];
     codon[2] = s[i+2];
     codon[3] = '\0';
     
     sAux[nAA] = getAADict(dAA,codon);
     nAA++;
     }
   
   /* Exon Remainder */
   for(i=p2 - rmd + 1; i <= p2; i++)
     {
       sAux[nAA] = s[i] + 32;
       nAA++;
     }
   sAux[nAA] = '\0';

   /* Corrections about aa not completely built */
   if (fra == 2)
     nAA--;
   if (rmd == 2)
     nAA--;
   
   return(nAA);
}

/* Translate a complete gene to protein */
void TranslateGen(exonGFF* e, char* s, dict* dAA, long nExons, 
		  int tAA[MAXEXONGEN][2], char* prot, long* nAA)
{  
  long i;
  short j;
  int totalAA;
  int currAA;
  char sAux[MAXAA];
  char* rs;
  long p1, p2;
  char codon[LENGTHCODON+1];
  char aa;
  short currFrame;
  short currRmd;
  int lAux;
  char rmdProt[LENGTHCODON+1];

  if (e->Strand == '+')
    {
      for(i=0, totalAA=0, currFrame=0; i<nExons; i++)
	{
	  p1 = (e->evidence)? 
	    e->Acceptor->Position + e->offset1 : 
	    e->Acceptor->Position + e->offset1 - COFFSET;
	  p2 = (e->evidence)? 
	    e->Donor->Position + e->offset2: 
	    e->Donor->Position + e->offset2 - COFFSET;

	  /* Remainder of the last exon from the gene */
	  if (!i)
	    {
	      switch (e->Remainder)
		{
		case 0:
		  rmdProt[0] = '\0'; 
		  break;
		case 1:
		  rmdProt[0] = s[p2] + 32; 
		  rmdProt[1] = '\0'; 
		  break;
		case 2:
		  /* curr RMD must be TWO */
		  rmdProt[0] = s[p2-1] + 32;
		  rmdProt[1] = s[p2] + 32;
		  rmdProt[2] = '\0'; 
		  break;
		}	
	      sprintf(prot,"%s",rmdProt);
	    }
	  
	  currAA = Translate(p1 + e->Frame,
			     p2 - (3 - e->Remainder)%3,
			     0,
			     0,
			     s, dAA, sAux);

	  /* Translating codon between this exon and its previous */
	  switch (currFrame)
	    {
	    case 1:
	      /* curr RMD must be TWO */
	      codon[0] = s[p2-1];
	      codon[1] = s[p2];
	      codon[3] = '\0';
	      /* translate this codon */
	      aa = getAADict(dAA,codon);
	      lAux = strlen(sAux);
	      sAux[lAux] = aa;
	      sAux[lAux+1] = '\0'; 
	      break;
	      
	    case 2:
	      /* curr RMD must be ONE */
	      codon[0] = s[p2];
	      codon[3] = '\0';
	      aa = getAADict(dAA,codon);
	      lAux = strlen(sAux);
	      sAux[lAux] = aa;
	      sAux[lAux+1] = '\0'; 
	      break;
	    }

	  /* Previous peptid after the translation of the current exon */
	  strcat(sAux,prot);
	  strcpy(prot,sAux);

	  currFrame = e->Frame;
	  /* Codon information for translating this shared codon next time */
	  switch (currFrame)
	    {
	    case 1:
	      codon[2] = s[p1];
	      break;
	      
	    case 2:
	      /* curr RMD must be TWO */
	      codon[1] = s[p1];
	      codon[2] = s[p1+1];
	      break;
	    }

	  /* Updating boundary aminoacids of this exon */
	  /* Not Shared codon */
	  tAA[i][0] = MAX(1,totalAA);
	  
	  if (!e->Remainder && i)
	    tAA[i][0]++;
	  	  
	  if (e->Frame)
	    totalAA++;
	  
	  totalAA = totalAA + currAA;
	  tAA[i][1] = totalAA;
	  
	  e = e->PreviousExon;
	} /* endfor */
      
      /* Frame of the first exon in sequence */
      switch (currFrame)
	    {
	    case 0:
	      rmdProt[0] = '\0'; 
	      break;
	    case 1:
	      rmdProt[0] = codon[2] + 32; 
	      rmdProt[1] = '\0'; 
	      break;
	    case 2:
	      /* curr RMD must be TWO */
	      rmdProt[0] = codon[1] + 32;
	      rmdProt[1] = codon[2] + 32;
	      rmdProt[2] = '\0'; 
	      break;
	    }
      
      sprintf(sAux,"%s%s",rmdProt,prot);
      strcpy(prot,sAux);
    }
  /* strand - */
  else
    {
      prot[0] = '\0'; 
      for(i=0, totalAA=0, currRmd=0; i<nExons; i++)
	{
	  p1 = (e->evidence)? 
	    e->Acceptor->Position + e->offset1: 
	    e->Acceptor->Position + e->offset1 - COFFSET;
	  p2 = (e->evidence)? 
	    e->Donor->Position + e->offset2: 
	    e->Donor->Position + e->offset2 - COFFSET;
	  
	  /* Reverse strand exon */
	  if ((rs = (char*) calloc(p2-p1+2,sizeof(char))) == NULL)
	    printError("Not enough space to translate reverse exon");
	 
	  ReverseSubSequence(p1, p2, s, rs);

	  /* Frame of the last exon in sequence */
	  if (!i)
	    {
	      switch (e->Frame)
		{
		case 0:
		  rmdProt[0] = '\0'; 
		  break;
		case 1:
		  rmdProt[0] = rs[0] + 32; 
		  rmdProt[1] = '\0'; 
		  break;
		case 2:
		  /* curr RMD must be TWO */
		  rmdProt[0] = rs[0] + 32;
		  rmdProt[1] = rs[1] + 32;
		  rmdProt[2] = '\0'; 
		  break;
		}	
	    }
	  
	  currAA = Translate(0 + e->Frame, p2-p1-(3 - e->Remainder)%3,
			  0,
			  0,
			  rs, dAA, sAux);
	  
	  /* Translating codon between this exon and its previous */
	  switch (currRmd)
	    {
	    case 1:
	      /* curr FRAME must be TWO */
	      codon[1] = rs[0];
	      codon[2] = rs[1];
	      codon[3] = '\0';
	      /* translate this codon */
	      aa = getAADict(dAA,codon);
	      lAux = strlen(prot);
	      prot[lAux] = aa;
	      prot[lAux+1] = '\0'; 
	      break;
	      
	    case 2:
	      /* curr FRAME must be ONE */
	      codon[2] = rs[0];
	      codon[3] = '\0';
	      aa = getAADict(dAA,codon);
	      lAux = strlen(prot);
	      prot[lAux] = aa;
	      prot[lAux+1] = '\0'; 
	      break;
	    }
	  
	  /* exon protein after remainder translated */
	  strcat(prot,sAux);
	  
	  /* Codon information for translating this shared codon the next time */
	  currRmd = (3 - e->Remainder)%3;
	  switch (currRmd)
	    {
	    case 1:
	      codon[0] = rs[p2-p1];
	      break;
	      
	    case 2:
	      /* curr RMD must be TWO */
	      codon[0] = rs[p2-p1-1];
	      codon[1] = rs[p2-p1];
	      break;
	    }
	  
	  /* Updating boundary aminoacids of this exon */

	  tAA[i][0] = MAX(1,totalAA);

  	  if (!e->Frame && i)
	    tAA[i][0]++;
	  
	  
	  if (currRmd)
	    totalAA++;
	  
	  totalAA = totalAA + currAA;
	  tAA[i][1] = totalAA;
	  
	  e = e->PreviousExon;	
	  free(rs); 
	}
      /* Frame of the last exon in sequence */
      sprintf(sAux,"%s%s",rmdProt,prot);
      strcpy(prot,sAux);
      
      /* Remainder of the first exon in sequence */
      switch (currRmd)
	{
	case 0:
	  rmdProt[0] = '\0'; 
	  break;
	case 1:
	  rmdProt[0] = codon[0] + 32; 
	  rmdProt[1] = '\0'; 
	  break;
	case 2:
	  /* curr RMD must be TWO */
	  rmdProt[0] = codon[0] + 32;
	  rmdProt[1] = codon[1] + 32;
	  rmdProt[2] = '\0'; 
	  break;
	}
      
      lAux = strlen(prot);
      for(j = 0; j < currRmd; j++)
	prot[lAux+j] = rmdProt[j];
      prot[lAux+j] = '\0'; 
    }
  *nAA = totalAA;  
}


void GetcDNA(exonGFF* e, char* s, long nExons, char* cDNA, long* nNN)
{  
  char* tmpDNA;
  long p1,p2;
  char* rs;
  int i;
  long j;

  if ((tmpDNA = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
    printError("Not enough space to store tmp cDNA");

  cDNA[0] = '\0';

  if (e->Strand == '+')
    {
      for(i=0, *nNN = 0; i<nExons; i++)
	{
	  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
	  p2 = e->Donor->Position + e->offset2 - COFFSET;
	  
	  /* Get the cDNA for this exon */
	  tmpDNA[0] = '\0';

	  for(j=p1; j <= p2; j++)
	    tmpDNA[j-p1] = s[j];
	  
	  tmpDNA[j-p1]='\0';
	  
	  /* Concat the current cDNA before the stored */
	  strcat(tmpDNA,cDNA);
	  strcpy(cDNA,tmpDNA);
	  
	  *nNN += p2-p1+1;
	  e = e->PreviousExon;
	}
    }

  /* Reverse and complement the sequence in this strand */
  else
    {
      if ((rs = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
	printError("Not enough space to reverse cDNA(-)");  

      for(i=0, *nNN = 0; i<nExons; i++)
	{
	  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
	  p2 = e->Donor->Position + e->offset2 - COFFSET;
	  
	  /* Get the cDNA for this exon */
	  tmpDNA[0] = '\0';
	  rs[0] = '\0';

	  for(j=p1; j <= p2; j++)
	    tmpDNA[j-p1] = s[j];
	  
	  tmpDNA[j-p1]='\0';
	  
	  ReverseSubSequence(0, j-p1-1, tmpDNA, rs);
	  rs[j-p1]='\0';

	  /* Concat the current cDNA after the stored */
	  strcat(cDNA,rs);

	  *nNN += p2-p1+1;
	  e = e->PreviousExon;
	}
    }
}

