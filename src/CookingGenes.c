/*************************************************************************
*                                                                        *
*   Module: CookingGenes                                                 *
*                                                                        *
*   Formatted output of geneid (GFF, default and extended)               *
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

/*  $Id: CookingGenes.c,v 1.7 2001-04-30 13:02:17 eblanco Exp $  */

#include "geneid.h"

extern int X10;
extern int GFF;
extern int XML;
extern int cDNA;

void PrintExonGFF (exonGFF *e, char Name[], char Source[])
{
  printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\t(%d)\n",
         Name,
         Source,
         e->Type,
         e->Acceptor->Position + e->offset1,
         e->Donor->Position + e->offset2,
         e->Score,
         e->Strand,
         e->Frame,
	 e->Group);
}

void PrintGeneGFF(exonGFF *e, char Name[], char Source[])
{
  if ((e-> PreviousExon)->Strand != '*') 
    PrintGeneGFF(e->PreviousExon,Name,Source);
  PrintExonGFF (e,Name, Source);
}

/* Data Structure to store stats about a gen */
typedef struct s_gen
{
  long nexons;
  double score;
  exonGFF *start;
  exonGFF *end;
} gen;

/* Printing a protein (Fasta format) */
void printProt(char* Name,long ngen, char* prot, long nAA, int mode)
{
  long j;
  char header[MAXLINE];

  if (mode == PROT)
    sprintf(header,"\n>%s|%s_predicted_protein_%ld|%ld_AA\n",
	    Name,
	    VERSION,
	    ngen,
	    nAA);
  else
    sprintf(header,"\n>%s|%s_predicted_cDNA_%ld|%ld_NN\n",Name,
	    VERSION,
	    ngen,
	    nAA);

  if (!XML)
    printf("%s",header);
  else
     printf("\t");
	  
  for(j=0; j < strlen(prot); j++)
    {
      printf("%c",prot[j]);
      if (!((j+1)%60))
	{
	   printf("\n");
	   if (XML)
	      printf("\t");
	}
    }
  printf("\n");
  if (!XML && (mode == PROT))
     printf("\n");  	
}

/* Select splice site profiles according to the type of a exon */
void selectFeatures(char* exonType, char exonStrand,
                    profile** p1, profile** p2,
                    int* type1, int* type2,
                    int* strand, gparam* gp)
{
  if (exonStrand == '+')
   {
     *strand = FORWARD;
     if (!strcmp(exonType,"First"))
       {
         *type1 = STA;
         *type2 = DON;
         *p1 = gp->StartProfile;
         *p2 = gp->DonorProfile;
       }
     if (!strcmp(exonType,"Internal"))
       {
         *type1 = ACC;
         *type2 = DON;
         *p1 = gp->AcceptorProfile;
         *p2 = gp->DonorProfile;    
       }
     if (!strcmp(exonType,"Terminal"))
       {
         *type1 = ACC;
         *type2 = STO;
         *p1 = gp->AcceptorProfile;
         *p2 = gp->StopProfile;    
       }
     if (!strcmp(exonType,"Single"))
       {
         *type1 = STA;
         *type2 = STO;
         *p1 = gp->StartProfile;
         *p2 = gp->StopProfile;    
      }
   }
  /* Reverse strand */
  else
   {
     *strand = REVERSE;
     if (!strcmp(exonType,"First"))
       {
         *type2 = STA;
         *type1 = DON;
         *p2 = gp->StartProfile;
         *p1 = gp->DonorProfile;    
      }
     if (!strcmp(exonType,"Internal"))
       {
         *type2 = ACC;
         *type1 = DON;
         *p2 = gp->AcceptorProfile;
         *p1 = gp->DonorProfile;    
      }
     if (!strcmp(exonType,"Terminal"))
       {
         *type2 = ACC;
         *type1 = STO;
         *p2 = gp->AcceptorProfile;
         *p1 = gp->StopProfile;    
       }
     if (!strcmp(exonType,"Single"))
       {
         *type2 = STA;
         *type1 = STO;
         *p2 = gp->StartProfile;
         *p1 = gp->StopProfile;    
       }
   }
}

/* Processing of genes. After this, cooked printing */
long CookingInfo(exonGFF *eorig, gen info[], long* nvExons)
{
  long igen=0;
  int stop,stop1,stop2;
  exonGFF *e;
  
  /* Pointer jumping back travelling over the exones of multiple genes */
  *nvExons = 0;
  e = eorig;
  for(igen=0; igen < MAXGENE; igen++)
    {
      info[igen].nexons = 0;
      info[igen].score = 0;
    }
  igen = 0;
  
  stop = (e->Strand == '*');
  while (!stop)
    {
      /* Single Genes */
      if (!strcmp(e->Type,"Single"))
	{
	  info[igen].start = e;
	  info[igen].end = e;
	  info[igen].nexons = 1;
	  if (e->Score==MAXSCORE)
	    (*nvExons)++;
	  else
	    info[igen].score = e->Score;
	  e = (e-> PreviousExon);
	  stop = (e->Strand == '*');
	}
      else
	{
	  /* Building a reverse-gene */
	  if (e->Strand == '-')
	    {
	      info[igen].start = e;
	      info[igen].end = e;
	      info[igen].nexons++;
	      if (e->Score==MAXSCORE)
		(*nvExons)++;
	      else
		info[igen].score += e->Score;
	      e = (e-> PreviousExon);
	      
	      /* stop means end of gene-sequence */
	      stop = (e->Strand == '*');
	      /* stop1 means change of gene */
	      stop1 = (!strcmp(e->Type,"First") || 
		       !strcmp(e->Type,"Single") ||  
		       e->Strand == '+'); 
	      while( !stop && !stop1 )
		{  
		  if (strcmp(e->Type,"Promoter"))  
		    {
		      info[igen].nexons++;
		      if (e->Score==MAXSCORE)
			(*nvExons)++;
		      else
			info[igen].score += e->Score;
		      info[igen].end = e;
		    }
		  e = (e-> PreviousExon);
		  stop = (e->Strand == '*');
		  stop1 = (!strcmp(e->Type,"First") ||  
			   !strcmp(e->Type,"Single") ||  
			   e->Strand == '+'); 
		}	
	    }
	  else
	    /* Building a forward-gene */
	    if (e->Strand == '+')
	      {
		info[igen].start = e;
		info[igen].end = e;
		info[igen].nexons++;
		if (e->Score==MAXSCORE)
		  (*nvExons)++;
		else
		  info[igen].score += e->Score;
		e = (e-> PreviousExon);
		
		stop = (e->Strand == '*');
		/* stop1 means change of gene */
		stop2 = (!strcmp(e->Type,"Terminal") ||  
			 !strcmp(e->Type,"Single") ||
			 e->Strand == '-'); 
		while( !stop && !stop2 )
		  { if (strcmp(e->Type,"Promoter"))  
		    {
		      info[igen].nexons++;
		      if (e->Score==MAXSCORE)
			(*nvExons)++;
		      else
			info[igen].score += e->Score;
		      info[igen].end = e;
		    }
		  e = (e-> PreviousExon);
		  stop = (e->Strand == '*');
		  stop2 = (!strcmp(e->Type,"Terminal") ||
			   !strcmp(e->Type,"Single") ||
			   e->Strand == '-'); 
		  }	
	      }
	}
      igen++;
    } 
  return (igen);
}

/* Print a gene according to formatted output selected */
void PrintGene(exonGFF* start, exonGFF* end, char Name[],
		char* s, gparam* gp, dict* dAA, long igen, long nAA,
		int** tAA, int nExon, int nExons) 
{
  exonGFF* eaux;
  profile* p1;
  profile* p2;
  int type1;
  int type2;
  int strand;

  if (start != end)
    {
      eaux = start -> PreviousExon;
      PrintGene(eaux,end,Name,s,gp,dAA,igen,nAA,tAA,nExon+1,nExons);
      if (XML)
        {
	   selectFeatures(start->Type,start->Strand,
			  &p1,&p2,&type1,&type2,&strand,gp);
	   PrintXMLExon(start,Name,igen,nExon+1,type1,type2,nExons);
	}
      else 
	if (X10)
	  {
	    /* Print both sites of exon Start */
	    selectFeatures(start->Type,start->Strand,
			   &p1,&p2,&type1,&type2,&strand,gp);
	    PrintSite(start->Acceptor,type1,Name,strand,s,p1);
	    PrintGExon(start,Name,s,dAA,igen,tAA[nExon][0],tAA[nExon][1],nAA);
	    PrintSite(start->Donor,type2,Name,strand,s,p2);
	  }
	else
	  PrintGExon(start,Name,s,dAA,igen,tAA[nExon][0],tAA[nExon][1],nAA);
    }
  else
    { 
      if (XML)
        {
	  selectFeatures(end->Type,end->Strand,
			 &p1,&p2,&type1,&type2,&strand,gp);
	  PrintXMLExon(end,Name,igen,nExon+1,type1,type2,nExons);
	}
      else 
	if (X10)
	  {
	    /* Print both sites of exon End */
	    selectFeatures(end->Type,end->Strand,
			   &p1,&p2,&type1,&type2,&strand,gp);
	    PrintSite(end->Acceptor,type1,Name,strand,s,p1);
	    PrintGExon(end,Name,s,dAA,igen,tAA[nExon][0],tAA[nExon][1],nAA);
	    PrintSite(end->Donor,type2,Name,strand,s,p2);
	  }
	else
	  PrintGExon(end,Name,s,dAA,igen,tAA[nExon][0],tAA[nExon][1],nAA);
    }
} 

void CookingGenes(exonGFF *e, char Name[], char* s,
                  gparam* gp, dict* dAA)
{
  long igen;
  long ngen;
  gen* info;
  char* prot;
  char* tmpDNA;
  long nAA, nNN;
  int** tAA;
  long nvExons;
  long i;

  /* Get info about each gen */
  if ((info = (gen*) calloc(MAXGENE,sizeof(gen))) == NULL)
    printError("Not enough space to store info_gen structure");
  
  if ((tAA = (int**) calloc(MAXEXONGENE,sizeof(int*))) == NULL)
    printError("Not enough space to store tAA general structure");
  
  for(i=0; i<MAXEXONGENE; i++)
    if ((tAA[i] = (int*) calloc(2,sizeof(int))) == NULL)
      printError("Not enough space to store tAA[] structure");
  
  ngen = CookingInfo(e,info,&nvExons);

  /* Protein */
  if ((prot = (char*) calloc(MAXAA,sizeof(char))) == NULL)
    printError("Not enough space to store protein");

  /* cDNA memory */
  if (cDNA)
    if ((tmpDNA = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
      printError("Not enough space to store cDNA");

  /* Header multiple gene */
    /* Header for  */
  if (XML)
    printf(" genes=\"%ld\" score =\"%.2f\">\n", 
	   ngen,e -> GeneScore - nvExons*MAXSCORE); 
  else
    printf("# Optimal Gene Structure. %ld genes. Score = %.6f \n", 
	   ngen,e -> GeneScore - nvExons*MAXSCORE); 

  for(igen=ngen-1; igen>=0; igen--)
    {
      /* Translate gen to protein */
      TranslateGen(info[igen].start,s,dAA,info[igen].nexons,tAA,prot,&nAA);
 
      if (cDNA)
	GetcDNA(info[igen].start,s,info[igen].nexons, tmpDNA, &nNN);
      
      /* Header gene */
      if (XML)
	printf("   <gene idGene=\"%s.G%ld\" strand =\"%s\" exons=\"%ld\" score=\"%.2f\">\n",
	       Name,
	       ngen-igen,
	       (info[igen].start->Strand == '+')? xmlFORWARD : xmlREVERSE, 
	       info[igen].nexons,
	       info[igen].score);
      else     
	printf("# Gene %ld (%s). %ld exons. %ld aa. Score = %f \n",
	       ngen-igen,
	       (info[igen].start->Strand == '+')? sFORWARD : sREVERSE,
	       info[igen].nexons,
	       nAA,
	       info[igen].score);

      PrintGene(info[igen].start, info[igen].end, Name, s, gp, dAA, ngen-igen,
		nAA,tAA,0,info[igen].nexons);

      /* [cDNA] and translated protein */
      if (XML)
      {
	if (cDNA)
	  {
	    printf("      <cDNA length=\"%ld\">\n",nNN);
	    /* cDNA in FASTA format */
	    printProt(Name,ngen-igen,tmpDNA,nNN,cDNA);
	    printf("      </cDNA>\n");
	  }
	
	printf("      <protein length=\"%ld\">\n",nAA);
	/* Protein in FASTA format */
	printProt(Name,ngen-igen,prot,nAA,PROT);
	printf("      </protein>\n");
	printf("   </gene>\n");
      }
      else
	if (!(GFF))
	  {
	    /* cDNA */
	    if (cDNA)
	      printProt(Name,ngen-igen,tmpDNA,nNN,cDNA);
	    
	    /* Protein in FASTA format */
	    printProt(Name,ngen-igen,prot,nAA,PROT);
	  }
    }

  /* Free operations */
  free(prot);
  free(info);
  for(i=0; i<MAXEXONGENE; i++)
    free(tAA[i]);
  if (cDNA)
    free(tmpDNA);
}
