/*************************************************************************
*                                                                        *
*   Module: PrintExons                                                   *
*                                                                        *
*   Formatted printing of predicted exons                                *
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

/*  $Id: PrintExons.c,v 1.18 2007-04-04 12:40:12 talioto Exp $  */

#include "geneid.h"

extern int GFF;
extern int GFF3;

/* Print a predicted exon from a list of exons */
void PrintExon(exonGFF *e, char Name[], char* s, dict* dAA, char* GenePrefix)
{
  char sAux[MAXAA];
  char* rs;
  long p1, p2;
  int nAA = 0;
/*   char saux[MAXTYPE]; */
/*   char saux2[MAXTYPE]; */
  
  /* Acquiring real positions of exon */
  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
  p2 = e->Donor->Position + e->offset2 - COFFSET;
  
  /* Translation of exon nucleotides into amino acids */
  strcpy(sAux,"\0");
  if(e->Donor->Position >= e->Acceptor->Position){
    if (e->Strand == '+')
      /* Translate codons to amino acids */
      nAA = Translate(p1,p2,
		      e->Frame,
		      (3 - e->Remainder)%3,
		      s, dAA, sAux);
    else
      {
	/* Reverse strand exon */
	if ((rs = (char*) calloc(p2-p1+2,sizeof(char))) == NULL)
	  printError("Not enough memory: translating reverse exon");
      
	/* Reversing and complementing the exon sequence */
	ReverseSubSequence(p1, p2, s, rs);
      
	/* Translate codons to aminoacids */
	nAA = Translate(0,p2-p1,
			e->Frame,
			(3 - e->Remainder)%3,
			rs, dAA, sAux);
	free(rs);
      }
  }
  /* According to selected options formatted output */
  if (GFF)
    {
      /* GFF format */
      printf ("%s\t%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s_%s_%ld_%ld\n",
			  /* correct stop codon position, Terminal- & Terminal+ */ 
			  Name,
			  (e->evidence)? EVIDENCE : EXONS,
			  e->Type,
			  (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
			  (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
			  (e->Score==MAXSCORE)? 0.0: e->Score,
			  e->Strand,
			  e->Frame,
			  Name,e->Type,(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
			  (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2 );
    }
  else
    {
      /* Default format */
      printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%d\t%s\n",
			  /* correct stop codon position, Terminal- & Terminal+ */ 
			  /* saux2,  */
	      e->Type,
			  (e->evidence)? e->Acceptor->Position: e->Acceptor->Position + e->offset1,
			  (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
			  (e->Score==MAXSCORE)? 0.0:e->Score,
			  e->Strand,
			  e->Frame,
			  (3 - e->Remainder)%3,
			  e->Acceptor->Score,
			  e->Donor->Score,
			  e->PartialScore,
			  e->HSPScore,
			  nAA,
			  sAux);
    }
}

/* Print a list of exons of the same type */
void PrintExons (exonGFF *e,
                 long ne,
                 int type,
                 char Name[],
                 long l1, long l2,
                 char* Sequence,
                 dict* dAA,
		 char* GenePrefix)
{
  long i;
  char Type[MAXTYPE];
  char strand;
  
  strand = (ne>0)? e->Strand : 'x'; 
  
  Type[0] = '\0';
  switch(type)
    {
    case FIRST: strcpy(Type,sFIRST);
	  break;
    case INTERNAL:strcpy(Type,sINTERNAL);
      break;
    case ZEROLENGTH:strcpy(Type,sZEROLENGTH);
      break;
    case TERMINAL:strcpy(Type,sTERMINAL);
      break;
    case SINGLE:strcpy(Type,sSINGLE);
      break;  
    case ORF:strcpy(Type,sORF);
      break;  
    default: strcpy(Type,sEXON);
      strand = 'x'; 
      break;
    }
  
  /* Output header for the following output */
  printf("# %ss(%c) predicted in sequence %s: [%ld,%ld]\n",
		 Type,
		 strand,
		 Name, l1, l2);
  
  /* Print list of exons (type) predicted in the current fragment */
  for (i=0;i<ne;i++)
    PrintExon((e+i), Name, Sequence, dAA,GenePrefix);   
}

/* Print a predicted exon from a assembled gene: gff/geneid format */
void PrintGExon(exonGFF *e,
                char Name[],
                char* s,
                dict* dAA,
                long ngen,
                int AA1,
                int AA2,
                int nAA,
		int nExon,
		char* GenePrefix
		)
{
  if (e->Donor->Position == e->Acceptor->Position - 1){return;}
  char attribute[MAXSTRING] = "";
  if (GFF3)
    {
      /* GFF3 format 5_prime_partial=true ???*/
      if (e->five_prime_partial) { strcpy(attribute,";5_prime_partial=true");}
      if (e->three_prime_partial) { strcpy(attribute,";3_prime_partial=true");}

      /* It's really not necessary to print out exons unless we can predict UTR exons (leave this commented for now) */
/*       printf ("%s\t%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=exon_%s%s_%ld.%i;Parent=mRNA_%s%s_%ld;type=%s\n", */
/* 	      /\* correct stop codon position, Terminal- & Terminal+ *\/  */
/* 	      Name, */
/* 	      (e->evidence)? EVIDENCE : VERSION,      */
/* 	      "exon", */
/* 	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1, */
/* 	      (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2, */
/* 	      (e->Score==MAXSCORE)? 0.0:e->Score, */
/* 	      e->Strand, */
/* 	      e->Frame, */
/* 	      GenePrefix, */
/* 	      Name, */
/* 	      ngen, */
/* 	      nExon, */
/* 	      GenePrefix, */
/* 	      Name, */
/* 	      ngen, */
/* 	      e->Type); */
      printf ("%s\t%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=cds_%s%s_%ld.%i;Parent=mRNA_%s%s_%ld;Target=%s_predicted_protein_%s%s_%ld %d %d%s\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,     
	      "CDS",
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
	      (e->Score==MAXSCORE)? 0.0:e->Score,
	      e->Strand,
	      e->Frame,
	      GenePrefix,
	      Name,
	      ngen,
	      nExon,
	      GenePrefix,
	      Name,
	      ngen,
	      VERSION,
	      GenePrefix,
	      Name,
	      ngen,
	      (e->Strand=='+')? nAA-AA2+COFFSET : AA1,
	      (e->Strand=='+')? nAA-AA1+COFFSET : AA2,
	      attribute);
     
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,     
		e->Type,
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		GenePrefix,
		Name,
		ngen);
      } else {
	/* Default format for genes */
	printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\tAA %3d:%3d %s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		e->Type,
		(e->evidence)? e->Acceptor->Position :e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		(3 - e->Remainder)%3,
		e->Acceptor->Score,
		e->Donor->Score,
		e->PartialScore,
		e->HSPScore,
		(e->Strand=='+')? nAA-AA2+COFFSET : AA1,
		(e->Strand=='+')? nAA-AA1+COFFSET : AA2,
		GenePrefix,
		Name,
		ngen);
      }
    }
}

/* Print a predicted gene from a assembled gene: gff/geneid format */
void PrintGGene(exonGFF *s,
		exonGFF *e,
                char Name[],
                long ngen,
                float score,
		char* GenePrefix)
{
  if (GFF3)
    {
      /* GFF3 format */
      printf ("%s\t%s\tgene\t%ld\t%ld\t%.2f\t%c\t.\tID=%s%s_%ld;Name=%s%s_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,     
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? s->Donor->Position : s->Donor->Position + s->offset2,
	      score,
	      e->Strand,
	      GenePrefix,
	      Name,
	      ngen,
	      GenePrefix,
	      Name,
	      ngen); 
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s\tgene\t%ld\t%ld\t%.2f\t%c\t.\t%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,     
		(e->evidence)? s->Acceptor->Position : s->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		score,
		e->Strand,
		GenePrefix,
		Name,
		ngen);
      }
    }
}
/* Print a predicted mRNA from a assembled gene: gff/geneid format */
void PrintGmRNA(exonGFF *s,
		exonGFF *e,
                char Name[],
                long ngen,
                float score,
		char* GenePrefix)
{
  if (GFF3)
    {
		
      /* GFF3 format */
      printf ("%s\t%s\tmRNA\t%ld\t%ld\t%1.2f\t%c\t.\tID=mRNA_%s%s_%ld;Parent=%s%s_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,     
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? s->Donor->Position : s->Donor->Position + s->offset2,
	      score,
	      e->Strand,
	      GenePrefix,
	      Name,
	      ngen,
	      GenePrefix,
	      Name,
	      ngen);
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s\tgene\t%ld\t%ld\t%1.2f\t%c\t.\t%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,     
		(e->evidence)? s->Acceptor->Position : s->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		score,
		e->Strand,
		GenePrefix,
		Name,
		ngen);
      }
    }
}
/* Print a predicted intron from a assembled gene: gff/geneid format */
void PrintGIntron(exonGFF *d,
		  exonGFF *a,
		  char Name[],
		  long ngen,
		  int numInt,
		  char* GenePrefix,
		  int evidence)
{
  float score = 0.0;
  char intronType[MAXTYPE]; 
  strcpy(intronType,"U2");
  char intronSubtype[MAXTYPE]; 
  strcpy(intronSubtype,"GT-AG");
  short phase = (3 - d->Remainder)%3;
  /* short phase = MIN(0, 3 - a->Frame); */
  long start = (a->evidence)? d->Donor->Position + 1: d->Donor->Position + 1 + d->offset2;
  long end = (a->evidence)? a->Acceptor->Position - 1: a->Acceptor->Position -1 + a->offset1;
  if (evidence){
    score = a->Score;
  }else{
    score = (d->Donor->Score + a->Acceptor->Score)/2;
  }
  if (!(strcmp(d->Donor->subtype,sU12gtag))||!(strcmp(d->Donor->subtype,sU12atac))){
    strcpy(intronType,"U12");
  }
  if (!(strcmp(d->Donor->subtype,sU12gtag))){
    strcpy(intronSubtype,"GT-AG");
  }
  if (!(strcmp(d->Donor->subtype,sU12atac))){
    strcpy(intronSubtype,"AT-AC");
  }
  if (d->Strand == '-'){
    phase = (3 - a->Remainder)%3;
  }
  if (GFF3) {
    /* GFF3 format */
    printf ("%s\t%s\tintron\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=intron_%s%s_%ld.%i;Parent=%s%s_%ld;type=%s;subtype=%s\n",
	    /* correct stop codon position, Terminal- & Terminal+ */ 
	    Name,
	    (a->evidence)? EVIDENCE : VERSION,     
	    start,
	    end,
	    (score==MAXSCORE)? 0.0:score,
	    a->Strand,
	    phase,
	    GenePrefix,
	    Name,
	    ngen,
	    numInt,
	    GenePrefix,
	    Name,
	    ngen,
	    intronType,
	    intronSubtype
	    );
  } else {
    if (GFF){		 
      /* GFF format */
      printf ("%s\t%s\t%s_intron\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s%s_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (a->evidence)? EVIDENCE : VERSION,     
	      intronType,
	      start,
	      end,
	      (score==MAXSCORE)? 0.0:score,
	      a->Strand,
	      phase,
	      GenePrefix,
	      Name,
	      ngen);	  	
    }
  }

}

/* Print a predicted exon from a assembled gene: XML format */
void PrintXMLExon(exonGFF *e,
                  char Name[],
                  long ngen,
                  long nExon,
                  int type1,
                  int type2,
                  int nExons,
		  char* GenePrefix) 
{
  char Type[MAXTYPE];
  
  /* XML format */
  /* exon*/
  printf("      <exon idExon=\"%s%s.G%ldE%ld\" type=\"%s\" frame=\"%hd\" score=\"%.2f\">\n",
		 GenePrefix,Name,ngen,nExons-nExon+1,e->Type,e->Frame,e->Score);
  
  /* Both sites */
  Type[0] = '\0';
  switch(type1)
    {
    case ACC: strcpy(Type,sACC); break;
    case STA: strcpy(Type,sSTA); break;
    case DON: strcpy(Type,sDON); break;
    case STO: strcpy(Type,sSTO); break;
    }
  
  printf("         <site idSite=\"%s%s.G%ldE%ldS1\" type=\"%s\" position=\"%ld\" score=\"%.2f\" />\n",
		 GenePrefix,Name,ngen,nExons-nExon+1,Type,e->Acceptor->Position + e->offset1,e->Acceptor->Score);
  
  Type[0] = '\0';
  switch(type2)
    {
    case ACC: strcpy(Type,sACC); break;
    case STA:strcpy(Type,sSTA); break;
    case DON:strcpy(Type,sDON); break;
    case STO:strcpy(Type,sSTO); break;
    }
  
  printf("         <site idSite=\"%s%s.G%ldE%ldS2\" type=\"%s\" position=\"%ld\" score=\"%.2f\" />\n",
		 GenePrefix,Name,ngen,nExons-nExon+1,Type,e->Donor->Position + e->offset2,e->Donor->Score);
  
  printf("      </exon>\n");
}
