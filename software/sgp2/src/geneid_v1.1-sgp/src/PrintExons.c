/*************************************************************************
*                                                                        *
*   Module: PrintExons                                                   *
*                                                                        *
*   Formatted printing of predicted exons                                *
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

/*  $Id: PrintExons.c,v 1.1 2003-09-10 14:53:34 gparra Exp $  */

#include "geneid.h"

extern int GFF;

/* Print a predicted exon from a list of exons */
void PrintExon(exonGFF *e, char Name[], char* s, dict* dAA)
{
  char sAux[MAXAA];
  char* rs;
  long p1, p2;
  int nAA;
  
  /* Acquiring real positions of exon */
  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
  p2 = e->Donor->Position + e->offset2 - COFFSET;
  
  /* Translation of exon nucleotides into amino acids */
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
  
  /* According to selected options formatted output */
  if (GFF)
    {
      /* GFF format */
      printf ("%s\t%s\t%s\t%ld\t%ld\t%5.2f\t%c\t%hd\n",
			  /* correct stop codon position, Terminal- & Terminal+ */ 
			  Name,
			  (e->evidence)? EVIDENCE : EXONS,
			  e->Type,
			  (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
			  (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
			  (e->Score==MAXSCORE)? 0.0: e->Score,
			  e->Strand,
			  e->Frame);
    }
  else
    {
      /* Default format */
      printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%d\t%s\n",
			  /* correct stop codon position, Terminal- & Terminal+ */ 
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
			  e->SRScore,
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
                 dict* dAA)
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
    PrintExon((e+i), Name, Sequence, dAA);   
}

/* Print a predicted exon from a assembled gene: gff/geneid format */
void PrintGExon(exonGFF *e,
                char Name[],
                char* s,
                dict* dAA,
                long ngen,
                int AA1,
                int AA2,
                int nAA)
{
  if (GFF)
    {
      /* GFF format */
      printf ("%s\t%s\t%s\t%ld\t%ld\t%5.2f\t%c\t%hd\tgene_%ld\n",
			  /* correct stop codon position, Terminal- & Terminal+ */ 
			  Name,
			  (e->evidence)? EVIDENCE : VERSION,     
			  e->Type,
			  (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
			  (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
			  (e->Score==MAXSCORE)? 0.0:e->Score,
			  e->Strand,
			  e->Frame,
			  ngen);
    }
  else
    {
      /* Default format for genes */
      printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\tAA %3d:%3d gene_%ld\n",
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
			  e->SRScore,
			  (e->Strand=='+')? nAA-AA2+COFFSET : AA1,
			  (e->Strand=='+')? nAA-AA1+COFFSET : AA2,
			  ngen);
    }
}

/* Print a predicted exon from a assembled gene: XML format */
void PrintXMLExon(exonGFF *e,
                  char Name[],
                  long ngen,
                  long nExon,
                  int type1,
                  int type2,
                  int nExons) 
{
  char Type[MAXTYPE];
  
  /* XML format */
  /* exon*/
  printf("      <exon idExon=\"%s.G%ldE%ld\" type=\"%s\" frame=\"%hd\" score=\"%.2f\">\n",
		 Name,ngen,nExons-nExon+1,e->Type,e->Frame,e->Score);
  
  /* Both sites */
  Type[0] = '\0';
  switch(type1)
    {
    case ACC: strcpy(Type,sACC); break;
    case STA: strcpy(Type,sSTA); break;
    case DON: strcpy(Type,sDON); break;
    case STO: strcpy(Type,sSTO); break;
    }
  
  printf("         <site idSite=\"%s.G%ldE%ldS1\" type=\"%s\" position=\"%ld\" score=\"%.2f\" />\n",
		 Name,ngen,nExons-nExon+1,Type,e->Acceptor->Position + e->offset1,e->Acceptor->Score);
  
  Type[0] = '\0';
  switch(type2)
    {
    case ACC: strcpy(Type,sACC); break;
    case STA:strcpy(Type,sSTA); break;
    case DON:strcpy(Type,sDON); break;
    case STO:strcpy(Type,sSTO); break;
    }
  
  printf("         <site idSite=\"%s.G%ldE%ldS2\" type=\"%s\" position=\"%ld\" score=\"%.2f\" />\n",
		 Name,ngen,nExons-nExon+1,Type,e->Donor->Position + e->offset2,e->Donor->Score);
  
  printf("      </exon>\n");
}
