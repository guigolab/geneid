/*************************************************************************
*                                                                        *
*   Module: PrintExons                                                   *
*                                                                        *
*   Formatted printing of predicted exons.                               *
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

/*  $Id: PrintExons.c,v 1.3 2000-10-03 10:53:11 jabril Exp $  */

#include "geneid.h"

extern int GFF;

/* Print a predicted exon without gene prediction */
void PrintExon(exonGFF *e, char Name[], char* s, dict* dAA)
{
  char sAux[MAXAA];
  char* rs;
  long p1, p2;
  int nAA;
  
  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
  p2 = e->Donor->Position + e->offset2 - COFFSET;
  
  
  /* Translation of nucleotids of exon to aminoacids */
  if (e->Strand == '+')
    /* Translate exon codons to aminoacids */
    nAA = Translate(p1,p2,
		    e->Frame,
		    (3 - e->Remainder)%3,
		    s, dAA, sAux);
  else
    {
      /* Reverse strand exon */
      if ((rs = (char*) calloc(p2-p1+2,sizeof(char))) == NULL)
	printError("Not enough space to translate reverse exon");
      
      ReverseSubSequence(p1, p2, s, rs);
      
      nAA = Translate(0,p2-p1,
		      e->Frame,
		      (3 - e->Remainder)%3,
		      rs, dAA, sAux);
      free(rs);
    }
  
  /* According to selected options print formatted */
  if (GFF)
    {
      /* GFF format */
      printf ("%s\t%s\t%8s\t%ld\t%ld\t%5.3f\t%c\t%hd\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : EXONS,
	      e->Type,
	      e->Acceptor->Position + e->offset1,
	      e->Donor->Position + e->offset2,
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
	      e->Acceptor->Position + e->offset1,
	      e->Donor->Position + e->offset2,
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

/* Print a group of exons of the same type */
void PrintExons (exonGFF *e, long ne, int type, char Name[],
		 long l1, long l2, char* Sequence,  dict* dAA) 
{
  long i;
  char Type[MAXTYPE];
  char strand;
  
  strand = (ne>0)? e->Strand : 'x'; 
  
  Type[0] = '\0';
  switch(type)
    {
    case FIRST: strcpy(Type,"First exon");
       break;
    case INTERNAL:strcpy(Type,"Internal exon");
      break;
    case TERMINAL:strcpy(Type,"Terminal exon");
      break;
    default: strcpy(Type,"Exon");
      strand = 'x'; 
      break;
    }

  printf("# %ss(%c) predicted in sequence %s: [%ld,%ld]\n", 
	 Type,
	 strand,
	 Name, l1, l2);

  for (i=0;i<ne;i++)
    PrintExon((e+i), Name, Sequence, dAA);   
}


/* Print the exon of a gene */
void PrintGExon(exonGFF *e, char Name[], char* s, dict* dAA, 
	       long ngen, int AA1, int AA2, int nAA) 
{
  if (GFF)
    {
      /* GFF format */
      printf ("%s\t%s\t%8s\t%ld\t%ld\t%5.3f\t%c\t%hd\tgene_%ld\t# AA-%d\tAA-%d\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,     
	      e->Type,
	      e->Acceptor->Position + e->offset1,
	      e->Donor->Position + e->offset2,
	      (e->Score==MAXSCORE)? 0.0:e->Score,
	      e->Strand,
	      e->Frame,
	      ngen,
	      (e->Strand=='+')? nAA-AA2+COFFSET : AA1,
	      (e->Strand=='+')? nAA-AA1+COFFSET : AA2);
    }
  else
    {
      /* Default format for genes */
      printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\tAA %3d:%3d gene_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      e->Type,
	      e->Acceptor->Position + e->offset1,
	      e->Donor->Position + e->offset2,
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
