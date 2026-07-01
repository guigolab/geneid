/*************************************************************************
*                                                                        *
*   Module: genamic                                                      *
*                                                                        *
*   Assembling genes from the input set of exons                         *
*                                                                        *
*   This file is part of the geneid 1.4 Distribution                     *
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

/*  $Id: genamic.c,v 1.20 2011-01-13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* This is the dynamic-programming core that chains scored exons into genes.
 *
 * Exons arrive in E[0..nExons) sorted by acceptor (left) position (done by
 * SortExons before this is called). genamic scans them once, left to right,
 * and for each exon looks up the best-scoring compatible partial gene that
 * could precede it, using two pieces of state carried in packGenes (pg)
 * across the whole run (all fragments/splits of one locus, reset only in
 * cleanGenes at the start of the next locus):
 *
 *   pg->d[etype]  - all exons seen so far that may end (via their Donor
 *                   site) a partial gene of gene-model class `etype`,
 *                   kept sorted by Donor position (built once per call by
 *                   BuildSort's insertion sort).
 *   pg->km[etype] - how many exons are currently in pg->d[etype].
 *   pg->je[etype] - how far pg->d[etype] has already been scanned into
 *                   pg->Ga (a forward-only cursor; see step 2 below).
 *   pg->Ga[etype][frame][spliceclass]
 *                 - the single best partial gene of class `etype` ending in
 *                   a given reading frame and acceptor splice class, i.e.
 *                   the DP table cell this exon will try to extend.
 *
 * "gene-model class" (etype, indexed 0..gp->nclass-1) is NOT the same as
 * exon type (First/Internal/Terminal/Single/...): it is a row of the .gm
 * gene-model file (read by ReadGeneModel) describing one legal exon-exon
 * transition, with its own min/max intron distance (gp->md/gp->Md) and
 * group-blocking rule (gp->block). gp->D maps an exon's "Type+Strand"
 * string to a dictionary key (the same numbering used by CookingGenes and
 * SortExons); gp->nc/gp->UC and gp->ne/gp->DE, indexed by that key, list
 * which gene-model classes an exon of this type may (a) contribute itself
 * into as a future predecessor (UC, "upstream compatible" -- this is what
 * BuildSort consults) and (b) look up as a predecessor for itself (DE,
 * "downstream equivalent" -- this is what the h-loop below iterates over).
 */

/* Complete gene prediction (sites and exons) or only assembling */
extern int GENEID;
extern int RSS;
extern float U12_SPLICE_SCORE_THRESH;
extern float U12_EXON_SCORE_THRESH;

/* E        exons to assemble, sorted by acceptor position (in/out: filled
 *          with GeneScore/PreviousExon chains on return)
 * nExons   number of exons in E
 * pg       DP state (Ga/d/km/je) persisted across the whole locus
 * gp       current isochore's parameters, including the gene model
 *          (D/nc/ne/UC/DE/md/Md/block/nclass) used to drive the DP */
void genamic(exonGFF* E, long nExons, packGenes* pg, gparam* gp)
{
  long i,j,j2;                 /* i: current exon in E; j/j2: cursors into pg->d[etype] */
  short frame,remainder,spliceclass,dclass; /* current exon's / candidate predecessor's DP-cell coordinates */
  int h;                       /* index over the current exon type's compatible predecessor classes (gp->DE) */
  char saux[MAXTYPE];          /* scratch "Type+Strand" string, looked up in the gene-model dictionary gp->D */
  int type,etype;              /* type: current exon's dictionary key; etype: gene-model class being tried */
  long MaxDist;                /* this class's max allowed intron/gap distance (gp->Md[etype]) */
  long MinDist;                /* this class's min allowed intron/gap distance (gp->md[etype]) */
  char mess[MAXSTRING];
  int current_exon_is_u12 = 0; /* is the current exon's acceptor a U12 (minor spliceosome) site? */
  int thresholdmet = 1;        /* does the candidate join clear the U2/U12 score threshold (see step 2c)? */

  /* 0. Starting process ... */
  printMess("-- Running gene assembling (genamic) --");

  /* geneid sends this set of exonsGFF together with the Ga and d-array*/
  sprintf(mess,"%ld exons GFF received from geneid",nExons);
  printMess(mess);

  /* Frame/remainder are stored for the forward (+) reading direction; on the
     reverse strand the DP below needs them read the other way round, so
     swap them here and swap back in step 3. Splits across fragments carry
     state in pg (Da = the previous fragment's leftover Ga/d, done once per
     call; GENEID off, i.e. -O, means genamic runs standalone on given
     exons/sites with no fragment history to carry). */
  if (GENEID)
    {
      /* Exons from the current fragment of sequence */
      SwitchFrames(E,nExons);

      /* Exons from the last fragment of sequence */
      SwitchFramesDa(pg,gp->nclass);
    }
  else
    {
      /* GENAMIC processing only: option -O */
      /* Exons from the current fragment of sequence */
      SwitchFrames(E,nExons);
    }
  
  /* 1. Add every exon from this fragment into pg->d[]: BuildSort inserts
     each exon (by insertion sort, by Donor position) into d[class] for
     every gene-model class it is UC-compatible with (see file header). */
  printMess("Sorting exons by donor...");

  /* Build a set of sorting exons by donor functions */
  BuildSort(gp->D, gp->nc, gp->ne,
	    gp->UC, gp->DE, gp->nclass,
	    pg->km, pg->dcap, pg->d, E, nExons);

  /* 2. Genamic Algorithm in linear time (size of input) */
  printMess("Assembling Genes...");

  /* For every exon verify non-blocking and distances well defined */
  /* and after this, it might be assemble with the best gene computed */
  for(i=0; i< nExons; i++)
    {
      /* What is the type of this exon? (skip the ghost/placeholder entries
         that carry an empty Type -- they are never assembled themselves) */
      if (strcmp((E+i)->Type,"")){
	saux[0]='\0';
	strcpy (saux, (E+i)->Type);
	strcat (saux, &((E+i)->Strand));
	type = getkeyDict(gp->D,saux);
	frame = (E+i)->Frame;
	/* Default: a lone exon, scoring just itself with no predecessor;
	   the h-loop below may improve on this by attaching a predecessor. */
	(E+i)->PreviousExon = pg->Ghost;
	(E+i)->GeneScore = (E+i)->Score;
	current_exon_is_u12 = 0;
	spliceclass = (E+i)->Acceptor->class;
	/* Classify by the current exon's acceptor signal: does joining it to a
	   predecessor cross a minor (U12) splice site, which is held to a
	   stricter score threshold below (step 2c)? */
	if (!strcmp((E+i)->Type,sEND) || !strcmp((E+i)->Type,sBEGIN)|| !strcmp((E+i)->Type,sSINGLE)){
	  current_exon_is_u12 = 0;
	}else{
	  if ((E+i)->Acceptor->class == U2){
	    current_exon_is_u12 = 0;
	  } else {
	    if (((E+i)->Acceptor->class == U12gtag)||((E+i)->Acceptor->class == U12atac)){
	      current_exon_is_u12 = 1;
	    }
	  }
	}
	if (type != NOTFOUND)
	  {
	    /* For every gene-model class this exon type may extend (gp->DE),
	       try to attach the best compatible predecessor found so far. */
	    for(h=0; h < gp->ne[type]; h++)
	      {
		etype = gp->DE[type][h];
		j = pg->je[etype];       /* where the forward scan (2b) left off last time */
		MaxDist = gp->Md[etype];
		MinDist = gp->md[etype];
		thresholdmet = 1;

		/* 2a. Has the current best-so-far predecessor for this DP cell
		   (pg->Ga[etype][frame][spliceclass]) fallen out of range, i.e.
		   is it now further back than MaxDist allows from this exon's
		   acceptor? (The "+ evidence - evidence" term nudges the
		   distance bound by +/-1 when exactly one side is an
		   evidence/annotation exon, to match their coordinate
		   convention against predicted signal coordinates.) If so it
		   can no longer be trusted, so re-derive the best predecessor
		   for every (remainder,dclass) cell by walking pg->d[etype]
		   backward from the current scan position (j) down to the
		   MaxDist boundary, keeping the best-scoring candidate seen
		   in each cell along the way. */
		if ((MaxDist != INFI) &&
		    (pg->Ga[etype][frame][spliceclass]->Strand !='*') &&
		    (pg->Ga[etype][frame][spliceclass]->Donor->Position
		     +
		     pg->Ga[etype][frame][spliceclass]->offset2)
		    <
		    ((E+i)->Acceptor->Position + (E+i)->offset1)
		    + ((E+i)->evidence - pg->Ga[etype][frame][spliceclass]->evidence)
		    - MaxDist)
		  {
		    /* loop backward searching another best gene matching MAX distance */
		    pg->Ga[etype][frame][spliceclass] = pg->Ghost;
		    j2 = j-1;
		    while (j2>=0 && j2 < pg->km[etype]  &&
			   ((pg->d[etype][j2]->Donor->Position
			     + pg->d[etype][j2]->offset2)
			    >=
			    ((E+i)->Acceptor->Position + (E+i)->offset1)
			    + ((E+i)->evidence - pg->d[etype][j2]->evidence)
			    - MaxDist))
		      {
			remainder = pg->d[etype][j2]->Remainder;
			dclass = pg->d[etype][j2]->Donor->class;

			if ((pg->d[etype][j2]->GeneScore >
			     pg->Ga[etype][remainder][dclass] -> GeneScore)
			    ){
			  pg->Ga[etype][remainder][dclass] = pg->d[etype][j2];
			}
			j2--;
		      }
		  }

		/* 2b. Forward scan: pg->d[etype] is sorted by Donor position, and
		   je[etype] (loaded into j above) remembers where the LAST call
		   left off, so this only ever advances -- each exon in
		   pg->d[etype] is folded into Ga exactly once, keeping the
		   whole outer loop linear in the total number of exons. Any
		   newly-reachable exon (Donor position now within MinDist of
		   this exon's acceptor) updates the best-partial-gene table
		   cell for its own (remainder,dclass), which is how later
		   exons will find it as a predecessor. */
		while(j < pg->km[etype] &&
		      ((pg->d[etype][j]->Donor->Position
			+ pg->d[etype][j]->offset2)
		       <=
		       ((E+i)->Acceptor->Position + (E+i)->offset1)
		       + ((E+i)->evidence - pg->d[etype][j]->evidence)
		       - MinDist)
		      )
		  {
		    remainder = pg->d[etype][j]->Remainder;
		    dclass = pg->d[etype][j]->Donor->class;
		    if ((frame == remainder && spliceclass == dclass &&
			 ((pg->d[etype][j]->Donor->Position 
			   + pg->d[etype][j]->offset2)
			  < 
			  ((E+i)->Acceptor->Position + (E+i)->offset1)
			  + ((E+i)->evidence - pg->d[etype][j]->evidence) 
			  - MaxDist))
			)
		      {
			/* Skip this exon because max distance not ok */
		      }
		    else
		      {
			if (pg->d[etype][j]->GeneScore > pg->Ga[etype][remainder][dclass]->GeneScore) 
			  {
			    pg->Ga[etype][remainder][dclass] = pg->d[etype][j];
			  }
		      }
		    j++;
		  }
		pg->je[etype] = j; /* remember the cursor for the next call */

		/* 2c. Assembling the exon with the best compatible gene before it */
		/* Verify group rules if there are evidence exons (annotations) */
		/* If this exon's acceptor is a minor (U12) site being joined to a
		   major (U2) donor, hold the join to a stricter combined score
		   threshold (unless either side is evidence-backed, which is
		   trusted outright); U12-to-U12 joins skip this check entirely
		   (the else branch below, thresholdmet stays 1). */
		if (current_exon_is_u12){
		  if (pg->Ga[etype][frame][spliceclass]->Donor->class != U2){				  
		    if((((pg->Ga[etype][frame][spliceclass]->Donor->Score + (E+i)->Acceptor->Score) > U12_SPLICE_SCORE_THRESH)
			&&
			((pg->Ga[etype][frame][spliceclass]->Score + (E+i)->Score) > U12_EXON_SCORE_THRESH)
			)
		       ||
		       (
			(E+i)->evidence || pg->Ga[etype][frame][spliceclass]->evidence
			)
		       ){
		      thresholdmet = 1;
		    } else {
		      thresholdmet = 0;
		    }
		  } else {
		    /* A minor-spliceosome acceptor can never follow a major (U2)
		       donor without clearing the threshold check above. */
		    thresholdmet = 0;
		  }
		} else {
		  /* Current exon is a normal (U2) acceptor: reject only if the
		     candidate predecessor's donor is itself minor-spliceosome
		     (a U12 donor must be matched with a U12 acceptor). */
		  if (pg->Ga[etype][frame][spliceclass]->Donor->class != U2){
		    thresholdmet = 0;
		  }
		}

		/* Group rule: if this gene-model class is BLOCKing, the two
		   exons must share the same evidence Group (or both be
		   NOGROUP) to be joined -- keeps predictions from mixing
		   incompatible annotation sources. Then require that joining
		   actually improves this exon's score, and that the U2/U12
		   threshold above was met. */
		if ((!(strcmp(pg->Ga[etype][frame][spliceclass]->Group,(E+i)->Group))
		     || gp->block[etype] == NONBLOCK)
		    &&
		    ((pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score) > (E+i)->GeneScore)
		    &&
		    (thresholdmet) 
		    )
		  {
			  
		    /* Introns are pure connectors (no sequence of their own to
		       join across a codon boundary): just inherit the predecessor's
		       split-codon flags unchanged rather than merging them. */
		    if (!strcmp((E+i)->Type,sINTRON)||!strcmp((E+i)->Type,sUTRINTRON)||!strcmp((E+i)->Type,sUTR5INTRON)||!strcmp((E+i)->Type,sUTR3INTRON)){			    
		      (E+i)->GeneScore = pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score;
		      (E+i)->PreviousExon = pg->Ga[etype][frame][spliceclass];
		      (E+i)->lValue = pg->Ga[etype][frame][spliceclass]->lValue;
		      (E+i)->rValue = pg->Ga[etype][frame][spliceclass]->rValue;
		    }else
		      /* Recursive/zero-length splice site: acceptor and donor are
			 adjacent, so no real coding bases sit between them either --
			 join unconditionally like the intron case above, but also
			 carry over Frame/Remainder since this exon added none of
			 its own. */
		      {if (RSS && ((E+i)->Donor->Position == (E+i)->Acceptor->Position -1)){
			  (E+i)->GeneScore = pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score;
			  (E+i)->PreviousExon = pg->Ga[etype][frame][spliceclass];
			  (E+i)->Frame = pg->Ga[etype][frame][spliceclass]->Frame;
			  (E+i)->Remainder = pg->Ga[etype][frame][spliceclass]->Remainder;
			  (E+i)->lValue = pg->Ga[etype][frame][spliceclass]->lValue;
			  (E+i)->rValue = pg->Ga[etype][frame][spliceclass]->rValue;
			}else{

			 /* Ordinary coding exon: joining splits a codon across the
			    predecessor's Donor and this exon's Acceptor. lValue/rValue
			    (set by ComputeStopInfo) encode which partial codon each
			    side's split leaves behind; if the two halves would spell a
			    stop codon (TAA/TAG/TGA) once joined, skip the join here and
			    keep this exon's standalone score from the top of the loop
			    instead. */
			 if (((E+i)->Strand == '+') && 
			     ((pg->Ga[etype][frame][spliceclass]->rValue == 1 && (E+i)->lValue == 1)
			      ||
			      (pg->Ga[etype][frame][spliceclass]->rValue == 2 && (E+i)->lValue == 2)
			      ||
			      (pg->Ga[etype][frame][spliceclass]->rValue == 3 && ((E+i)->lValue == 2 || (E+i)->lValue == 3))))
			   {
			     /* FWD: Avoiding building a stop codon */  
			   }
			 else
			   {
			     if (((E+i)->Strand == '-') && 
				 ((pg->Ga[etype][frame][spliceclass]->lValue == 1 && (E+i)->rValue == 1)
				  ||
				  (pg->Ga[etype][frame][spliceclass]->lValue == 2 && (E+i)->rValue == 2)
				  ||
				  ((pg->Ga[etype][frame][spliceclass]->lValue == 2 || pg->Ga[etype][frame][spliceclass]->lValue == 3) && (E+i)->rValue == 3)))
			       {
				 /* RVS: Avoiding building a stop codon */
			       }
			     else
			       {
				 (E+i)->GeneScore = pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score;
				 (E+i)->PreviousExon = pg->Ga[etype][frame][spliceclass];
			       }
			   }
		       }
		      }
		  }
	      }

	    /* Updating the best gene assembled so far for the whole locus:
	       track the single highest-scoring exon reached by any chain --
	       except a chain ending on a bare intron/UTR-intron (i.e. nothing
	       real was ever attached past the Ghost) can't stand alone as a
	       complete gene. */
	    if ((((E+i)->GeneScore) > (pg->GOptim -> GeneScore))){
	      if (((E+i)->PreviousExon->Strand == '*')&&(!strcmp((E+i)->Type,sINTRON) || !strcmp((E+i)->Type,sUTRINTRON) || !strcmp((E+i)->Type,sUTR5INTRON) || !strcmp((E+i)->Type,sUTR3INTRON))){
	      }else{
		pg->GOptim = (E+i);
	      }
	    }
	  }
      }
    }

  /* 3. Undo the frame/remainder swap from the prologue, on both this
     fragment's exons and the DP state carried over from the previous
     fragment (Db is SwitchFramesDa's matching undo), so downstream code
     (CookingGenes, Translate, ...) sees the normal forward-reading values. */
  if (GENEID)
    {
      /* Exons predicted in the current fragment of sequence */
      SwitchFrames(E,nExons);

      /* Only changing exons predicted in the last fragment of sequence */
      SwitchFramesDb(pg,gp->nclass);
    }
  else
    {
      /* GENAMIC processing only: option -O */
      UndoFrames(E,nExons);
    }
  
  /* Finishing process */
  printMess("-- Finishing gene assembling (genamic) --\n");
}

