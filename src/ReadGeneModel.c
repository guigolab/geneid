/*************************************************************************
*                                                                        *
*   Module: ReadGeneModel                                                *
*                                                                        *
*   Reading the rules used to assemble best predicted genes              *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
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

/*  $Id: ReadGeneModel.c,v 1.9 2011-01-13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* The gene model (isochores[0]->D/nc/ne/UC/DE/md/Md/block/nclass) is only  */
/* ever loaded once, from the .gm/param file, into the first isochore's     */
/* gparam struct (see ReadGeneModel below). Every other isochore uses the   */
/* exact same assembly rules, so this just points their pointers/arrays at  */
/* isochore 0's copy instead of re-parsing or duplicating the data.         */
void shareGeneModel(gparam** isochores, int nIsochores)
{ 
  int i,j,k;
  int nTypes;
  
  /* Original gene model is loaded in the first isochore */
  nTypes = isochores[0]->D->nextFree; 
  
  /* Sharing global parameters for any isochore: i */
  for(i=1; i < nIsochores; i++)
    {
      isochores[i]->D      = isochores[0]->D; 
      isochores[i]->nc     = isochores[0]->nc; 
      isochores[i]->ne     = isochores[0]->ne; 
      isochores[i]->md     = isochores[0]->md;
      isochores[i]->Md     = isochores[0]->Md;
      isochores[i]->nclass = isochores[0]->nclass;
      
      /* Copying info for every feature/exon type */
      for(j=0; j < nTypes; j++)
		{
		  /* Frequency as left part on a rule (id. rule) */
		  for(k=0; k < isochores[i]->nc[j]; k++)
			isochores[i]->UC[j][k] = isochores[0]->UC[j][k];
		  
		  /* Frequency as right part on a rule (id. rule) */
		  for(k=0; k < isochores[i]->ne[j]; k++)
			isochores[i]->DE[j][k] = isochores[0]->DE[j][k];
		}
	  
      /* Copy info for every class: rules requiring group checkpoint */
      for(j=0; j < isochores[0]->nclass; j++)
		isochores[i]->block[j]  = isochores[0]->block[j];
    }
}

/* Loading the gene model rules to build correct genes.
 *
 * The gene model section of the parameter file lists, one per (non-comment,
 * non-blank) line, the legal ways an already-built exon of one type may be
 * followed by an exon of another type to form part of a gene. Each such
 * line is a "rule", and the rule's line number (0-based, in file order) IS
 * its class id ("nlines" below) -- the same class id that genamic.c's DP
 * uses to index Ga[class]/d[class] and that this function's own md/Md/
 * block arrays are indexed by. There is no separate class-name column:
 * the class is simply "the n-th rule in the file".
 *
 * Rule syntax (columns separated by literal space characters):
 *
 *   UC1:UC2:(...):UCn   DE1:DE2:(...):DEm   dMin:dMax   [block]
 *
 * - Column 1 (UC, "upstream compatible"): colon-separated list of exon
 *   types that may end a partial gene being extended INTO this class, i.e.
 *   the predecessor exon's type must appear here for genamic to consider
 *   this rule at all (see UC[][] below, and its use as gp->UC in genamic.c).
 * - Column 2 (DE, "downstream equivalent"): colon-separated list of exon
 *   types that may be the new exon completing this class once the distance
 *   test passes (see DE[][] below, gp->DE in genamic.c).
 * - Column 3: the allowed gap between the two exons, as dMin:dMax
 *   (introns use this as the intron-length window). Use the literal string
 *   SINFI ("Infinity", case-sensitive) in place of dMax to mean unbounded.
 * - Column 4 (optional): if this column is present AT ALL, its content is
 *   never inspected -- any token here (traditionally "block") switches this
 *   rule's block[] entry to BLOCK, meaning genamic must not let the two
 *   exons belong to different evidence Groups (group checkpoint). Omitting
 *   the column entirely leaves it NONBLOCK.
 *
 * Common parameter-file pitfalls with this format:
 * - Columns must be separated by literal space characters, not tabs --
 *   strtok(line," ") only splits on ' ', so a tab between columns leaves
 *   the whole line stuck together as line1 and triggers the "Wrong format"
 *   error below (line2/line3 come back NULL).
 * - A line is only treated as a comment if '#' or '\n' is the very FIRST
 *   character (line[0]). Leading whitespace before a '#' is NOT recognized
 *   as a comment and will instead be parsed as a (malformed) rule.
 * - Feature type names in columns 1-2 (e.g. "First+", "Internal-",
 *   "Single+", "Terminal-") are matched against the shared type dictionary
 *   by exact, case-sensitive string comparison (setkeyDict/getkeyDict). A
 *   typo or wrong case does NOT raise an error -- it silently registers a
 *   brand-new, never-produced type, so the rule quietly becomes dead/
 *   unreachable instead of failing loudly. Copy existing type strings
 *   rather than retyping them.
 * - dMax must be either a plain integer or exactly the string SINFI; any
 *   other non-numeric token there is misread as 0 via atol() with no error.
 *
 * Returns how many rules (classes) were loaded. */
long ReadGeneModel (FILE* file, dict* d,
                    int nc[], int ne[],
                    int UC[][MAXENTRY],
                    int DE[][MAXENTRY],
                    long md[], long Md[],
                    int block[])
{
  char line[MAXLINE];
  char lineCopy[MAXLINE];
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  
  /* Identifier for feature (from dictionary) */
  int a;
  
  /* Identifier for class (assembling rule) */
  int nlines;
  
  char mess[MAXSTRING];
  char *t1;

  /* Input lines from parameter file */
  nlines=0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      /* For every line extracting features (upstream/downstream), */
      /* the minMax distances and the (optional) block */
      /* line number is the class/rule identifier */
      if (line[0] !='#' && line[0] !='\n')
	{
	  /* 0. Backup the line to display errors */
	  strcpy(lineCopy,line);
		  
	  /* 1. Splitting line into 4 parts: UC DE Dist block */
	  line1 = (char *) strtok(line," ");
	  line2 = (char *) strtok(NULL," ");
	  line3 = (char *) strtok(NULL," "); 
	  line4 = (char *) strtok(NULL," ");
	  	  
	    /* Three first columns are mandatory, last one is optional */
	    if (line1 == NULL || line2 == NULL || line3 == NULL)
	      {
		sprintf(mess,"Wrong format in gene model rule:\n%s",lineCopy);
		printError(mess);
	      }
		  
	    /* 2. Processing upstream compatible features list */
	    for ( t1 =(char *) strtok(line1,":");
		  t1 != NULL;
		  t1 = (char *) strtok(NULL, ":") )
	      { 
		/* Extracting and adding to the dictionary of types */
		a = setkeyDict(d,t1);
		/* Save: this feature appeared in the UC list of this rule */
		UC[a][nc[a]++] = nlines;
	      }
		  
	    /* 3. Processing downstream equivalent features list */
	    for ( t1 =(char *) strtok(line2,":");
		  t1 != NULL;
		  t1 = (char *) strtok(NULL, ":") )
	      { 
		/* Extracting and adding to the dictionary of types */
		a = setkeyDict(d,t1);
		/* Save: this feature appeared in the DE list of this rule */
		DE[a][ne[a]++]=nlines;
	      } 
		  
	    /* 4. Read the distances xx:yy [block] */
	    /* a) minimum distance */
	    t1 =(char *) strtok(line3,":");      
	    if (t1 == NULL)
	      {
		sprintf(mess,"Wrong distance range (min) in gene model rule:\n%s",lineCopy);
		printError(mess);
	      }     
	    md[nlines] = atol(t1); 
		  
	    /* b) maximum distance */
	    t1 = (char *) strtok(NULL, ":");      
	    if (t1 == NULL)
	      {
		sprintf(mess,"Wrong distance range (max) in gene model rule:\n%s",lineCopy);
		printError(mess);
	      }          
	    /* To forget the DMAX requirement use the string SINFI */
	    /* Extracting \n in case there aren't block word behind */
	    if (t1[strlen(t1)-1] == '\n') 
	      t1[strlen(t1)-1] = '\0';      
	    if (!(strcmp(t1,SINFI)))
	      Md[nlines] = INFI;
	    else
	      Md[nlines] = atol(t1); 
		  
	    /* 5. Read the block record (to preserve group)... if exists */
	    if (line4 != NULL) 
	      block[nlines] = BLOCK;
	    else 
	      block[nlines] = NONBLOCK;
		  
	    nlines++;
	  
	} /* End of if-comment */
    } /* Next rule to read */
  
  return(nlines);
}

/* Used instead of ReadGeneModel when forcing single-gene prediction (-S):
 * hardcodes the same 4 rules (UC/DE/md/Md/block, same field meanings as
 * documented above ReadGeneModel) that a "normal" one-gene-per-model .gm
 * file would express, without needing one on disk. Rules 1-2 chain
 * First->Internal->Terminal exons (forward/reverse) with the usual
 * 20:40000 intron window and group blocking; rules 3-4 close the gene off
 * at the fragment boundaries (sBEGINFWD/sBEGINRVS, sENDFWD/sENDRVS) with
 * an unbounded (0:Infinity) gap and no group check.
 * Returns how many rules (always 4) were loaded. */
long ForceGeneModel (dict* d,
                    int nc[], int ne[],
                    int UC[][MAXENTRY],
                    int DE[][MAXENTRY],
                    long md[], long Md[],
                    int block[])
{
  /* Identifier for feature (from dictionary) */
  int a;
  
  /* Identifier for class (assembling rule) */
  int nlines;
  
  printMess("Force one gene prediction");

  nlines = 0;
 
  /* 1. First+:Internal+     Internal+:Terminal+       20:40000 block */
  a = setkeyDict(d,"First+");
  UC[a][nc[a]++] = nlines;
  a = setkeyDict(d,"Internal+");
  UC[a][nc[a]++] = nlines;

  a = setkeyDict(d,"Internal+");
  DE[a][ne[a]++]=nlines;
  a = setkeyDict(d,"Terminal+");
  DE[a][ne[a]++]=nlines;

  md[nlines] = 20;
  Md[nlines] = 40000;
  block[nlines] = BLOCK;
  nlines++;

  /* 2. Terminal-:Internal-  First-:Internal-          20:40000 blockr */
  a = setkeyDict(d,"Terminal-");
  UC[a][nc[a]++] = nlines;
  a = setkeyDict(d,"Internal-");
  UC[a][nc[a]++] = nlines;

  a = setkeyDict(d,"First-");
  DE[a][ne[a]++]=nlines;
  a = setkeyDict(d,"Internal-");
  DE[a][ne[a]++]=nlines;

  md[nlines] = 20;
  Md[nlines] = 40000;
  block[nlines] = BLOCK;
  nlines++;
  
  /* 3. BEGIN+:BEGIN-      First+:Terminal-:Single+:Single-     0:Infinity */
  a = setkeyDict(d,sBEGINFWD);
  UC[a][nc[a]++] = nlines;
  a = setkeyDict(d,sBEGINRVS);
  UC[a][nc[a]++] = nlines;
  
  a = setkeyDict(d,"First+");
  DE[a][ne[a]++]=nlines;
  a = setkeyDict(d,"Terminal-");
  DE[a][ne[a]++]=nlines;
  a = setkeyDict(d,"Single+");
  DE[a][ne[a]++]=nlines;
  a = setkeyDict(d,"Single-");
  DE[a][ne[a]++]=nlines;
  
  md[nlines] = 0;
  Md[nlines] = INFI;
  nlines++;
  
  /* 4. Terminal+:First-:Single+:Single-     END+:END-      0:Infinity*/
  a = setkeyDict(d,"Terminal+");
  UC[a][nc[a]++] = nlines;
  a = setkeyDict(d,"First-");
  UC[a][nc[a]++] = nlines;
  a = setkeyDict(d,"Single+");
  UC[a][nc[a]++] = nlines;
  a = setkeyDict(d,"Single-");
  UC[a][nc[a]++] = nlines;
  
  a = setkeyDict(d,sENDFWD);
  DE[a][ne[a]++]=nlines;
  a = setkeyDict(d,sENDRVS);
  DE[a][ne[a]++]=nlines;
  
  md[nlines] = 0;
  Md[nlines] = INFI;
  nlines++;
 
  return(nlines);
}


