/*************************************************************************
*                                                                        *
*   Module: ReadExonsGFF                                                 *
*                                                                        *
*   Reading exons (GFF format) from file                                 *
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

/*  $Id: ReadExonsGFF.c,v 1.1 2003-09-10 14:53:34 gparra Exp $  */

#include "geneid.h"
#define MAXSITESEVIDENCES 3*MAXEVIDENCES

/* Read annotations (exons) to improve or fixed some gene prediction */
/* GFF format: tab "\t" is the field separator and # for comments */
/* Name  Source  Type  Begin  End  Score  Strand  Frame  [group] */
long ReadExonsGFF (char *FileName, packEvidence* pv, dict* d)
{
  /* File handle */
  FILE *file;
  
  /* Final number of exons loaded from file (including copies) */
  long i;
  
  /* Split every input line into several tokens (gff records) */
  char line[MAXLINE];
  char lineCopy[MAXLINE];
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;
  char *line7;
  char *line8;
  char *line9;
  
  /* Control of good sorting property: starting position, increasing */
  long lastAcceptor;
  long currAcceptor;

  /* If frame = '.' then make three copies of current exon (3 frames) */
  int three;

  char saux[MAXTYPE];
  char c;
  int slen;
  char mess[MAXSTRING];
  

  /* 0. Open exons file to read the information */
  if ((file=fopen(FileName, "r"))==NULL)
    printError("The exonsGFF file can not be opened to read");
  
  /* 1. Reset counters */
  i = 0;
  three = 0; 
  lastAcceptor = -INFI;
  
  /* 2. Read while there are exons left in the input file */
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      /* 2.a. Comment or empty line: forget it */
      if (line[0]=='#' || line[0]=='\n')
		printMess("Skipping comment line in evidences file");
      else
		{
		  /* 2.b. Processing exon (annotation) */
		  /* Backup copy of the line to show error messages */
      strcpy(lineCopy,line);
	  
      /* Extracting GFF features */
	  line1 = (char *) strtok(line,"\t");
	  line2 = (char *) strtok(NULL,"\t");
	  line3 = (char *) strtok(NULL,"\t"); 
	  line4 = (char *) strtok(NULL,"\t");
	  line5 = (char *) strtok(NULL,"\t");
	  line6 = (char *) strtok(NULL,"\t");
	  line7 = (char *) strtok(NULL,"\t");
	  line8 = (char *) strtok(NULL,"\t");
	  line9 = (char *) strtok(NULL,"\n");
	  
	  /* There are 8 mandatory columns and the last one is optional */
	  if (line1 == NULL || line2 == NULL || line3 == NULL ||
	      line4 == NULL || line5 == NULL || line6 == NULL ||
	      line7 == NULL || line8 == NULL)
	    {
		  sprintf(mess,"Wrong GFF format in annotations (number of records):\n-->%s\n",lineCopy);
	      printError(mess);
	    }
	  
	  /* 1/2. Sequence and Source records not used anymore (redundant) */
	  
	  /* 3. Exon feature: Single, First, Internal, Terminal, ... */
	  if (sscanf(line3,"%s",(pv->vExons+i)->Type) != 1)
	    {
		  sprintf(mess,"Wrong GFF format in annotations (feature):\n-->%s\n",lineCopy);
	      printError(mess);
	    }
	  
	  /* 4. Starting position */
	  if (sscanf(line4,"%ld",&((pv->vSites + pv->nvSites)->Position)) != 1)
	    {
		  sprintf(mess,"Wrong GFF format in annotations (starting position):\n-->%s\n",lineCopy);
		  printError(mess);
	    }
	  
	  /* 5. Finishing position */
	  if (sscanf(line5,"%ld",&((pv->vSites + pv->nvSites + 1)->Position)) != 1)
	    {
		  sprintf(mess,"Wrong GFF format in annotations (finishing position):\n-->%s\n",lineCopy);
	      printError(mess);
	    }
	  
	  /* 6. Score = float value or '.'(infinitum) */
	  if (sscanf(line6,"%lf",&((pv->vExons+i)->Score)) != 1)
	    {
		  if ((sscanf(line6,"%c",&c) != 1) || (c != '.'))
			{
			  sprintf(mess,"Wrong GFF format in annotations (score):\n-->%s\n",lineCopy);
			  printError(mess);
			}
	      (pv->vExons+i)->Score = MAXSCORE;
	    }
	  
	  /* 7. Strand (reading sense) [+|-] */
	  if ((sscanf(line7,"%c",&((pv->vExons+i)->Strand))!= 1) ||
		  (((pv->vExons+i)->Strand != '+') && ((pv->vExons+i)->Strand != '-')))
	    {       
		  sprintf(mess,"Wrong GFF format in annotations (strand):\n-->%s\n",lineCopy);
	      printError(mess);
	    }
	  
	  /* 8. Frame = integer or '.' */
	  if (sscanf(line8,"%hd",&((pv->vExons+i)->Frame)) != 1)
	    {
		  /* Is it a dot? */
		  if ((sscanf(line8,"%c",&c)!= 1) || (c!='.'))
			{
			  sprintf(mess,"Wrong GFF format in annotations (frame):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  /* make three copies */
		  three = 1;
	    }
	  else
		{
		  /* Checking input frame between 0..2 */
		  if (((pv->vExons+i)->Frame < 0) || ((pv->vExons+i)->Frame > 2))
			{
			  sprintf(mess,"Wrong GFF value in annotations (frame not between 0..2):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		}
	  
	  /* 9. Group: optional, string */
	  if (line9 != NULL)
	    {
		  if (sscanf(line9,"%s",(pv->vExons+i)->Group) != 1)
			{
			  sprintf(mess,"Wrong GFF value in annotations (group):\n-->%s\n",lineCopy);
			  printError(mess);
			}
	    }
	  else
        {
		  /* This exon will be allowed to join to ab initio predictions */
		  strcpy((pv->vExons+i)->Group, NOGROUP);
        }
	  
	  /* 2.c. Process current exon */
	  /* (A). Checking exon feature (gene model) */
	  saux[0]='\0';
	  strcpy (saux, (pv->vExons+i)->Type);
	  slen = strlen(saux);
	  saux[slen++] = (pv->vExons+i)->Strand;
	  saux[slen] = '\0';
	  
	  if (getkeyDict(d,saux) == NOTFOUND)
	    {
         /* Forget it: exon type is not in current gene model */
		  sprintf(mess,"Wrong GFF feature in annotations (unknown):\n-->%s\n",lineCopy);
	      printMess(mess); 
	    }
	  else
	    {
		  /* (B). Well-sorted (by start position) list of read exons */
	      currAcceptor = (pv->vSites + pv->nvSites)->Position;
	      if (lastAcceptor > currAcceptor)
			{
			  sprintf(mess,"Wrong order in annotations (starting position %ld):\n-->%s\n",
					  lastAcceptor,
					  lineCopy);  
			  printError(mess);  
			}
	      else
			{
			  lastAcceptor = (pv->vSites + pv->nvSites)->Position;
			  
			  /* (C). Setting dummy sites to this exon */
			  (pv->vExons+i)->Acceptor = (pv->vSites + pv->nvSites);
			  (pv->vExons+i)->Donor = (pv->vSites + pv->nvSites + 1); 
			  
			  /* Updating information about sites to the range 0..L-1 */
			  (pv->vExons+i)->offset1 = -COFFSET;
			  (pv->vExons+i)->offset2 = -COFFSET;
			  
			  /* (D). Making three (two more) copies if needed */
			  if (three)
				{
				  /* Creating three exons (3 frames): sharing sites */
				  (pv->vExons+i+1)->Acceptor = (pv->vSites + pv->nvSites);
				  (pv->vExons+i+1)->Donor = (pv->vSites + pv->nvSites + 1); 
				  (pv->vExons+i+2)->Acceptor = (pv->vSites + pv->nvSites);
				  (pv->vExons+i+2)->Donor = (pv->vSites + pv->nvSites + 1); 
				  
				  /* Updating information about sites to the range 0..L-1 */
				  (pv->vExons+i+1)->offset1 = -COFFSET;
				  (pv->vExons+i+1)->offset2 = -COFFSET;
				  (pv->vExons+i+2)->offset1 = -COFFSET;
				  (pv->vExons+i+2)->offset2 = -COFFSET;
				  
				  /* Setting frame values */
				  (pv->vExons+i)->Frame = 0;
				  (pv->vExons+i+1)->Frame = 1;
				  (pv->vExons+i+2)->Frame = 2;
				  
				  /* Copy some exon attributes */
				  strcpy((pv->vExons+i+1)->Type,(pv->vExons+i)->Type);
				  strcpy((pv->vExons+i+2)->Type,(pv->vExons+i)->Type);
				  (pv->vExons+i+1)->Score = (pv->vExons+i)->Score;
				  (pv->vExons+i+2)->Score = (pv->vExons+i)->Score;
				  (pv->vExons+i+1)->Strand = (pv->vExons+i)->Strand;
				  (pv->vExons+i+2)->Strand = (pv->vExons+i)->Strand;
		      
				  /* Computing remainder from frame value for copies 1,2 */
				  (pv->vExons+i+1)->Remainder = 
					((3 - (((pv->vExons+i+1)->Donor->Position - 
							(pv->vExons+i+1)->Acceptor->Position - 
							(pv->vExons+i+1)->Frame + 1)%3)) %3);
				  
				  (pv->vExons+i+2)->Remainder = 
					((3 - (((pv->vExons+i+2)->Donor->Position - 
							(pv->vExons+i+2)->Acceptor->Position - 
							(pv->vExons+i+2)->Frame + 1)%3)) %3);
				  
				  /* The same group */
				  strcpy((pv->vExons+i+1)->Group,(pv->vExons+i)->Group);
				  strcpy((pv->vExons+i+2)->Group,(pv->vExons+i)->Group);
				  
				  /* Evidence flag activated */
				  (pv->vExons+i+1)->evidence = 1;
				  (pv->vExons+i+2)->evidence = 1;
				  
				} /* End of if(three) */      
			  
			  /* (E). Doing the same for the original exon: */
			  /* Computing remainder from frame value */
			  (pv->vExons+i)->Remainder = 
				((3 - (((pv->vExons+i)->Donor->Position - 
						(pv->vExons+i)->Acceptor->Position - 
						(pv->vExons+i)->Frame + 1)%3)) %3);
			  
			  /* Evidence flag activated */
			  (pv->vExons+i)->evidence = 1;
			  
			  /* Updating and checking loop counters */
			  i = (three)? i+3 : i+1;
			  pv->nvSites = pv->nvSites + 2;
			  three = 0;
			  
			  if ((i + FRAMES) > MAXEVIDENCES)
				printError("Too many evidences: increase MAXEVIDENCES definition");
			  
			  if ((pv->nvSites + (2*FRAMES)) > MAXSITESEVIDENCES)
				printError("Too many site evidences: increase MAXEVIDENCES definition");
			  
			} /* End of sorting checkpoint */
		} /* End of feature_in_gene_model checkpoint */
		} /* End of checkpoint for comments (#) */
    } /* End of while (input exons) */
  fclose(file);
  
  /* Return the number of created exons (including replications) */
  return(i);
}

