/*************************************************************************
*                                                                        *
*   Module: ReadExonsGFF                                                 *
*                                                                        *
*   GenAmic alone: It reads exons in GFF format.                         *
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

/*  $Id: ReadExonsGFF.c,v 1.3 2000-08-08 14:19:47 eblanco Exp $  */

#include "geneid.h"

long ReadExonsGFF (char *FileName, packEvidence* pv, dict* d)
{
  long i;
  FILE *file;
  char line[MAXLINE];
  char saux[MAXTYPE+1];
  short fraux;
  char c;
  int three;
  long lastAcceptor;
  long currAcceptor;
  char mess[MAXSTRING];
 
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;
  char *line7;
  char *line8;
  char *line9;

  if ((file=fopen(FileName, "r"))==NULL)
    printError("The exonsGFF file cannot be open for read");
  
  /* Coments: line begins with # */
  /* gff format = "Name  Source  Type  Begin  End  Score  Strand  Frame [group] */
  i = 0; 
  three = 0; 
  lastAcceptor = -INFI;

  pv->i1vExons = 0;
  pv->i2vExons = 0;
  pv->nvSites = 0;

  while(fgets(line,MAXLINE,file)!=NULL)
    {
      if(line[0]=='#' || line[0]=='\n')
	{
	  /* Skip this line */
	  printMess("Skipping comment line");
	}
      else
	{
	  /* For each line extract the features (GFF format) */
	  /* Split line in four parts: UC DE Dist block */
          line1 = (char *) strtok(line,"\t");
          line2 = (char *) strtok(NULL,"\t");
          line3 = (char *) strtok(NULL,"\t"); 
          line4 = (char *) strtok(NULL,"\t");
	  line5 = (char *) strtok(NULL,"\t");
	  line6 = (char *) strtok(NULL,"\t");
          line7 = (char *) strtok(NULL,"\t");
          line8 = (char *) strtok(NULL,"\t");
          line9 = (char *) strtok(NULL,"\n");

	  if (line1 == NULL || line2 == NULL || line3 == NULL ||
	      line4 == NULL || line5 == NULL || line6 == NULL ||
	      line7 == NULL || line8 == NULL)
	    {
	      sprintf(mess, "Bad format: Exon GFF %ld\n",i);
	      printError(mess);
	    }

	  /* 1/2. Sequence and Source not used */
	  
	  /* 3. Exon Type */
	  if (sscanf(line3,"%s",(pv->vExons+i)->Type) != 1)
	    {
	      sprintf(mess, "Bad format Type: Exon %ld\n",i);
	      printError(mess);
	    }
	  
	  /* 4. Left position */
	  if (sscanf(line4,"%ld",&((pv->vSites + pv->nvSites)->Position)) != 1)
	    {
	      sprintf(mess, "Bad format Acceptor: Exon %ld\n",i);
	      printError(mess);
	    }

	  /* 5. Right position */
	  if (sscanf(line5,"%ld",&((pv->vSites + pv->nvSites + 1)->Position)) != 1)
	    {
	      sprintf(mess, "Bad format Donor: Exon %ld\n",i);
	      printError(mess);
	    }

	  /* 6. Score = '.' or float */
	  if (sscanf(line6,"%lf",&((pv->vExons+i)->Score)) != 1)
	    {
	      if ((sscanf(line6,"%c",&c)!= 1) || (c!='.'))
		{
		  sprintf(mess, "Bad format Score: Exon %ld\n",i);
		  printError(mess);
		}
	      (pv->vExons+i)->Score = MAXSCORE;
	    }

	  /* 7. Strand */
	  if (sscanf(line7,"%c",&((pv->vExons+i)->Strand))!= 1)
	    {
	      sprintf(mess, "Bad format Strand: Exon %ld\n",i);
	      printError(mess);
	    }
	  
	  /* 8. Frame = '.' or integer */
	  if (sscanf(line8,"%hd",&((pv->vExons+i)->Frame)) != 1)
	    {
	      if ((sscanf(line8,"%c",&c)!= 1) || (c!='.'))
		{
		  sprintf(mess, "Bad format Frame: Exon %ld\n",i);
		  printError(mess);
		}
	      three = 1; 
	    }	  

	  /* 9. Group */
	  if (line9 != NULL)
	    {
	      if (sscanf(line9,"%d",&((pv->vExons+i)->Group)) != 1)
		{
		  sprintf(mess, "Bad format Group: Exon %ld\n",i);
		  printError(mess);
		}
	    }
	  else
	      (pv->vExons+i)->Group = NOGROUP; 

	  /* Process features from current exon */

	  /* What is the type of this exon? */
	  saux[0]='\0';
	  strcpy (saux, (pv->vExons+i)->Type);
	  strcat (saux, &((pv->vExons+i)->Strand));

	  if (getkeyDict(d,saux) == NOTFOUND)
	    {
	      /* Type of this exon does not appear in GeneModel */
	      /* Skip this exon */
	      sprintf(mess, "Skipping exon %ld(unknown type)\n",i);
	      printMess(mess); 
	    }
	  else
	    {
	      currAcceptor = (pv->vSites + pv->nvSites)->Position;
	      /* File must be ordered by acceptor positions */
	      if (lastAcceptor > currAcceptor)
		{
		  sprintf(mess,"Exons file bad ordered by acceptor position(line %ld)",i);
		  printError(mess);  
		}
	      else
		{
		  lastAcceptor = (pv->vSites + pv->nvSites)->Position;
		  
		  /* Assign fool sites to this exon */
		  (pv->vExons+i)->Acceptor = (pv->vSites + pv->nvSites);
		  (pv->vExons+i)->Donor = (pv->vSites + pv->nvSites + 1); 
	      
		  if (three)
		    {
		      /* Creating three exons(3 frames) */
		      (pv->vExons+i+1)->Acceptor = (pv->vSites + pv->nvSites);
		      (pv->vExons+i+1)->Donor = (pv->vSites + pv->nvSites + 1); 
		      (pv->vExons+i+2)->Acceptor = (pv->vSites + pv->nvSites);
		      (pv->vExons+i+2)->Donor = (pv->vSites + pv->nvSites + 1); 
		      
		      (pv->vExons+i)->Frame = 0;
		      (pv->vExons+i+1)->Frame = 1;
		      (pv->vExons+i+2)->Frame = 2;
		      
		      strcpy((pv->vExons+i+1)->Type,(pv->vExons+i)->Type);
		      strcpy((pv->vExons+i+2)->Type,(pv->vExons+i)->Type);
		      (pv->vExons+i+1)->Score = (pv->vExons+i)->Score;
		      (pv->vExons+i+2)->Score = (pv->vExons+i)->Score;
		      (pv->vExons+i+1)->Strand = (pv->vExons+i)->Strand;
		      (pv->vExons+i+2)->Strand = (pv->vExons+i)->Strand;
		      
		      /* Remainder is ... */
		      (pv->vExons+i+1)->Remainder = 
			((3 - (((pv->vExons+i+1)->Donor->Position - 
				(pv->vExons+i+1)->Acceptor->Position - 
				(pv->vExons+i+1)->Frame + 1)%3)) %3);
		      
		      /* Remainder is ... */
		      (pv->vExons+i+2)->Remainder = 
			((3 - (((pv->vExons+i+2)->Donor->Position - 
				(pv->vExons+i+2)->Acceptor->Position - 
				(pv->vExons+i+2)->Frame + 1)%3)) %3);
		      
		      /* The same group */
		      (pv->vExons+i+1)->Group = (pv->vExons+i)->Group;
		      (pv->vExons+i+2)->Group = (pv->vExons+i)->Group;  
		      
		      /* Evidence exons */
		      (pv->vExons+i+1)->evidence = 1;
		      (pv->vExons+i+2)->evidence = 1;
		      
		      /* If strand (-), frame and remainder must be exchanged */
		      if ((pv->vExons+i)->Strand =='-')  
			{ 
			  fraux = (pv->vExons+i+1)->Frame; 
			  (pv->vExons+i+1)->Frame = (pv->vExons+i+1)->Remainder;
			  (pv->vExons+i+1)->Remainder = fraux;
			  fraux = (pv->vExons+i+2)->Frame; 
			  (pv->vExons+i+2)->Frame = (pv->vExons+i+2)->Remainder;
			  (pv->vExons+i+2)->Remainder = fraux;
			}
		    }      
		  
		  /* Remainder is ... */
		  (pv->vExons+i)->Remainder = 
		    ((3 - (((pv->vExons+i)->Donor->Position - 
			    (pv->vExons+i)->Acceptor->Position - 
			    (pv->vExons+i)->Frame + 1)%3)) %3);
		  
		  /* Evidence exon */
		  (pv->vExons+i)->evidence = 1;
		  
		  /* If strand (-), frame and remainder must be exchanged */
		  if ((pv->vExons+i)->Strand =='-')  
		    { 
		      fraux = (pv->vExons+i)->Frame; 
		      (pv->vExons+i)->Frame = (pv->vExons+i)->Remainder;
		      (pv->vExons+i)->Remainder = fraux;
		    }
		  
		  i = (three)? i+3 : i+1;
		  three = 0;
		  
		  pv->nvSites = pv->nvSites + 2;
	    
		  if (i > NUMEEVIDENCES)
		    printError("Too many evidence exons: Change NUMEEVIDENCES parameter");
		  
		  if (pv->nvSites > NUMSEVIDENCES)
		    printError("Too many evidence sites: Change NUMSEVIDENCES parameter");
		}
	    }
	}
    }
  fclose(file);
  return(i);
}
