/*************************************************************************
*                                                                        *
*   Module: readparam                                                    *
*                                                                        *
*   Reading of prediction parameters of geneid and Gene Model.           *
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

/*  $Id: readparam.c,v 1.3 2000-08-16 07:45:00 eblanco Exp $  */

#include "geneid.h"

void readLine(FILE *File, char* line)
{
  char* res;
 
  res = fgets(line,MAXLINE,File);
  while((res !=NULL) && ((line[0]=='#')|| (line[0]=='\n')))
    res = fgets(line,MAXLINE,File);
}

void readHeader(FILE *File, char* line)
{
  char* res;

  res = fgets(line,MAXLINE,File);
  while((res !=NULL) && ((line[0]=='#') || (line[0]=='\n')))
      res = fgets(line,MAXLINE,File);
  
  if (res == NULL)
    printError("Parameter file not complete");
}

void ReadProfile(FILE *RootFile, profile* p, char* signal)
{
  int i,j;
  char line[MAXLINE];
  
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%d %d %f %d", 
	      &(p->dimension),
	      &(p->offset), 
	      &(p->cutoff),
	      &(p->order)))!=4)
    printError("Bad format: Profile parameters");  
  
  RequestMemoryProfile(p);
  
  /* Transition probabilities */
  for(i=0; i < p->dimension; i++)
    for(j=0; j < p->dimensionTrans; j++)
      {
	readLine(RootFile,line);
	if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j])))!=1)
	  printError("Bad format. Transition profile values"); 
      }
  
  PrintProfile(p,signal);
}

void ReadIsochore(FILE *RootFile, gparam* gp)
{
  float lscore;       
  int OligoLength_1;
  int i,j,f;
  char line[MAXLINE];
  char mess[MAXSTRING];

  /* 0. read boundaries of isochores */
  readHeader(RootFile,line);
  readLine(RootFile,line);

  if ((sscanf(line,"%d %d\n", 
	      &(gp->leftValue),
	      &(gp->rightValue)))!=2)
    printError("Bad format: Boundaries of isochore");

  sprintf(mess,"Percent Boundaries: %d,%d", 
	 gp->leftValue,
	 gp->rightValue);
  printMess(mess); 

  /* 1. read cutoffs over final score of exons */
  readHeader(RootFile,line);
  readLine(RootFile,line);

  if ((sscanf(line,"%f %f %f %f\n",
	      &(gp->Initial->ExonCutoff),
	      &(gp->Internal->ExonCutoff),
	      &(gp->Terminal->ExonCutoff),
	      &(gp->Single->ExonCutoff)))!=4)
    printError("Bad format: ExonCutoffs");  

  sprintf(mess,"Exon cutoffs: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
	  gp->Initial->ExonCutoff,
	  gp->Internal->ExonCutoff,
	  gp->Terminal->ExonCutoff,
	  gp->Single->ExonCutoff);
  printMess(mess); 

  /* 2. read cutoffs over Oligonucleotids statistic score of exons */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",     
	      &(gp->Initial->OligoCutoff),
	      &(gp->Internal->OligoCutoff),
	      &(gp->Terminal->OligoCutoff),
	      &(gp->Single->OligoCutoff)))!=4)
    printError("Bad format: OligoCutoffs");  
  
  sprintf(mess,"Oligo cutoffs: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
	  gp->Initial->OligoCutoff,
	  gp->Internal->OligoCutoff,
	  gp->Terminal->OligoCutoff,
	  gp->Single->OligoCutoff);
  printMess(mess); 

  /* 3. read weigths of oligos statistic score */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",
	      &(gp->Initial->OligoWeight),
	      &(gp->Internal->OligoWeight),
	      &(gp->Terminal->OligoWeight),
	      &(gp->Single->OligoWeight)))!=4)
        printError("Bad format: OligoWeights");  
  
  sprintf(mess,"Oligo weights: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
	  gp->Initial->OligoWeight,
	  gp->Internal->OligoWeight,
	  gp->Terminal->OligoWeight,
	  gp->Single->OligoWeight);
  printMess(mess); 

  /* 4. read weigths to correct the score of exons after cutoff */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",
	      &(gp->Initial->ExonWeight),
	      &(gp->Internal->ExonWeight),
	      &(gp->Terminal->ExonWeight),
	      &(gp->Single->ExonWeight)))!=4)
    printError("Bad format: ExonWeight");  

  sprintf(mess,"Exon weights: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
	  gp->Initial->ExonWeight,
	  gp->Internal->ExonWeight,
	  gp->Terminal->ExonWeight,
	  gp->Single->ExonWeight);
  printMess(mess);

  /* 5. Read splice site profiles */
  /* (a).start codon profile */
  ReadProfile(RootFile, gp->StartProfile , sSTA);

  /* (b).acceptor site profile */
  ReadProfile(RootFile, gp->AcceptorProfile , sACC);

  /* (c).donor site profile */
  ReadProfile(RootFile, gp->DonorProfile , sDON);

  /* (d).stop codon profile */
  ReadProfile(RootFile, gp->StopProfile , sSTO);
  
  /* 6. read Markov log-likelyhood file */
  readHeader(RootFile,line);
  readLine(RootFile,line); 

  if ((sscanf(line,"%d", &(gp->OligoLength)))!=1)
    printError("Bad format: Oligo Lenght");

  sprintf(mess,"Oligo word length: %d",gp->OligoLength);
  printMess(mess);

  /* (a). Initial probability matrix */
  printMess("Reading Markov Initial probability matrix");

  gp->OligoDim=(int)pow((float)4,(float)gp->OligoLength);

  sprintf(mess,"Used oligo array size: %ld",gp->OligoDim * 3);
  printMess(mess);

  readHeader(RootFile,line);
  for(j = 0; j < gp->OligoDim * 3; j++)
    { 
      readLine(RootFile,line); 
      if ((sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3)
	printError("Bad format: Initial Markov values");

      gp->OligoLogsIni[f][i]=lscore;
    }

  /* (b). Transition probability matrix */
  printMess("Reading Markov transition probability matrix");
  
  OligoLength_1= gp->OligoLength + 1;
  gp->OligoDim_1=(int)pow((float)4,(float)OligoLength_1);

  sprintf(mess,"Used oligo array size: %ld",gp->OligoDim_1 * 3);
  printMess(mess);

  readHeader(RootFile,line);
  for(j = 0; j < gp->OligoDim_1 * 3; j++)
    { 
      readLine(RootFile,line); 
      if ((sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3)
	printError("Bad format: Transition Markov values");

      gp->OligoLogsTran[f][i]=lscore;
    }

  /* 7. read maximum number of donors per acceptor site */
  readHeader(RootFile,line);
  readLine(RootFile,line); 
  
  if ((sscanf(line,"%d", &(gp->MaxDonors)))!=1)
    printError("Bad format: MaxDonors value");   
    
  sprintf(mess,"Maximum donors by acceptor = %d\n", gp->MaxDonors);
  printMess(mess);
}

int readparam (char *name, gparam** isochores) 
{
  FILE *RootFile;
  char *Geneid;
  char ExternalFileName[FILENAMELENGTH];

  /* oligo sequence and oligo score in HexaLogFile */
  int i;
  char line[MAXLINE];
  char mess[MAXSTRING];
  int nIsochores;

  printRes("\n\n\t\t\t** Executing geneid 1.0 2000 geneid@imim.es **\n\n");
  Geneid=getenv("GENEID");

  /* 0. Select parameters filename for reading it */
  /* Filename can be: environment option P, env.var GENEID or default */
  /* Using -P option */
  if (strcmp(PARAMETERFILE,name)) 
      {
	sprintf(mess,"Getting param file using -P option: %s",name);
	if ((RootFile = fopen(name,"rb"))==NULL)
	  printError("The parameter file (-P) can not be opened");
      }
  /* Using GENEID environment var */
  else
    if (Geneid)
      {
	sprintf(mess,"Getting param file using GENEID(env var): %s",Geneid);
	sprintf(ExternalFileName,"%s", Geneid);
	if ((RootFile = fopen(ExternalFileName,"rb"))==NULL) 
	  printError("The parameter file (GENEID env.var) can not be opened ");
      }
  /* Using default parameter file */
    else
      {
	sprintf(mess,"Getting param file default");
	if ((RootFile = fopen(name,"rb"))==NULL)
	  printError("The default parameter file can not be opened");
      }
  printMess(mess);

  /* How many isochores? */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  
  if ((sscanf(line,"%d\n", &(nIsochores)))!=1)
      printError("Bad format: Number of isochores");
    
  sprintf(mess,"Number of isochores: %d", nIsochores);
  printMess(mess);

  if (nIsochores > MAXISOCHORES || !(nIsochores>0))
    printError("Bad format: number of isochores(MAXISOCHORES)");

  /* Reading every one of the isochores */
  for(i=0; i<nIsochores; i++)
    {
      sprintf(mess,"Reading isochore number %d", i+1);
      printMess(mess);
      ReadIsochore(RootFile,isochores[i]);
    }

  /* 8. read the GeneModel */
  readHeader(RootFile,line);
  
  /* Ready to update dictionary of exon types */
  resetDict(isochores[0]->D);
  printMess("Dictionary ready");
  
  printMess("Reading GeneModel");
  isochores[0]->nclass = ReadGeneModel(RootFile,
				       isochores[0]->D,
				       isochores[0]->nc,
				       isochores[0]->ne,
				       isochores[0]->UC,
				       isochores[0]->DE,
				       isochores[0]->md,
				       isochores[0]->Md,
				       isochores[0]->block);
  sprintf(mess,"%d GM-lines read and stored", isochores[0]->nclass);
  printMess(mess);   

  shareGeneModel(isochores, nIsochores);
  return(nIsochores);
}



