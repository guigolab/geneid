/*************************************************************************
*                                                                        *
*   Module: readparam                                                    *
*                                                                        *
*   Reading statistic parameters and gene construction model             *
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

/*  $Id: readparam.c,v 1.1 2003-09-10 14:53:34 gparra Exp $  */

#include "geneid.h"

/* Numeric values: read one line skipping comments and empty lines */
void readLine(FILE *File, char* line)
{
  char* res;
 
  res = fgets(line,MAXLINE,File);
  while((res !=NULL) && ((line[0]=='#')|| (line[0]=='\n')))
    res = fgets(line,MAXLINE,File);

  if (res == NULL)
    printError("Parameter file: unexpected end of reading");
}

/* Headers of values: read one line skipping comments and empty lines */
void readHeader(FILE* File, char* line)
{
  char* res;

  res = fgets(line,MAXLINE,File);
  while((res !=NULL) && ((line[0]=='#') || (line[0]=='\n')))
	res = fgets(line,MAXLINE,File);
  
  if (res == NULL)
    printError("Parameter file: unexpected end of reading");
}

void SetProfile(profile* p, FILE* RootFile, char* signal)
{
  int i,j,x,y;
  char line[MAXLINE];
  char mess[MAXSTRING];

  /* According to the order of Markov chain, select a different method */
  /* Position weight array: transition probabilities in every position */
  switch(p->order)
	{
	case 0:
	  /* 5 combinations / pos */
	  for(i=0; i < p->dimension; i++)  
		{
		  j = 0;

		  /* Reading A,C,G and T: creating N transition */
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+1])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+2])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+3])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  /* Generating the N value */
		  p->transitionValues[i][j+4] = 
			(p->transitionValues[i][j] +
			p->transitionValues[i][j+1] +
			p->transitionValues[i][j+2] +
			p->transitionValues[i][j+3]) / 4;
		}
	  break;
	case 1:
	  /* 25 combinations / pos */
	  for(i=0; i < p->dimension; i++)  
		{
		  /* Reading AX,CX,GX and TX: creating XN */
		  for(j=0; j < p->dimensionTrans-5; j=j+5)
			{
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+1])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+2])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+3])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  /* Generating the XN value every iteration: AN,CN,GN,TN */
			  p->transitionValues[i][j+4] = 
				(p->transitionValues[i][j] +
				 p->transitionValues[i][j+1] +
				 p->transitionValues[i][j+2] +
				 p->transitionValues[i][j+3]) / 4;
			}
		  
		  /* Creating NA,NC,NG,NT transition values */
		  for(x=0; x < 4; x++)
			{
			   p->transitionValues[i][j+x] =  
				 (p->transitionValues[i][x] +
				  p->transitionValues[i][5+x] +
				  p->transitionValues[i][10+x] +
				  p->transitionValues[i][15+x]) / 4;
			}

		  /* Creating the value NN */
		  p->transitionValues[i][j+4] =
			(p->transitionValues[i][j] +
			 p->transitionValues[i][j+1] +
			 p->transitionValues[i][j+2] +
			 p->transitionValues[i][j+3]) / 4;
		}
	  break;
	case 2:
	  /* 125 combinations / pos */
	  for(i=0; i < p->dimension; i++) 
		{
		  /* Reading AXX,CXX,GXX,TXX and creating ANX,CNX,GNX,TNX */
		  for(j=0; j < p->dimensionTrans-25; j=j+25)
			{
			  /* 20 = 5^order - 5 */
			  for(y=0; y < 20; y=y+5)
				{
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+1+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+2+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+3+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  /* Generating the XXN value (5th value): XAN,XCN,XGN,XTN */
				  p->transitionValues[i][j+4+y] = 
					(p->transitionValues[i][j+y] +
					 p->transitionValues[i][j+1+y] +
					 p->transitionValues[i][j+2+y] +
					 p->transitionValues[i][j+3+y]) / 4;
				}
			  
			  /* Creating XNA,XNC,XNG,XNT */
			  for(x=0; x < 4; x++)
				{
				  p->transitionValues[i][j+y+x] =  
					(p->transitionValues[i][j+x] +
					 p->transitionValues[i][j+5+x] +
					 p->transitionValues[i][j+10+x] +
					 p->transitionValues[i][j+15+x]) / 4;
				}
			  
			  /* Creating the value XNN */
			  p->transitionValues[i][j+4+y] =
				(p->transitionValues[i][j+y] +
				 p->transitionValues[i][j+1+y] +
				 p->transitionValues[i][j+2+y] +
				 p->transitionValues[i][j+3+y]) / 4;
			}
		  
		  /* Creating NAX,NCX,NGX,NTX (j=100)*/
		  for(y=0; y < 20; y=y+5)
			{
			  /* Creating NAX... */
			  for(x=0; x < 4; x++)
				{
				  p->transitionValues[i][j+y+x] =  
					(p->transitionValues[i][y+x] +
					 p->transitionValues[i][25+y+x] +
					 p->transitionValues[i][50+y+x] +
					 p->transitionValues[i][75+y+x]) / 4;
				}
			  /* Computing NAN (x=4)*/
			  p->transitionValues[i][j+y+x] = 
				(p->transitionValues[i][j+y] + 
				 p->transitionValues[i][j+y+1] + 
				 p->transitionValues[i][j+y+2] +
				 p->transitionValues[i][j+y+3]) / 4;
			}
		  /* Creating NNX (j=100,y=20)*/
		  for(x=0; x < 4; x++)
			{
			  p->transitionValues[i][j+y+x] = 
				(p->transitionValues[i][j+x] +
				 p->transitionValues[i][j+5+x] +
				 p->transitionValues[i][j+10+x] +
				 p->transitionValues[i][j+15+x]) /4;
			}
		  /* Finally, creating NNN (j=100,y=20,x=4) */
		  p->transitionValues[i][j+y+x] = 
			(p->transitionValues[i][j+y] +
			 p->transitionValues[i][j+y+1] +
			 p->transitionValues[i][j+y+2] +
			 p->transitionValues[i][j+y+3]) /4;
		}
	  break;
	default:
	  printError("Sorry, Markov order higher than 2 not implemented yet");
	  break;
	}
}

/* Read information useful to predict signals */
void ReadProfile(FILE* RootFile, profile* p, char* signal)
{
  char line[MAXLINE];
  char mess[MAXSTRING];
  
  /* Definition parameters: Length, offset, cutoff and order (Markov chain) */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%d %d %f %d", 
			  &(p->dimension),
			  &(p->offset), 
			  &(p->cutoff),
			  &(p->order)))!=4)
	{
	  sprintf(mess,"Wrong format: Definition parameters in %s profile",signal);
	  printError(mess);
	}
  
  /* Memory to allocate the data with these fixed dimensions */
  RequestMemoryProfile(p);

  /* Useful to check everything is OK */
  PrintProfile(p,signal);

  /* Prepare and read values */
  SetProfile(p,RootFile,signal);
}

/* Read information about signal and exon prediction in one isochore */
/* - isochores are specific DNA regions according to the G+C content - */
void ReadIsochore(FILE* RootFile, gparam* gp)
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
    printError("Wrong format: isochore boundaries (G+C percent)");

  sprintf(mess,"Isochores boundaries(min/max percentage): %d,%d", 
		  gp->leftValue,
		  gp->rightValue);
  printMess(mess); 

  /* 1. read cutoff (final score) to accept one predicted exon */
  readHeader(RootFile,line);
  readLine(RootFile,line);

  if ((sscanf(line,"%f %f %f %f\n",
			  &(gp->Initial->ExonCutoff),
			  &(gp->Internal->ExonCutoff),
			  &(gp->Terminal->ExonCutoff),
			  &(gp->Single->ExonCutoff)))!=4)
    printError("Wrong format: exon score cutoffs (number/type)");  

  sprintf(mess,"Exon cutoffs: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
		  gp->Initial->ExonCutoff,
		  gp->Internal->ExonCutoff,
		  gp->Terminal->ExonCutoff,
		  gp->Single->ExonCutoff);
  printMess(mess); 

  /* 2. read cutoff (potential coding score) to accept one predicted exon */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",     
			  &(gp->Initial->OligoCutoff),
			  &(gp->Internal->OligoCutoff),
			  &(gp->Terminal->OligoCutoff),
			  &(gp->Single->OligoCutoff)))!=4)
    printError("Wrong format: potential coding score cutoffs (number/type)");  
  
  sprintf(mess,"Oligo cutoffs: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
		  gp->Initial->OligoCutoff,
		  gp->Internal->OligoCutoff,
		  gp->Terminal->OligoCutoff,
		  gp->Single->OligoCutoff);
  printMess(mess); 

  /* 3. Weight of coding potential score in exon score */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",
			  &(gp->Initial->OligoWeight),
			  &(gp->Internal->OligoWeight),
			  &(gp->Terminal->OligoWeight),
			  &(gp->Single->OligoWeight)))!=4)
	printError("Wrong format: weight of coding potential score (number/type)");  
  
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
    printError("Wrong format: exon weight values (number/type)");  

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
  
  /* 6. read coding potential log-likelihood values (Markov chains) */
  readHeader(RootFile,line);
  readLine(RootFile,line); 

  if ((sscanf(line,"%d", &(gp->OligoLength)))!=1)
    printError("Wrong format: oligonucleotide length");

  sprintf(mess,"Oligonucleotide (word) length: %d",gp->OligoLength);
  printMess(mess);

  /* (a). Initial probability matrix */
  printMess("Reading Markov Initial likelihood matrix");

  /* Computing the right number of initial values to read */      
  gp->OligoDim=(int)pow((float)4,(float)gp->OligoLength);

  sprintf(mess,"Used oligo array size: %ld",gp->OligoDim * 3);
  printMess(mess);

  readHeader(RootFile,line);
  for(j = 0; j < gp->OligoDim * 3; j++)
    { 
      readLine(RootFile,line); 
      if ((sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3)
        {
          sprintf(mess, "Wrong format/nunber (%s): Initial Markov value", line);
          printError(mess);
        }
      gp->OligoLogsIni[f][i]=lscore;
    }

  /* (b). Transition probability matrix */
  printMess("Reading Markov Transition likelihood matrix");

  OligoLength_1= gp->OligoLength + 1;
  gp->OligoDim_1=(int)pow((float)4,(float)OligoLength_1);

  sprintf(mess,"Used oligo array size: %ld",gp->OligoDim_1 * 3);
  printMess(mess);

  readHeader(RootFile,line);
  for(j = 0; j < gp->OligoDim_1 * 3; j++)
    { 
      readLine(RootFile,line); 
      if ((sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3)
        {
          sprintf(mess, "Wrong format/number (%s): Transition Markov value", line);
          printError(mess);
        }
      gp->OligoLogsTran[f][i]=lscore;
    }

  /* 7. read maximum number of donors per acceptor site (BuildExons) */
  readHeader(RootFile,line);
  readLine(RootFile,line); 
  
  if ((sscanf(line,"%d", &(gp->MaxDonors)))!=1)
    printError("Bad format: MaxDonors value");   
    
  sprintf(mess,"Maximum donors by acceptor = %d\n", gp->MaxDonors);
  printMess(mess);
}

/* Read the input of statistics data model */
int readparam (char* name, gparam** isochores)
{
  FILE* RootFile;
  char* Geneid;
  char ExternalFileName[FILENAMELENGTH];

  int i;
  char line[MAXLINE];
  char mess[MAXSTRING];
  int nIsochores;

  /* 1. Select parameters filename for reading it */
  /* Filename must be: option P, env.var GENEID or default (none) */
  Geneid=getenv("GENEID");

  /* (a) Using -P option */
  if (strcmp(PARAMETERFILE,name))
	{
	  sprintf(mess,"Loading parameter file by using -P option: %s",name);
	  if ((RootFile = fopen(name,"rb"))==NULL)
		printError("Parameter file (-P) can not be open to read");
	}
  /* (b) Using GENEID environment var */
  else
    if (Geneid)
      {
        sprintf(mess,"Loading parameter file by using GENEID (env. var): %s",Geneid);
		sprintf(ExternalFileName,"%s", Geneid);
		if ((RootFile = fopen(ExternalFileName,"rb"))==NULL) 
          printError("Parameter file (GENEID env.var) can not be open to read");
      }
  /* (c) Using default parameter file */
    else
      {
        sprintf(mess,"Loading parameter file default");
		if ((RootFile = fopen(name,"rb"))==NULL)
          printError("Parameter file (default) can not be open to read");
      }

  /* rootfile will be the parameter file handle descriptor */
  printMess(mess);

  /* 2. Read the number of isochores */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  
  if ((sscanf(line,"%d\n", &(nIsochores)))!=1)
	printError("Wrong format: Number of isochores");
    
  sprintf(mess,"Number of isochores: %d", nIsochores);
  printMess(mess);

  if (nIsochores > MAXISOCHORES || !(nIsochores>0))
    printError("Wrong value: number of isochores(MAXISOCHORES)");

  /* 3. Reading every one of the isochores */
  for(i=0; i<nIsochores; i++)
    {
      sprintf(mess,"Reading isochore %d", i+1);
      printMess(mess);
      ReadIsochore(RootFile,isochores[i]);
    }

  /* 4. Reading the GeneModel */
  readHeader(RootFile,line);
  
  /* Ready to update dictionary of exon types */
  resetDict(isochores[0]->D);
  printMess("Dictionary ready to acquire information");
  
  printMess("Reading Gene Model rules");
  isochores[0]->nclass = ReadGeneModel(RootFile,
									   isochores[0]->D,
									   isochores[0]->nc,
									   isochores[0]->ne,
									   isochores[0]->UC,
									   isochores[0]->DE,
									   isochores[0]->md,
									   isochores[0]->Md,
									   isochores[0]->block);

  sprintf(mess,"%d Gene Model rules have been read and saved",
		  isochores[0]->nclass);
  printMess(mess);   

  /* Replication of gene model information for each isochore */
  shareGeneModel(isochores, nIsochores);

  return(nIsochores);
}



