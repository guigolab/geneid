/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: ReadHSP                                            *
*                                                              *
*   Reading of blast hsps                                      *
***************************************************************/   

#include "blast2gff.h"

extern int GFFIN;

hsp* allocateNewHSP()
{
  hsp* p;
  
  /* New HSP */
  if ((p = (hsp *) malloc(sizeof(hsp))) == NULL)
    printError("Not enough space to hold one new HSP");
  
  return(p);
}

long ReadHSP_RAW(packHSP* h,char *FileName)
{
  int i;
  FILE *file;
  char line[MAXLINE];
  char mess[MAXSTRING];
  long pos1, pos2;
  int score;
  char strand;
  short frame;
  float bits;
  double prob, exp;
  int id,sim;
  long pos3,pos4;
  int indexFrame;
  long L;
  char Subject[LOCUSLENGTH];

  if ((file=fopen(FileName, "r"))==NULL)
    printError("The HSP file file cannot be open for read");
  
  i = 0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      if(line[0]=='#')
	{
	  /* Skip this line */
	}
      else
	{
	  if ((sscanf(line, "%*s %*s %ld %ld %c %hd %d %f %lf %lf %d %d %s %ld %ld %*c %*d %ld",
		      &pos1,
		      &pos2,
		      &strand,
		      &frame,
		      &score,
		      &bits,
		      &prob,
		      &exp,
		      &id,
		      &sim,
		      Subject,
		      &pos3,
		      &pos4,
		      &L) != 14) || (frame<1) || (frame>3))
	    {
	      sprintf(mess, "Error reading HPS: line %d\n",i);
	      printError(mess);
	    }
	  else
	    {
	      /* Allocating HSPs in packHSP */
	      if (strand == '+')
		indexFrame = (frame % 3);
	      else
		indexFrame = (frame % 3) + FRAMES;

	      /* New item HSP */
	      h->hsps[indexFrame][h->nHsps[indexFrame]] = allocateNewHSP();

	      /* Setting attributes */
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StartQ = pos1;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->EndQ = pos2;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StrandQ = strand;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->FrameQ = frame;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Score = score;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Bits = bits;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Probability = prob;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Expect = exp;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Similarity = sim;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Identity = id;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StartS = pos3;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->EndS = pos4;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Lenght = L;
	      strcpy(h->hsps[indexFrame][h->nHsps[indexFrame]]->Subject,Subject);
	      h->nHsps[indexFrame]++;
	      h->nTotalHsps++;
	    }
	}
      i++;
    }
  fclose(file);
  
  return(h->nTotalHsps);
}

long ReadHSP_GFF (packHSP* h,char *FileName, char Query[LOCUSLENGTH])
{
  int i;
  FILE *file;
  char line[MAXLINE];
  char mess[MAXSTRING];
  long pos1, pos2;
  char strand;
  short frame;
  float score;
  int indexFrame;
  char Subject[LOCUSLENGTH];


  if ((file=fopen(FileName, "r"))==NULL)
    printError("The HSP file file cannot be open for read");
  
  i = 0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      if(line[0]=='#')
	{
	  /* Skip this line */
	}
      else
	{
	  if ((sscanf(line, "%s %*s %*s %ld %ld %f %c %hd %s",
		      Query,
		      &pos1,
		      &pos2,
		      &score,
		      &strand,
		      &frame,
		      Subject) != 7) || (frame<1) || (frame>3))
	    {
	      sprintf(mess, "Error reading HPS: line %d\n",i);
	      printError(mess);
	    }
	  else
	    {
	      /* Allocating HSPs in packHSP */
	      if (strand == '+')
		indexFrame = (frame % 3);
	      else
		indexFrame = (frame % 3) + FRAMES;

	      /* New item HSP */
	      h->hsps[indexFrame][h->nHsps[indexFrame]] = allocateNewHSP();
	      
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StartQ = pos1;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->EndQ = pos2;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StrandQ = strand;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->FrameQ = frame;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Score = score;
	      strcpy(h->hsps[indexFrame][h->nHsps[indexFrame]]->Subject,Subject);

	      h->nHsps[indexFrame]++;
	      h->nTotalHsps++;
	    }
	}
      i++;
    }
  fclose(file);
  
  return(h->nTotalHsps);
}

void ReadHSP (packHSP* allHsp,char* HSPFile, char Query[LOCUSLENGTH])
{
  if (GFFIN)
    ReadHSP_GFF(allHsp,HSPFile,Query); 
  else
    {
      ReadHSP_RAW(allHsp,HSPFile); 
      strcpy(Query,NONAME);
    }
}
