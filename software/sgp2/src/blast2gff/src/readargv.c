/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: readargv                                           *
*                                                              *
*   Get setup options and filenames of user input              *
***************************************************************/ 

#include "blast2gff.h"

extern int VRB, GFFIN, GFFOUT, PSR;
                
extern char* optarg;
extern int optind;

char *USAGE="Incorrect usage:\nNAME\n\tblast2gff - A program to build similarity regions from blastx outputs\nSYNOPSIS\n\tblast2gff <HSP_file>\n\n";

void printHelp()
{
  printf("Short Manual for blast2gff:\n");
  printf("------------------------\n\n");
  printf("Setup Options:\n\n");
 
  printf("\t-g: Input HSP in GFF format\n");
  printf("\t-G: Print HSP in GFF format\n");
  printf("\t-o: Not SR construction\n");
  printf("\t-v: Verbose\n");
  printf("\t-h: Show this Short Manual\n");
}

void readargv (int argc,char* argv[],
	       char* HSPFile) 
{
  int c;
  int error=0;
  char mess[MAXSTRING];

  /* Reading setup options */
  while ((c = getopt(argc,argv,"Ggvoh")) != -1)
    switch(c)
      { 
      case 'g': GFFIN++; 
	break;
      case 'G': GFFOUT++; 
	break;
      case 'o': PSR--; 
	break;
      case 'v': VRB++; 
	break;
      case 'h': printHelp();
	exit(1);
	break;
      }

  if (error)
    printError(USAGE);

  /* Setup Errors: Wrong number of filenames */
  /* Get the name of the file *.hsp */
  if (optind < argc)
    {
      strcpy(HSPFile,argv[optind]);
      optind++;
      if (optind < argc)
	{ 
	  sprintf(mess,"Too many files. Only one filename needed\n%s",USAGE);
	  printError(mess);
	}
    }
  else
    {
      sprintf(mess,"Where is the HSP file?\n%s",USAGE);
      printError(mess);
    }
}










