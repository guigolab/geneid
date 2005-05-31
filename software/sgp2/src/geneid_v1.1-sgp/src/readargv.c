/*************************************************************************
*                                                                        *
*   Module: readargv                                                     *
*                                                                        *
*   Read set up options and filenames from user input                    *
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

/*  $Id: readargv.c,v 1.2 2005-05-31 13:21:11 arnau Exp $  */

#include "geneid.h"

/* geneid.c external vars */
extern int 	SFP,SDP,SAP,STP,
            EFP,EIP,ETP,EXP,ESP,EOP,
            VRB,
            FWD,RVS,
            GENEID, GENAMIC,
            GFF, X10,
            EVD, SRP, BEG,
            scanORF, XML, cDNA;
extern float EW;

/* required by getopts */
extern char* optarg;
extern int optind;

char* USAGE="NAME\n\tgeneid - a program to predict genes in genomic sequences\nSYNOPSIS\n\tgeneid\t[-bdaefitxsz] [-D] [-Z]\n\t\t[-G] [-X] [-M] [-m]\n\t\t[-WC] [-o]\n\t\t[-O gff_exons_file]\n\t\t[-R gff_evidences-file]\n\t\t[-S gff_similarity_regions-file]\n\t\t[-E exonweight] [-P ParameterFile]\n\t\t[-Bv] [-h]\n\t\t<locus_seq_in_fasta_format>\nRELEASE\n\tgeneid v 1.1a\n";

void printHelp()
{
  printf(USAGE);

  printf ("OPTIONS\n");
  
  printf("\t-b: Output Start codons\n");
  printf("\t-d: Output Donor splice sites\n");
  printf("\t-a: Output Acceptor splice sites\n");
  printf("\t-e: Output Stop codons\n");
  
  printf("\t-f: Output Initial exons\n");
  printf("\t-i: Output Internal exons\n");
  printf("\t-t: Output Terminal exons\n");
  printf("\t-s: Output Single genes\n");
  printf("\t-x: Output all predicted exons\n");
  printf("\t-z: Output Open Reading Frames\n\n");
  
  printf("\t-D: Output genomic sequence of exons in predicted genes\n\n");
  
  printf("\t-G: Use GFF format to print predictions\n");
  printf("\t-X: Use extended-format to print gene predictions\n");
  printf("\t-M: Use XML format to print gene predictions\n");
  printf("\t-m: Show DTD for XML-format output \n\n");
  
  printf("\t-W: Only Forward sense prediction (Watson)\n");
  printf("\t-C: Only Reverse sense prediction (Crick)\n");
  printf("\t-o: Only running exon prediction (disable gene prediction)\n");
  printf("\t-O exons_filename: Only running gene prediction (not exon prediction)\n");
  printf("\t-Z: Activate Open Reading Frames searching\n\n");
  
  printf("\t-R exons_filename: Provide annotations to improve predictions\n");
  printf("\t-S SR_filename: Using information about sequence homology to improve predictions\n\n");
  
  printf("\t-E: Adding this value to the exon weight parameter (see parameter file)\n");
  printf("\t-P parameter_file: Use other than default parameter file (human)\n\n");
  
  printf("\t-B: Display memory required to execute geneid given a sequence\n");
  printf("\t-v: Verbose. Display info messages\n");
  printf("\t-h: Show this help\n");

  printf ("AUTHORS\n");
  printf("\tgeneid v 1.1 has been developed by Enrique Blanco and Roderic Guigo.\n\tParameter files have been computed by G. Parra. Any bug or suggestion\n\tcan be reported to geneid@imim.es\n");

  printf("\n\n\n");
}

void printDTD()
{

  printf("<?xml version=\"1.0\" ?>");
  printf("<!-- DTD for XML format in geneid output -->");
  
  printf("<!-- Element declarations -->");
  printf("<!ELEMENT prediction (gene*)>");
  printf("<!ELEMENT gene ((exon+),cDNA,protein)>");
  printf("<!ELEMENT exon (site,site)>");
  printf("<!ELEMENT cDNA (#PCDATA)>");
  printf("<!ELEMENT protein (#PCDATA)>");
  
  printf("<!-- Attribute declarations -->");
  printf("<!ATTLIST prediction");
  printf("\tlocus    CDATA   #REQUIRED");
  printf("\tlength   CDATA   #IMPLIED");
  printf("\tsource   CDATA   #IMPLIED");
  printf("\tdate     CDATA   #IMPLIED");
  printf("\tgenes    CDATA   #REQUIRED");
  printf("\tscore    CDATA   #REQUIRED>");
  
  printf("<!ATTLIST gene");
  printf("\tidGene   ID      #REQUIRED");
  printf("\tstrand   (fwd|rvs)   #IMPLIED");
  printf("\tnExons   CDATA   #IMPLIED");
  printf("\tscore    CDATA   #REQUIRED>");
  
  printf("<!ATTLIST exon");  
  printf("\tidExon   ID      #REQUIRED");
  printf("\ttype     (First | Internal | Terminal | Single) #REQUIRED");
  printf("frame    (0|1|2) #REQUIRED");
  printf("\tscore    CDATA   #REQUIRED>");
  
  printf("<!ATTLIST site");
  printf("\tidSite   ID      #REQUIRED");
  printf("\ttype     (Acceptor | Donor | Start | Stop) #REQUIRED");
  printf("position CDATA   #REQUIRED");
  printf("\tscore    CDATA   #REQUIRED>\n\n");
  
}

void readargv (int argc,char* argv[],
			   char* ParamFile, char* SequenceFile,
			   char* ExonsFile, char* SRFile) 
{
  int c;
  int error=0;
  int geneidOpts = 0;
  int genamicOpts = 0;
  int printOptions =0;
  char mess[MAXSTRING];
  
  /* Reading setup options */
  while ((c = getopt(argc,argv,"oO:bdaefitsxDzZXmMGBvE:R:S:WCP:h")) != -1)
    switch(c)
      {
      case 'Z': scanORF++;
		geneidOpts++;
		break;
      case 'o': GENAMIC--;
		break;
      case 'O': GENEID--;
		strcpy (ExonsFile,optarg);
		break;
      case 'R': EVD++;
		strcpy (ExonsFile,optarg);
		geneidOpts++;
		break;
      case 'S': SRP++;
		strcpy (SRFile,optarg);
		geneidOpts++;
		break;
      case 'E': EW = atof(optarg);
		geneidOpts++;
		break;	
      case 'b': SFP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'd': SDP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'a': SAP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'e': STP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'f': EFP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'i': EIP++;
		geneidOpts++;
		printOptions++;
		break;
      case 't': ETP++;
		printOptions++;
		geneidOpts++;
		break;
      case 'x': EXP++;
		geneidOpts++;
		printOptions++;
		break;
      case 's': ESP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'z': EOP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'B': BEG++; 
		break;
      case 'v': VRB++;
		break;
      case 'W': RVS--;
		geneidOpts++;
		break; 
      case 'C': FWD--;
		geneidOpts++;
		break;
      case 'P': strcpy (ParamFile,optarg); 
		break;
      case 'G': GFF++;
		break;
      case 'D': cDNA++;
		genamicOpts++;
		break;
      case 'X': X10++;
		genamicOpts++;
		break;
      case 'M': XML++;
        genamicOpts++;
		break;
      case '?':error++;
		break; 
      case 'm': printDTD();
		exit(1);
		break;
      case 'h': printHelp();
		exit(1);
		break;
      }
  
  /* Setup Errors (a): Incompatible options selected */
  if (!GENEID && geneidOpts)
    printError("Incompatible options (with -O)");
  
  if (!GENAMIC && genamicOpts)
    printError("Incompatible options (with -o)");
  
  if (!GENAMIC && !GENEID)
    printError("Incompatible options (-o | -O)");
 
  if (XML && printOptions)
    printError("Incompatible options (-M | print gene features)"); 

  if (XML && (GFF || X10))
    printError("Incompatible options (XML and other output formats)");

  if (cDNA && GFF)
    printError("Incompatible options( -D | -G)");
  
  if (error)
	{
	  sprintf(mess,"Wrong usage of options\n%s",USAGE);
	  printError(mess);
	}

  /* Setup Errors (b): Wrong number of filenames */
  /* Read the name of the input fasta file */
  if (optind < argc)
    {
      strcpy(SequenceFile,argv[optind]);
      optind++;
      if (optind < argc)
		{ 
          sprintf(mess,"Only one filename required but more than one presented\n%s",USAGE);
		  printError(mess);
		}
    }
  else
    if (GENEID)
      {
        sprintf(mess,"One filename is needed (DNA sequence, fasta format)\n%s",USAGE);
		printError(mess);
      }
  
  /* Default parameter file selected if option -P not used */
  if (!strcmp(ParamFile,""))
    strcpy(ParamFile,PARAMETERFILE);
}
