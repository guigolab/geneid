/*************************************************************************
*                                                                        *
*   Module: readargv                                                     *
*                                                                        *
*   Get setup options and filenames of user input.                       *
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

/*  $Id: readargv.c,v 1.8 2001-04-23 13:56:56 eblanco Exp $  */

#include "geneid.h"

extern int 	SFP,SDP,SAP,STP,
                EFP,EIP,ETP,EXP,ESP,EOP,
                VRB,
                FWD,RVS,
                GENEID, GENAMIC,
                GFF, X10,
                EVD, SRP,
                scanORF, XML, cDNA;

extern float EW;

extern char* optarg;
extern int optind;

char *USAGE="Incorrect usage:\nNAME\n\tgeneid - a program to predict genes\nSYNOPSIS\n\tgeneid\t[-bdaefitxsz]\n\t\t[-v] [-Z] [-W] [-C] [-D] [-M] [-m]\n\t\t[-P ParameterFile]\n\t\t[-h] [-X] [-G] [-o] \n\t\t[-O gff_exons_file]\n\t\t[-R gff_evidences-file]\n\t\t[-S gff_similarity_regions-file]\n\t\t[-E exonweight]\n\t\t<locus_seq_in_fasta_format>\n\n";

void printHelp()
{
  printf("Short Manual for geneid:\n");
  printf("------------------------\n\n");
  printf("Setup Options:\n\n");
 
  printf("\t-b: Print Start Codons\n");
  printf("\t-d: Print Donor Sites\n ");
  printf("\t-a: Print Acceptor Sites\n");
  printf("\t-e: Print Stop Codons\n");
  printf("\t-f: Print Initial Exons\n");
  printf("\t-i: Print Internal Exons\n");
  printf("\t-t: Print Terminal Exons\n");
  printf("\t-s: Print Single Genes\n");
  printf("\t-z: Print ORFs\n");
  printf("\t-x: Print All exons\n\n");

  printf("\t-D: Print cDNA in predicted genes\n\n");
  printf("\t-Z: Activate ORFs searching\n\n");

  printf("\t-G: Print output predictions in GFF-format\n");
  printf("\t-X: Print extended-format Output Gene Predicted \n\n");
  printf("\t-M: Print XML-format Output for Gene Prediction \n\n");
  printf("\t-m: Print XML-DTD for XML-format Output \n\n");

  printf("\t-W: Only Forward Prediction(Watson)\n");
  printf("\t-C: Only Reverse Prediction(Crick)\n");
  printf("\t-o: Only runs exon prediction(not gen prediction)\n");
  printf("\t-O: Only runs gen prediction(not exon prediction)\n\n");
  printf("\t-R: Runs geneid predictions with evidences-file\n\n");
  printf("\t-S: Runs geneid predictions with similarity regions file\n\n");

  printf("\t-E: Exon Weight parameter\n\n");

  printf("\t-P: Change the name of the Parameter File\n");
  printf("\t-v: Verbose. Print all messages\n");
  printf("\t-h: Show this Short Manual\n");
}

void printDTD()
{
  printf("<?xml version=\"1.0\" ?>
   <!-- DTD for XML format in geneid output -->

   <!-- Element declarations -->
   <!ELEMENT prediction (gene*)>
   <!ELEMENT gene ((exon+),cDNA,protein)>
   <!ELEMENT exon (site,site)>
   <!ELEMENT cDNA (#PCDATA)>
   <!ELEMENT protein (#PCDATA)>

   <!-- Attribute declarations -->
   <!ATTLIST prediction
        locus    CDATA   #REQUIRED
        length   CDATA   #IMPLIED
        source   CDATA   #IMPLIED
        date     CDATA   #IMPLIED
        genes    CDATA   #REQUIRED
	score    CDATA   #REQUIRED>

   <!ATTLIST gene
        idGene   ID      #REQUIRED
        strand   (fwd|rvs)   #IMPLIED
        nExons   CDATA   #IMPLIED
        score    CDATA   #REQUIRED>

   <!ATTLIST exon  
        idExon   ID      #REQUIRED
        type     (First | Internal | Terminal | Single) #REQUIRED
	frame    (0|1|2) #REQUIRED
        score    CDATA   #REQUIRED>

   <!ATTLIST site
        idSite   ID      #REQUIRED
        type     (Acceptor | Donor | Start | Stop) #REQUIRED
	position CDATA   #REQUIRED
        score    CDATA   #REQUIRED>\n\n");
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
  while ((c = getopt(argc,argv,"oO:bdaefitsxDzZXmMGvE:R:S:WCP:h")) != -1)
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

  /* Setup Errors: Incompatible options selected */
  if (!GENEID && geneidOpts)
    printError("Incompatible options(-O)");
  
  if (!GENAMIC && genamicOpts)
    printError("Incompatible options(-o)");

  if (!GENAMIC && !GENEID)
    printError("Incompatible options(-o|-O)");
 
  if (XML && printOptions)
    printError("Incompatible options(-M| print genic features)"); 

  if (XML && (GFF || X10))
    printError("Incompatible options(-M|-X,-G)");

  if (cDNA && GFF)
    printError("Incompatible options(-D|-G)");
  
  if (error)
    printError(USAGE);

  /* Setup Errors: Wrong number of filenames */
  /* Get the name of the file *.fasta */
  if (optind < argc)
    {
      strcpy(SequenceFile,argv[optind]);
      optind++;
      if (optind < argc)
	{ 
	  sprintf(mess,"Too many files. Only one filename needed\n%s",USAGE);
	  printError(mess);
	}
    }
  else
    if (GENEID)
      {
	sprintf(mess,"Where is the fasta file(DNA Sequence)?\n%s",USAGE);
	printError(mess);
      }
  
  /* Default parameters file if option P not used */
  if (!strcmp(ParamFile,""))
    strcpy(ParamFile,PARAMETERFILE);
}










