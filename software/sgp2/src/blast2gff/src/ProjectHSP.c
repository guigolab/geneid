/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: ProjectHSP                                         *
*                                                              *
*   Project the HSP set against a SR-regions line              *
***************************************************************/

#include "blast2gff.h"

void printSR(sr_t* t, long i)
{
  sr_t* tmp; 

  printf("Resultados iteracion %ld\n",i);
  
  tmp = t;
  while(tmp != NULL)
    {
      printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\n",
	     NONAME,
	     VERSION,
	     SR,
	     tmp->Start,
	     tmp->End,
	     tmp->Score,
	     tmp->Strand,
	     tmp->Frame); 
      tmp = tmp->next;
    }
}

int FreeItems(sr_t* q)
{
  int numElems;
  
  if (q == NULL)
    numElems = 0;
  else
    {
      numElems = FreeItems(q->next);
      numElems++;
      free(q);
    }
  return(numElems);
}

/* Returns the number of free SR */
int CutAndPaste(sr_t** Left, sr_t** Right,
                sr_t** insLeft, sr_t** insRight,
                sr_t** sr)
{
  sr_t* tmp;
  int numElems;
  
/*    printf("nr0: %ld %ld %d\n",Left->Start,Left->End,Left->Score); */
/*    printf("sr: %ld %ld %d\n",sr->Start,sr->End,sr->Score); */
/*    printf("firstSR: %ld %ld %d\n",Right->Start,Right->End,Right->Score); */
/*    printf("fnSR: %ld %ld %d\n",insLeft->Start,insLeft->End,insLeft->Score); */
/*    printf("newSR: %ld %ld %d\n",insRight->Start,insRight->End,insRight->Score); */
  
  /* Cut and Paste at the beginning of the list */
  if (*Left == NULL)
    {
      /* tmp used for freeing the useless nodes*/
      tmp = *sr;
      *sr = *insLeft;
    }
  else
    {
      /* tmp used for freeing the useless nodes*/
      tmp = (*Left)->next;
      (*Left)->next = *insLeft;
    }

  
  (*insRight)->next = (*Right)->next;
  
  /* Free nodes between [Left->next | sr] and Right */
  (*Right)->next = NULL;

  numElems = FreeItems(tmp);
  
  return(numElems);
}

void Chain(sr_t** newSR, sr_t* auxSR, sr_t** firstNewSR)
{
  /* Chaining the new SR to the chain */
  if (*newSR != NULL)
    {
      (*newSR)->next = auxSR;
      (*newSR) = (*newSR)->next;
    }
  else
    {
      /*        printMess("INT: iniciando la nueva cadena"); */
      
      *newSR = auxSR;
      *firstNewSR = auxSR;
    }
}

sr_t* RequestMemoryNewSR()
{
  sr_t* s;
  
  if ((s = (sr_t *) malloc(sizeof(sr_t))) == NULL)
    printError("Not enough space to hold one new SR");
  
  return(s);
}

void project1HSP(hsp** hsps, long nHSP, packSR *allSr)
{
  int index;
  long i;
  sr_t* beforeFirstSR;
  sr_t* firstSR;
  sr_t* auxSR;
  sr_t* newSR;
  sr_t* firstNewSR;
  int numElems;
  
  /* For each HSP do... */
  for(i=0; i<nHSP; i++)
    {
      /*    printf("\nTratando HSP %ld: %ld %ld %d %hd\n", */
      /*  	     i,hsps[i]->StartQ,hsps[i]->EndQ,hsps[i]->Score,hsps[i]->FrameQ); */
      
      index = (hsps[i]->StrandQ == '+')?
      hsps[i]->FrameQ % 3 :  
	(hsps[i]->FrameQ % 3) + FRAMES;
      
      /* Searching into the list, the first SR after this HP */
      /* beforeFirstSR: the last before the first SR after this HP */
      /* Possible values for firstSR and beforeFirstSR */
      /* a. Reading from the start (1 element) */
      if (allSr->nr[index] > 0 && allSr->nr0[index]==NULL)
	{
	  firstSR = allSr->sr[index];
	  beforeFirstSR = NULL;
	}
      else
	{
	  /* b. Empty list */
	  if (allSr->nr[index] == 0 && allSr->nr0[index] == NULL)         
	    firstSR = NULL;         
	  /* c. Normal situation: inside the list */
	  else
	    firstSR = allSr->nr0[index];
	  
	  beforeFirstSR = firstSR;
	}
      
      /* scanning the list of SR */
      while(firstSR != NULL && firstSR->End < hsps[i]->StartQ)
	{
	  beforeFirstSR = firstSR;
	  firstSR = firstSR->next;
	}
      
      /* last SR before first selected SR */
      allSr->nr0[index] = beforeFirstSR;
      
      /* Empty list: Insert first SR (the same HP) */
      if(firstSR == NULL && beforeFirstSR == NULL)
	{
	  /*  printMess("Lista vacia: creando primer SR\n"); */
	  
	  /* Creating a new SR placed after end of this list */
	  auxSR = RequestMemoryNewSR();
	  
	  /* Setting HSP values */
	  auxSR->Start = hsps[i]->StartQ;
	  auxSR->End = hsps[i]->EndQ;
	  auxSR->Strand = hsps[i]->StrandQ;
	  auxSR->Frame = hsps[i]->FrameQ;
	  auxSR->Score = hsps[i]->Score;
	  strcpy(auxSR->Locus,hsps[i]->Subject);
	  auxSR->next = NULL;
	  
	  /* Insert at the beginning */
	  allSr->sr[index] = auxSR;
	  allSr->nr[index]  = 1;
	}
      else
	/* Insert SR at the end of the list */
	if(firstSR == NULL && beforeFirstSR != NULL)
	  {
	    /*  printMess("Recorrida lista completa: insertar al final\n"); */
	    
	    /* Creating a new SR placed after end of this list */
	    auxSR = RequestMemoryNewSR();
	    
	    /* Setting HSP values */
	    auxSR->Start = hsps[i]->StartQ;
	    auxSR->End = hsps[i]->EndQ;
	    auxSR->Strand = hsps[i]->StrandQ;
	    auxSR->Frame = hsps[i]->FrameQ;
	    auxSR->Score = hsps[i]->Score;
	    strcpy(auxSR->Locus,hsps[i]->Subject);
	    auxSR->next = NULL;
	    
	    /* Insert at the end */
	    beforeFirstSR->next = auxSR;
	    allSr->nr[index]++;
	  }
	else
	  {
	    /* We work with the sublist starting at firstSR */
	    /* newSR is the current created node */
	    /* firstNewSR is the beginning of the new SR chain */
	    newSR = NULL;
	    firstNewSR = NULL;
	    
	    /*  printMess("Intersecciones"); */
	    /*  	    printf("INT: SR1: %ld %ld\n",firstSR->Start,firstSR->End); */
	    
	    /* 1. Creating a new SR: START(firstSR) - START(HSP) */
	    if (firstSR->Start < hsps[i]->StartQ)
	      {
		newSR = RequestMemoryNewSR();
		firstNewSR = newSR;
		
		/* Setting HSP values */
		newSR->Start = firstSR->Start;
		newSR->End = hsps[i]->StartQ-1;
		newSR->Strand = firstSR->Strand;
		newSR->Frame = firstSR->Frame;
		newSR->Score = firstSR->Score;
		strcpy(newSR->Locus,hsps[i]->Subject);
		newSR->next = NULL; 
		
		allSr->nr[index]++;
		
		/*  printf("INT: Generado primer bloque\n"); */
		/*  		printf("%ld %ld\n",newSR->Start,newSR->End); */
	      }
	    
	    /* 2. Creating the rest of intersections until last*/
	    auxSR = RequestMemoryNewSR();
	    auxSR->Start = hsps[i]->StartQ;
	    
	    while((hsps[i]->EndQ > firstSR->End) && (firstSR->next != NULL))
	      {
		/*  printf("INT: Dentro del while\n"); */
		auxSR->End = firstSR->End;
		auxSR->Score = MAX(firstSR->Score,hsps[i]->Score);
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		/*  auxSR->Locus = strcat(firstSR->Locus,hsps[i]->Subject); */
		strcpy(auxSR->Locus,hsps[i]->Subject);
		auxSR->next = NULL;
		
		allSr->nr[index]++;
		
		/*  printf("INT: Creado nuevo SR(while)\n"); */
		
		Chain(&newSR,auxSR,&firstNewSR);
		
		/* jumping to the next SR on the chain */
		firstSR = firstSR->next;
		
		auxSR = RequestMemoryNewSR();
		auxSR->Start = firstSR->Start;
	      }/* endwhile */
	    
	    /* 3. Creating the last intersections */
	    /* NO Sobresale por detras el SR */
	    if (hsps[i]->EndQ >= firstSR->End)
	      {
		auxSR->End = firstSR->End;
		auxSR->Score = MAX(firstSR->Score,hsps[i]->Score);
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		/*  auxSR->Locus = strcat(firstSR->Locus,hsps[i]->Subject); */
		strcpy(auxSR->Locus,hsps[i]->Subject);
		auxSR->next = NULL;
		
		/*  printf("auxSR %ld %ld %d\n",auxSR->Start,auxSR->End,auxSR->Score); */

		allSr->nr[index]++;
		
		Chain(&newSR, auxSR,&firstNewSR); 

		/* Coinciden los finales */
		if (hsps[i]->EndQ == firstSR->End)
		  {}/*  printf("INT2: igual final. No se crean mas\n"); */
		else
		  {
		    /* Creating a new SR: END(firstSR) - END(HSP) */
 		    auxSR = RequestMemoryNewSR();
		    auxSR->Start = firstSR->End+1;
		    auxSR->End = hsps[i]->EndQ;
		    auxSR->Strand = firstSR->Strand;
		    auxSR->Frame = firstSR->Frame;
		    auxSR->Score = hsps[i]->Score;
		    strcpy(auxSR->Locus,hsps[i]->Subject);
		    auxSR->next = NULL;
		    
		    allSr->nr[index]++;
		    
		    Chain(&newSR,auxSR,&firstNewSR);
		  }
	      }
	    /* Sobresale por detras el SR */
	    else
	      {
		auxSR->End = hsps[i]->EndQ;
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		auxSR->Score = MAX(firstSR->Score,hsps[i]->Score);
		/*  auxSR->Locus = strcat(firstSR->Locus,hsps[i]->Subject); */
		strcpy(hsps[i]->Subject,auxSR->Locus);
		auxSR->next = NULL;
		
		allSr->nr[index]++;
		
		Chain(&newSR,auxSR,&firstNewSR);
		
		auxSR = RequestMemoryNewSR();
		
		auxSR->Start = hsps[i]->EndQ +1;
		auxSR->End = firstSR->End;
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		auxSR->Score = firstSR->Score;
		strcpy(auxSR->Locus,firstSR->Locus);
		auxSR->next = NULL;

		allSr->nr[index]++;

		Chain(&newSR,auxSR,&firstNewSR);
	      }

	    /* Recompute SR chain using saved before backup pointers */
	    /*  printMess("Cut and Paste!"); */
	    numElems = CutAndPaste(&allSr->nr0[index], 
				   &firstSR,
				   &firstNewSR, 
				   &newSR,
				   &allSr->sr[index]);
	    
	    /*  printf("%d erased elems\n",numElems); */

	    allSr->nr[index] = allSr->nr[index] - numElems;
	    
	  }/* end working with sublist*/
      
      /*  printSR(allSr->sr[index],i); */
    }/* endfor */
}

void ProjectHSP(packHSP *allHsp, packSR *allSr)
{
  int i;
  
  /* projecting all the strands and frames */ 
  for(i=0; i<STRANDS*FRAMES; i++)
    project1HSP(allHsp->hsps[i],allHsp->nHsps[i],allSr);
}
