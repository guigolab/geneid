/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: account                                            *
*                                                              *
*   Accounting of results stats and time of execution          *
***************************************************************/ 

#include "blast2gff.h"

/* Init accounting */
account* InitAcc()
{
  account* m; 
  
  m = (account *) malloc(sizeof(account));

  /* Reset counters */

  printMess("Reset accounting");
  /* What time is it? */ 
  clock();
  (void) time(&m->tStart);

  return(m);
}
