# 
#   this file is part of the geneid-a0 distribution
#   (c) Enrique Blanco, 2000
#
#   makefile for blast2gff
#

INCLUDE= ./include
CDIR= ./src
OBJ = ./objects
BIN = ./bin
HEADERS = $(INCLUDE)/blast2gff.h 
PROGRAM= blast2gff
PRODUCT= $(BIN)/$(PROGRAM)
CC= gcc
OPTS=-I$(INCLUDE) -Wall

#######

OBJECTS = $(OBJ)/readargv.o $(OBJ)/account.o $(OBJ)/Output.o $(OBJ)/SortHSP.o\
	$(OBJ)/RequestMemory.o $(OBJ)/ReadHSP.o $(OBJ)/JoinSR.o $(OBJ)/ProjectHSP.o\
	$(OBJ)/SortSR.o\
#######

$(PRODUCT): $(BIN) $(OBJ) $(OBJ)/$(PROGRAM).o $(OBJECTS) $(HEADERS)
	$(CC) $(OPTS) -o $(PRODUCT) $(OBJ)/$(PROGRAM).o $(OBJECTS) -lm

$(BIN) :
	mkdir $(BIN); 

$(OBJ) :
	mkdir $(OBJ); 

$(OBJ)/$(PROGRAM).o : $(CDIR)/$(PROGRAM).c  $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/$(PROGRAM).c -o $(OBJ)/$(PROGRAM).o

$(OBJ)/readargv.o :  $(CDIR)/readargv.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/readargv.c -o $(OBJ)/readargv.o

$(OBJ)/RequestMemory.o :  $(CDIR)/RequestMemory.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/RequestMemory.c -o $(OBJ)/RequestMemory.o

$(OBJ)/ReadHSP.o : $(CDIR)/ReadHSP.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadHSP.c -o $(OBJ)/ReadHSP.o

$(OBJ)/SortHSP.o : $(CDIR)/SortHSP.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SortHSP.c -o $(OBJ)/SortHSP.o

$(OBJ)/SortSR.o : $(CDIR)/SortSR.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SortSR.c -o $(OBJ)/SortSR.o

$(OBJ)/ProjectHSP.o : $(CDIR)/ProjectHSP.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ProjectHSP.c -o $(OBJ)/ProjectHSP.o

$(OBJ)/JoinSR.o : $(CDIR)/JoinSR.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/JoinSR.c -o $(OBJ)/JoinSR.o

$(OBJ)/Output.o: $(CDIR)/Output.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Output.c -o $(OBJ)/Output.o

$(OBJ)/account.o :  $(CDIR)/account.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/account.c -o $(OBJ)/account.o

clean:
	rm -f $(OBJ)/*.o $(BIN)/$(PROGRAM) *~ $(INCLUDE)/*~ $(CDIR)/*~ core 
	rmdir $(BIN) $(OBJ);


