# 
#   this file is part of the geneid distribution
#   (c) by Roderic Guigo and Enrique Blanco, 1999
#
#   makefile for geneid
#

INCLUDE= ./include
CDIR= ./src
OBJ = ./objects
HEADERS = $(INCLUDE)/geneid.h 
PROGRAM= geneid
PRODUCT= bin/$(PROGRAM)
CC= gcc
OPTS=-I$(INCLUDE) -Wall

#######

OBJECTS = $(OBJ)/readparam.o $(OBJ)/readargv.o\
	$(OBJ)/PrintSites.o $(OBJ)/PrintExons.o \
	$(OBJ)/ReadExonsGFF.o\
	$(OBJ)/GetStopCodons.o  $(OBJ)/GetSitesWithProfile.o \
	$(OBJ)/BuildInitialExons.o $(OBJ)/BuildInternalExons.o $(OBJ)/BuildTerminalExons.o \
	$(OBJ)/SortExons.o $(OBJ)/ScoreExons.o \
	$(OBJ)/FetchSequence.o $(OBJ)/RecomputePositions.o $(OBJ)/SwitchPositions.o\
	$(OBJ)/genamic.o $(OBJ)/BuildSort.o $(OBJ)/ReadGeneModel.o \
	$(OBJ)/CookingGenes.o $(OBJ)/Dictionary.o \
	$(OBJ)/CorrectExon.o $(OBJ)/SwitchFrames.o $(OBJ)/account.o \
	$(OBJ)/Output.o $(OBJ)/RequestMemory.o \
	$(OBJ)/BackupGenes.o $(OBJ)/Translate.o $(OBJ)/DumpHash.o \
	$(OBJ)/BuildSingles.o $(OBJ)/ReadSequence.o $(OBJ)/SearchEvidenceExons.o \
	$(OBJ)/ReadSR.o

#######

$(PRODUCT): $(OBJ)/$(PROGRAM).o $(OBJECTS) $(HEADERS)
	$(CC) $(OPTS) -o $(PRODUCT) $(OBJ)/$(PROGRAM).o $(OBJECTS) -lm

$(OBJ)/$(PROGRAM).o : $(CDIR)/$(PROGRAM).c  $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/$(PROGRAM).c -o $(OBJ)/$(PROGRAM).o 

$(OBJ)/ReadSequence.o : $(CDIR)/ReadSequence.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadSequence.c -o $(OBJ)/ReadSequence.o

$(OBJ)/readparam.o : $(CDIR)/readparam.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/readparam.c -o $(OBJ)/readparam.o

$(OBJ)/readargv.o :  $(CDIR)/readargv.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/readargv.c -o $(OBJ)/readargv.o

$(OBJ)/account.o :  $(CDIR)/account.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/account.c -o $(OBJ)/account.o

$(OBJ)/PrintSites.o :  $(CDIR)/PrintSites.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/PrintSites.c -o $(OBJ)/PrintSites.o

$(OBJ)/Translate.o :  $(CDIR)/Translate.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Translate.c -o $(OBJ)/Translate.o

$(OBJ)/PrintExons.o :  $(CDIR)/PrintExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/PrintExons.c -o $(OBJ)/PrintExons.o

$(OBJ)/RequestMemory.o :  $(CDIR)/RequestMemory.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/RequestMemory.c -o $(OBJ)/RequestMemory.o

$(OBJ)/ReadExonsGFF.o : $(CDIR)/ReadExonsGFF.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadExonsGFF.c -o $(OBJ)/ReadExonsGFF.o

$(OBJ)/ReadSR.o : $(CDIR)/ReadSR.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadSR.c -o $(OBJ)/ReadSR.o

$(OBJ)/GetStopCodons.o : $(CDIR)/GetStopCodons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/GetStopCodons.c -o $(OBJ)/GetStopCodons.o

$(OBJ)/GetSitesWithProfile.o : $(CDIR)/GetSitesWithProfile.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/GetSitesWithProfile.c -o $(OBJ)/GetSitesWithProfile.o

$(OBJ)/BuildInitialExons.o : $(CDIR)/BuildInitialExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildInitialExons.c -o $(OBJ)/BuildInitialExons.o

$(OBJ)/BuildInternalExons.o : $(CDIR)/BuildInternalExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildInternalExons.c -o $(OBJ)/BuildInternalExons.o

$(OBJ)/BuildTerminalExons.o : $(CDIR)/BuildTerminalExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildTerminalExons.c -o $(OBJ)/BuildTerminalExons.o

$(OBJ)/BuildSingles.o : $(CDIR)/BuildSingles.c $(HEADERS)	
	$(CC) -c $(OPTS) $(CDIR)/BuildSingles.c -o $(OBJ)/BuildSingles.o

$(OBJ)/SearchEvidenceExons.o : $(CDIR)/SearchEvidenceExons.c $(HEADERS)	
	$(CC) -c $(OPTS) $(CDIR)/SearchEvidenceExons.c -o $(OBJ)/SearchEvidenceExons.o

$(OBJ)/DumpHash.o: $(CDIR)/DumpHash.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/DumpHash.c -o $(OBJ)/DumpHash.o

$(OBJ)/BackupGenes.o : $(CDIR)/BackupGenes.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BackupGenes.c -o $(OBJ)/BackupGenes.o

$(OBJ)/SortExons.o : $(CDIR)/SortExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SortExons.c -o $(OBJ)/SortExons.o

$(OBJ)/ScoreExons.o : $(CDIR)/ScoreExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ScoreExons.c -o $(OBJ)/ScoreExons.o

$(OBJ)/FetchSequence.o : $(CDIR)/FetchSequence.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/FetchSequence.c -o $(OBJ)/FetchSequence.o

$(OBJ)/RecomputePositions.o : $(CDIR)/RecomputePositions.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/RecomputePositions.c -o $(OBJ)/RecomputePositions.o

$(OBJ)/SwitchPositions.o : $(CDIR)/SwitchPositions.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SwitchPositions.c -o $(OBJ)/SwitchPositions.o

$(OBJ)/SwitchFrames.o : $(CDIR)/SwitchFrames.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SwitchFrames.c -o $(OBJ)/SwitchFrames.o

$(OBJ)/CorrectExon.o : $(CDIR)/CorrectExon.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/CorrectExon.c -o $(OBJ)/CorrectExon.o

$(OBJ)/genamic.o : $(CDIR)/genamic.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/genamic.c -o $(OBJ)/genamic.o 

$(OBJ)/CookingGenes.o: $(CDIR)/CookingGenes.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/CookingGenes.c -o $(OBJ)/CookingGenes.o

$(OBJ)/Output.o: $(CDIR)/Output.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Output.c -o $(OBJ)/Output.o

$(OBJ)/Dictionary.o: $(CDIR)/Dictionary.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Dictionary.c -o $(OBJ)/Dictionary.o

$(OBJ)/ReadGeneModel.o: $(CDIR)/ReadGeneModel.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadGeneModel.c -o $(OBJ)/ReadGeneModel.o

$(OBJ)/BuildSort.o: $(CDIR)/BuildSort.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildSort.c -o $(OBJ)/BuildSort.o

clean:
	rm -f $(OBJ)/*.o $(PRODUCT) *~ $(INCLUDE)/*~ $(CDIR)/*~ core


