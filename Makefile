# 
#   This file is part of the geneid distribution
#   (c) by Roderic Guigo and Enrique Blanco, 2003
#
#   Makefile for geneid
#

INCLUDE= ./include
CDIR= ./src
OBJ = ./objects
BIN = ./bin
HEADERS = $(INCLUDE)/geneid.h 
PROGRAM= geneid
PRODUCT= $(BIN)/$(PROGRAM)
CC=gcc

#######
# Memory profile: pick the array-sizing constants for the target machine.
#   make             -> MEM=auto: detect THIS machine's RAM and pick a profile
#   make MEM=low      ~5 GB reserved  (small machines; tight on dense genomes)
#   make MEM=medium   ~13 GB reserved (16-24 GB machines)
#   make MEM=high     ~45 GB reserved (big servers; lazy alloc, RSS stays ~4 GB)
# Profiles override the #ifndef'd defaults in include/geneid.h via -D.
#
# *** auto detects the BUILD machine's RAM. *** If you build and run on
# different machines (e.g. a cluster login node vs a fat compute node),
# auto can guess wrong -- pass an explicit MEM=low|medium|high matching the
# RUN machine. The resolved profile is printed on every build, and
# `make print-mem` shows the detected RAM. Always `make clean` when
# switching profiles (the constants are compile-time).
#######
MEM ?= auto

# auto -> resolve to low/medium/high from physical RAM (macOS sysctl / Linux /proc)
ifeq ($(MEM),auto)
  MEM_BYTES := $(shell sysctl -n hw.memsize 2>/dev/null || awk '/MemTotal/{print $$2*1024}' /proc/meminfo 2>/dev/null || echo 0)
  MEM_RESOLVED := $(shell b="$(MEM_BYTES)"; \
    if   [ "$$b" -ge 51539607552 ]; then echo high;   \
    elif [ "$$b" -ge 17179869184 ]; then echo medium; \
    elif [ "$$b" -gt 0 ];          then echo low;     \
    else echo medium; fi)
else
  MEM_BYTES := 0
  MEM_RESOLVED := $(MEM)
endif

ifeq ($(MEM_RESOLVED),low)
  MEMFLAGS = -DLENGTHSi=220000 -DREXONS=2 -DFSORT=12 -DFDARRAY=5 -DRBSITES=150 -DRBEXONS=250
else ifeq ($(MEM_RESOLVED),medium)
  MEMFLAGS =
else ifeq ($(MEM_RESOLVED),high)
  MEMFLAGS = -DLENGTHSi=500000 -DREXONS=0.25 -DFSORT=20 -DFDARRAY=10 -DRBSITES=75 -DRBEXONS=125
else
  $(error Unknown MEM='$(MEM)' (resolved '$(MEM_RESOLVED)'). Use one of: auto, low, medium, high)
endif

$(info geneid memory profile: MEM=$(MEM) -> '$(MEM_RESOLVED)')

OPTS=-I$(INCLUDE) -Wall -O3 $(MEMFLAGS)
#OPTS=-I$(INCLUDE) -Wall -g $(MEMFLAGS)

# Keep the program the default goal (print-mem is defined after it below)
.DEFAULT_GOAL := $(PRODUCT)
#######

OBJECTS = $(OBJ)/BackupGenes.o $(OBJ)/PeakEdgeScore.o $(OBJ)/GetTranscriptTermini-usingslopes.o $(OBJ)/BuildAcceptors.o $(OBJ)/BuildU12Acceptors.o $(OBJ)/BuildDonors.o \
	$(OBJ)/BuildInitialExons.o $(OBJ)/BuildInternalExons.o $(OBJ)/BuildZeroLengthExons.o $(OBJ)/BuildUTRExons.o\
	$(OBJ)/BuildORFs.o $(OBJ)/BuildSingles.o $(OBJ)/BuildSort.o $(OBJ)/BuildTerminalExons.o \
	$(OBJ)/ComputeStopInfo.o $(OBJ)/CookingGenes.o $(OBJ)/CorrectExon.o \
	$(OBJ)/Dictionary.o $(OBJ)/DumpHash.o $(OBJ)/FetchSequence.o \
	$(OBJ)/GetSitesWithProfile.o $(OBJ)/GetStopCodons.o $(OBJ)/Output.o \
	$(OBJ)/PrintExons.o $(OBJ)/PrintSites.o $(OBJ)/ReadExonsGFF.o \
	$(OBJ)/ReadGeneModel.o $(OBJ)/ReadSequence.o $(OBJ)/ReadHSP.o $(OBJ)/RecomputePositions.o \
	$(OBJ)/RequestMemory.o $(OBJ)/ProcessHSPs.o $(OBJ)/ScoreExons.o $(OBJ)/SearchEvidenceExons.o \
	$(OBJ)/SetRatios.o $(OBJ)/SortExons.o $(OBJ)/SortSites.o $(OBJ)/SortHSPs.o $(OBJ)/SwitchFrames.o \
	$(OBJ)/SwitchPositions.o $(OBJ)/Translate.o \
	$(OBJ)/account.o $(OBJ)/beggar.o $(OBJ)/genamic.o $(OBJ)/manager.o \
	$(OBJ)/readparam.o $(OBJ)/readargv.o \

#######

$(PRODUCT): $(BIN) $(OBJ) $(OBJ)/$(PROGRAM).o $(OBJECTS) $(HEADERS)
	$(CC) $(OPTS) -o $(PRODUCT) $(OBJ)/$(PROGRAM).o $(OBJECTS) -lm

$(BIN) :
	mkdir $(BIN); 

$(OBJ) :
	mkdir $(OBJ); 

$(OBJ)/$(PROGRAM).o : $(CDIR)/$(PROGRAM).c  $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/$(PROGRAM).c -o $(OBJ)/$(PROGRAM).o 

$(OBJ)/account.o :  $(CDIR)/account.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/account.c -o $(OBJ)/account.o

$(OBJ)/beggar.o :  $(CDIR)/beggar.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/beggar.c -o $(OBJ)/beggar.o

$(OBJ)/genamic.o : $(CDIR)/genamic.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/genamic.c -o $(OBJ)/genamic.o

$(OBJ)/manager.o :  $(CDIR)/manager.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/manager.c -o $(OBJ)/manager.o

$(OBJ)/readargv.o :  $(CDIR)/readargv.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/readargv.c -o $(OBJ)/readargv.o

$(OBJ)/readparam.o : $(CDIR)/readparam.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/readparam.c -o $(OBJ)/readparam.o

$(OBJ)/BackupGenes.o : $(CDIR)/BackupGenes.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BackupGenes.c -o $(OBJ)/BackupGenes.o

$(OBJ)/GetTranscriptTermini-usingslopes.o : $(CDIR)/GetTranscriptTermini-usingslopes.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/GetTranscriptTermini-usingslopes.c -o $(OBJ)/GetTranscriptTermini-usingslopes.o

$(OBJ)/BuildAcceptors.o : $(CDIR)/BuildAcceptors.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildAcceptors.c -o $(OBJ)/BuildAcceptors.o

$(OBJ)/BuildU12Acceptors.o : $(CDIR)/BuildU12Acceptors.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildU12Acceptors.c -o $(OBJ)/BuildU12Acceptors.o

$(OBJ)/BuildDonors.o : $(CDIR)/BuildDonors.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildDonors.c -o $(OBJ)/BuildDonors.o

#$(OBJ)/BuildU12Donors.o : $(CDIR)/BuildU12Donors.c $(HEADERS)
#	$(CC) -c $(OPTS) $(CDIR)/BuildU12Donors.c -o $(OBJ)/BuildU12Donors.o

$(OBJ)/BuildInitialExons.o : $(CDIR)/BuildInitialExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildInitialExons.c -o $(OBJ)/BuildInitialExons.o

$(OBJ)/BuildInternalExons.o : $(CDIR)/BuildInternalExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildInternalExons.c -o $(OBJ)/BuildInternalExons.o

$(OBJ)/BuildZeroLengthExons.o : $(CDIR)/BuildZeroLengthExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildZeroLengthExons.c -o $(OBJ)/BuildZeroLengthExons.o

$(OBJ)/BuildUTRExons.o : $(CDIR)/BuildUTRExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildUTRExons.c -o $(OBJ)/BuildUTRExons.o

$(OBJ)/BuildORFs.o : $(CDIR)/BuildORFs.c $(HEADERS)	
	$(CC) -c $(OPTS) $(CDIR)/BuildORFs.c -o $(OBJ)/BuildORFs.o

$(OBJ)/BuildSingles.o : $(CDIR)/BuildSingles.c $(HEADERS)	
	$(CC) -c $(OPTS) $(CDIR)/BuildSingles.c -o $(OBJ)/BuildSingles.o

$(OBJ)/BuildSort.o: $(CDIR)/BuildSort.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildSort.c -o $(OBJ)/BuildSort.o

$(OBJ)/BuildTerminalExons.o : $(CDIR)/BuildTerminalExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/BuildTerminalExons.c -o $(OBJ)/BuildTerminalExons.o

$(OBJ)/ComputeStopInfo.o : $(CDIR)/ComputeStopInfo.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ComputeStopInfo.c -o $(OBJ)/ComputeStopInfo.o

$(OBJ)/CookingGenes.o: $(CDIR)/CookingGenes.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/CookingGenes.c -o $(OBJ)/CookingGenes.o

$(OBJ)/CorrectExon.o : $(CDIR)/CorrectExon.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/CorrectExon.c -o $(OBJ)/CorrectExon.o

$(OBJ)/Dictionary.o: $(CDIR)/Dictionary.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Dictionary.c -o $(OBJ)/Dictionary.o

$(OBJ)/DumpHash.o: $(CDIR)/DumpHash.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/DumpHash.c -o $(OBJ)/DumpHash.o

$(OBJ)/FetchSequence.o : $(CDIR)/FetchSequence.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/FetchSequence.c -o $(OBJ)/FetchSequence.o

$(OBJ)/GetSitesWithProfile.o : $(CDIR)/GetSitesWithProfile.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/GetSitesWithProfile.c -o $(OBJ)/GetSitesWithProfile.o

$(OBJ)/GetStopCodons.o : $(CDIR)/GetStopCodons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/GetStopCodons.c -o $(OBJ)/GetStopCodons.o

$(OBJ)/Output.o: $(CDIR)/Output.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Output.c -o $(OBJ)/Output.o

$(OBJ)/PeakEdgeScore.o :  $(CDIR)/PeakEdgeScore.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/PeakEdgeScore.c -o $(OBJ)/PeakEdgeScore.o

$(OBJ)/PrintExons.o :  $(CDIR)/PrintExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/PrintExons.c -o $(OBJ)/PrintExons.o

$(OBJ)/PrintSites.o :  $(CDIR)/PrintSites.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/PrintSites.c -o $(OBJ)/PrintSites.o

$(OBJ)/ProcessHSPs.o :  $(CDIR)/ProcessHSPs.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ProcessHSPs.c -o $(OBJ)/ProcessHSPs.o

$(OBJ)/ReadGeneModel.o: $(CDIR)/ReadGeneModel.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadGeneModel.c -o $(OBJ)/ReadGeneModel.o

$(OBJ)/ReadExonsGFF.o : $(CDIR)/ReadExonsGFF.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadExonsGFF.c -o $(OBJ)/ReadExonsGFF.o

$(OBJ)/ReadSequence.o : $(CDIR)/ReadSequence.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadSequence.c -o $(OBJ)/ReadSequence.o

$(OBJ)/ReadHSP.o : $(CDIR)/ReadHSP.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ReadHSP.c -o $(OBJ)/ReadHSP.o

$(OBJ)/RecomputePositions.o : $(CDIR)/RecomputePositions.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/RecomputePositions.c -o $(OBJ)/RecomputePositions.o

$(OBJ)/RequestMemory.o :  $(CDIR)/RequestMemory.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/RequestMemory.c -o $(OBJ)/RequestMemory.o

$(OBJ)/ScoreExons.o : $(CDIR)/ScoreExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/ScoreExons.c -o $(OBJ)/ScoreExons.o

$(OBJ)/SearchEvidenceExons.o : $(CDIR)/SearchEvidenceExons.c $(HEADERS)	
	$(CC) -c $(OPTS) $(CDIR)/SearchEvidenceExons.c -o $(OBJ)/SearchEvidenceExons.o

$(OBJ)/SetRatios.o : $(CDIR)/SetRatios.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SetRatios.c -o $(OBJ)/SetRatios.o

$(OBJ)/SortExons.o : $(CDIR)/SortExons.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SortExons.c -o $(OBJ)/SortExons.o

$(OBJ)/SortSites.o : $(CDIR)/SortSites.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SortSites.c -o $(OBJ)/SortSites.o

$(OBJ)/SortHSPs.o : $(CDIR)/SortHSPs.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SortHSPs.c -o $(OBJ)/SortHSPs.o

$(OBJ)/SwitchFrames.o : $(CDIR)/SwitchFrames.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SwitchFrames.c -o $(OBJ)/SwitchFrames.o

$(OBJ)/SwitchPositions.o : $(CDIR)/SwitchPositions.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/SwitchPositions.c -o $(OBJ)/SwitchPositions.o

$(OBJ)/Translate.o :  $(CDIR)/Translate.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Translate.c -o $(OBJ)/Translate.o

clean:
	rm -f $(OBJ)/*.o $(PRODUCT) *~ $(INCLUDE)/*~ $(CDIR)/*~ core;
	rmdir $(BIN) $(OBJ);

# Show the active memory profile, detected RAM, and the -D flags it injects
print-mem:
	@echo "MEM=$(MEM)  resolved=$(MEM_RESOLVED)  detected_RAM_bytes=$(MEM_BYTES)"
	@echo "MEMFLAGS=$(MEMFLAGS)"
.PHONY: print-mem clean

