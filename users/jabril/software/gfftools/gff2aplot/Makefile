################################################################################
#
# WEAVING AND TANGLING FROM NOWEB DOCUMENTS
#
################################################################################
#
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#        1         2         3         4         5         6         7         8
#
# $Id: Makefile,v 1.2 2003-06-10 10:36:29 jabril Exp $
#
# RUNNING on GNU Make version 3.79.1
#
################################################################################
#
### MAKE VARS

.PHONY    : clean cleanall forceall
.DEFAULT  : all
.PRECIOUS : %.sh %.awk %.pl %.pm %.sql \
            %.nw %.tex %.bib %.ps %.tgz %.gz \
            %.rc %.gff %.gb %.fa \
            .bash_VARS README
.SECONDARY : %.shs %.awks %.pls %.pms %.dvi %.s.ps
.SUFFIXES :

#
### PATHS

    BBIN = /bin
    UBIN = /usr/bin
    LBIN = /usr/local/bin

#
### INTERPRETERS

    DATE = $(shell date)
    USER = $(shell whoami)
    GAWK = $(firstword $(shell which gawk))
   GAWKV = $(shell $(GAWK) --version | \
             $(GAWK) '$$0 ~ /^GNU Awk/ { print $$0; }')
    BASH = $(firstword $(shell which bash))
   BASHV = $(shell $(BASH) --version | \
             $(GAWK) '$$0 ~ /^GNU bash/ { print $$0; }')
    PERL = $(firstword $(shell which perl))
   PERLV = $(shell $(PERL) --version | \
             $(GAWK) '$$0 ~ /^This/ { \
                        gsub(/^This is perl,[ \t]+/,"",$$0); \
                        print $$0; \
                      }') #'
#
### COMMANDS

     V = -v

   RM = $(firstword $(shell which rm   )) $(V)
   RD = $(firstword $(shell which rmdir)) $(V)
   CP = $(firstword $(shell which cp   )) $(V)
   MV = $(firstword $(shell which mv   )) $(V)
   LN = $(firstword $(shell which ln   )) $(V)
   MD = $(firstword $(shell which mkdir)) $(V)
  CHM = $(firstword $(shell which chmod)) $(V) -c
  CHO = $(firstword $(shell which chown)) $(V) -c

 TTAR = $(firstword $(shell which tar  )) $(V) -ztf
 CTAR = $(firstword $(shell which tar  )) $(V) -zcf
 XTAR = $(firstword $(shell which tar  )) $(V) -zxf
 GZIP = $(firstword $(shell which gzip )) $(V)

  SSH = $(firstword $(shell which ssh  )) -oProtocol=2
  SCP = $(firstword $(shell which scp  )) -oProtocol=2

 LATEX = $(UBIN)/latex
BIBTEX = $(UBIN)/bibtex
 DVIPS = $(UBIN)/dvips -D 600 -t a4
 PSNUP = $(UBIN)/psnup -pa4 -Pa4 -d -2 -m0.25cm
TANGLE = $(LBIN)/no tangle
 WEAVE = $(LBIN)/no weave
  CPIF = $(LBIN)/cpif

#
### LOCAL PATHS

    WORK = .
   NOWEB = $(WORK)/noweb
    DOCS = $(WORK)/docs
     SRC = $(WORK)/src
     BIN = $(WORK)/bin
     MAN = $(BIN)/man# cwd already set to $BIN
    DATA = $(WORK)/data
    TEST = $(WORK)/tests

 MAINDIR = $(NOWEB) $(DOCS) $(SRC) $(BIN) $(DATA) $(TEST)

#
### LOCAL VARS

   PROJECT = gff2aplot.nw
   TEXDEFS = defs.tex
   TEXTODO = todo.tex

   TEXDOCS = $(addprefix $(DOCS)/, $(PROJECT:%.nw=%.tex))
    RMDOCS = $(addprefix $(DOCS)/, $(TEXDEFS)  \
                                   $(TEXTODO)) \
             $(TEXDOCS) \
             $(wildcard $(BIN)/$(MAN)/*.3)

    BIBTEX = $(addprefix $(DOCS)/, $(PROJECT:%.nw=%.bib))

#
### SCRIPT FILES

     SRCEXT = .shs .awks .pls .pms
     BINEXT = .sh .awk .pl .pm

#
### DATA FILES

   BASHVARS = $(addprefix $(WORK)/, .bash_VARS)
 BASHVARS_R = "BASH: Environment Variables"

#
### DEFINES

define echo_weave
	echo "### Weaving from \"$<\" -> \"$@\"" ;
endef

define echo_tangle
	echo "### Tangling from \"$(<)\" -> \"$(@)\"" ;
endef

#
### MAKE RULES

all : header main trailer

man : header mkman trailer

docs : header weaving trailer

code : header tangling trailer

data : header datasets trailer

forceall : header cleansrc main trailer

clean : header cleansrc trailer

cleanall : header cleanbin cleansrc cleandoc trailer

cleanbin :
	-@$(RM) -rf $(BINCODE) $(SQLCODE) ;

cleansrc :
	-@$(RM) -rf $(RMCODE) ;

cleandoc :
	-@$(RM) -rf $(RMDOCS) ;

main : dirs weaving tangling datasets

#
### MAIN RULES #####################

dirs : dirsmsg mkdirs

dirsmsg :
	@echo "###";
	@echo "### CHECKING DIRS...";
	@echo "###";

mkdirs : $(MAINDIR)
	@for D in $(MAINDIR); \
	   do { \
	     echo "#---> $$D"; \
	     [ -d "$$D" ] || $(MD) -p $$D; \
	     }; \
	   done;

#
### WEAVING DOCS #####################

weaving : weavemsg report

weavemsg :
	@echo "###";
	@echo "### WEAVING DOCS...";
	@echo "###";

report : $(addprefix $(DOCS)/, $(PROJECT:%.nw=%.s.ps))

$(TEXDOCS) : $(addprefix $(DOCS)/, $(TEXDEFS))

$(addprefix $(DOCS)/, $(REPORT:%.nw=%.tex)) : $(addprefix $(DOCS)/, $(TEXTODO))

$(addprefix $(DOCS)/, $(REPORT:%.nw=%.tex)) : \
                      $(addprefix $(NOWEB)/, $(REPORT)) $(BIBTEX)
	@$(echo_weave)
	@$(WEAVE) -t4 -delay -index -filter 'elide "HIDE: *"' $< | $(CPIF) $@ ;

$(addprefix $(DOCS)/, $(TEXDEFS)) : $(addprefix $(NOWEB)/, $(REPORT))
	@$(echo_tangle)
	@$(TANGLE) -R"HIDE: LaTeX new definitions" $< | \
       $(CPIF) $(addprefix $(DOCS)/, $(TEXDEFS)) ;

$(addprefix $(DOCS)/, $(TEXTODO)) : $(addprefix $(NOWEB)/, $(REPORT))
	@$(echo_tangle)
	@$(TANGLE) -R"HIDE: TODO" $< | \
       $(CPIF) $(addprefix $(DOCS)/, $(TEXTODO)) ;


#
### COMPILING DOCS

define whole_latex
	@echo "### WHOLE LaTeXing on \"$(<)\" -> \"$(@)\"" ;
	@cd $(DOCS); pwd; \
	 texfile=$(<F); \
	 logfile=$(<F:.tex=.log); \
	 bibfile=$(<F:.tex=); \
	 echo "### ...LaTeXing" ; \
	 $(LATEX) $$texfile ; \
	 if [ "`grep undefined $$logfile`" != "" ]; \
	   then \
	     echo "### ...BibTeXing + LaTeXing" ; \
	     $(BIBTEX) $$bibfile ; \
	     $(LATEX)  $$texfile ; \
	   fi ; \
	 if [ "`grep Rerun $$logfile`" != "" ]; \
	   then \
	     echo "### ...LaTeXing (when Rerun found)" ; \
	     $(LATEX) $$texfile ; \
	   fi ;
endef

$(addprefix $(DOCS)/, $(REPORT:%.nw=%.dvi)) : $(addprefix $(DOCS)/, $(REPORT:%.nw=%.tex))
	$(whole_latex)

$(addprefix $(DOCS)/, %.ps) : $(addprefix $(DOCS)/, %.dvi)
	@echo "### Running DVIPS on \"$<\" -> \"$@\"" ;
	@cd $(DOCS); pwd; \
	 dvifile=$(<F); \
	 psfile=$(@F); \
	 $(DVIPS) $$dvifile -o $$psfile ; # -Z # -E -O 0cm,0cm -T 21cm,29.7cm

$(addprefix $(DOCS)/, %.s.ps) : $(addprefix $(DOCS)/, %.ps)
	@echo "### Running PSNUP on \"$<\" -> \"$@\"" ;
	@cd $(DOCS); pwd; \
	 psfile=$(<F); \
	 pssfile=$(@F); \
	 $(PSNUP) $$psfile $$pssfile ;

#
### TANGLING CODE #####################

tangling : tanglemsg mainprog auxprogs isexec

tanglemsg :
	@echo "###";
	@echo "### TANGLING CODE...";
	@echo "###";

mainprog:

auxprogs:

#
### FINISHING CODE

$(addprefix $(BIN)/, %.sh) : $(addprefix $(SRC)/, %.shs)
	@echo "### Finishing BASH script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(BASH)"; cat $< ) | $(CPIF) $@;

$(addprefix $(BIN)/, %.awk) : $(addprefix $(SRC)/, %.awks)
	@echo "### Finishing GAWK script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(GAWK) -f"; cat $< ) | $(CPIF) $@;

$(addprefix $(BIN)/, %.pl) : $(addprefix $(SRC)/, %.pls)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL)"; cat $< ) | $(CPIF) $@;
	@$(PERL) -c $@;
    # -w is deprecated as 'use warnings' has been included into use strict

$(addprefix $(BIN)/, %.pm) : $(addprefix $(SRC)/, %.pms)
	@echo "### Finishing PERL module from \"$<\" -> \"$@\"" ;
	@$(CP) -u $< $@ ;
#	@( echo "#!$(PERL)"; cat $< ) | $(CPIF) $@;

isexec : $(BINCODE)
	@$(CHM) 755 $^ ;

chkmodperl :
	@perl -e 'print STDERR "###\n### CHECKING PERL MODULES...\n###\n"; \
              foreach $$module (@ARGV) { \
	              eval "require $$module" || do { \
	                  print STDERR "\t### MODULE $$module NOT FOUND !!!\n"; \
	                  next; \
	              }; \
	              ($$mod = $$module) =~ s%::%/%og; \
	              $$mod .= ".pm"; \
	              open(R,"locate $$mod |"); \
                  @mod = <R>; \
                  close(R); \
	              print STDERR "\t### MODULE $$module FOUND AT:\n\t". \
                               join("\t", @mod)."\n"; \
	          }; \
	          exit(0);' $(REQUIRED_PERL_MODULES)

#'
### Making man pages

mkman : $(PERLMODULES)
	@echo "### Manifying \"$<\"...";
	cd $(BIN); \
	 pod2man -quotes "none" $(notdir $<) | \
	   $(CPIF) $(MAN)/$(addsuffix .3, $(basename $(notdir $<)));

#
### TANGLING DATA #####################

datasets : datasetsmsg getdataset

datasetsmsg :
	@echo "###";
	@echo "### TANGLING DATASETS...";
	@echo "###";

getdataset : $(DATASETS)

$(BASHVARS) : $(addprefix $(NOWEB)/, $(REPORT))
	@$(echo_tangle)
	( $(TANGLE) -R$(BASHVARS_R) $< ; \
	  echo 'BASH=$(BASH) ;'; \
	  echo 'GAWK=$(GAWK) ;'; \
	  echo 'PERL=$(PERL) ;'; \
	  echo 'export BASH GAWK PERL ;'; \
	  echo '#' ) | $(CPIF) $@ ;

#
### INFO RULES

info : header getinfo trailer

getinfo :
	@echo "### BASH ###  $(BASH)  ###  $(BASHV)";
	@echo "### GAWK ###  $(GAWK)  ###  $(GAWKV)";
	@echo "### PERL ###  $(PERL)  ###  $(PERLV)";

header :
	@echo "###";
	@echo "### RUNNING MAKEFILE";
	@echo "###";
	@echo "### $(DATE) -- $(USER)";
	@echo "###";

trailer :
	@echo "###";
	@echo "### MAKEFILE DONE...";
	@echo "###";

################################################################################
# MAKEFILE = $(WORK)/Makefile
# #
# all: $(NOWEB)
#
# $(MAKEFILE) : $(NOWEB)
#         no tangle -R'Makefile' $< > $@~; \
#     DIF=`diff -bq $@ $@~`; \
#     if [ -z "$DIF" ]; \
#       then \
#         echo "### NO DIFERENCES IN MAKE !!!"; \
#       else \
#         echo "### UPDATING MAKEFILE !!!"; \
#         mv -v $@~ $@; \
#         touch $(NOWEB); \
#         make all; \
#         exit; \
#           fi
