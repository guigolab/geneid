#
PROG = gff2aplot
WORK = /home/ug/jabril/development/softjabril/gfftools/gff2aplot
MAKEFILE = $(WORK)/Makefile
NOWEB = $(WORK)/$(PROG).nw
#
all: $(NOWEB)

$(MAKEFILE) : $(NOWEB)
        no tangle -R'Makefile' $< > $@~; \
    DIF=`diff -bq $@ $@~`; \
    if [ -z "$DIF" ]; \
      then \
        echo "### NO DIFERENCES IN MAKE !!!"; \
      else \
        echo "### UPDATING MAKEFILE !!!"; \
        mv -v $@~ $@; \
        touch $(NOWEB); \
        make all; \
        exit; \
          fi
