# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                   INSTALLING SGP2 FILES
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#               The following Makefile was developed
#     on a Red Hat Linux box running on GNU Make version 3.79.1 
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     Copyright (C) 1999-2003 - Josep Francesc ABRIL FERRANDO
#                                        Genis PARRA FARRE
#                                      Roderic GUIGO SERRA
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  INSTALLDIR is the only variable that should be modified
#    (unless you were not using GNUmake, of course, because in such
#     case you must fix all the GNUmake specific commands... ;^D).
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
### PATHS

       BBIN = /bin
       UBIN = /usr/bin
       LBIN = /usr/local/bin
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #  REMEMBER: INSTALLDIR is the only variable that should be modified
 INSTALLDIR = $(LBIN)
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#
### MAKE VARS

.PHONY    : clean cleanall forceall
.DEFAULT  : all
.PRECIOUS : %.sh %.awk %.pl %.pm %.ps %.tgz %.gz \
            %.rc %.gff %.gb %.fa README
.SECONDARY :
.SUFFIXES :

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
# You can update the following paths to your system settings

         CP = $(firstword $(shell which cp   )) 
         CHM = $(firstword $(shell which chmod))
#
### LOCAL PATHS

       WORK = .
        SRC = $(WORK)/src
        BIN = $(WORK)/bin

#
### LOCAL VARS

     BINEXT = .sh .awk .pl .pm

       SGP2 = sgp2
 PARSEBLAST = parseblast
     GENEID = geneid
     HSP2SR = blast2gff


     GENEIDD = $(SRC)/geneid_v1.1-sgp
     HSP2SRD = $(SRC)/blast2gff

   GENEIDSRC = $(wildcard $(addprefix $(GENEIDD)/src/, *.c)) \
               $(wildcard $(addprefix $(GENEIDD)/include/, *.h))
   HSP2SRSRC = $(wildcard $(addprefix $(HSP2SRD)/src/, *.c)) \
               $(wildcard $(addprefix $(HSP2SRD)/include/, *.h))


    SRCCODE = $(GENEID)-sgp  $(PARSEBLAST) \
              $(HSP2SR) $(SGP2)
    BINCODE = $(addprefix $(BIN)/, $(SRCCODE))

#
### MAKE RULES

all: header main trailer

test :
	@echo "### SORRY !!! Tests were not implemented yet...";

install : header installbin trailer

forceall : header cleanbin main trailer

clean : header cleanbin trailer

cleanall : header cleanbin trailer

cleanbin :
	-@$(RM) -f $(wildcard $(BIN)/*);
	rmdir $(BIN);
	@echo "### Cleaning GENEID..."
	@cd $(GENEIDD); \
	 $(MAKE) clean;
	@echo "### Cleaning HSP2SR..."
	@cd $(HSP2SRD); \
	 $(MAKE) clean;


main : $(BIN) $(BINCODE)
	@$(CHM) 755 $^;

installbin : $(BINCODE)
	@echo "### COPYING BIN FILES TO $(INSTALLDIR)...";
	@$(CP) -p $^ $(INSTALLDIR)/;

#
### FINISHING CODE

$(BIN) :
	mkdir $(BIN);

$(addprefix $(BIN)/, $(GENEID)-sgp) : $(addprefix $(GENEIDD)/bin/, $(GENEID))
	@${CP} $(addprefix $(GENEIDD)/bin/, $(GENEID)) \
           $(addprefix $(BIN)/, $(GENEID)-sgp)

$(addprefix $(GENEIDD)/bin/, $(GENEID)) : $(GENEIDSRC)
	@echo "### Making GENEID..."
	@cd $(GENEIDD); \
	 $(MAKE) -f ./Makefile;

$(addprefix $(BIN)/, $(HSP2SR)) : $(addprefix $(HSP2SRD)/bin/, $(HSP2SR))
	@${CP} $(addprefix $(HSP2SRD)/bin/, $(HSP2SR)) \
           $(addprefix $(BIN)/, $(HSP2SR))

$(addprefix $(HSP2SRD)/bin/, $(HSP2SR)) : $(HSP2SRSRC)
	@echo "### Making HSP2SR..."
	@cd $(HSP2SRD); \
	 $(MAKE) -f ./Makefile;

$(addprefix $(BIN)/, $(SGP2)) : $(addprefix $(SRC)/, $(SGP2).pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, $(PARSEBLAST)) : $(addprefix $(SRC)/, $(PARSEBLAST).pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

isexec : $(BINCODE)
	@$(CHM) 755 $<;

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
