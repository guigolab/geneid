#
# PREPARING FILES for GFF2APLOT DOWNLOADS
#
# $Id: Makefile,v 1.2 2003-03-17 18:02:11 jabril Exp $
#
# RUNNING on GNU Make version 3.79.1
#
### PATHS

  BINDIR = /bin
  USRBIN = /usr/bin
LOCALBIN = /usr/local/bin
    WORK = .
    ODIR = $(WORK)/ftp

#
### PROGRAMS

    V = -v
   RM = $(BINDIR)/rm $(V)
   RD = $(BINDIR)/rmdir $(V)
   CP = $(BINDIR)/cp $(V)
   LN = $(BINDIR)/ln $(V)
   MD = $(BINDIR)/mkdir $(V)
  CHM = $(BINDIR)/chmod $(V)
  CHO = $(BINDIR)/chown $(V)
 PERL = $(USRBIN)/perl
 TTAR = $(BINDIR)/tar $(V) -ztf
 CTAR = $(BINDIR)/tar $(V) -zcf
 XTAR = $(BINDIR)/tar $(V) -zxf
 GZIP = $(BINDIR)/gzip $(V)
  SSH = $(LOCALBIN)/ssh2 -oProtocol=2
  SCP = $(LOCALBIN)/scp2 -oProtocol=2

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
### MAKE VARS

.PHONY    : clean forceall
.DEFAULT  : all
.PRECIOUS : .pl .ps .tgz .gz \
            .rc .gff .gb .fa
.SUFFIXES :

#
### VARS

  APLOT = gff2aplot.pl
 PBLAST = parseblast.pl
   PSIM = sim2gff.pl
   PMUM = ali2gff
  PBLAT = blat2gff.pl
  PROGS = $(APLOT) $(PBLAST) $(PSIM) $(PMUM) $(PBLAT)
   DOCS = GFF2APLOT_MANUAL.ps 
 DOCSGZ = $(addsuffix .gz, $(DOCS))

LICENSE = GNU-GPL
  RDOCS = INSTALL README
 MKFILE = dist_Makefile

#
# PATHS for ORI FILES

 GFFTOOLS = ../..
 APLOTDIR = $(GFFTOOLS)/gff2aplot
PBLASTDIR = $(GFFTOOLS)/parseblast/src
  PSIMDIR = $(GFFTOOLS)/sim2gff/src
  PMUMDIR = $(GFFTOOLS)/sim2gff/old
 PBLATDIR = $(GFFTOOLS)/blat2gff/src

 APLOTBIN = $(APLOTDIR)/src
  DOCSDIR = $(APLOTDIR)/docs

 VERSFILE = $(APLOTDIR)/current_version
  TARBALL = $(shell sed 's/ /\-/; s/\./\_/;' $(VERSFILE))
  TARFILE = $(addsuffix .tgz, $(TARBALL))
     LDIR = $(ODIR)/$(TARBALL)
     SDIR = $(LDIR)/src
     DDIR = $(LDIR)/doc

   EXPATH = $(DOCSDIR)/html/gff2aplot_docs
      EX1 = wublast
      EX2 = mhcregion
      EX3 = splice_sites
   EXDIRS = $(EX1) $(EX2) $(EX3)
  EXFILES = README *.rc *.gff *.gb *.fa
    EXDIR = examples
     XDIR = $(LDIR)/$(EXDIR)
   FILES1 = $(foreach file, \
              $(addprefix $(addprefix $(EXPATH)/, $(EX1))/,$(EXFILES)), \
              $(wildcard $(file)))
   FILES2 = $(foreach file, \
              $(addprefix $(addprefix $(EXPATH)/, $(EX2))/,$(EXFILES)), \
              $(wildcard $(file)))
   FILES3 = $(foreach file, \
              $(addprefix $(addprefix $(EXPATH)/, $(EX3))/,$(EXFILES)), \
              $(wildcard $(file)))

 AUXFILES = GNU-GPL README INSTALL

  FTPUSER = gmaster
  FTPHOST = monstre1.imim.es
   FTPDIR = /home/ftp/pub/software/gff_tools/gff2aplot
  FTPLDIR = $(FTPDIR)/last
  PUBFTP := $(FTPUSER)@$(FTPHOST):$(FTPDIR)
   FTPCMD = "( [ -d $(FTPLDIR) ] && rm -rf $(FTPLDIR)/* ) || mkdir $(FTPLDIR)"

#
### DEPENDENCIES

all : header main trailer

header :
	@echo "###" ;
	@echo "### BUILDING $(TARBALL) DISTRIBUTION FILES..." ;
	@echo "###" ;
	@echo "### $(DATE) -- $(USER)";
	@echo "###";

trailer :
	@echo "###";
	@echo "### MAKEFILE DONE...";
	@echo "###";

main : tarball ftp

forceall : clean

clean : 
	-$(RM) -rf $(LDIR) $(addprefix $(ODIR)/, $(TARFILE)) ;
#	-$(RM) -rf $(addprefix $(LDIR)/, $(PROGS)) \
#               $(addprefix $(LDIR)/, $(DOCS) $(DOCSGZ)) \
#               $(addprefix $(ODIR)/, $(TARFILE)) ;
#	-$(RD)     $(LDIR) ;

#
#########################

tarball : $(ODIR)/$(TARFILE)

ftp : ftp_tarball ftp_last

ftp_tarball : tarballmsg $(ODIR)/$(TARFILE)
	@echo "###";
	@echo "### COPYING TARBALL TO REMOTE FTP SERVER...";
	@echo "###";
	@$(SCP) $(ODIR)/$(TARFILE) $(PUBFTP)/ ;

tarballmsg : 
	@echo "###";
	@echo "### BUILDING DISTRIBUTION TARBALL...";
	@echo "###";

ftp_last : $(addprefix $(SDIR)/, $(PROGS)) $(addprefix $(DDIR)/, $(DOCSGZ))
	@echo "###";
	@echo "### COPYING LATEST VERSIONS TO REMOTE FTP SERVER...";
	@echo "###";
	$(SSH) -l $(FTPUSER) $(FTPHOST) $(FTPCMD) ;
	@$(SCP) -r $^ $(PUBFTP)/last/ ;

$(ODIR)/$(TARFILE) : mkdirs \
                     $(addprefix $(SDIR)/, $(PROGS)) \
                     $(addprefix $(DDIR)/, $(DOCSGZ)) \
                     $(addprefix $(LDIR)/, $(EXDIR)) \
                     auxfiles examples ownership
	@cd $(ODIR) ; pwd ; \
	 $(CTAR) $(TARFILE) $(addprefix $(TARBALL)/, \
	                      $(AUXFILES) Makefile bin \
                          $(addprefix src/, $(PROGS)) \
                          $(addprefix doc/, $(DOCSGZ)) $(EXDIR)); \
	 $(CHM) 664 $(TARFILE) ; \
	 $(CHO) 500:20 $(TARFILE) ; \
	 cd .. ;

auxfiles : $(addprefix $(LDIR)/, $(AUXFILES) Makefile)

$(addprefix $(LDIR)/, GNU-GPL) : $(addprefix $(WORK)/, GNU-GPL)
	@$(CP) $< $@;

define UpdateVersion
	$(PERL) -ne 'BEGIN{ $$f = shift @ARGV; \
	                    open(T,"< $$f"); \
	                    $$t = <T>; chomp($$t); \
	                    ($$prog,$$ver) = split /\s+/o,$$t; \
	                    close(T); \
	                    $$Prog = uc($$prog); \
	              }; \
	              s/\#\#\#PROGRAM\#\#\#/$$Prog/; \
	              s/\#\#\#program\#\#\#/$$prog/; \
	              s/\#\#\#VERSION\#\#\#/$$ver/; \
	              print STDOUT $$_; \
	            ' $(VERSFILE) $(<) > $(@);
endef

$(addprefix $(LDIR)/, README INSTALL) : $(VERSFILE)

$(addprefix $(LDIR)/, README) : $(addprefix $(WORK)/, README)
	@$(UpdateVersion)

$(addprefix $(LDIR)/, INSTALL) : $(addprefix $(WORK)/, INSTALL)
	@$(UpdateVersion)

$(addprefix $(LDIR)/, Makefile) : $(addprefix $(WORK)/, dist_Makefile)
	@$(CP) $< $@;

ownership :
	@$(CHO) -R 500:20 $(LDIR) ;

mkdirs :
	-@$(MD) $(LDIR) ;
	-@$(MD) $(addprefix $(LDIR)/, src bin doc);
	@$(CHM) -R 775 $(LDIR) ;

examples: $(XDIR) $(addprefix $(XDIR)/, $(EXDIRS))

$(XDIR) :
	-@$(MD) $(XDIR); \
	  $(CHM) -R 775 $(XDIR) ;

$(addprefix $(XDIR)/, $(EX1)) : $(FILES1)
	-@$(MD) $@ ; $(CHM) -R 775 $@ ;
	@$(CP) $(FILES1) $@/; $(CHM) 664 $@/* ;

$(addprefix $(XDIR)/, $(EX2)) : $(FILES2)
	-@$(MD) $@ ; $(CHM) -R 775 $@ ;
	@$(CP) $(FILES2) $@/; $(CHM) 664 $@/* ;

$(addprefix $(XDIR)/, $(EX3)) : $(FILES3)
	-@$(MD) $@ ; $(CHM) 775 $@ ;
	@$(CP) $(FILES3) $@/; $(CHM) 664 $@/* ;

# $(XDIR): $(addprefix $(addsuffix /$(EXDIRS)/, $(EXPATH)), $(EXFILES))
# 	@$(MD) $(XDIR) ; \
# 	 $(CHM) 775 $(XDIR) ; \
# 	 $(CHO) 500:20 $(XDIR) ;

$(addprefix $(SDIR)/,$(APLOT)) : $(addprefix $(APLOTBIN)/,$(APLOT))
	@$(CP) $< $@ ; \
	 $(PERL) -c $@ ;
	-@$@ --version ;

$(addprefix $(SDIR)/,$(PBLAST)) : $(addprefix $(PBLASTDIR)/,$(PBLAST))
	@$(CP) $< $@ ; \
	 $(PERL) -c $@ ;
	-@$@ --version ;

$(addprefix $(SDIR)/,$(PSIM)) : $(addprefix $(PSIMDIR)/,$(PSIM))
	@$(CP) $< $@ ; \
	 $(PERL) -c $@ ;
	# -@$@ --version ;

$(addprefix $(SDIR)/, $(PMUM)) : $(addprefix $(PMUMDIR)/,$(PMUM))
	@$(CP) -r $< $@ ; \
	 $(CHM) 777 $@ ; \
	 $(CHM) 775 $@/* ; # \
#	 $(TTAR) $@ ;

$(addprefix $(SDIR)/, $(PBLAT)) : $(addprefix $(PBLATDIR)/, $(PBLAT))
	@$(CP) $< $@ ; \
	 $(PERL) -c $@ ;
	# -@$@ --version ;

$(addprefix $(DDIR)/, $(DOCS)) : $(addprefix $(DOCSDIR)/, $(DOCS))
	@$(CP) $< $@ ;
$(addprefix $(DDIR)/, $(DOCSGZ)) : $(addprefix $(DDIR)/, $(DOCS))
	@$(GZIP) -c $< > $@ ; $(CHM) 664 $@ ;

### EOM ###