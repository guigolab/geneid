#
# PREPARING FILES for GFF2APLOT DOWNLOADS
#
# $Id: Makefile,v 1.1 2003-03-07 20:39:06 jabril Exp $
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
 LTAR = $(BINDIR)/tar $(V) -ztf
 MTAR = $(BINDIR)/tar $(V) -zcf
 GZIP = $(BINDIR)/gzip $(V)
  SSH = $(LOCALBIN)/ssh2 -oProtocol=2
  SCP = $(LOCALBIN)/scp2 -oProtocol=2

#
### MAKE VARS

.PHONY    : clean forceall
.DEFAULT  : all
.PRECIOUS : .pl .ps .tgz .gz
.SUFFIXES :

#
### VARS

 APLOT = gff2aplot.pl
PBLAST = parseblast.pl
  PSIM = sim2gff.pl
  PMUM = ali2gff.tgz
 PBLAT = blat2gff.pl
 PROGS = $(APLOT) $(PBLAST) $(PSIM) $(PMUM) $(PBLAT)
  DOCS = GFF2APLOT_MANUAL.ps 
DOCSGZ = $(addsuffix .gz, $(DOCS))

#
# PATHS for ORI FILES

 GFFTOOLS = ../..
 APLOTDIR = $(GFFTOOLS)/gff2aplot
PBLASTDIR = $(GFFTOOLS)/parseblast/bin
  PSIMDIR = $(GFFTOOLS)/sim2gff/bin
  PMUMDIR = $(GFFTOOLS)/sim2gff/old
 PBLATDIR = $(GFFTOOLS)/blat2gff/bin

 APLOTBIN = $(APLOTDIR)/bin
  DOCSDIR = $(APLOTDIR)/docs

 VERSFILE = $(APLOTDIR)/current_version
  TARBALL = $(shell sed 's/ /\-/; s/\./\_/;' $(VERSFILE))
  TARFILE = $(addsuffix .tgz, $(TARBALL))
     LDIR = $(ODIR)/$(TARBALL)

  FTPUSER = gmaster
  FTPHOST = monstre1.imim.es
   FTPDIR = /home/ftp/pub/software/gff_tools/gff2aplot
  FTPLDIR = $(FTPDIR)/last
  PUBFTP := $(FTPUSER)@$(FTPHOST):$(FTPDIR)
   FTPCMD = "( [ -d $(FTPLDIR) ] && rm -f $(FTPLDIR)/* ) || mkdir $(FTPLDIR)"
#
### DEPENDENCIES

all :
	@echo "###" ;
	@echo "### BUILDING $(TARBALL) DISTRIBUTION FILES..." ;
	@echo "###" ;
	@$(MAKE) tarball ftp ;

forceall : clean

clean : 
	-$(RM) -f $(addprefix $(LDIR)/, $(PROGS)) \
              $(addprefix $(LDIR)/, $(DOCS) $(DOCSGZ)) \
              $(addprefix $(ODIR)/, $(TARFILE)) ;
	-$(RD)    $(LDIR) ;

#
#########################

tarball : $(ODIR)/$(TARFILE)

ftp : ftp_tarball ftp_last

ftp_tarball : $(ODIR)/$(TARFILE)
	@$(SCP) $(ODIR)/$(TARFILE) $(PUBFTP)/ ;

ftp_last : $(addprefix $(LDIR)/, $(PROGS) $(DOCSGZ))
	$(SSH) -l $(FTPUSER) $(FTPHOST) $(FTPCMD) ;
	@$(SCP) $(addprefix $(LDIR)/, $(PROGS) $(DOCSGZ)) $(PUBFTP)/last/ ;

$(ODIR)/$(TARFILE) : $(LDIR) $(addprefix $(LDIR)/, $(PROGS) $(DOCSGZ))
	@cd $(ODIR) ; pwd ; \
	 $(MTAR) $(TARFILE) $(addprefix $(TARBALL)/, $(PROGS) $(DOCSGZ)); \
	 $(CHM) 664 $(TARFILE) ; \
	 $(CHO) 500:20 $(TARFILE) ; \
	 cd .. ;

$(LDIR) :
	@$(MD) $(LDIR) ; \
	 $(CHM) 775 $(LDIR) ; \
	 $(CHO) 500:20 $(LDIR) ;

$(addprefix $(LDIR)/,$(APLOT)) : $(addprefix $(APLOTBIN)/,$(APLOT))
	@$(CP) $< $@ ; \
	 $(CHM) 775 $@ ; \
	 $(CHO) 500:20 $@ ; \
	 $(PERL) -c $@ ;
	-@$@ --version ;

$(addprefix $(LDIR)/,$(PBLAST)) : $(addprefix $(PBLASTDIR)/,$(PBLAST))
	@$(CP) $< $@ ; \
	 $(CHM) 775 $@ ; \
	 $(CHO) 500:20 $@ ; \
	 $(PERL) -c $@ ;
	-@$@ --version ;

$(addprefix $(LDIR)/,$(PSIM)) : $(addprefix $(PSIMDIR)/,$(PSIM))
	@$(CP) $< $@ ; \
	 $(CHM) 775 $@ ; \
	 $(CHO) 500:20 $@ ; \
	 $(PERL) -c $@ ;
	# -@$@ --version ;

$(addprefix $(LDIR)/,$(PMUM)) : $(addprefix $(PMUMDIR)/,$(PMUM))
	@$(CP) $< $@ ; \
	 $(CHM) 775 $@ ; \
	 $(CHO) 500:20 $@ ; \
	 $(LTAR) $@ ;

$(addprefix $(LDIR)/,$(PBLAT)) : $(addprefix $(PBLATDIR)/,$(PBLAT))
	@$(CP) $< $@ ; \
	 $(CHM) 775 $@ ; \
	 $(CHO) 500:20 $@ ; \
	 $(PERL) -c $@ ;
	# -@$@ --version ;

$(addprefix $(LDIR)/,$(DOCS)) : $(addprefix $(DOCSDIR)/,$(DOCS))
	@$(CP) $< $@ ; \
	 $(CHM) 664 $@ ; \
	 $(CHO) 500:20 $@ ;

$(addprefix $(LDIR)/,$(DOCSGZ)) : $(addprefix $(LDIR)/,$(DOCS))
	@$(GZIP) -c $< > $@ ; \
	 $(CHM) 664 $@ ; \
	 $(CHO) 500:20 $@ ;


### EOM ###