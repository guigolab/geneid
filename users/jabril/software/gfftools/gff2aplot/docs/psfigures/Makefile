#
# BUILDING PSFIGURES for GFF2APLOT DOCUMENTATION
#
# $Id: Makefile,v 1.2 2003-02-05 16:15:01 jabril Exp $
#
# RUNNING on GNU Make version 3.79.1
#
### PATHS

BINDIR = /bin
USRBIN = /usr/bin
LOCALBIN = /usr/local/bin

#
### PROGRAMS

RM = ${BINDIR}/rm -v
LN = ${BINDIR}/ln -vfs
SLEEP = ${BINDIR}/sleep
TOUCH = ${BINDIR}/touch
LATEX = ${USRBIN}/latex
BIBTEX = ${USRBIN}/bibtex
DVIPS = ${USRBIN}/dvips
GS = ${LOCALBIN}/gs
EPSEPS = ${LOCALBIN}/eps2eps

#
### VARS

DEBUG =  # debug flag is off by default

V_fc_GFF = flowchart_GFFrecordparsing.ps
V_fc_VAR = flowchart_settingvariables.ps
V_flowcharts = ${V_fc_GFF} ${V_fc_VAR}
V_bd_flt = blockdiagram_filter.ps
V_bd_int = blockdiagram_interactive.ps
V_blocks = ${V_bd_flt} ${V_bd_int}
V_makepsfigs = ${V_flowcharts} ${V_blocks}

#
### MAKE VARS

.PHONY : clean forceall
.DEFAULT: all
.SECONDARY: ${V_makepsfigs:.ps=.dvi}    \
            ${V_makepsfigs:.ps=.tmp.ps} \
            ${V_makepsfigs:.ps=.eps}
# Now not using SUFFIX RULES to build PS files from LATEX 
.SUFFIXES :
#.SUFFIXES : .epsi .ps .eps .tmp.ps .dvi .tex

#
### DEPENDENCIES

# ${V_makepsfigs:.ps=.tex} : showgrid.tex
#	${SLEEP} 1s; ${TOUCH} $@

%.dvi : %.tex showgrid.tex
	${LATEX} $* ;                           \
    if [ "`grep undefined $*.log`" != "" ]; \
      then                                  \
        ${BIBTEX} $* ;                      \
        ${LATEX}  $* ;                      \
      fi ;                                  \
    if [ "`grep Rerun $*.log`" != "" ];     \
      then                                  \
        ${LATEX}  $* ;                      \
      fi
#	${TOUCH} $*.tex

%.tmp.ps : %.dvi
	${DVIPS} $< -o $*.tmp.ps -D 600 -Z # -E -O 0cm,0cm -T 21cm,29.7cm

%.eps : %.tmp.ps
	${GS} -dBATCH -dNOPAUSE -r600 -sPAPERSIZE=a4 \
          -sDEVICE=epswrite -sOutputFile=$*.eps $<

%.ps : %.eps
	${EPSEPS} $< $*.ps

# %.epsi : %.tmp.ps
#	ps2epsi $< $*.epsi

#
### Building PostScript files

all : makepsfigs

forceall :
	${TOUCH} ${V_makepsfigs:.ps=.tex}
	${MAKE} all

debug :
	if [ -z "${DEBUG}" ];                              \
      then                                             \
        if [ "`diff -q ./showgrid_off.tex ./showgrid.tex`" != "" ]; \
          then                                         \
            ${LN} showgrid_off.tex showgrid.tex ;      \
          else                                         \
            echo "## showgrid.tex unchanged...";       \
          fi;                                          \
      else                                             \
        if [ "`diff -q ./showgrid_on.tex ./showgrid.tex`" != "" ]; \
          then                                         \
            ${LN} showgrid_on.tex  showgrid.tex ;      \
           else                                         \
            echo "## showgrid.tex unchanged...";       \
          fi;                                          \
      fi
# "diff -q" reports whether the files differ, not the details of the diffs

fc_GFF : debug ${V_fc_GFF}

fc_VAR : debug ${V_fc_VAR}

flowcharts : debug ${V_flowcharts}

bd_flt : debug ${V_bd_flt}

bd_int : debug ${V_bd_int}

blocks : debug ${V_blocks}

makepsfigs : debug ${V_makepsfigs}

#
# Other targets

clean :
	${RM} *~ ${V_makepsfigs:.ps=.dvi} \
             ${V_makepsfigs:.ps=.aux} \
             ${V_makepsfigs:.ps=.log} ;

#### EOM ####