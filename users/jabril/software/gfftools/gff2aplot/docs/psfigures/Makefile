#
# BUILDING PSFIGURES for GFF2APLOT DOCUMENTATION
#
# $Id: Makefile,v 1.1 2002-12-24 13:07:37 jabril Exp $
#

#
# Using SUFFIX RULES to build PS files from LATEX 
#
.SUFFIXES: .epsi .eps .ps .tmp.ps .dvi .tex

# T='flowchart_GFFrecordparsing';
.tex.dvi:
	latex $<

.dvi.tmp.ps:
	dvips $< -E -o $*.tmp.ps

.tmp.ps.ps:
	ps2ps   $< $*.ps
.tmp.ps.eps:
	eps2eps $< $*.eps
.tmp.ps.epsi:
	ps2epsi $< $*.epsi

#
# Building PostScript files

all: fc_GFF fc_VAR

fc_GFF:	flowchart_GFFrecordparsing.ps
fc_VAR:	flowchart_settingvariables.ps

clean:
	/bin/rm -v *.dvi *.aux *.log
