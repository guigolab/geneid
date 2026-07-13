For more information on how to install, train and run geneid check the [Wiki](https://github.com/guigolab/geneid/wiki)

******************** geneid v1.6 README File ********************

Summary:
A. What's geneid ?
B. Installing geneid
C. File Listing
D. Compiling geneid
E. To run geneid
F. Authors and help

***************************************

A. What's geneid ?
------------------

geneid is a program to predict genic elements as splice sites, exons 
or genes, along eukaryotic DNA sequences. geneid offers also some 
rudimentary support to integrate predictions for multiple sources. 
The program is written in ANSI C, and runs on UNIX operating systems.

Installation, setup and usage of geneid is very easy, and there is a 
range of options to configure output predictions and program behaviour.

geneid source code, compiled binaries and documentation are available 
under the GNU GENERAL PUBLIC LICENSE.

Comments and questions are welcome.      

***************************************

B. Installing geneid
--------------------

The geneid distribution contains several directories and files. Source 
code, compiled binaries for a number of architectures and operating 
systems, and documentation files are included in the distribution.

The distribution is archived and compressed in a single file using the
command tar -zcvf. The compressed file name is geneid.tar.gz (or something
similar depending on compiled binaries included). The geneid files can 
be extracted following these instructions:

** UNIX:

Type:
>gzip -d geneid.tar.gz
>tar -xvf geneid.tar

** LINUX:

Type: 
>tar -zxvf geneid.tar.gz

After executing these commands, the directory geneid will be created 
in your working directory. 

***************************************

C. File Listing
---------------

The geneid distribution contains the following files and directories:

** bin/
compiled binaries

** docs/
documentation: a short handbook.

** include/
geneid.h: The geneid header file.

** objects/
object files after compiling the source code.

** param/
Parameter files for several species.

** samples/
Test sequences.

** src/
source code of geneid program.

** GNULicense
This software is registered under GNU license.

** Makefile
This file is required to build geneid binary file.

** README
This file.

***************************************

D. Compiling geneid
-------------------

Move into the geneid directory.

Type:
>make
to compile geneid.

This will generate the geneid executable file within the bin/ subdirectory.

geneid's working arrays now grow on demand, so memory follows the data and
there is nothing to size at compile time -- just "make" (the older
MEM=low|medium|high build profiles are gone). Resident memory scales with
the input and stays modest (a few GB even on a large, dense genome).

Dependencies:

  * zlib (libz) -- REQUIRED. Used by the indexed bigBed (-R) and bigWig (-S)
    evidence readers. It is preinstalled on virtually every UNIX/Linux/macOS
    system (package "zlib1g-dev" / "zlib-devel"); the Makefile links it with
    -lz automatically, so a plain "make" just works.

  * htslib -- OPTIONAL, only for reading indexed BAM evidence directly
    (-R / -S / -Y with a .bam; see section E). To enable it, build with:

    >make WITH_HTSLIB=1

    If htslib is not under /usr/local, point the build at its install prefix:

    >make WITH_HTSLIB=1 HTSLIB_PREFIX=/opt/homebrew        # e.g. macOS/Homebrew
    >make WITH_HTSLIB=1 HTSLIB_PREFIX=/path/to/htslib

    The build bakes HTSLIB_PREFIX/lib into the binary's runtime search path
    (-Wl,-rpath), so the resulting geneid finds libhts.so.* on its own -- no
    "module load" or LD_LIBRARY_PATH needed at run time. (If you ever move the
    htslib install, either rebuild or set
    LD_LIBRARY_PATH=<prefix>/lib:$LD_LIBRARY_PATH.)

    The default build (plain "make") never references htslib; BAM input just
    isn't available in that binary. Install htslib from
    https://github.com/samtools/htslib (or "brew install htslib",
    "apt-get install libhts-dev", a cluster module, etc.).

Type:
>geneid -h

to test the binary file has been correctly created.

***************************************

E. To run geneid
-----------------

To run geneid type:
>geneid -P <parameter_filename> <Sequence_filename>.

Alternatively you can set the parameter file using the environment
variable GENEID.

For example,
>bin/geneid -vP param/human3iso.param samples/example1.fa

Run "geneid -h" for the full option list. The main options for integrating
external evidence (annotations, protein homology and RNA-seq) are:

  -R <file>   Annotations that guide the prediction (e.g. exon/intron features).
              Accepts a GFF file, a bigBed of the same records, or -- in a
              WITH_HTSLIB build -- an indexed BAM, whose spliced-read (CIGAR N)
              junctions become intron evidence. For all three, only the evidence
              overlapping each split is read (per-fragment range queries).
              NOTE: a text GFF must be sorted by acceptor position; the bigBed
              and BAM readers sort internally.

  -S <file>   RNA-seq / homology signal that scores exons. Accepts:
                - a text GFF (protein-alignment HSPs, or expression coverage);
                - bigWig coverage: "-S plus.bw,minus.bw" (stranded) or
                  "-S cov.bw" (unstranded);
                - (WITH_HTSLIB) an indexed BAM: "-S reads.bam" -- per-base
                  read depth is computed per split.
              Coverage scores exons with or without -u.

  -u          Turn on UTR prediction (uses the -S coverage to place UTR ends).

  -y <mode>   Library strandedness for BAM -S coverage: "rf" (dUTP / reverse,
              the common Illumina protocol), "fr" (forward), or "none"
              (unstranded, the default).

  -Y <reads.bam>
              (WITH_HTSLIB) Convenience: use ONE indexed BAM as BOTH intron
              evidence (-R junctions) AND RNA-seq coverage (-S), instead of
              passing the same file to -R and -S separately.

BAM inputs must be coordinate-sorted and indexed (a .bai/.csi alongside).
For BAM introns, the junction strand is taken from the read's XS tag, else the
minimap2 ts tag, else the GT-AG/CT-AC splice motif.

For example, RNA-seq-guided prediction from a single STAR/minimap2 BAM
(htslib build):
>bin/geneid -3U -u -Y aligned.bam -P param/human.rnaseq.param genome.fa

***************************************

F. Authors and help
-------------------

geneid has been written by Enrique Blanco, Tyler Alioto and Roderic Guigó (CRG).
Parameter files have been generated by Genis Parra, Francisco Camara and Tyler Alioto.

geneid home page is "http://genome.crg.es/software/geneid" and 
geneid distributions can be obtained via github at "https://github.com/guigolab/geneid".

If you need help, send a message to "geneid@crg.es".

