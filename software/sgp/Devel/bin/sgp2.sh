#! /bin/sh 
#
# SGP-2, syntenic gene prediction
#      gene prediction by phylogenetic footprinting
#
# (c) roderic guigo, imim,  and thomas wiehe, max plank, 
#
# 1999-2000
#
# modified by jabril - 2000/08/10
#
# $Id: sgp2.sh,v 1.1 2000-10-04 18:07:06 jabril Exp $
#

# --------------------------------------------------
# -. defaults
# --------------------------------------------------
NF=0;                # number of input files
SGP2=$HOME/research/humus/SGP2-4/sgp2dir/bin       # binaries, scripts and params in SGP2
TMP=/tmp              # temporal files
SGPTMP=$TMP/sgp2_$$  # prefix for temporal files
trap 'rm -f $SGPTMP; exit 1' 0 1 2 3 9 15 


# tblastx defaults
BLASTPROGRAM=tblastx
PRESSDB=pressdb 
BLMATRIX=blosum62mod 
BLASTOPTIONS="-matrix $BLMATRIX -hspmax=10000 -nogap"


# geneid defaults
#GENEID=$SGP2/geneid.sun
GENEID=$SGP2/geneid

PARAM=$SGP2/human3iso.param

# parsing tblastx output defaults
S_CUTOFF=50
SCF=12 # substract to tblastx scores S_CUTOFF - SCF;
## ATTENTION!!! CHANGE THIS IN EXTRACT
DSC=`expr $S_CUTOFF - $SCF`
#SHSP=6         # shrink hsp by $SHSP
SHSP=0


# --------------------------------------------------
# 0. pre-process
# --------------------------------------------------


USAGE="sgp [-hv] [-o 'options'] [-g 'options] [-P filename] [-p filename] [-k filename] [-c value] [-s value] -1 seqfile_1 -2 seqfile_2\n\n 
-h \t\t produces this message\n 
-v \t\t verbose mode\n 
-g \t\t geneid options\n  
-o \t\t tblastx options\n  
-c value \t\t tblastx score cuttof\n  
-s value \t\t shrink hsp's by value\n  
-t filename \t\t read tblastx file\n  
    -f prefix \t\t read hsp gff files with in directory prefix and extension .hsp-rs\n  
-k prefix \t\t keep intermediate files with prefix\n  
-p filename \t ps output in filename file \n
-P filename\t\t geneid parameter file\n"

# process parameters
while getopts 1:2:a:g:f:i:k:p:t:fhlstv opt
do
    case $opt in
        1) SEQ1=$OPTARG;NF=`expr $NF + 1`;;     # seq file 1
        2) SEQ2=$OPTARG;NF=`expr $NF + 1`;;     # seq file 2
        o) BLASTOPTIONS=$OPTARG;;               # tblastx options
        g) GENEIDOPTIONS=$OPTARG;;              # geneid options
	P) GENEIDPARAM=$OPTARG;;                # geneid parameter file
        c) GENEIDOPTIONS=$OPTARG;;              # tblastx score cutoff
	s) GENEIDPARAM=$OPTARG;;                # shrink hsp's by
	k) OFN=$OPTARG;;                        # intermediate filename
	t) TBX=$OPTARG;;                        # read tblastx from file
	f) HSP=$OPTARG;;                        # read HSP files in directory
        p) PSO=$OPTARG;;                        # postscript output
        v) VRB=1;;                              # verbose
        \? | h) echo $USAGE; exit 2;;
    esac
done
shift `expr $OPTIND - 1`;


# check if enough parameters
if [ $NF -ne 2 ] 
then   
    echo "two input files required. USAGE:\n"  1>&2
    echo $USAGE  1>&2 ; exit 2;
fi

# check if input files exit
if [ ! -f $SEQ1 ]
then 
   echo $SEQ1: $ER_01 1>&2 ; exit 2;
fi


if [ ! -f $SEQ2 ]
then 
   echo $SEQ2: $ER_01 1>&2 ; exit 2;
fi


# preprocess input files
# base name
SEQ1_NAME=`basename $SEQ1`
SEQ2_NAME=`basename $SEQ2`

# tmpfiles prefix
SGPTMP1=$SGPTMP.$SEQ1_NAME
SGPTMP2=$SGPTMP.$SEQ2_NAME

# check if files in the correct format.
# to be done

# get locus names. Assuming input files in fasta format
LOC1=`gawk "NR==1{print substr(\\\$1,2); exit}"  $SEQ1`
LOC2=`gawk "NR==1{print substr(\\\$1,2); exit}"  $SEQ2`
  #####################
  ## Added by jabril ##
SGPTMPG="${SGPTMP}.${LOC1}_${LOC2}"  
  #####################

if [ $LOC1 = $LOC2 ]
then
  echo ERROR: sequence $SEQ1 and sequence $SEQ2 must have different locus names 1>&2
fi


# --------------------------------------------------
# PROGRAM
# --------------------------------------------------

# goto geneid if hsps provided

if [ $HSP ]
 then
  #######################
  ## Removed by jabril ##
#     cp $HSP$LOC1.hsp-rs $SGPTMP1.hsp-rs 
#     cp $HSP$LOC2.hsp-rs $SGPTMP2.hsp-rs 
#     cp $HSP$LOC1"_"$LOC2.aln $SGPTMP.aln
  #####################
  ## Added by jabril ##
cp $HSP${LOC1}_${LOC2}.srQ $SGPTMPG.srQ
cp $HSP${LOC1}_${LOC2}.srS $SGPTMPG.srS
cp $HSP${LOC1}_${LOC2}.aln $SGPTMPG.aln
  #####################
 else 

# --------------------------------------------------
# 1. run tblastx 
# --------------------------------------------------

if [ $TBX ]
then 
    # reading tblastx from file. 
    # WARINING, it assumes tblastx SEQ1 SEQ2
    # 
    cp $TBX  $SGPTMP.tbx
else  
    if [ $VRB ] 
    then
	echo running tblastx $SEQ1 $SEQ2 1>&2
    fi

    # temporal hack. not know how to do it better
    cp $SGP2/$BLMATRIX . # very dangerous
    # convert the first fasta file to database
    SEQ1DB=$SGPTMP.$SEQ1_NAME  #temporal name for database
    cp $SEQ1 $SEQ1DB
    $PRESSDB $SEQ1DB >/dev/null
    $BLASTPROGRAM $SEQ1DB $SEQ2 $BLASTOPTIONS > $SGPTMP.tbx

    rm $BLMATRIX
fi

if [ $OFN ]
then
    cp $SGPTMP.tbx $OFN$LOC1"_"$LOC2.tbx 1>&2
fi


# --------------------------------------------------
# 2. extract and process hsp's -- in development
# --------------------------------------------------

  #######################
  ## Removed by jabril ##
  #######################
# 
# if [ $VRB ] 
# then
#  echo processing tblastx output  1>&2
# fi
# 
# # to aplot
# $SGP2/blast2aplot scropt="bits" $SGPTMP.tbx > $SGPTMP.aln
# if [ $VRB ] 
# then
#  echo `wc -l $SGPTMP.aln | gawk '{print $1}'` hsp alignments   1>&2
# fi
# 
# 
# if [ $OFN ]
# then
#     cp $SGPTMP.aln $OFN$LOC1"_"$LOC2.aln
# fi
# 
# 
# # get the projected hsp's
# touch  $SGPTMP2.hsp  $SGPTMP1.hsp
# $SGP2/aln2gff $SGPTMP2.hsp $SGPTMP1.hsp $S_CUTOFF $SGPTMP.aln
# 
# if [ $VRB ] 
# then
#     echo `wc $SGPTMP1.hsp | gawk '{print $1}'` projected hsps over $S_CUTOFF in $SEQ1>&2
#     echo `wc $SGPTMP2.hsp | gawk '{print $1}'` projected hsps over $S_CUTOFF in $SEQ2>&2
# fi
# 
# if [ $OFN ]
# then
#     cp  $SGPTMP1.hsp $OFN$LOC1.hsp
#     cp  $SGPTMP2.hsp $OFN$LOC2.hsp
# fi
# 
# 
# # sort projected hsps by end position, cluster and shrink
# sort +4n $SGPTMP1.hsp | $SGP2/getregs | sort +3n | awk "{\$4+=$SHSP;\$5-=$SHSP; \$6+=-($S_CUTOFF - $SCF); print \$0;}" > $SGPTMP1.hsp-rs
# 
# sort +4n $SGPTMP2.hsp | $SGP2/getregs | sort +3n | awk "{\$4+=$SHSP;\$5-=$SHSP; print \$0;}" > $SGPTMP2.hsp-rs 
# 
# if [ $VRB ] 
# then
#     echo `wc $SGPTMP1.hsp-rs | gawk '{print $1}'` clustered hsps in $SEQ1 1>&2
#     echo `wc $SGPTMP2.hsp-rs | gawk '{print $1}'` clustered hsps in $SEQ2 1>&2
# fi
# 
# if [ $OFN ]
# then
#     cp  $SGPTMP1.hsp-rs $OFN$LOC1.hsp-rs
#     cp  $SGPTMP2.hsp-rs $OFN$LOC2.hsp-rs
# fi
#
  #####################
  ## Added by jabril, modified by rguigo ##
EXTRACT () {
  gawk 'BEGIN{
    DSC=38;
    SHSP=0;
    s1=ARGV[1]; ARGV[1]="";
    if (ARGC>3) { s2=ARGV[2]; ARGV[2]=""; scnd=1; }
    else scnd=0;
  }
  { gsub(/\"/,"",$10);
              print  $1, $2, $3,  $4+SHSP,  $5-SHSP, $6-DSC,  $7,  $8 > s1;
    if (scnd) print $10, $2, $3, $11+SHSP, $12-SHSP, $6-DSC, $16, $18 > s2;
    }' "$@"
}

MYVRB=""
if [ $VRB ] 
 then
    MYVRB="-v"
 fi

$SGP2/GetSRsAln $MYVRB -bC $S_CUTOFF -H $SGPTMPG -W $SGPTMPG -- $SGPTMP.tbx

if [ $VRB ] 
 then
   echo "Processing HSP/SR files from BLAST" 1>&2
 fi

EXTRACT $SGPTMP1.hsp-rs $SGPTMPG.srS   # I suposse that $SGPTMP1 is Subject
EXTRACT $SGPTMP2.hsp-rs $SGPTMPG.srQ   # and $SGPTMP2 is query, else modify this lines.
EXTRACT $SGPTMP2.hsp $SGPTMP1.hsp $SGPTMPG.hsp

gawk '{ 
  gsub(/\"/,"",$10); 
  print $1":"$10, $2, "hsp:hsp", $4":"$11, $5":"$12, $6, $7":"$16, $8":"$18 ;
}' $SGPTMPG.hsp > $SGPTMP.aln # Remember that GetSRsAln also removes HSP lower than $S_CUTOFF.
#
if [ $OFN ]
 then
     cp $SGPTMPG.srQ    $OFN$LOC1"_"$LOC2.srQ  # Those files have the relationships between
     cp $SGPTMPG.srS    $OFN$LOC1"_"$LOC2.srS  # Query and Subject coordinates so you can
     cp $SGPTMPG.hsp    $OFN$LOC1"_"$LOC2.hsp  # recover the sequence substring for both.
     cp $SGPTMP1.hsp-rs $OFN$LOC1.hsp-rs
     cp $SGPTMP2.hsp-rs $OFN$LOC2.hsp-rs
     cp $SGPTMP1.hsp    $OFN$LOC1.hsp
     cp $SGPTMP2.hsp    $OFN$LOC2.hsp
     cp $SGPTMP.aln     $OFN$LOC1"_"$LOC2.aln
 fi
  ## Added by jabril ##
  #####################

 fi

# --------------------------------------------------
# 3. run geneid
# --------------------------------------------------

if [ $VRB ] 
then
 echo running geneid in $SEQ1 with tblastx output 1>&2
fi
#$GENEID -xXP $PARAM -S  $SGPTMP1.hsp-rs $SEQ1 >  $SGPTMP1.sgp
$GENEID -GP $PARAM -S  $SGPTMP1.hsp-rs $SEQ1 >  $SGPTMP1.sgp

if [ $VRB ] 
then
 echo running geneid in $SEQ2 with tblastx output 1>&2
fi
#$GENEID -xXP $PARAM -S  $SGPTMP2.hsp-rs $SEQ2 >  $SGPTMP2.sgp
$GENEID -GP $PARAM -S  $SGPTMP2.hsp-rs $SEQ2 >  $SGPTMP2.sgp 

if [ $OFN ]
then
    cp $SGPTMP1.sgp $OFN$LOC1.sgp
    cp $SGPTMP2.sgp $OFN$LOC2.sgp
fi

cat $SGPTMP1.sgp $SGPTMP2.sgp


# --------------------------------------------------
# 4. graphical output
# --------------------------------------------------


if [ $PSO ]
then
 echo creating postscript output  1>&2
# 9.1 predicted genes separately -- gff2ps
$SGP2/gff2ps $SGPTMP1.sgp $SGPTMP1.hsp-rs > $PSO$LOC1.ps  2> /dev/null
$SGP2/gff2ps $SGPTMP2.sgp $SGPTMP2.hsp-rs > $PSO$LOC2.ps  2> /dev/null

# 9.2 predicted genes in aplot

echo $LOC1:$LOC2 SOURCE seqbounds 1:1 $LSEQ1:$LSEQ2 0  . 0 0 >  $SGPTMP.aplot
gawk 'BEGIN{strand["+"]=0;strand["-"]=1;}substr($1,1,1)!="#" && NF>0 {print $1, $2, "exon", $4, $5, 1.00, strand[$7],$8, "gene" $9 ":"}'  $SGPTMP1.sgp >> $SGPTMP.aplot

gawk 'BEGIN{strand["+"]=0;strand["-"]=1;}substr($1,1,1)!="#" && NF>0 {print $1, $2, "exon", $4, $5, 1.00, strand[$7],$8, "gene" $9 ":"}'  $SGPTMP2.sgp >> $SGPTMP.aplot
# 9.2 monotonic alignment
    gawk 'BEGIN{
  # defaults ----------------------------------------------------------
  # input file data structure
  # gff data structure
    seqname  = 1;
    source  = 2;
    feature  = 3;
    start    = 4;
    end      = 5;
    score    = 6;
    strand   = 7;
    frame    = 8;
    group    = 9;
}

substr($1,1,1)!="#" && NF>0 {
       split($start,start_,":");
       split($end,end_,":");
       split($strand,strand_,":");

       if (strand_[1]=="-") {
          aux=start_[1];
          start_[1]=end_[1];
          end_[1]=aux;
        }


       if (strand_[2]=="-") {
          aux=start_[2];
          start_[2]=end_[2];
          end_[2]=aux;
        }
# here score is divided by 1000 to use the %box in gff2aplot. 
# ideally, similarity instead of score should be use. 
#    print $seqname, $source, "align", start_[1] ":" start_[2],  end_[1] ":" end_[2], $score, ".", ".";
    print $seqname, $source, "align", start_[1] ":" start_[2],  end_[1] ":" end_[2], $score/2, ".", ".";
}' $SGPTMP.aln >> $SGPTMP.aplot


if [ $OFN ]
then
      cp $SGPTMP.aplot  $OFN$LOC1"_"$LOC2.aplot
fi

# 9.3 creating postscript
echo "$SGP2/gff2aplot -vP -O $SGP2/gff2aplot.sgp.rc $OFN$LOC1"_"$LOC2.aplot > $PSO$LOC1"_"$LOC2.ps" 1>&2;
    $SGP2/gff2aplot -vP -O $SGP2/gff2aplot.sgp.rc $SGPTMP.aplot > $PSO$LOC1"_"$LOC2.ps
fi

# removing files
rm -f $SGPTMP
exit 1









