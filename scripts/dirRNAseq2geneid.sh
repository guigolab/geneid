#!/usr/bin/env bash
b=`basename $1 .bam`;
GENOME=$2
BEDTOANN=`which bed12toAnnotation.awk`;
echo "Converting to bed...";
bedtools bamtobed -bed12 -i $1 | perl -ane '$F[3]=~/(\/| )([12]).*$/;if($2==1){$F[5]=$F[5]eq"-"?"+":"-";}print join("\t",@F),"\n";' | bedtools sort > $b.stranded.bed
echo "Extracting introns...";
cat $b.stranded.bed | awk -v type=intron -f $BEDTOANN | cut -f1-3,6 | sort | uniq -c | perl -ane 'print join("\t",($F[1],".","Intron",$F[2]+1,$F[3],$F[0],$F[4],".",".")),"\n";' > $b.introns.gff
echo "Creating bedgraph per strand and converting to gff...";
bedtools genomecov -g $GENOME -i $b.stranded.bed -split -bg -strand - | perl -ane 'print join("\t",($F[0],".",".",$F[1]+1,$F[2],$F[3],"-",".", ".")),"\n"' > $b.stranded.minus.gff
bedtools genomecov -g $GENOME -i $b.stranded.bed -split -bg -strand + | perl -ane  'print join("\t",($F[0],".",".",$F[1]+1,$F[2], $F[3],"+",".",".")), "\n"' > $b.stranded.plus.gff
echo "Combining results...";
if [ -x "$(command -v shuf)" ]; then
    cat $b.stranded.plus.gff $b.stranded.minus.gff | shuf >  $b.stranded.expression.shuffled.gff
elif [ -x "$(command -v gshuf)" ]; then
    cat $b.stranded.plus.gff $b.stranded.minus.gff | gshuf >  $b.stranded.expression.shuffled.gff
else
    cat $b.stranded.plus.gff $b.stranded.minus.gff >  $b.stranded.expression.gff
    echo "$b.stranded.expression.gff is not well shuffled. For efficient processing by geneid, it is recommended to use shuf (UNIX/LINUX) or gshuf (on MacOSX, run 'brew install coreutils' to install gshuf.)"
fi

