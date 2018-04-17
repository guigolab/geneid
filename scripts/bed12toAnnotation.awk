# ===============================================================
# awk script to extract intron/5UTR/3UTR/CDS from bed12 annotation
# Author: Xianjun Dong
# Date: 29 Jan 2016
# the input file is a bed12 file, which can be converted from GTF using gtf2bed
# the output file is a bed12 file for the required part of gene structure (cds, intron, 3utr, 5utr)
# Usage: cat input.bed12 | awk -v type=5utr -f bed12toAnnotation.awk
# Note: for CDS, it's recommended to run "grep protein_coding.protein_coding" before hand.e.g.
# fgrep protein_coding___protein_coding gencode.v19.annotation.bed12 | awk -v type=cds -f ~/pipeline/bin/bed12toAnnotation.awk > gencode.v19.annotation.gtf.cds.bed12
#
#
# Copyright (c) 2017 Xianjun Dong (xdong@rics.bwh.harvard.edu)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
# ALSO, IT WOULD BE NICE IF YOU LET ME KNOW YOU USED IT.
# ===============================================================
#!/bin/awk -f

{
  OFS="\t";
  split($11,blockSizes,","); 
  split($12,blockStarts,","); 
  blockCount=$10;
  A=""; B="";
  N=0;
  
  ## 5' UTR
  if(type=="5utr" || type=="5UTR" || type=="utr5" || type=="UTR5")
  {
    if($7==$8) next; 
    if($6=="+" && $2<$7) {
      start=$2;end=$7;
      for(i=1;i<=blockCount;i++) if(($2+blockStarts[i]+blockSizes[i])<=$7) {
        A=A""blockSizes[i]",";B=B""blockStarts[i]","; end=($2+blockStarts[i]+blockSizes[i]); N++;} 
        else { 
          if(($2+blockStarts[i])<$7) {A=A""($7-$2-blockStarts[i])",";B=B""blockStarts[i]","; N++; end=$7;} 
          break; 
        } 
      print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;
    } 
    if($6=="-" && $8<$3) {
      start=$8;end=$3;
      for(i=1;i<=blockCount;i++) if(($2+blockStarts[i])>=$8) 
      {
        if(start==0) {A=blockSizes[i];B=0; start=$2+blockStarts[i];}
        else {A=A","blockSizes[i];B=B","($2+blockStarts[i]-start);}
        N++;
      } else 
      {
        if(($2+blockStarts[i]+blockSizes[i])>$8) {
            A=($2+blockStarts[i]+blockSizes[i]-$8);B=0; N++; start=$8;
          }
        if(($2+blockStarts[i]+blockSizes[i])==$8) start=0; 
      }
      print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;
    }  
  }
  
  ## 3' UTR
  if(type=="3utr" || type=="3UTR" || type=="utr3" || type=="UTR3")
  {
    if($7==$8) next; 
    if($6=="-" && $2<$7) {
      start=$2;end=$7;
      for(i=1;i<=blockCount;i++) if(($2+blockStarts[i]+blockSizes[i])<=$7) {A=A""blockSizes[i]",";B=B""blockStarts[i]","; end=($2+blockStarts[i]+blockSizes[i]); N++;} else { if(($2+blockStarts[i])<$7) {A=A""($7-$2-blockStarts[i])",";B=B""blockStarts[i]","; N++; end=$7;} break; } 
      print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;
    } 
    if($6=="+" && $8<$3) {
      start=$8;end=$3;
      for(i=1;i<=blockCount;i++) if(($2+blockStarts[i])>$8) 
      {
        if(start==0) {A=blockSizes[i];B=0; start=$2+blockStarts[i];}
        else {A=A","blockSizes[i];B=B","($2+blockStarts[i]-start);}
        N++;
      } else 
      {
        if(($2+blockStarts[i]+blockSizes[i])>$8) {
          A=($2+blockStarts[i]+blockSizes[i]-$8);B=0; N++; start=$8;
        }
        if(($2+blockStarts[i]+blockSizes[i])==$8) start=0;  # for case where the thickEnd is at the last nt of an exon
      }
      print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;
    }  
  }
  
  ## CDS
  if(type=="cds" || type=="CDS")
  {
    for(i=1;i<=blockCount;i++) if(($2+blockStarts[i]+blockSizes[i])>$7 && ($2+blockStarts[i])<$8) {
      N++; 
      start=$2+blockStarts[i]-$7; size=blockSizes[i]; 
      if(($2+blockStarts[i])<=$7) {start=0;size=size-($7-($2+blockStarts[i]));} 
      if(($2+blockSizes[i]+blockStarts[i])>=$8) {size=size-($2+blockSizes[i]+blockStarts[i]-$8);} 
      A=A""size",";B=B""start",";
    } 
    print $1,$7,$8,$4,$5,$6,$7,$8,$9,N,A,B;
  }
  
  ## intron
  if(type=="intron" || type=="INTRON")
  {
    if(blockCount>1) {
      for(i=1;i<blockCount;i++) {A=A""(blockStarts[i+1]-blockStarts[i]-blockSizes[i])",";B=B""(blockStarts[i]+blockSizes[i]-(blockStarts[1]+blockSizes[1]))",";} 
      print $1,$2+blockSizes[1], $3-blockSizes[blockCount], $4,$5,$6,$2+blockSizes[1], $3-blockSizes[blockCount],$9,blockCount-1,A,B;
    }
  }
  
}