BEGIN{
  # codon table
  codon["GCA"]= "A";
  codon["GCC"]= "A";
  codon["GCG"]= "A";
  codon["GCT"]= "A";
  codon["TGC"]= "C";
  codon["TGT"]= "C";
  codon["GAC"]= "D";
  codon["GAT"]= "D";
  codon["GAA"]= "E";
  codon["GAG"]= "E";
  codon["TTC"]= "F";
  codon["TTT"]= "F";
  codon["GGA"]= "G";
  codon["GGC"]= "G";
  codon["GGG"]= "G";
  codon["GGT"]= "G";
  codon["CAC"]= "H";
  codon["CAT"]= "H";
  codon["ATA"]= "I";
  codon["ATC"]= "I";
  codon["ATT"]= "I";
  codon["AAA"]= "K";
  codon["AAG"]= "K";
  codon["TTA"]= "L";
  codon["TTG"]= "L";
  codon["CTA"]= "L";
  codon["CTC"]= "L";
  codon["CTG"]= "L";
  codon["CTT"]= "L";
  codon["ATG"]= "M";
  codon["AAC"]= "N";
  codon["AAT"]= "N";
  codon["CCA"]= "P";
  codon["CCC"]= "P";
  codon["CCG"]= "P";
  codon["CCT"]= "P";
  codon["CAA"]= "Q";
  codon["CAG"]= "Q";
  codon["AGA"]= "R";
  codon["AGG"]= "R";
  codon["CGA"]= "R";
  codon["CGC"]= "R";
  codon["CGG"]= "R";
  codon["CGT"]= "R";
  codon["AGC"]= "S";
  codon["AGT"]= "S";
  codon["TCA"]= "S";
  codon["TCC"]= "S";
  codon["TCG"]= "S";
  codon["TCT"]= "S";
  codon["ACA"]= "T";
  codon["ACC"]= "T";
  codon["ACG"]= "T";
  codon["ACT"]= "T";
  codon["GTA"]= "V";
  codon["GTC"]= "V";
  codon["GTG"]= "V";
  codon["GTT"]= "V";
  codon["TGG"]= "W";
  codon["TAC"]= "Y";
  codon["TAT"]= "Y";
  codon["TAA"]= "*";
  codon["TAG"]= "*";
  codon["TGA"]= "*";

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

  # gff formatted output
  gff_format="%s\t%s\t%s\t%d\t%d\t%.3f\t%s\t%d\t%s\n";


  # defaults
  FRAMES = 3;

  DSC=38  # substract this from hsp score


  # read sequences
  # sequence 1
  getline<ARGV[1];
  locus1=$1;
  seq[locus1]=toupper($2);

  getline<ARGV[2];
  if ($1!=locus1) 
    print "WARNING: locus names do not match"
  seq_r[$1]=toupper($2);

  lenseq[locus1]=length(seq[locus1]);

  # sequence 2
  getline<ARGV[3];
  locus2=$1;
  seq[locus2]=toupper($2);

  getline<ARGV[4];
  if ($1!=locus2) 
    print "WARNING: locus names do not match"
  seq_r[locus2]=toupper($2);

  lenseq[locus2]=length(seq[locus2]);

  # read substitution matrix
  while (getline<ARGV[5]>0)
    SMAT[$1,$2]=$3;

  # read SRs regions
  while (getline<ARGV[6]>0) {
    if (substr($0,1,1) != "\#"){
      nsr++;

      if ($1 != locus1) {
	print "ERROR: not matching locus names in sequence 1";
	exit;
      }

      if (substr($10,2,length($10)-2) != locus2) {
	print "ERROR: not matching locus names in sequence 2";
	exit;
      }

      # get hsp info
      hsp[nsr,"start" ]  = $4;
      hsp[nsr,"end"   ]  = $5;
      hsp[nsr,"length"]  = $5-$4+1;
      hsp[nsr,"frame" ]  = $8;
      hsp[nsr,"strand"]  = $7;
      
      hsp[nsr,"score"]   = $6-DSC;

      # hsp sequence
      if (hsp[nsr,"strand"] == "\+") {
	rstart=hsp[nsr,"start"]; # required later
	hsp[nsr,"sequence"]=substr(seq[locus1],rstart,hsp[nsr,"length"]);
      }
      else {
	rstart=lenseq[locus1] - hsp[nsr,"end"] + 1;
	hsp[nsr,"sequence"]=substr(seq_r[locus1],rstart,hsp[nsr,"length"]);
      }
      
      # hsp aaseq
      fra = (hsp[nsr,"frame"] + 3 - (rstart % 3)) % 3; #absolute frame
      hsp[nsr,"absframe"]  = fra;
      hsp[nsr,"aaseq"]     = translate(codon,hsp[nsr,"sequence"],hsp[nsr,"length"],fra);

      # get hsp alseq info

      hsp[nsr,"alstart" ]  = $11;
      hsp[nsr,"alend"   ]  = $12;
      hsp[nsr,"allength"]  = $12-$11+1;
      hsp[nsr,"alframe" ]  = $18;
      hsp[nsr,"alstrand"]  = $16;

      # hsp al. sequence
      if (hsp[nsr,"alstrand"] == "\+") {
	rstart=hsp[nsr,"alstart"]; # required later
	hsp[nsr,"alsequence"]=substr(seq[locus2],rstart,hsp[nsr,"allength"]);
      }
      else {
	rstart=lenseq[locus2] - hsp[nsr,"alend"] + 1;
	hsp[nsr,"alsequence"]=substr(seq_r[locus2],rstart,hsp[nsr,"allength"]);
      }
      
      # hsp al. aaseq
      fra = (hsp[nsr,"alframe"] + 3 - (rstart % 3)) % 3; #absolute frame

      hsp[nsr,"alabsframe"]  = fra;
      hsp[nsr,"alaaseq"]   = translate(codon,hsp[nsr,"alsequence"],hsp[nsr,"allength"],fra);

      if (hsp[nsr,"length"] != hsp[nsr,"allength"]) 
      print "WARNING: hsp differ in sequences length"
      
      #recompute blosum scores
      hsp[nsr,"blscore"]=sim(hsp[nsr,"aaseq"],hsp[nsr,"alaaseq"],SMAT);

    }
  }
#  print nsr;
#  for (i=1;i<=nsr;i++) {

#    print locus1, hsp[i,"start"], hsp[i,"end"], hsp[i,"strand"], hsp[i,"frame"], hsp[i,"score"], hsp[i,"blscore"];
#      print locus2, hsp[i,"alstart"], hsp[i,"alend"], hsp[i,"alstrand"], hsp[i,"alframe"];
#  print hsp[i,"aaseq"];
#  print hsp[i,"alaaseq"];
#  print "\n";
#  }

  ARGV[1]=ARGV[2]=ARGV[3]=ARGV[4]=ARGV[5]=ARGV[6]="";
}
{
  # for each exon, walk over all srs
  SRscore=0;

  # get exon frame compatible with hsps
  if ($strand == "+")
      ExonFrame=(($start+$frame) % 3);
  else 
      ExonFrame=(($end-$frame) % 3);

 
  for (i=1; i<=nsr;i++) {
  # get hsp frame compatible with exon
  ind = hsp[i,"strand"] == "+" ? 0 : FRAMES;
  HSPFrame = (hsp[i,"frame"] % 3) + ind;

#  if (HSPFrame==5)
#    HSPFrame=4;
#  else if (HSPFrame==4)
#      HSPFrame=5;

    print ExonFrame, HSPFrame , $strand , hsp[i,"strand"], hsp[i,"start"], hsp[i,"end"], $start, $end;
    if (ExonFrame == HSPFrame && $strand == hsp[i,"strand"]) { 
      print $start, $end, hsp[i,"start"], hsp[i,"end"];

      # compatible frames and strands
      a = min($end,hsp[i,"end"]);
      b = max($start,hsp[i,"start"]);
      if (b <= a) {
	print "overlap"
      #overlap, rescore exons
	lIntersection = a - b + 1;
	SRscore+=hsp[i,"score"] * (lIntersection/hsp[i,"length"]);
	print hsp[i,"score"],SRscore, lIntersection, hsp[i,"length"];
      }
    }
  }
  print  $seqname, $source, $feature, $start, $end, $score+SRscore, $strand, $frame, SRscore, $score;
#  printf   gff_format, $seqname, $source, $feature, $start, $end, $score+SRscore, $strand, $frame;
}

function max(x,y) {
  return x > y ? x : y;
}

function min(x,y) {
  return x < y ? x : y;
}

function translate(c,s,l,f,   i) {
  aas="";

  for (i=1+f;i<l-1;i+=3) {
    triplet = substr(s,i,3);
    if (triplet in c)
      aas = aas c[triplet];
    else
      aas = aas "X";
  }

  return aas;
}


function sim(s1,s2,wmat) {
  
  d=0;
  l1=length(s1);
  l2=length(s2);

  l = l1 < l2 ? l1 : l2;

  for (k=1; k<=l;k++)
    d+=wmat[substr(s1,k,1),substr(s2,k,1)];
  return d;
}
