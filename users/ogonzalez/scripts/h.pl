#!/usr/bin/perl

# QUAN FIQUIS NUMEROS DE PROVA, INTENTA QUE A FALTA DE 4 ANYS TINGUIS
# UN DEUTE DE 33000 EUROS APROX. ES LA MANERA DE MAXIMITZAR ELS
# "BENEFICIS": L'ESTAT ET TORNARA 1800x4 = 7200 APROX I ELS INTERESSOS
# QUE PAGARAS SERAN UNS 3000 ==> EN TREURAS UNS 4000!!!
#
# HE DE REVISAR LA DECLARACIO PERQUE AIXO NO ES AIXI 100%, PERO VAJA...


# calcula la mensualitat que correspon a un nombre donat de mesos per
# pagar, un deute inicial donat i uns interesos donats
# sembla tenir un error aproximat, possiblement per problemes de precisio
# de +-1 pts / +-1 centim
sub mensualitat {

   $nf=$_[0]; # mesos fins al final
   $d0=$_[1]; # deute inicial
   $i=$_[2];  # interes

   $numerador=$d0*exp($nf*log(1+$i/1200)); # faig les potencies amb exp i log
                                           # perque no tinc funcio per fer-ho
   $denominador=0;
   for ($j=0;$j<=$nf-1;$j++) {
      $denominador+=exp($j*log(1+$i/1200)); # idem
   }

   return ($numerador/$denominador);

}


sub amortitzacio {

   $m=$_[0]; # mensualitat
   $do=$_[1]; # deute inicial
   $i=$_[2]; # interes
   $n=$_[3]; # mes actual

   $primer=($m-$i*$d0/1200);
   $segon=exp(($n-1)*log(1+$i/1200));

   return ($primer*$segon);

}


# retorna l'interes
# a partir de $#interes, suposo que l'interes no varia
sub interes {
   if ( $indexany > $#interes ) {
      return $interes[$#interes];
   } else {
      return $interes[$indexany];
   }
}


# retorna les despeses de correu
# a partir de $#correu, suposo que es mantindra con l'any anterior
sub correu {
   if ( $indexany > $#correu ) {
      return $correu[$#correu];
   } else {
      return $correu[$indexany];
   }
}



# valors inicials - comenco a partir del canvi a euros de la hipoteca!!!
$prestec=96675.28;
$duracio=340;        # mesos que queden des del mes de deute $prestec
$mesact=11;          # mes en que abans de pagar dec $prestec
$anyact=2001;        # any de idem

# primer any: 11-2001 a 01-2002 (recorda: NOMES des del canvi a euros)
$interes[1]=5.25;    # interes de l'any zero
$correu[1]=0.25;     # despeses de correu

# segon any: 02-2002 a 01-2003
$interes[2]=4.75;    # interes de l'any 1
$extra{"0204"}=3000; # amortitzo a efectes del mes 4 (abril)
$correu[2]=0.25;

# tercer any: 02-2003 a 01-2004
$interes[3]=4.25;
$extra{"0303"}=6000;
$reduc[3]=48;        # redueixo 4 anys ==> 48 mesos
$correu[3]=0.26;

# quart any: 02-2004 a 01-2005
$interes[4]=3.25;    # passo a euribor+0.75
$extra{"0403"}=5000;
$reduc[4]=84;

# cinque any: 02-2005 a 01-2006
$interes[5]=3.50;
$extra{"0503"}=3300;
$reduc[5]=0;

# sise any: 02-2006 a 01-2007
$interes[6]=4.00;
$extra{"0603"}=3300;
$reduc[6]=12;

# sete any: 02-2007 a 01-2008
$interes[7]=4.50;
$extra{"0703"}=3600;
$reduc[7]=12;

# vuite any: 02-2008 a 01-2009
$interes[8]=5.00;
$extra{"0803"}=3600;
$reduc[8]=12;

# nove any: 02-2009 a 01-2010
$interes[9]=5.25;
$extra{"0903"}=3900;
$reduc[9]=12;

# dese any: 02-2010 a 01-2011
$interes[10]=5.25;
$extra{"1003"}=3900;
$reduc[10]=0;

# onze any: 02-2011 a 01-2012
$interes[11]=5.25;
$extra{"1103"}=4200;
$reduc[11]=12;

# ja no queden mes anys. final de la hipoteca

$men=mensualitat($duracio,$prestec,$interes[1]);
$total=1;

$indexany=1; # any d'hipoteca en que em trobo

$totint=0; # interes total de l'any
$acuint=0; # interes acumulat durant el prestec
$totamo=0; # amortitzacio total de l'any
$acuamo=0; # amortitzacio acumulada durant el prestec
while ( $duracio>0 ) {
   if ( $mesact == 3 ) { # reviso la hipoteca el mes 3: març
      $indexany++;
      if ( $indexany <= $#interes ) {
         $total=1;
         $duracio-=$reduc[$indexany];
         $men=mensualitat($duracio,$prestec,interes);
      }
   }
   $anyambcero=substr("00".$indexany,length($indexany),2);
   $mesambcero=substr("00".$mesact,length($mesact),2);
   $anyimes=$anyambcero.$mesambcero;
   if ( exists($extra{$anyimes})) {
      $prestec=$prestec-$extra{$anyimes};
      $men=mensualitat($duracio,$prestec,interes);
      $totamo+=$extra{$anyimes};
      printf ("%2i-%4i        -",$mesact,$anyact);
      print " "x(5-length($extra{$anyimes}));
      printf ("%7.2f",$extra{$anyimes});
      printf ("%21s --> %8.2f\n","",$prestec);
      $total=1;
   }
   $amo=amortitzacio($men,$prestec,interes,$total);
   $prestec=$prestec-$amo;
   $totint+=($men-$amo);
   $totamo+=$amo;

   printf ("%2i-%4i (%4.2f) -  ",$mesact,$anyact,&interes);
   printf ("%6.2f %6.2f %4.2f = ",$amo,($men-$amo),&correu);
   printf ("%6.2f --> %8.2f\n",($men+&correu),$prestec);
   if ( $mesact == 2 ) {
      $acuint+=$totint;
      $acuamo+=$totamo;
      printf ("%40s interesos   : anuals    %9.2f\n","",$totint);
      printf ("%54s acumulats %9.2f\n","",$acuint);
      printf ("%40s amortitzacio: anual     %9.2f\n","",$totamo);
      printf ("%54s acumulada %9.2f\n","",$acuamo);
      printf ("%40s total pagat : anual     %9.2f\n","",($totint+$totamo));
      printf ("%54s acumulat  %9.2f\n","",($acuint+$acuamo));
      if ( $reduc[$indexany] != 0 ) {
         print " "x41;
         print "vaig reduir " . $reduc[$indexany]/12 . " any";
         print "s" if ( $reduc[$indexany]/12 > 1 );
         print "\n";
      }
      print " "x41;
      print "ara falten " . ($duracio-1)/12 . " anys\n\n";
      $totint=0;
      $totamo=0;
   }

   $mesact++;
   if ( $mesact > 12 ) {
      $mesact=1;
      $anyact++;
   }

   $total++;
   $duracio--;

}
