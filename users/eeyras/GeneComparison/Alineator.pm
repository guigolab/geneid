##########################################################################
#
#  ALINEATOR
#                                                                        #
#  written by Jesus Sanchez y Patricia Resa                              #
#  
#  modified by Eduardo Eyras
#                                                                        #
#  This program is free software; you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation; either version 2 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program; if not, write to the Free Software           #
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             #
##########################################################################
 
############################################################
# Este programa le va a permitir alinear exon 
# fingerprint de genes ortologos, proporcionando 
# una puntuación basada en el macht o mismatch de los 
# elementos del triplete que conforma el 
# exon fingerprint, penalizando los gaps.                                         
# Este programa utiliza programación dinámica.
############################################################

package GeneComparison::Alineator;

@ISA = qw(ClusterMerge::Root);

#########################################################################
 
sub new{
    my ($class, @args)=@_;
     
    if (ref($class)){
        $class = ref($class);
    }
    my $self = {};
    bless($self,$class);
     
    return $self;
}

############################################################

sub compare_fingerprints{
    my ($self,$ID1,$nex1,$finger1,$ID2,$nex2,$finger2) = @_;
    
    my @score;    #objetivo ultimo final e inconmensurable. 
    #Contiene el score del mejor alineamiento de los 2 fingerprints.
    
    # Para cada par de genes vamos a querer:
    # - Una variable singular que contenga el número de exones
    # - otra variable singular que guarde el ID.
    # - Un vector que contenga los fingerprints.
    
    # $nex1;     #contiene el número de exones del gen de la primera especie
    # $nex2;     #contiene el número de exones del gen de la segunda especie
    # $ID1;      #contiene el identificador del gen de la primera especie
    # $ID2;      #contiene el identificador del gen de la segunda especie
    
    # a fingerprint $finger is an arrayref
    # with  [$id,$nex,$exon1_finger,...,$exonN_finger]
    #
    my @finger1 = @$finger1;  #fingerprint del gen de la primera especie
    my @finger2 = @$finger2;  #fingerprint del gen de la segunda especie
    
    my @vectemp;  #vector temporal
    
    ## dado q necesitaremos el tamaño real de los genes(sólo el CDS) 
    # para compararlos, vamos a calcularlo ahora  
    ## tomaremos el tamanyo de cada exon y los iremos sumando, 
    # hasta tener el tamanyo total del CDS.
    
    my $tam1;   #tamanyo del 1er gen.
    my $tam2;   #tamanyo del 2o gen.
    
    $tam1 = $self->tamanyo(\@finger1);
    $tam2 = $self->tamanyo(\@finger2);
    
    
    ############################################################
    ## CONDICIONES INICIALES ##
    # estos dos genes pueden hallarse en una serie de situaciones distintas:
    # 1) mismo numero de exones, mismo tamanyo de CDS, los exones tienen tamanyos parecidos.
    # 2) mismo numero de exones, mismo tamanyo de CDS, diferente tamanyo de exones entre si.
    # 3) mismo numero de exones, diferente tamanyo del CDS.
    # 4) diferente numero de exones, mismo tamanyo del CDS.
    # 5) diferente numero de exones, diferente tamanyo del CDS.
    
    # definimos que el tamanyo de dos exones es parecido 
    # cuando difieren en menos de 6 nucleotidos (2 codones) por exon
    # para comparar el tamaño del CDS completo, multiplicamos 
    # por 6 el numero medio de exones entre los dos genes.
    # segun en que situacion se encuentre, recibira o no una 
    # penalizacion en el score final del alineamiento de los dos fingerprints.

    my $bonus;                #sera la penalizacion segun el caso.
    
    ############################################################
    # same number of exons
    if ($nex1==$nex2) {
	
	############################################################
	# if they differ by approximately 6bp per exon
	if ( abs($tam1-$tam2)<=(6*($nex1+$nex2)/2) ){
	    
	    #determina si los exones de cada fingerprint son parecidos en tamanyo o no.
	    my $difexon= $self->diferenciaexones(\@finger1, \@finger2);

	    ############################################################
	    # CASE 1
	    #si los exones tienen tamanyos parecidos (no son diferentes)
	    #si se cumplen estas condiciones estamos en el caso 1).
	    if ($difexon==0) {                  
		$bonus=1;
	    }
	    ############################################################
	    # CASE 2
	    #si los exones tienen tamanyos diferentes
	    else {                            
		$bonus=0.95;
	    }
	}
	############################################################
	# CASE 3
	else { 
	    $bonus=1-(abs($tam1-$tam2)/($tam1+$tam2));
	}
    }
    ############################################################
    # CASE 4
    # same size CDS (approx) but different number of exons
    else {
	if ( abs($tam1-$tam2) <=(6*($nex1+$nex2)/2) ){
	    $bonus=1-(abs($nex1-$nex2)/($nex1+$nex2));
	}
	############################################################
	# CASE 5
	else { 
	    $bonus=1 - ( (abs($nex1-$nex2)/($nex1+$nex2))*(abs($tam1-$tam2)/($tam1+$tam2)));
	}
    }
    ############################################################
    # ALINEAMIENTO DE FASES 
    ############################################################
    
    #################
    ##### PASO 1 ####
    my @score1;   #resultado de alinear los dos fingerprints sin realizarles ningun cambio.
    my @score2;   #resultado de alinear los dos fingerprints intentando combinar exones.
    
    @score1=$self->alinea(\@finger1,\@finger2);   #calculamos @score1.
    
    #################
    ##### PASO 2 ####
    my @refinger=$self->refingerprint(\@finger1,\@finger2);
    my @nuevfinger1=@{$refinger[0]};
    my @nuevfinger2=@{$refinger[1]};
    
    @score2=$self->alinea(\@nuevfinger1,\@nuevfinger2);   #calculamos @score2.

    ######################
    ## SCORE + BONUS #####
    ### NORMALIZACION ####
    
    
    my $normal = 0;
    if ( scalar(@score1)-1 ){
	$normal = 1;
	$score1[0]= ( ($score1[0]/( (scalar(@score1)-1)* 2 ) )*$bonus)* 100;
    }
    else{
	print "problem: ".$self->print_fingerprint($finger1);
    }
    if ( scalar(@score2)-1 ){
	$normal = 1;
	$score2[0]= ( ($score2[0]/( (scalar(@score2)-1)* 2 ) )*$bonus* 0.9)* 100;
    }
    else{
	print "problem: ".$self->print_fingerprint($finger2);
    }

    #print "Normalized\n" if ($normal);
    #print "score1[0] = $score1[0]\n";
    #print "score2[0] = $score2[0]\n";
    
    ######################
    
    ##### SCORE OPTIMO:###
    my @f1;                        #fingerprint utilizado para el alineamiento del gen 1
    my @f2;                        #fingerprint utilizado para el alineamiento del gen 2
    
    if ($score1[0]>=$score2[0]) {
	@score=@score1;
	@f1=@finger1;
	@f2=@finger2;
    }
    else {
	@score=@score2;
	@f1=@nuevfinger1;
	@f2=@nuevfinger2;
    }

    ######################
    ### SALIDA DE DATOS ##
   
    #print "#\t$vectemp[0]\n#$vectemp[1]\n";  

    #@vectemp es un vector temporal donde inicialmente guardamos cada linia, 
    #es decir almacena los dos fingerprints que alineamos.
    
    print "#ALINEATOR\n";
    $self->print_fingerprint($finger1);
    $self->print_fingerprint($finger2);
    
    $score[0] = 0.1 unless ($score[0]>0);
    print "SCORE\t$score[0]\n";
    
    my $i;             #contador de posiciones del vector con el alineamiento
    my $x;             #contador de posiciones del vector con el fingerprint
    
    # vectores que contienen las fases de los exones con los 
    # GAPs correspondientes intercalados:
    my @alineafases1;      
    my @alineafases2;      

    my $gaps1 = 0;
    my $gaps2 = 0;

    ############################################################
    #primer gen
    $i=1;
    $x=0;
    while ($i<scalar(@score)) {  
	$score[$i]=~ m/([\*\w\:]+)\s\|\s[\*\w\:]+/;
	if($1 eq '**GAP**') {
	    $gaps1++;
	    $alineafases1[$i-1]=$1;
	    $i=$i+1;
	}
	else {
	    $alineafases1[$i-1]=$f1[$x];
	    $i=$i+1;
	    $x=$x+1;
	}
    }

    ############################################################
    #segundo gen
    $i=1;
    $x=0;
    while ($i<scalar(@score)) {
  
	$score[$i]=~ m/[\w\*\:]+\s\|\s([\w\*\:]+)/;
	if($1 eq '**GAP**') {
	    $gaps2++;
	    $alineafases2[$i-1]=$1;
	    $i=$i+1;
	}
	else {
	    $alineafases2[$i-1]=$f2[$x];
	    $i=$i+1;
	    $x=$x+1;
	}
    }

    
    my $alignment = "@alineafases1\n@alineafases2";
    print "@alineafases1\n@alineafases2\n--\n";

    return ($score[0],$finger1,$finger2,\@alineafases1,\@alineafases2,$gaps1,$gaps2);

}


###############################################################################################################
###########################################  FUNCIONES  #######################################################
###############################################################################################################


###############################################################################################################

#NOMBRE: tamanyo
#OBJETIVO: calcular el tamanyo del CDS sumando los tamanyos de los exones codificantes.
#RETORNA: el tamanyo del CDS ($tamanyoCDS).

sub tamanyo{
    my ($self,$fingerprint) = @_;
    
    #print "tamanyo:\n";
    #$self->print_fingerprint($fingerprint);
    my $i=0;                         #contador de exones.
    
    my ($id, $nex, @fingerprint) = @$fingerprint;
    my $tamanyoCDS=0;                #suma de los tamanyos de los exones (tamanyo del CDS).
    
    while ($i<scalar(@fingerprint)){
	$fingerprint[$i]=~ m/\d\:\d\:(\d+)/;
	$tamanyoCDS = $tamanyoCDS+$1;
	$i=$i+1;
    }
    return $tamanyoCDS;
}
    
################################################################################################################

###############################################################################################################

#NOMBRE: diferenciaexones
#OBJETIVO: determinar si entre si los exones de un par de genes 
#con el mismo numero de exones difieren o no en tamanyo.
#RETORNA: the number of exons with differences which are greater than $max_diff

sub diferenciaexones{
    my ($self, $f1, $f2 ) = @_;
    
    $max_diff = 6;

    #print "difex:\n";
    #$self->print_fingerprint($f1);
    #$self->print_fingerprint($f2);


    my ($id1,$nex1,@f1) = @$f1;
    my ($id2,$nex2,@f2) = @$f2;
    
    my $i=0;          #contador de exones.
    my $difex = 0;    # counts the number of exons with 
                      # difference in length greater than max_diff
    my $a;            #almacena la resta en valor absoluto
    my $t1;           #tamanyo exon primera secuencia
    my $t2;           #tamanyo exon segunda secuencia
    
    while ($i<scalar(@f1)) {
	
	$f1[$i]=~ m/\d\:\d\:(\d+)/;
	$t1=$1;
	
	$f2[$i]=~ m/\d\:\d\:(\d+)/;
	$t2=$1;
	
	$a = abs($t1-$t2);

	unless ($a <= $max_diff) {
	    $difex++;
	}
    }
    
    return $difex;
}
###############################################################################################################


############################################################
#NOMBRE: alinea
#OBJETIVO: calcula el score del mejor alineamiento posible entre dos fingerprints.
#RETORNA: retorna el score del mejor alineamiento ($score);

sub alinea{
    my ($self,$f1,$f2) = @_;

    #print "alinea:\n";
    #$self->print_fingerprint($f1);
    #$self->print_fingerprint($f2);

    my ($id1,$nex1,@f1) = @$f1;
    my ($id2,$nex2,@f2) = @$f2;
    
    my $i=0;          #contador de exones.
    my @fases1;       #fases 5':3' de todos los exones del primer gen
    my @fases2;       #fases 5':3' de todos los exones del segundo gen
    my @matpunt;      #matriz de puntuaciones, resultado de alinear los dos fingerprints.
    my @score;        #puntuacion del alineamiento optimo y los emparejamientos.
    my $penGAP=-1;    #penalizacion por alinear con un GAP.
    
    #capturem les fases del 1er gen.
    while ($i<scalar(@f1)) {
	$f1[$i]=~ m/(\d\:\d)\:\d+/;   
	$fases1[$i]=$1;
	$i=$i+1;
    }
    #capturem les fases del 2on gen.
    $i=0;
    while ($i<scalar(@f2)) {
	$f2[$i]=~ m/(\d\:\d)\:\d+/;   
	$fases2[$i]=$1;
	$i=$i+1;
    }

    #contador de posiciones en la matriz (filas), 
    # la aprovecharemos para contabilizar los exones
    my $n1=1;       
    
    #contador de posiciones en la matriz (columnas), 
    # la aprovecharemos para contabilizar los exones
    my $n2=1;       

    $matpunt[0][0]=0;                  #definimos la posicion (0,0) de la matriz

    while ($n2<=scalar(@f2)){           #definimos la primera fila de la matriz como -1.
	$matpunt[0][$n2]=$penGAP;
	$n2=$n2+1;
    }
    while ($n1<=scalar(@f1)){           #definimos la primera columna de la matriz como -1;
	$matpunt[$n1][0]=$penGAP;
	$n1=$n1+1;
    }

    #reinicializamos los contadores de exones.
    $n1=1; 
    $n2=1;

    # backtracking:
    my $a; # go up-left (diagonal)
    my $b; # go up
    my $c; # go left

    ## para la recuperacion del alineamento.
    
    # matriz que contiene un codigo que define el mejor movimiento desde cualquier posicion.
    my @retroceso;     
    
    #indica el valor maxiomo de a,b y c y a cual corresponde.
    my @abcmax;        
    
    #definimos la posicion (0,0) de la matriz de retroceso.
    $retroceso[0][0]='x';       

    #definimos la primera fila de la matriz de retroceso como c.
    while ($n2<=scalar(@f2)){           
	$retroceso[0][$n2]='c';
	$n2=$n2+1;
    }
    #definimos la primera columna de la matriz como b;
    while ($n1<=scalar(@f1)){           
	$retroceso[$n1][0]='b';
	$n1=$n1+1;
    }
    $n1=1; 
    $n2=1;
    while ($n1<=scalar(@f1)){
	
	#reinicializamos los valores de las columnas.
	$n2=1;                     
	
	while ($n2<=scalar(@f2)){
	    
	    # $a es el resultado de alinear las fases de dos exones en las posiciones n1 y n2.
	    $a= $matpunt[$n1-1][$n2-1] + $self->puntuacion ($fases1[$n1-1],$fases2[$n2-1]);
	    
	    #$b es el resultado de alinear el exon en posicion n1 del primer gen con un gap.
	    $b= ($matpunt[$n1-1][$n2]) + $penGAP;   
	    
	    #$c es el resultado de alinear el exon en posicion n2 del segundo gen con un gap.
	    $c= ($matpunt[$n1][$n2-1]) + $penGAP;   
	    
	    #print "a=$a b=$b c=$c\n";
	    #calculamos el valor maximo de a,b y c.
	    @abcmax = $self->maxabc($a,$b,$c); 

	    #escogemos el alineamiento que obtiene la puntuacion maxima
	    $matpunt[$n1][$n2]  =$abcmax[0];   
	    $retroceso[$n1][$n2]=$abcmax[1];
	    
	    $n2=$n2+1;                    #avanzamos en las columnas.
	}

	$n1=$n1+1;                        #avanzamos en las filas.
    }
    
    #tomamos el valor final del alineamiento como puntuacion maxima.
    $score[0]= $matpunt[$n1-1][$n2-1];        

    $n1 = $n1 -1;
    $n2 = $n2 -1;
    #print "n1=$n1 n2=$n2\n";
    
    ############################################################
    # backtracking:
    my @temp= $self->retroceso(\@retroceso,$n1,$n2,\@fases1,\@fases2);
    push(@score,@temp);

    return @score;
}

############################################################
#NOMBRE: puntuacion
#OBJETIVO: especificar la puntuacion de alinear las fases de dos exones.
#RETORNA: el resultado de alinear estas fases.($result)
############################################################

sub puntuacion{
    my ($self, $fase1, $fase2) = @_;
    #fases del primer exon
    #fases del segundo exon

    my $result;       #puntuacion del alineamiento.
    
    $fase1=~ m/(\d)\:(\d)/;
    my $x=$1;         #fase 5' del primer exon
    my $y=$2;         #fase 3' del primer exon
    
    $fase2=~ m/(\d)\:(\d)/;
    my $w=$1;         #fase 5' del segundo exon
    my $z=$2;         #fase 3' del segundo exon

    ############################################################
    # phases-score == nuber of coincident phases
    if   ( $x==$w && $y==$z ){
	$result=2;
    }
    elsif( $x!=$w && $y!=$z ){
	$result=0;
    }
    else{
	$result=1;
    }
    return $result;
}


############################################################
#NOMBRE: maxabc
#OBJETIVO: calcula el valor de $a,$b y $c, e indica cual es.
#RETORNA: un hash de una posicion indicando la letra con puntuacion maxima y que valor tiene.
############################################################

sub maxabc{
    my ($self,@y) = @_;
    
    my $x=-99999999;   #variable que almacenara el valor maximo.
    my @abcmax;        #indica el valor maximo y en la segunda posicion, la letra a la que corresponde.
    
    if ($y[0]>$x){
	$x=$y[0];
	$abcmax[0]=$x;
	$abcmax[1]='a';
    }
    if ($y[1]>$x){
	$x=$y[1];
	$abcmax[0]=$x;
	$abcmax[1]='b';
    }
    if ($y[2]>$x){
	$x=$y[2];
	$abcmax[0]=$x;
	$abcmax[1]='c';
    }
    return @abcmax;
}



############################################################
#NOMBRE: retroceso
#OBJETIVO: retornar los emparejamientos de exones.
#RETORNA: un vector con los exones emparejados con otros exones o con GAPs(@apareamiento)
#input:  matriz de punteros, longitud fingerprint1 + 2, long fingerprint2 +2, fases1, fases2
#uso:     my @temp= retroceso(\@retroceso,$n1, $n2,\@fases1, \@fases2);
############################################################

sub retroceso{
    my ($self, $r, $n1, $n2, $f1, $f2,) = @_;
    
    #matriz con los movimientos optimos.
    my @recorrido=@{$r};             
    my @f1 = @$f1; #vector con las fases del 1er gen.
    my @f2 = @$f2; #vector con las fases del 2o gen.
    
    # $n1 indica en que exon de 1 nos encontramos.
    # $n2 indica en que exon de 2 nos encontramos.

    my @apareamiento;                   #contiene todos los apareamientos.
    
    #print "fases1: ".scalar(@f1).": @f1\n";
    #print "fases2: ".scalar(@f2).": @f2\n";
    #print "n1 = $n1 n2 = $n2\n";
    #print "rec = @recorrido\n";
    #print scalar(@recorrido)."x".scalar(@{$recorrido[0]})."\n";
    #print "re[n1=$n1][n2=$n2] = $r[$n1][$n2]\n";
    
    while ($recorrido[$n1][$n2] ne 'x') {
        
	#movimiento en diagonal.
	if ($recorrido[$n1][$n2] eq 'a'){
	    unshift(@apareamiento,"$f1[$n1-1] | $f2[$n2-1]");
	    $n1=$n1-1;
	    $n2=$n2-1;
	}
	#movimiento hacia la izquierda (alineamos n2 con un gap);
	elsif ($recorrido[$n1][$n2] eq 'c'){
	    unshift(@apareamiento,"**GAP** | $f2[$n2-1]");
	    #corremos hacia la izquierda.
	    $n2=$n2-1;    
	}
	#movimiento hacia arriba (alineamos n1 con un gap);
	elsif ($recorrido[$n1][$n2] eq 'b'){
	    unshift(@apareamiento,"$f1[$n1-1] | **GAP**");
	    #corremos hacia arriba.
	    $n1=$n1-1;    
	}
    }

    return @apareamiento;
}

############################################################
#NOMBRE: refingerprint
#OBJETIVO:toma los fingerprint de dos genes,compara 
#         el tamaño de los exones, y si es posible los fusiona en 
#         un solo exon, cambiando el fingerprint del gen. Fusionara todos los exones posibles.
#         se consideraran tamaños distintos si es superior a 6.
#RETORNA: dos vectores, cada uno con el nuevo fingerprint y además el nuevo numero de exones.
############################################################

sub refingerprint{
    my ($self,$f1,$f2) = @_;

    my $max_diff = 6;
 
    my ($id1,$nex1,@finger1) = @$f1;
    my ($id2,$nex2,@finger2) = @$f1;

    my @f1;                  
    my @f2;
  
    my $i1=0;          #contador de exones de primer gen
    my $i2=0;          #contador de exones del segundo gen

    my @tamanyos1;     #vector q guardo los tamanyos de los exones originales del primer gen
    my @tamanyos2;     #vector q guardo los tamanyos de los exones originales del segundo gen

    my $tammod=0;
    my $tammod2=0;
    my $tammodb=0;
    my $tammod2b=0;
    
    my $refing;
    my $refingb;

    while ($i1<scalar(@finger1)) {
	$finger1[$i1]=~ m/\d\:\d\:(\d+)/;
	$tamanyos1[$i1]=$1;
	$i1=$i1+1;
    }
    while ($i2<scalar(@finger2)){
	$finger2[$i2]=~ m/\d\:\d\:(\d+)/;
	$tamanyos2[$i2]=$1;
	$i2=$i2+1;
    }

    $i1=0;
    $i2=0;
    
  GRAN_LOOP:
    while ( ($i1<scalar(@finger1)) && ($i2<scalar(@finger2)) )  {
	
	if(($self->compara($tamanyos1[$i1],$tamanyos2[$i2]))==1) {    
	    
	    #los dos exones tienen el mismo tamaño
	    #add el fingerprint de ese exon al vector definitivo
	    
	    push(@f1,$finger1[$i1]);                           
	    
	    #add el fingerprint de ese exon al vector definitivo
	    push(@f2,$finger2[$i2]);                           

	    $i1=$i1+1;
	    $i2=$i2+1;
	}
	elsif(($tamanyos1[$i1]>$tamanyos2[$i2])&&($tamanyos2[$i2+1])) { 

	    #los exones tienen distinto tamanyo	   
	    $tammod=(($tamanyos2[$i2])+($tamanyos2[$i2+1])); 
	    
	    #el nuevo exon tiene tamanyo similar al q era +grande
	    if($self->compara($tamanyos1[$i1],$tammod)==1) {            
		
		#calcula las nuevas fases e incluye el nuevo tamanyo
		$refingb= $self->refing($finger2[$i2],$finger2[$i2+1]); 

		$refing="$refingb:$tammod";
       
		push(@f1,$finger1[$i1]);
		push(@f2,$refing);
		
		$i1=$i1+1;
		#tiene q avanzar dos posiciones en el vector original
		$i2=$i2+2;                                     
	    }
	    elsif($tamanyos1[$i1]>$tammod && ($tamanyos2[$i2+2]) ){   

		#el nuevo tamanyo sigue siendo menor, tratamos de meter un tercer exon
		$tammod2=($tammod+$tamanyos2[$i2+2]);
		
		if($tamanyos1[$i1]+ $max_diff < $tammod2) {                  

		    #no entra un tercer exon
		    $refingb=$self->refing($finger2[$i2],$finger2[$i2+1]);
		    $refing="$refingb:$tammod";
		    
		    push(@f1,$finger1[$i1]);
		    push(@f2,$refing);
		    
		    $i1=$i1+1;
		    $i2=$i2+2;
                }
		else{
		    $refingb=$self->refing($finger2[$i2],$finger2[$i2+2]);
		    $refing="$refingb:$tammod2";
		    
		    push(@f1,$finger1[$i1]);
		    push(@f2,$refing);
		    
		    $i1=$i1+1;
		    $i2=$i2+3;
		}
		
	    }
	    else {  
		#el nuevo tamanyo es más grande asi q lo descartamos
		#add el fingerprint de ese exon al vector definitivo
		push(@f1,$finger1[$i1]);                           
		
		#add el fingerprint de ese exon al vector definitivo
		push(@f2,$finger2[$i2]);                           
	                                              
		$i1=$i1+1;
		$i2=$i2+1;
		
	    }
	}
	elsif($tamanyos1[$i1+1]) {
	    
	    $tammodb=(($tamanyos1[$i1])+($tamanyos1[$i1+1])); 
	    
	    #el nuevo exon tiene tamanyo similar al q era +grande
	    if($self->compara($tamanyos2[$i2],$tammodb)==1) {            
		
		#calcula las nuevas fases e incluye el nuevo tamanyo
		$refingb=$self->refing($finger1[$i1],$finger1[$i1+1]); 
		$refing="$refingb:$tammodb";
	
		push(@f2,$finger2[$i2]);
		push(@f1,$refing);

		#tiene q avanzar dos posiciones en el vector original
		$i1=$i1+2;
		$i2=$i2+1;                                     
	    }

	    elsif( ($tamanyos2[$i2]>$tammodb) &&($tamanyos1[$i1+2]) ){     
	
		#el nuevo tamanyo sigue siendo menor, tratamos de meter un tercer exon
		$tammod2b=($tammodb+$tamanyos1[$i1+2]);
		
		if($tamanyos2[$i2]+6<$tammod2b) {                  
		    
		    #no entra un tercer exon
		    $refingb=$self->refing($finger1[$i1],$finger1[$i1+1]);
		    $refing="$refingb:$tammodb";
		   		    
		    push(@f2,$finger2[$i2]);
		    push(@f1,$refing);

		    $i1=$i1+2;
		    $i2=$i2+1;
                }
		else{
		    $refingb=$self->refing($finger1[$i1],$finger1[$i1+2]);
		    $refing="$refingb:$tammod2b";
		    
		    push(@f2,$finger2[$i2]);
		    push(@f1,$refing);

		    $i1=$i1+3;
		    $i2=$i2+1;
		}
		
	    }

	    else {  
		#el nuevo tamanyo es más grande asi q lo descartamos
		#add el fingerprint de ese exon al vector definitivo
		push(@f1,$finger1[$i1]);
		
		#add el fingerprint de ese exon al vector definitivo
		push(@f2,$finger2[$i2]);                           	                                              
		$i1=$i1+1;
		$i2=$i2+1;
	    }
	}
	else {
	    #el nuevo tamanyo es más grande asi q lo descartamos
	    #add el fingerprint de ese exon al vector definitivo
	    push(@f1,$finger1[$i1]);                         
	    
	    #add el fingerprint de ese exon al vector definitivo
	    push(@f2,$finger2[$i2]);                           
	    
	    $i1=$i1+1;
	    $i2=$i2+1;
	}
	
    } #end of GRAN_LOOP
    ############################################################
    if ($i1<scalar(@finger1)) {
	while ($i1<scalar(@finger1)) {
	    push(@f1,$finger1[$i1]);
	    $i1=$i1+1;
	}
    }
    if ($i2<scalar(@finger2)) {
	while ($i2<scalar(@finger2)) {
	    push(@f2,$finger2[$i2]);
	    $i2=$i2+1;
	}
    }    
    return (\@f1,\@f2);
}


############################################################
#NOMBRE:compara
#OBJETIVO:dados dos numeros 
#         (correspondientes a los tamanyos de dos exones) compara si son iguales
#         (igual significa que no difieran en mas de $max_diff residuos).
#RETORNA: son iguales ($x=1) o diferentes ($x=0)
############################################################

sub compara{
    my ($self,$a,$b) = @_;
    
    my $max_diff = 6;
    my $x;
    if (abs($a-$b)>$max_diff){
	$x=0;
    }
    else{
	$x=1;
    }
    return $x;
}

############################################################
#NOMBRE:   refing
#OBJETIVO: merge the fingerprint from two exons
#RETORNA:  the new intron phases for the merged exon ($x)
############################################################

sub refing{
    my ($self,$exon1,$exon2) = @_;
    my $cinco;                 #fase en 5'
    my $tres;                  #fase en 3'
    my $x;
    $exon1=~m/(\d)\:\d\:\d+/;
    $cinco=$1;
    $exon2=~m/\d\:(\d)\:\d+/;
    $tres=$1;
    $x="$cinco:$tres";
    return $x;
}

############################################################

sub print_fingerprint{
    my ($self,$f) = @_;
    my $s = join "\t", @$f;
    print $s."\n";
}

1;

