#!/bin/bash

# script ultra-guarro que compara les quotes de n fitxers donats
#
# compara la quota dels usuaris mostrats al fitxer n amb les que tenien al
# fitxer n-1. els usuaris a mostrar els trec del darrer fitxer (hauria de
# ser el mes recent), aixi que alguns apareixeran una vegada, altres no
# apareixeran i altres apareixeran 2 vegades.
# aixo es aixi perque de vegades nomes examinem la quota de ug o de um
#
# amb el flag -g o -m nomes comparem els usuaris de /home/ug o /home/um,
# independentment de que hi hagi mes usuaris al fitxer de quota.
#
# amb el flag -b nomes revisem l'ocupacio de backup; amb -n la de no_backup;
# amb -t la total. son mutuament excloents i tenen preferencia -t > -b > -n.
# per omisio faig -t
#
# amb -h l'output es en format html (per enviar per mail o jo que se), pero
# es ultra-guarro
#
# NO ORDENA ELS FITXERS PROPORCIONATS PER DATA!!! t'has d'assegurar que
# els passes ordenats de mes antic a mes nou


# extreu la data del nom del fitxer
data() {
   datanoquota=${1#*quota-}
   datanotxt=${datanoquota%-????.txt}
   any=${datanotxt%????}
   mesdia=${datanotxt#$any}
   echo -n $any-${mesdia%??}-${mesdia#??}
}


pattern=""
camp=0
output=txt
nomesflags=0
while [ $# != 0 -a $nomesflags = 0 ]
do
   case $1 in
      -g | -ug ) pattern=g${pattern#g}
                 shift
                 ;;
      -m | -um ) pattern=${pattern%m}m
                 shift
                 ;;
            -t ) camp=4 # nomes revisem ocupacio total
                 shift
                 ;;
            -b ) if [ $camp != 4 ]
                 then
                    camp=1 # nomes revisem ocupacio de backup
                 fi
                 shift
                 ;;
            -n ) if [ $camp = 0 ]
                 then
                    camp=2 # nomes revisem ocupacio de no_backup
                 fi
                 shift
                 ;;
            -h ) output=html
                 shift
                 ;;
            -* ) echo Flag $1 no reconegut
                 exit
                 ;;
             * ) nomesflags=1
                 ;;
   esac
done
if [ "$pattern" = "" ]
then
   pattern="gm"
fi
if [ $camp = 0 ]
then
   camp=4
fi

typeset -i index
index=0
while [ $# != 0 ]
do
   if [ ! -s $1 ]
   then
      echo $1 no existeix
   else
      let index=index+1
      fit[$index]=$1
      total[$index]=0
   fi
   shift
done

# capceleres: data (agafada del nom) de cada fitxer i subratllat
if [ $output = html ]
then
   echo "<html><head><style type=\"text/css\">td { font-size: 12 }</style></head>"
   echo "<body><table>"
   echo -n "<tr><td></td><td align=center><b>"
else
   echo
   echo -n "            "
fi
i=1
while (( i <= index ))
do
   data ${fit[$i]}
   if (( i != index ))
   then
      if [ $output = html ]
      then
         echo -n "</b></td><td align=center><b>"
      else
         echo -n " -> "
      fi
   fi
   let i=i+1
done
if [ $output = html ]
then
   echo "</b></td></tr>"
else
   echo
   echo -n "            "
   i=1
   while (( i <= index ))
   do
      echo -n "----------"
      if (( i != index ))
      then
         echo -n "    "
      fi
      let i=i+1
   done
   echo
fi

# vull mantenir l'ordenacio per ocupacio, pero "grep -f fitxer $1 $2" mostra
# primer tots els "match" de $1 i despres tots els de $2, un sort per nom es
# carrega l'ordenacio per ocupacio, i un sort per ocupacio es carrega la del
# nom. per tant, cerco linia per linia guarrament

egrep -h " u[$pattern]   " ${fit[$index]} | awk '{ print $6 }' >/tmp/quota.1
for usuari in `cat /tmp/quota.1`
do
   if [ $output = html ]
   then
      echo -n "<tr><td align=left><b>$usuari</b></td>"
   else
      printf "%-12s" $usuari
   fi
   i=1
   while (( i <= index ))
   do
      ocupacio=`grep -h " $usuari$" ${fit[$i]} | awk -v camp=$camp '{ print $camp }'`
      let tmp=`echo $ocupacio | tr -d "."`+${total[$i]}
      total[$i]=$tmp # ho he de fer aixi perque no m'agafa arrays a l'esquerra d'un let
      if [ $output = html ]
      then
         echo -n "<td align=right>&nbsp;&nbsp;$ocupacio&nbsp;&nbsp;</td>"
      else
         printf "%10s" $ocupacio
      fi
      if (( i != index ))
      then
         if [ $output != html ]
         then
            if [ "$ocupacio" != "" ]
            then
               echo -n " -> "
            else
               echo -n "    "
            fi
         fi
      fi
      let i=i+1
   done
   if [ $output = html ]
   then
      echo "</td></tr>"
   else
      echo
   fi
done

if [ $output = html ]
then
   echo "<tr><td></td>"
else
   echo -n "            "
   i=1
   while (( i <= index ))
   do
      echo -n "----------"
      if (( i != index ))
      then
         echo -n "    "
      fi
      let i=i+1
   done
   echo
   echo -n "            "
fi

i=1
while (( i <= index ))
do
   tmp=${total[$i]}
   let tmp=tmp/1024/1024
   if [ $output = html ]
   then
      echo "<td align=right><b>$tmp Gb&nbsp;&nbsp;</b></td>"
   else
      printf "%7s Gb" $tmp
      if (( i != index ))
      then
         echo -n "    "
      fi
   fi
   let i=i+1
done

if [ $output = html ]
then
   echo "</tr></table></body></html>"
else
   echo
fi

rm /tmp/quota.1
