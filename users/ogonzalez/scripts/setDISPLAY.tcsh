#!/bin/tcsh
#
# per definir facilment el DISPLAY mentre no trobem una alternativa millor
#
# sense parametres, defineix la nostra maquina de sempre com a DISPLAY
# amb parametres, agafa la maquina de $1 i la defineix com a DISPLAY
#
# per tal que quedi definit, s'ha de fer un source d'aquest script, aixi
# que defineixo un alias a profile per fer-ho. es diu 'd':
#
#	alias d 'source <path>/setDISPLAY'
#

if ( $# == 0 ) then
   set jo=`whoami`
   grep $jo - <<EOF | awk '{ print $1 }' >/tmp/maquina
      apolo     mburset
      korax     ogonzalez oscar
      cel       rguigo roderic
      icarus    gparra
      tmatik    scaste
      basileia  eblanco
      xinxan    jabril
      zapiron   rcastelo
      nautilus  fcamara
      kalamos   malba
EOF
   set maquina=`cat /tmp/maquina`
   set monitor=:0.0
else
   set maquina=`echo $1 | sed -e "s%:.*%%"`
   set monitor=`echo $1 | sed -e "s%$maquina%%"`
   if ( xx$monitor == xx ) then
      set monitor=:0.0
   endif
endif


set DISPLAY=$maquina$monitor
unset maquina monitor

