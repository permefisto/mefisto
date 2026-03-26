#!/bin/sh
   echo
   echo "Modifications X11T => xvue"
   echo  ==========================
   echo nom de la bibliotheque derriere nef?
   read dir
   labibl=$mefisto/$dir 
   echo "Modifications X11T => xvue dans $labibl"
   cd  $labibl
#
#  destruction des fichiers .bak du repertoire
   rm *.bak   
#
#  destruction eventuelle de MODIFS
   rm  -R  MODIFS
#
   echo  
   echo les fichiers initiaux de $labibl
   echo =================================================================  
   ls $labibl
#
#  creation de MODIFS
   mkdir MODIFS
#
#  la liste des fichiers presentes un par ligne
   rm LISTFIC
   ls -1 >LISTFIC
#
#  modifications des fichiers .f 
   $mefisto/bin/ajax
#
#  destruction de la liste des fichiers
   rm  LISTFIC
#
#  destruction des .f modifies dans la bibliotheque initiale
#  rm *.f
#
#  copie de MODIFS dans la bibliotheque initiale
   cp MODIFS/*.f .
#
#  destruction de MODIFS 
   rm -R MODIFS

   echo  
   echo les fichiers finaux de $labibl
   echo ====================================================================
   ls $labibl
