#!/bin/bash
#  lecture et mise sous le directory 
   echo
   echo MODIFICATION DE TEXTES DANS LES LIBRAIRIES MEFISTO
   echo -n nom de la bibliotheque?
   read dir
   cd  $mefisto/$dir/f
#
#  destruction des fichiers .bak
   rm    *.bak  
#
#  destruction eventuelle de f1
   rm    ../f1/*
   rmdir ../f1
#
#  creation de f1
   mkdir ../f1
#
#  la liste des fichiers 1 par ligne
   ls -1 >../bibftn 
#
#  traitement 
   ~/bin/modiftexte
#
#  destruction de la liste des fichiers
   rm ../bibftn 
#
#  changement des noms f -> f0  et f1 -> f
   chn $mefisto/$dir/f  $mefisto/$dir/f0
   chn $mefisto/$dir/f1 $mefisto/$dir/f
#
#  le resultat   
   echo  
   echo les fichiers de $dir/f0 initial
   echo ===============================  
   ls $mefisto/$dir/f0     
   echo  
   echo les fichiers de $dir/f   final
   echo ==============================  
   ls $mefisto/$dir/f

