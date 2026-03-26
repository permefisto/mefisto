#!/bin/bash
#
#  Alain Perronnet                                            janvier 2022
echo
echo ===============================================================================================
echo MEFISTO INSTALLATION with 8 steps
echo ===============================================================================================
echo
INITDIR=${PWD#}
echo "1) instalsource.bash install mefistosource.tgz in $INITDIR"
echo ===============================================================================================
if test -f mefistosource.tgz
then
    echo "The file mefistosource.tgz is in $INITDIR"
    ls -l mefistosource.tgz
    echo
    echo "tar -xzf $INITDIR/mefistosource.tgz"
    echo
    tar -xzf mefistosource.tgz
else
    echo
    echo "ls -l $INITDIR"
    ls -l $INIDIR
    echo
    echo "Attention: The file mefistosource.tgz IS NOT in $INITDIR"
    echo "Type  cd mefistosource.tgzDirectoryName and execute again instalsource.bash"
    exit
fi

echo
echo ===============================================================================================
echo "2) Choice of the Mefisto MENU LANGUAGE: English or French ?"
echo "   Type E or F"
read mot
if [ $mot = E ] || [ $mot = e ]
then
  LANGAGE=1
  echo "ENGLISH is the language of Mefisto"
  echo "If problem during installation, send me the listing of this execution"
  echo "at perronnet@ann.jussieu.fr"
else
  LANGAGE=0
  echo "FRANCAIS est le langage de Mefisto"
  echo "Si probleme durant l installation, envoyer le listage de cette execution"
  echo "a perronnet@ann.jussieu.fr"
fi
echo
echo  "MEFISTO=$INITDIR/mefistosource"
export MEFISTO=$INITDIR/mefistosource
if test -d $MEFISTO
then
# Le repertoire $MEFISTO EXISTE
  if [ $LANGAGE = 0 ]
  then
    echo "$MEFISTO contient initialement"
  else
    echo "$MEFISTO contains at the starting"
  fi
  ls $MEFISTO
else
# Le repertoire $MEFISTO N'EXISTE PAS
  if [ $LANGAGE = 0 ]
  then
    echo "$MEFISTO NE CONTIENT PAS LE SOURCE de Mefisto"
    echo "Taper  cd NomRepertoire   contenant mefistosource.tgz"
    echo "puis relancer instalsource.bash a partir de ce repertoire"
  else
    echo "$MEFISTO DOES NOT CONTAIN THE SOURCE of Mefisto"
    echo "Type  cd NameDirectory   which contains mefistosource.tgz"
    echo "Then execute again instalsource.bash from this directory"
  fi
  exit
fi

#calcul du nombre NBCHAR de caracteres de $MEFISTO
cat <<@ >ba
$MEFISTO
@
#nombre de caracteres de $MEFISTO suivi de ' ba'
wc -m ba >nbc
#suppression de ' ba' dans nbc
ed nbc <<@
s/ ba//
w
q
@
#NBCHAR=variable contenant le nombre de caracteres de $MEFISTO
NBCHAR=`cat nbc`
NBCHAR=` expr $NBCHAR - 1 `
#construction du fichier $MEFISTO/incl/homdir.inc
cat <<@ >$MEFISTO/incl/homdir.inc
C     $MEFISTO/incl/homdir.inc
      CHARACTER*$NBCHAR HOMDIR
      PARAMETER (HOMDIR=
     %'$MEFISTO'
     %)
@
echo Creation $MEFISTO/incl/homdir.inc
cat $MEFISTO/incl/homdir.inc
#destruction des fichiers inutiles
rm ba nbc

echo
echo ===============================================================================================
NomMachine=lnx64
if [ $Machine = Lnx64 ] || [ $Machine = lnx64 ]
then
  NomMachine=lnx64
  if [ $LANGAGE = 0 ]
  then
    echo "Installation de Mefisto sous LINUX 64bits avec les compilateurs gcc gfortran"
  else
    echo "Installation of Mefisto under 64bits LINUX with gcc gfortran compilers"
  fi
fi

if [ NomMachine = 0 ]
then
  if [ $LANGAGE = 0 ]
  then
      echo "Desole: SYSTEME INCONNU pour installer les sources de Mefisto"
  else
      echo "Sorry: UNKNOWN SYSTEM to install the sources of Mefisto"
  fi
  exit
fi

echo
echo ===============================================================================================
if [ $LANGAGE = 0 ]
then
  echo "4) Installation des LIBRAIRIES version en FRANCAIS"
  cd $MEFISTO/td
  rm -rf d;  rm i;  rm -rf m;
  cp -pr df d
  cp -p  if i
  cp -pr mf m
  cd $MEFISTO
  rm -rf doc
  cp -pr docf  doc
else
  echo "4) Installation of LIBRARIES ENGLISH version"
  cd $MEFISTO/td
  rm -rf d;  rm -rf i;  rm -rf m;
  cp -pr da d
  cp -p  ia i
  cp -pr ma m
  cd $MEFISTO
  rm -rf doc
  cp -r  doca  doc
fi
echo "cp -pr  bin.$NomMachine  bin"
cp -pr bin.$NomMachine  bin

echo
if [ $LANGAGE = 0 ]
then
  echo "$MEFISTO contient maintenant"
else
  echo "$MEFISTO contains now"
fi
ls $MEFISTO

echo
echo ===============================================================================================
if [ $LANGAGE = 0 ]
then
  echo "5) ATTENDEZ PATIEMMENT la fin des compilations, editions de liens, creations des librairies"
else
  echo "5) WAIT PATIENTLY the end of compilations, bindings and libraries constructions"
fi
./bin/cbl_tout >ListMefisto 2>>WarningMefisto

rm  $MEFISTO/pp/pppoba

echo "ls -l $MEFISTO/pp"
ls -l $MEFISTO/pp

if test -f $MEFISTO/pp/ppflui
then
  if [ $LANGAGE = 0 ]
  then
    echo "Le repertoire $MEFISTO/pp contient"
    echo "ppelas  ppflui  ppinit  ppmail  ppnlse  ppther  pxyz  => Mefisto est CORRECT"
  else
    echo "The directory $MEFISTO/pp contains"
    echo "ppelas  ppflui  ppinit  ppmail  ppnlse  ppther  pxyz  => Mefisto is OK"
  fi
else
  if [ $LANGAGE = 0 ]
  then
    echo "PROBLEME a la COMPILATION: Mefisto est INCORRECT"
    echo "EDITER les fichiers $MEFISTO/ListMefisto and $MEFISTO/WarningMefisto"
  else
    echo "Sorry, PROBLEM of COMPILATION:  Mefisto is INCORRECT"
    echo "EDIT the files $MEFISTO/ListMefisto and $MEFISTO/WarningMefisto"
  fi
  exit
fi

echo
echo ===============================================================================================
#Le repertoire d'execution des projets Mefisto est cree
export MEFISTOX=$HOME/mefistox
if test ! -d $HOME/mefistox
then
  if [ $LANGAGE = 0 ]
  then
    echo "6) Creation du repertoire $MEFISTOX pour contenir les PROJETS UTILISATEUR"
  else
    echo "6) $MEFISTOX directory creation to contain the USER PROJECTS"
  fi
  mkdir $HOME/mefistox
else
  if [ $LANGAGE = 0 ]
  then
    echo "6) Le repertoire $MEFISTOX contiendra les PROJETS UTILISATEUR"
  else
    echo "6) The directory $MEFISTOX will contain the USER'S PROJECTS"
  fi
fi

echo
echo ===============================================================================================
if [ $LANGAGE = 0 ]
then
  echo "7) Ajout des variables d environnement MEFISTO MEFISTOX dans ~/.bashrc"
else
  echo "7) Environment Variables MEFISTO MEFISTOX are ADDED into ~/.bashrc"
fi

cat <<@ >>ba
# 
#added during the execution of instalsource.bash
export MEFISTO=$INITDIR/mefistosource
export MEFISTOX=$HOME/mefistox
export PATH=.:$PATH:$MEFISTO/bin
export CDPATH=.:$HOME:$MEFISTO:$MEFISTOX
alias ll='ls -l'
alias lrt='ls -rtl'
@

cat ba >>~/.bashrc
rm ba

if [ $LANGAGE = 0 ]
then
  echo
  echo "Votre fichier $HOME/.bashrc est maintenant"
  echo .............................................................................................
  cat $HOME/.bashrc
  echo .............................................................................................
else
  echo
  echo "Your file $HOME/.bashrc is now"
  echo .............................................................................................
  cat $HOME/.bashrc
  echo .............................................................................................
fi


echo
echo
echo ===============================================================================================
if test -f $MEFISTO/pp/pxyz
then
  rm ListMefisto
  rm WarningMefisto
  if [ $LANGAGE = 0 ]
  then
    echo
    echo ===============================================================================================
    echo "8) Fin de l'installation de MEFISTO Version $NomMachine dans $MEFISTO"
    echo Vous pouvez taper: INITIER et repondre aux questions ...
    echo Ensuite, vous pouvez taper: MAILLER et repondre aux questions ...
    echo Ensuite, vous pouvez taper: THERMICER or ELASTICER or FLUIDER et repondre aux questions ...
    echo
    echo BIENVENUE et BON TRAVAIL avec Mefisto!
    echo ===============================================================================================
  else
    echo
    echo ===============================================================================================
    echo "8) End of the installation of MEFISTO Version $NomMachine in $MEFISTO"
    echo You may type: INITIER and answer questions ...
    echo After, you may type: MESHER  and answer questions ...
    echo After, you may type: HEATER or ELASTICER or FLUIDER and answer questions ...
    echo
    echo WELCOME and GOOD WORK with Mefisto!
    echo ===============================================================================================
  fi
else
  if [ $LANGAGE = 0 ]
  then
    echo
    echo ===============================================================================================
    echo "FIN INCORRECTE de l installation de Mefisto"
    echo "Le fichier $MEFISTO/pp/pxyz N EXISTE PAS"
    echo "EDITER les fichiers $MEFISTO/ListMefisto and $MEFISTO/WarningMefisto"
    echo ===============================================================================================
  else
    echo ===============================================================================================
    echo "INCORRECT EXIT of Mefisto installation"
    echo "The file $MEFISTO/pp/pxyz DOES NOT EXIST"
    echo "EDIT the files $MEFISTO/ListMefisto and $MEFISTO/WarningMefisto"
    echo ===============================================================================================
  fi
fi
echo
echo
