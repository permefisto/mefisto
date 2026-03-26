#!/bin/bash
#
#  Alain Perronnet                       Version Fevrier 2022
echo
echo ============================================================================
echo "Execution du fichier install.bash des EXECUTABLES de Mefisto dans"
pwd
echo ============================================================================
echo
echo "1) Choix du langage pour les menus de Mefisto"
echo "English ou Francais (E/F)"
read bidon
if [ $bidon = E ] || [ $bidon = e ]
then
  LANGAGE=1
  echo "ENGLISH is the language of Mefisto"
  echo "If problem during installation, send me the listing of this execution"
  echo "at perronnet@ann.jussieu.fr"
  echo
  echo "Attention: To install and compile the SOURCES of Mefisto"
  echo "see http://www.ljll.math.upmc.fr/sourceload.php"
else
  LANGAGE=0
  echo "Le FRANCAIS est le langage de Mefisto"
  echo "Si probleme durant l installation, envoyer le listage de cette execution"
  echo "a perronnet@ann.jussieu.fr"
  echo
  echo "Attention: Pour installer et compiler les SOURCES de Mefisto"
  echo "voir http://www.ljll.math.upmc.fr/sourceload.php"
fi

echo
if [ $LANGAGE = 0 ]
then
  echo "Installation des modules executables du logiciel MEFISTO"
  echo "Nom du repertoire contenant le fichier Mefisto.tgz"
else
  echo "Installation of MEFISTO load-modules"
  echo "Name of directory where is the file Mefisto.tgz"
fi
pwd
REPINIT=${PWD#}
echo "Repertoire courant actuel: " $REPINIT
ls -l

NomMachine=0
NomBin=0

if test -f mefisto.lnx64.tgz
then
  if [ $LANGAGE = 0 ]
  then
    echo "Installation de la VERSION EXECUTABLE LINUX Processeur 64bits"
  else
    echo "Installation of LINUX EXECUTABLE VERSION Processor 64bits"
  fi
  cp -p mefisto.lnx64.tgz  mefisto.lnx.tgz
  NomMachine=lnx64
  NomBin=lnx64
fi

if [ $NomMachine = 0 ]
then
  if [ $LANGAGE = 0 ]
  then
    echo
    echo "MEFISTO: Version INCONNUE"
    echo "Dans ce repertoire, PAS de fichier mefisto.NomMachine.tgz"
    echo "Cf https://www.ljll.math.upmc.fr/perronnet/mefisto.charger.php"
  else
    echo
    echo "MEFISTO: UNKNOWN Version"
    echo "In this directory, NO file mefisto.MachineName.tgz"
    echo "Cf https://www.ljll.math.upmc.fr/perronnet/mefistoa.charger.php"
  fi
  exit
fi

echo
echo
if [ $LANGAGE = 0 ]
then
  echo "Installation du fichier mefisto.$NomMachine.tgz"
else
  echo "Installation of file mefisto.$NomMachine.tgz"
fi


if test ! -f mefisto.$NomMachine.tgz
then
  if [ $LANGAGE = 0 ]
  then
    ls
    echo "Le fichier mefisto.$NomMachine.tgz n'est pas disponible dans ce repertoire"
    echo "Changer de repertoire ou venir l'y mettre et relancer"
    exit
  else
    ls
    echo "UNKNOWN file mefisto.$NomMachine.tgz  NOT in this directory"
    echo "Change the directory to this one where mefisto.$NomMachine.tgz is"
    exit
  fi
fi

REPMEFISTO=$REPINIT/mefisto
mkdir mefisto
cd mefisto
tar xzvf ../mefisto.$NomMachine.tgz

echo
if [ $LANGAGE = 0 ]
then
  echo "Le nom du repertoire mefisto est $REPMEFISTO qui contient"
else
  echo "The name of MEFISTO directory is $REPMEFISTO which contains"
fi
ls -l $REPMEFISTO

echo
echo
echo ============================================================================
#Le repertoire de Mefisto devient /usr/local/mefisto  au lien pres ...
export MEFISTO=/usr/local/mefisto
if [ $LANGAGE = 0 ]
then
  echo "2) Attention: Tentative de lien entre $REPMEFISTO et $MEFISTO"
  echo "Si un probleme survient ensuite, c'est sans doute du au fait que"
  echo "les droits d'acces de /usr/local ne vous donnent pas la permission en ecriture"
  echo "et que vous n etes pas le SUPER-UTILISATEUR root seul habilite a le faire!"
else
  echo "2) If a problem appears when linking between $REPMEFISTO and $MEFISTO"
  echo "May be, you are NOT the super-user or modify the WRITE ACCESS RIGHT of $MEFISTO"
fi

if test -d $MEFISTO
then
# Le repertoire $MEFISTO existe . Est il $mefisto?
  ls -l $MEFISTO
  if [ $LANGAGE = 0 ]
  then
    echo "En lisant la ligne ci-dessus verifiez que"
    echo "/usr/local/mefisto -> $REPMEFISTO"
    echo "est le bon lien?"
    echo
    echo "Taper un caractere pour continuer"
  else
    echo "Read above and verify"
    echo "/usr/local/mefisto -> $REPMEFISTO"
    echo " is OK?"
    echo
    echo "Type a character to continue"
  fi
else
# Le repertoire $MEFISTO n'existe pas. Le lien est cree
  echo "sudo ln -s $REPMEFISTO /usr/local/mefisto"
  sudo ln -s $REPMEFISTO /usr/local/mefisto
# Le repertoire de Mefisto devient /usr/local/mefisto  au lien pres ...
  if test ! -d $MEFISTO
  then
    if [ $LANGAGE = 0 ]
    then
      echo "Lien symbolique NON REALISE de $mefisto vers /usr/local/mefisto"
      echo "Revoyez les droits d'acces ou travaillez en tant que SUPER-UTILISATEUR root"
      echo "seul habilite a realiser ce lien symbolique"
    else
      echo "Problem: You are NOT the SUPER-USER or modify the WRITE ACCESS RIGHT of $MEFISTO"
    fi 
    read bidon
    exit
  fi
fi

echo ============================================================================
#Le repertoire d'execution des projets Mefisto
export MEFISTOX=$HOME/mefistox
if test ! -d $HOME/mefistox
then
  if [ $LANGAGE = 0 ]
  then
    echo "3) Creation du repertoire $MEFISTOX pour contenir les projets utilisateur"
  else
    echo "3) Creation of directory $MEFISTOX to contain the user's projects"
  fi
  mkdir $HOME/mefistox
else
  if [ $LANGAGE = 0 ]
  then
    echo "3) Le repertoire $MEFISTOX pour contenir les projets utilisateur existe"
  else
    echo "3) The directory $MEFISTOX to contain the user's projects exists"
  fi
fi

echo
echo
echo ============================================================================
if [ $LANGAGE = 0 ]
then
  echo "4) Le langage de Mefisto est le FRANCAIS"
  cd $MEFISTO/td
  rm -rf d;  rm i;  rm -rf m;
  cp -r df d
  cp    if i
  cp -r mf m
  cd $MEFISTO
  rm -rf doc
  cp -r  docf doc
else
  echo "4) The language of Mefisto is ENGLISH"
  cd $MEFISTO/td
  rm -rf d;  rm i;  rm -rf m;
  cp -r da d
  cp    ia i
  cp -r ma m
  cd $MEFISTO
  rm -rf doc
  cp -r  doca doc
fi

echo
echo
echo ============================================================================
if [ $LANGAGE = 0 ]
then
  echo "5) Le repertoire bin.$NomBin devient bin"
else
  echo "5) bin.$NomBin becomes bin"
fi
cp -pr  bin.$NomBin  bin


echo
echo
echo ============================================================================
if [ $LANGAGE = 0 ]
then
  echo "6) Le repertoire pp.$NomMachine devient pp"
else
  echo "6) The directory pp.$NomMachine becomes pp"
fi
cp -pr  pp.$NomMachine  pp

echo
echo ===============================================================================================
if [ $LANGAGE = 0 ]
then
  echo "7) Ajout des variables d environnement MEFISTO MEFISTOX dans ~/.bashrc"
else
  echo "7) Environment Variables MEFISTO MEFISTOX are ADDED into ~/.bashrc"
fi

export MEFISTO=$REPMEFISTO
export MEFISTOX=$HOME/mefistox
export PATH=.:$PATH:$MEFISTO/bin
export CDPATH=.:$HOME:$MEFISTO:$MEFISTOX
alias ll='ls -l'
alias lrt='ls -rtl'


cat <<@ >>ba
# 
#added during the execution of instalsource.bash
export MEFISTO=$REPMEFISTO
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
echo ============================================================================
if [ $LANGAGE = 0 ]
then
  echo "8) Test de l installation de MEFISTO"
  echo "Contenu actuel du repertoire $MEFISTO"
else
  echo "8) Test of the MEFISTO installation"
  echo "Contents of the directory $MEFISTO"
fi
ls -l $MEFISTO

echo
if [ $LANGAGE = 0 ]
then
  echo "Taper un caractere pour continuer"
else
  echo "Type a character to continue"
fi
read bidon

echo
if [ $LANGAGE = 0 ]
then
  echo "DANS UNE AUTRE FENETRE xterm, TESTER l'INSTALLATION de Mefisto"
  echo "Pour cela dans cette nouvelle fenetre, taper /usr/local/mefisto/bin/INITIER"
  echo "Alors, il doit apparaitre:"
  echo "=> MEFISTO-INITIER: nom du projet ?"
  echo
  echo "Si ce texte n'apparait pas"
  echo "alors, verifier que \$PATH contient \$MEFISTO/bin"
  echo "sinon, taper 1 fois un nom (en minuscules) de projet"
  echo
  echo "Il apparait"
  echo 
  echo "=> Execution MEFISTO-INITIER dans \$MEFISTOX/NomProjet"
  echo 
  echo "---------------------------------------------------------"
  echo "| MM   MM  EEEEEE  FFFFFF  IIII  SSSSSS  TTTTTT OOOOOOO |"
  echo "| MMM MMM  EE      FF       II   SS        TT   OO   OO |"
  echo "| MM M MM  EEEEE   FFFF     II   SSSSSS    TT   OO   OO |  FAIT L'EF !"
  echo "| MM   MM  EE      FF       II       SS    TT   OO   OO |"
  echo "| MM   MM  EEEEEE  FF      IIII  SSSSSS    TT   OOOOOOO | Version Avril 2012"
  echo "---------------------------------------------------------"
  echo
  echo "=>  Tapez ? pour obtenir la DOCUMENTATION"
  echo "=>  Tapez @ pour ABANDONNER un menu"
  echo 
  echo "=>  NOM DU PROJET ?"
  echo
  echo "=> Taper encore 1 fois le nom du projet"
  echo
  echo
  echo "=> Des informations sur les fichiers numeriques du repertoire"
  echo "=> $MEFISTOX/NomProjet apparaissent"
  echo
  echo "=> Taper ls \$MEFISTOX/NomProjet pour vous en convaincre"
  echo
  echo
  echo "Taper un caractere pour continuer"
  read bidon
  echo
  echo "Taper /usr/local/mefisto/bin/MAILLER"
  echo "=> MEFISTO-MAILLAGE : nom du projet ?   apparait"
  echo "Taper 1 fois le nom du projet"
  echo
  echo "=> Une fenetre s'ouvre avec la banniere Mefisto"
  echo
  echo "Cliquer dans la fenetre"
  echo "ou taper un caractere au clavier"
  echo "Le menu Debut apparait"
  echo 
  echo "Choisissez votre option dans les menus"
  echo "BON TRAVAIL avec Mefisto ..."
  echo
  echo "Pour terminer et SAUVEGARDER le travail fait, revenir au menu debut"
  echo "POUR TERMINER taper 99;"
  echo "et SURTOUT PAS Control C car vous PERDRIEZ tout ce que vous venez de calculer!"
  echo
  echo "=================================================================================="
  echo "Fin de l'installation de MEFISTO Version $NomMachine dans $REPMEFISTO"
  echo "=================================================================================="
else
  echo "IN AN OTHER WINDOW xterm, Test the Mefisto installation"
  echo "Type /usr/local/mefisto/bin/INITIER"
  echo "Then, it appears:"
  echo "MEFISTO-INITIER: Project (low case) name ?"
  echo
  echo "If this text does not appear"
  echo "then, Verify \$PATH contains \$MEFISTO/bin"
  echo "else, Type once the project name"
  echo
  echo "Execution MEFISTO INITIER in directory $MEFISTOX/ProjectName/ "
  echo "---------------------------------------------------------"
  echo "| MM   MM  EEEEEE  FFFFFF  IIII  SSSSSS  TTTTTT OOOOOOO |"
  echo "| MMM MMM  EE      FF       II   SS        TT   OO   OO |"
  echo "| MM M MM  EEEEE   FFFF     II   SSSSSS    TT   OO   OO |  FAIT L EF !"
  echo "| MM   MM  EE      FF       II       SS    TT   OO   OO |"
  echo "| MM   MM  EEEEEE  FF      IIII  SSSSSS    TT   OOOOOOO | Version February 2022"
  echo "---------------------------------------------------------"
  echo
  echo "Project (low case) name ?"
  echo "Type an other time the same name of project"
  echo
  echo "=> Some informations on the numerical files of this project appear"
  echo "=> Type ls \$MEFISTOX/ProjectName to see them"
  echo
  echo "Type a character to continue"
  read bidon
  echo
  echo
  echo "Type  $MEFISTO/bin/MAILLER"
  echo "it appears"
  echo "==============================================================="
  echo "   MEFISTO-MESHER: MESH an OBJECT on LINUX PC"
  echo "==============================================================="
  echo "Project (low case) name ?"
  echo "Type once the project name"
  echo
  echo "=> A window opens with the logo Mefisto"
  echo
  echo "Click in the window"
  echo "or type a character on keyboard"
  echo "The menu Debut appears"
  echo 
  echo "Choose an option in the menus ..."
  echo "GOOD WORK with Mefisto ..."
  echo
  echo "To end and SAVE the work, go back to the Debut menu"
  echo "To FINISH type 99;  to keep all you have done!"
  echo
  echo ==================================================================================
  echo "End of the installation of MEFISTO Version $NomMachine in $REPMEFISTO"
  echo ==================================================================================
fi
echo
echo
