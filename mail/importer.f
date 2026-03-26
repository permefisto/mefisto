      SUBROUTINE IMPORTER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPORTER un MAILLAGE d'un PLSV DEFINI par UN FICHIER ASCII
C ----  ISSU de Mefisto ou d'un LOGICIEL de CREATION d'un FICHIER.obj
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            AVRIL 2020
C2345X7..............................................................012

C     LECTURE DU MOT-CLE OU OPTION A EXECUTER PARMI TOUTES LES OPTIONS
C     CONTENUES DANS LE MENU DE NOM 'importer'
C     ----------------------------------------------------------------
      CALL LIMTCL( 'importer', NMTCL )

C     TRAITEMENT DE L'OPTION NMTCL
      IF( NMTCL .LT.  0 ) GOTO 9999

      GOTO( 10, 20, 30, 40, 9999 ), NMTCL

C     un FICHIER xyznsef.plsv.nom d'UN PLSV -> Mefisto
 10   CALL IMPXYZNSEF
      GOTO 9999

C     un FICHIER .XyzNoEF issu de NEF -> Mefisto
 20   CALL NEF2MEF
      GOTO 9999

C     un FICHIER.obj d'UNE SURFACE -> Mefisto
 30   CALL FOBJ2MEF
      GOTO 9999

C     2 FICHIERS NomSurface.xyzlsd et NomSurface.nselsd -> Mefisto
 40   CALL IMPLSDYNA

C     Sortie
 9999 RETURN
      END


