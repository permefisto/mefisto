      SUBROUTINE EXPORTER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  EXPORTER un MAILLAGE Mefisto d'un PLSVO sous forme d'UN FICHIER
C -----  ADAPTE a DIFFERENTS LOGICIELS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C MODIFS : Alain PERRONNET  Saint PIERRE du PERRAY         Decembre 2024
C2345X7..............................................................012

C     LECTURE DU MOT-CLE OU OPTION A EXECUTER PARMI TOUTES LES OPTIONS
C     CONTENUES DANS LE MENU DE NOM 'exporter'
C     ----------------------------------------------------------------
 1    CALL LIMTCL( 'exporter', NMTCL )

C     TRAITEMENT DE L'OPTION NMTCL
      IF( NMTCL .LE. 0 ) GOTO 9999

      GOTO( 10, 20, 30, 40, 50, 60, 70, 80, 90, 9999 ), NMTCL

C     FICHIER MAILLAGE des TMS XYZSOMMET NSEF d'UN PLSV -> MEFISTO
 10   CALL EXPXYZNSEF
      GOTO 1

C     FICHIER MAILLAGE XYZNOEUD NPEF d'un OBJET -> MEFISTO
 20   CALL XYZNPE
      GOTO 1

C     FICHIER du MAILLAGE d'UN OBJET -> NEKTON
 30   CALL NEKTON
      GOTO 1

C     FICHIERS du MAILLAGE d'UN OBJET -> FIDAP
 40   CALL FIDAP
      GOTO 1

C     FICHIERS du MAILLAGE d'UN OBJET -> FREEFEM
 50   CALL FREEFEM
      GOTO 1

C     FICHIERS du MAILLAGE d'UN OBJET -> OPENFOAM
 60   CALL OPENFOAM
      GOTO 1

C     FICHIER MAILLAGE nopomodulef.NomObjet d'UN OBJET -> MODULEF
 70   CALL NOPOMODULEF
      GOTO 1

C     FICHIER MAILLAGE d'une SURFACE + un VOLUME -> LS-DYNA
 80   CALL LSDYNA
      GOTO 1

C     FICHIER.stl de la TRIANGULATION d'une SURFACE FERMEE -> Imprimante 3D
 90   CALL STLXYZNSEF
      GOTO 1

 9999 RETURN
      END
