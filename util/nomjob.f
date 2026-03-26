      SUBROUTINE NOMJOB( EXIST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RETROUVER LE NOM DU PROJET
C -----   EN LISANT LE SUPER-FICHIER 'ms10' DU REPERTOIRE COURANT
C
C SORTIES :
C ---------
C EXIST : TRUE SI LE FICHIER ms10 EST DANS LE REPERTOIRE COURANT
C              ET LE COMMON / MSNOM / CONTIENT ALORS LE NOM DU PROJET
C         FALSE SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  JUILLET 1989
C.......................................................................
      include"./incl/pp.inc"
      include"./incl/msvaau.inc"
      include"./incl/nmproj.inc"
      include"./incl/typnoobj.inc"
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      LOGICAL           EXIST
C
C     PAS D'OBJET EN TRAITEMENT
      NUMTYPOBJ = 0
      NUMOBJLX  = 0

C     NOMBRE DE MOTS D'UNE PAGE DU SUPER-FICHIER MOPASF
      MOPASF = 128
C
C     INITIALISATION DES VARIABLES DEPENDANT MACHINE
C     ==============================================
      CALL MSINIT
C
C     ESSAI D'OUVERTURE DU SUPER-FICHIER SF
C     =====================================
      CALL TRUNIT( NOFISF )
C
C     NOM DU FICHIER 'ms10' DEPENDANT MACHINE
C
      INQUIRE ( FILE='ms10' , EXIST=EXIST )
      IF( .NOT. EXIST ) GOTO 9000
C
      OPEN( UNIT=NOFISF , ERR=9000 , STATUS='OLD' ,
     %      FILE='ms10' , ACCESS='DIRECT' ,
     %      FORM='UNFORMATTED' , IOSTAT=IOERR ,
     %      RECL=MOPASF * NBCHMO )
C
C     LECTURE DU NOM DU PROJET DE LA MS
C     =================================
      NPA = 2
      READ(  UNIT=NOFISF , REC=NPA ) NUEXEC,NMPROJ
C
C     FERMETURE DU FICHIER
C     ====================
      CLOSE( UNIT=NOFISF , STATUS='KEEP' )
      EXIST = .TRUE.
      RETURN
C
C     ERREUR
 9000 EXIST = .FALSE.
      RETURN
      END
