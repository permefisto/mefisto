      SUBROUTINE INITIA( EXIST )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISER UN ENCHAINEMENT MEFISTO
C -----
C
C ENTREE :
C --------
C EXIST : TRUE LE NOM DU PROJET EST DANS COMMON /MSNOM/
C         FALSE SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1984
C ......................................................................
      PARAMETER         (NBLPDC=10,NBEPDC=201)
      include"./incl/lu.inc"
      include"./incl/pp.inc"
      include"./incl/ppmck.inc"
      COMMON             MCN(MOTMCN)
      include"./incl/nbcamo.inc"
      include"./incl/homdir.inc"
      include"./incl/nmproj.inc"
      include"./incl/pilect.inc"
      include"./incl/impres.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ctemps.inc"
      include"./incl/langue.inc"
      include"./incl/darete.inc"
C
      COMMON /DOCUME/    INDOCU(NBEPDC)
      COMMON /MSNMFI/    NMFIC
      DOUBLE PRECISION   CPUOLD
      COMMON /TEMCPU/    CPUOLD
      COMMON /EPSSSS/    EPZERO,EPSXYZ
      COMMON /TRAVA1/    NTITRE(20),NDATE(2),NOMCRE(6)
      COMMON /UNITES/LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NOAFTS,NUNIT(26)
      CHARACTER*160      KNOM
      CHARACTER*28       NMVERS
      CHARACTER*16       NMFIC
C
      CHARACTER*(NC1MCK) KBLANC
      CHARACTER          KDATE*80,KCREAT*80,KINFO*80,KHEURE*80
      LOGICAL            EXIST
C
      DATA               KBLANC
     %/'                                                '/
 1000 FORMAT(/
     &'---------------------------------------------------------'/
     &'| MM   MM  EEEEEE  FFFFFF  IIII  SSSSSS  TTTTTT OOOOOOO |'/
     &'| MMM MMM  EE      FF       II   SS        TT   OO   OO |'/
     &'| MM M MM  EEEEE   FFFF     II   SSSSSS    TT   OO   OO |',
     &'  FAIT L''EF !'/
     &'| MM   MM  EE      FF       II       SS    TT   OO   OO |'/
     &'| MM   MM  EEEEEE  FF      IIII  SSSSSS    TT   OOOOOOO | ',
     & A /
     &'---------------------------------------------------------'/)
C
C     RECONNAISSANCE DE LA LANGUE DES DONNEES DE MEFISTO
      CALL LANGUE
C
C     MISE A BLANC DU SUPER-TABLEAU DE CARACTERES MCK
      DO 5 I=1,MOTMCK
         MCK(I) = KBLANC
 5    CONTINUE
C
C     LECTEUR
      LECTEU = IINFO('LECTEUR INITIAL')
      LPLECT(1) = LECTEU
C
C     IMPRIMANTE
      IMPRIM = IINFO('IMPRIMANTE INITIALE')
C
C     NUMERO D'AFFICHAGE DES TMS : ICI IL N'EST PAS FORCE
      NOAFTS = 0
C
C     LA BANNIERE ET LE NOM DE LA VERSION
      CALL VRSION( NMVERS )
      I = NUDCNB( NMVERS )
      WRITE (IMPRIM,1000) NMVERS(1:I)
C
C     LE PARAMETRE D'IMPRESSION
      IMPRE  = 0
C
C     ENTREE DU TITRE DU TRAVAIL
      IF( .NOT. EXIST ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10005)
         ELSE
            WRITE(IMPRIM,20005)
         ENDIF
         READ (LECTEU,'(A)') NMPROJ
      ENDIF
10005 FORMAT('nom (en MINUSCULES) du projet?')
20005 FORMAT('Project (low case) name?')
C
C     INITIALISATION DU TEMPS CPU
      CPUOLD = 0D0
C
C     LE NOM PARAMETRE POUR LES FICHIERS MS
      J=0
      DO 6 I=1,72
         IF( NMPROJ(I:I) .NE. ' ' ) THEN
            J = J + 1
            NMFIC(J:J) = NMPROJ(I:I)
            IF( J .GE. 15 ) GOTO 8
         ENDIF
 6    CONTINUE
C
 8    J = J + 1
      NMFIC(J:J) = '.'
C
C     LA DOCUMENTATION ET LE JOURNAL
      NFDOCU = 8
      I      = NUDCNB( HOMDIR )
      KNOM   = HOMDIR(1:I) // '/doc/' // 'journal'
      OPEN( FILE=KNOM, UNIT=NFDOCU , STATUS='OLD' , IOSTAT=IOERR )
      IF( IOERR .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'PROBLEME POUR OUVRIR LE FICHIER ',KNOM
         ELSE
            WRITE(IMPRIM,*) 'PROBLEM TO OPEN THE FILE ',KNOM
         ENDIF
         STOP
      ENDIF
      DO 30 I = 1 , MXKLG
         READ( UNIT=NFDOCU, IOSTAT=IOERR ) KLG(I)
         IF( IOERR .NE. 0 ) GOTO 40
         WRITE(IMPRIM,*) KLG(I)
 30   CONTINUE
C
C     LES EPS POUR IDENTIFIER LES XYZ DES SOMMETS
 40   EPZERO = 1E-6
      EPSXYZ = 1E-4
C
C     DATE
      KDATE  = KINFO( 'DATE' )
      KHEURE = KINFO( 'HEURE' )
      IF( LANGAG .EQ. 0 ) THEN
         NDATE(1) = ICHARX(KDATE(10:11) // '/' // KDATE(13:13))
         NDATE(2) = ICHARX(KDATE(14:14) // '/' // KDATE(25:26))
      ELSE
         NDATE(1) = ICHARX(KDATE(13:14) // '/' // KDATE(10:10))
         NDATE(2) = ICHARX(KDATE(11:11) // '/' // KDATE(25:26))
      ENDIF
      WRITE (IMPRIM,*)'DATE         : ',KDATE(1:8),'   ',
     %                 KHEURE(1:2) // 'H ' // KHEURE(3:4) // 'M ' //
     %                 KHEURE(5:6) // 'S '
C
C     NOM DE L'UTILISATEUR
      KCREAT = KINFO('UTILISATEUR')
      DO 100 I = 1,6
         NOMCRE(I) = ICHARX(KCREAT(4*I-3:4*I))
100   CONTINUE
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,'('' AUTEUR       : '',A)') KCREAT(1:24)
         WRITE (IMPRIM,'('' PROJET       : '',A)') NMPROJ
      ELSE
         WRITE (IMPRIM,'('' AUTHOR       : '',A)') KCREAT(1:24)
         WRITE (IMPRIM,'('' PROJECT      : '',A)') NMPROJ
      ENDIF

C     LECTURE DU DICTIONNAIRE DES TABLEAUX DESCRIPTEURS
      CALL LEDICO
      IF( LANGAG .EQ. 0 ) THEN

         WRITE (IMPRIM,10900)
10900 FORMAT(/'DONNEES de MEFISTO en LANGUE FRANCAISE'/
     %       'Tapez ? pour obtenir la DOCUMENTATION'/
     %       'Tapez @ pour ABANDONNER l''entree d''une DONNEE')
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,10998) EPZERO
         WRITE(IMPRIM,10999) EPSXYZ
10998 FORMAT('PRECISION Autour de l''ORIGINE EPZERO =',G15.7)
10999 FORMAT('PRECISION Loin   de l''ORIGINE EPSXYZ =',G15.7)

      ELSEIF( LANGAG .EQ. 1 ) THEN

         WRITE (IMPRIM,20900)
20900 FORMAT(/'INPUT DATA of MEFISTO in ENGLISH LANGUAGE'/
     %       'Type ? to obtain the DOCUMENTATION'/
     %       'Type @ to ESCAPE a MENU or INPUT DATA')
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,20998) EPZERO
         WRITE(IMPRIM,20999) EPSXYZ
20998 FORMAT('PRECISION     AROUND ZERO EPZERO =',G15.7)
20999 FORMAT('PRECISION NOT AROUND ZERO EPSXYZ =',G15.7)

      ENDIF
      WRITE(IMPRIM,*)

C     IMPRESSION DU NOM DES TS
      IMNMTS = 1
C
C     LES TEMPS ACTUEL, INITIAL et FINAL
      TEMPS    = 0
      TEMPSINI = 0
      TEMPSFIN = 0
C
      RETURN
      END
