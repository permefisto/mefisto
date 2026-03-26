      SUBROUTINE FUSUFX( NOM , SUFFIX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AJOUTER LE SUFFIXE AU NOM A PARTIR DU 1-ER BLANC RENCONTRE
C ----- OU A LA PLACE DU DERNIER SUFFIXE DEJA EXISTANT DANS NOM
C
C ENTREES :
C ---------
C NOM    : LE NOM
C SUFFIX : LE SUFFIXE
C          SI LE SUFFIXE COMMENCE PAR ' ' IL N'EST PAS AJOUTE AU NOM
C
C SORTIES :
C ---------
C NOM    : NOM // SUFFIXE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    MARS 1989
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES /  LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL               RMCN(1)
      EQUIVALENCE       (RMCN(1),MCN(1))
      CHARACTER*(*)      NOM,SUFFIX
      CHARACTER*(MXSUFX) VASUFX
      CHARACTER*4        CHARX
      DOUBLE PRECISION   D
      INTEGER            ID(2)
      EQUIVALENCE       (ID(1) , D )
C
C     RECHERCHE DU DERNIER CARACTERE NON BLANC DU SUFFIXE
      NCFSUF = INDEX( SUFFIX , ' ' )
      IF( NCFSUF .EQ. 1 ) RETURN
C
C     LE DERNIER CARACTERE NON BLANC DU SUFFIXE
      IF( NCFSUF .LE. 0 ) THEN
         NCFSUF = LEN( SUFFIX )
      ELSE
         NCFSUF = NCFSUF - 1
      ENDIF
C
C     DECODAGE EVENTUEL DU SUFFIXE
      IF( SUFFIX(1:1) .NE. '^' ) THEN
C
C        LA VALEUR EST UNE CHAINE DE CARACTERES
C        ======================================
         VASUFX = SUFFIX(1:NCFSUF)
C
      ELSE
C
C        LA VALEUR EST A DECODER
C        =======================
         VASUFX = SUFFIX(2:NCFSUF)
C        RECHERCHE DE L'IDENTIFICATEUR A DECODER DANS VASUFX
         L      = NCFSUF - 1
         DO 100 I = NBIDEN(LHTMS) , 1 , -1
C           LE NUMERO DE LIGNE ET COLONNE DU 1-ER CARACTERE
C           DE L'IDENTIFICATEUR I
            NCD = IDENT(5,I)
            NCF = NCD + L - 1
            IF( VASUFX(1:L) .EQ. KTD(IDENT(4,I))(NCD:NCF) ) THEN
C              LE NUMERO DU TMS EN COURS DE TRAITEMENT
               LHA = LAPILE(0,LHPILE)
C              L'ADRESSE MCN DU DEBUT DU TABLEAU MS
               MN = MCTAMS( LHA )
C              L'ADRESSE DANS LE TABLEAU MCN DU 1-ER MOT DE L'ENTIER
               MN = MN + IDENT(3,I)
C              L'ADRESSE DE LA VARIABLE DU TABLEAU ACTUEL ( VARIABLE=>0 )
               MN = MN + NBVATA(LHTMS)
C              LE TYPE DE L'IDENTIFICATEUR
               NOTYPE = IDENT(1,I)
C
C              CONVERSION EN CARACTERES SELON LE TYPE . ICI MOT = ENTIER !
               GOTO( 90 , 20 , 90 , 40 , 50 , 60 , 90 , 90 ,
     %               90 , 40 , 40 , 90 ) , NOTYPE
C
C              CARACTERE4
 20            VASUFX = CHARX( MCN(MN) )
               GOTO 500
C
C              ENTIER
 40            WRITE( VASUFX , '(1X,I15)' ) MCN(MN)
               GOTO 500
C
C              REEL
 50            WRITE( VASUFX , '(1X,G15.7)' ) RMCN(MN)
               GOTO 500
C
C              REEL2
C              LE DOUBLE PRECISION EST COPIE DANS D  VIA L'EQUIVALENCE  SUR ID
 60            ID(1) = MCN( MN )
               ID(2) = MCN( MN + 1 )
               WRITE( VASUFX , '(1X,G15.7)' ) D
               GOTO 500
C
C              ERREUR
 90            NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NOTYPE
               KERR(1)='FUSUFX: TYPE '//KERR(MXLGER)(1:4)
     %                  //' NON FORMATTABLE'
               CALL LEREUR
            ENDIF
 100     CONTINUE
         NBLGRC(NRERR) = 1
         KERR(1) = 'FUSUFX: SUFFIXE INCORRECT ='//SUFFIX
         CALL LEREUR
         RETURN
      ENDIF
C
C     RECHERCHE DE L'EVENTUEL DERNIER SUFFIXE LEQUEL SERAIT REMPLACE
C     --------------------------------------------------------------
 500  LNM = LEN( NOM )
      DO 510 I=LNM,1,-1
         IF( NOM(I:I) .EQ. '>' ) GOTO 550
         IF( NOM(I:I) .EQ. KSUFIX ) THEN
            NCNM = I
            GOTO 570
         ENDIF
 510  CONTINUE
C
C     PAS DE DERNIER SUFFIXE : RECHERCHE DU 1-ER BLANC DANS NOM
C     ---------------------------------------------------------
 550  NCNM = INDEX( NOM , ' ' )
C
C     RECHERCHE DU PREMIER I ET DERNIER NCFSUF CARACTERE NON BLANC DE VASUFX
C     ----------------------------------------------------------------------
 570  DO 580 I=1,MXSUFX
         IF( VASUFX(I:I) .NE. ' ' ) GOTO 590
 580  CONTINUE
 590  NCFSUF = INDEX( VASUFX(I:MXSUFX) , ' ' )
      IF( NCFSUF .LE. 0 ) THEN
         NCFSUF = MXSUFX
      ELSE
         NCFSUF = I - 2 + NCFSUF
      ENDIF
C
      IF( NCNM+NCFSUF-I+1 .GT. LNM ) THEN
C        PAS ASSEZ DE PLACE POUR METTRE LE SUFFIXE
         NBLGRC(NRERR) = 3
         KERR(1) ='FUSUFX:NOM TROP LONG. PLUS DE PLACE POUR LE SUFFIXE'
         KERR(2) = NOM(1:NCNM)
         KERR(3) = SUFFIX(1:NCFSUF)
         CALL LEREUR
         RETURN
      ENDIF
C
C     CONCATENATION DU SUFFIXE
      NOM = NOM(1:NCNM-1) // KSUFIX // VASUFX(I:NCFSUF)
C
      RETURN
      END
