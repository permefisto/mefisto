      SUBROUTINE SDDEF6(NTLXSD,NTOBRE,NTEST,IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RECHERCHE DES POINTS DE L'INTERFACE CONTENUS DANS
C -----   LES OBJETS AUX LIMITES . RETRAIT DE CES POINTS
C
C ENTREE :
C --------
C NTLXSD : LE NUMERO DE LEXIQUE DE L'OBJET SOUS-DOMAINE
C NTOBRE : LE NUMERO DE LEXIQUE DE L'OBJET A RETROUVER
C
C SORTIE :
C --------
C NTEST  : VALEUR TEST : 1 SI L'OBJET EST INCLUS DANS L'INTERFACE
C                        0 SINON
C IERR   : CODE D'ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  FEVRIER 1990
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___interface.inc"
      include"./incl/a___xyzsommet.inc"
C
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE TABLEAU 'XYZSOMMET' DE L'OBJET A RETROUVER
      CALL LXTSOU( NTOBRE , 'XYZSOMMET' , NTSORE , MNSORE )
      IF( NTSORE .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : TABLEAU SOMMET NON RETROUVE'
         CALL LEREUR
         CALL LXIM( NTOBRE )
         IERR = 1
         RETURN
      ENDIF
      NBSORE = MCN(MNSORE+WNBSOM)
      NXYZRE = MNSORE+WYZSOM-1
C
C     LE TABLEAU 'INTERFACE' DE L'OBJET SOUS-DOMAINE
      CALL LXTSOU( NTLXSD , 'INTERFACE' , NTINTE , MNINTE )
      IF( NTINTE .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : TABLEAU INTERFACE NON RETROUVE'
         CALL LEREUR
         CALL LXIM( NTLXSD )
         IERR = 1
         RETURN
      ENDIF
      NBSOIN = MCN(MNINTE+WBNOIN)
      LONUNO = MCN(MNINTE+WONUNO)
      NXYZIN = MNINTE + WBOBNO + NBSOIN + LONUNO - 1
C
C     RECHERCHE DES SOMMETS DANS L'INTERFACE
C     --------------------------------------
      NTEST = 0
      NB    = 0
      DO 101 NBSRE = 1 , NBSORE
         NUM1 = NXYZRE + ( NBSRE - 1 ) * 3
         DO 102 NBSIN = 1 , NBSOIN
            NUM2 = NXYZIN + ( NBSIN - 1 ) * 3
            DO 103 J=1,3
               R1 = RMCN(NUM1+J)
               R2 = RMCN(NUM2+J)
               IF     ( ABS(R1) .LE. EPZERO ) THEN
                  IF  ( ABS(R2) .GT. EPZERO ) GOTO 102
               ELSE IF( ABS(R2) .LE. EPZERO ) THEN
                  GOTO 102
               ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
                  GOTO 102
               ENDIF
 103        CONTINUE
C           LE SOMMET EST RETROUVE
            NB = NB + 1
            GO TO 101
 102     CONTINUE
 101  CONTINUE
      IF ( NB .EQ. NBSORE ) THEN
C        TOUS LES SOMMETS SONT RETROUVES
         NTEST = 1
         RETURN
      ELSE
C        TOUS LES SOMMETS NE SONT PAS RETROUVES
         RETURN
      END IF
C
      END
