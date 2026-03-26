      SUBROUTINE POEX35( NTLXPO , LADEFI , NTSOPO , MNSOPO , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXTRAIRE DES SOMMETS D'UNE LIGNE
C -----    CREER LE TABLEAU 'XYZSOMMET' DE CE "POINT"
C          AUCUNE TANGENTE N'EST EXTRAITE
C
C ENTREES:
C --------
C NTLXPO : NUMERO DU TABLEAU TS DU LEXIQUE DU POINT
C LADEFI : TABLEAU ENTIER DE DEFINITION DU POINT
C          CF ~/td/d/a_point__definition
C
C SORTIES:
C --------
C NTSOPO : NUMERO      DU TMS 'XYZSOMMET' DU POINT
C MNSOPO : ADRESSE MCN DU TMS 'XYZSOMMET' DU POINT
C IERR   : 0 SI PAS D'ERREUR
C          1 SI LIGNE INCONNUE OU INCORRECTE
C          4 SI AUCUN SOMMET EXTRAIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1993
C2345X7..............................................................012
      IMPLICIT           INTEGER(W)
      include"./incl/gsmenu.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/ntmnlt.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
C
C     LA LIGNE INITIALE
C     =================
C     LE NOM DE CETTE LIGNE
      NULGSX = LADEFI( WULGSX )
C     LE TABLEAU LEXIQUE DE CETTE LIGNE
      CALL LXNLOU( NTLIGN , NULGSX , NTLXLG , MN )
      IF( NTLXLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE INITIALE INCONNUE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE NUMERO DU PREMIER ET DERNIER SOMMET A EXTRAIRE
      NU1STX = LADEFI( WU1STX )
      NUDSTX = LADEFI( WUDSTX )
      IF( NU1STX .LE. 0 .OR. NUDSTX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NO NEGATIF OU NUL D''ARETE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
      IF( NU1STX .GT. NUDSTX ) THEN
         NU1STX = NUDSTX
         NUDSTX = LADEFI( WU1STX )
      ENDIF
C
C     LE TABLEAU 'NSEF' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG , 'NSEF' , NTARLG , MNARLG )
      IF( NTARLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS NSEF'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG , 'XYZSOMMET' , NTSOLG , MNSOLG )
      IF( NTSOLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS SOMMETS'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      NBSOM = MCN( MNSOLG + WNBSOM )
C
      IF( NUDSTX .GT. NBSOM ) THEN
C        SI TROP DEMANDE LE MAXIMUM EST IMPOSE
         NUDSTX = NBSOM
      ENDIF
C
C     RESERVATION DE LA PLACE NECESSAIRE AUX SOMMETS EXTRAITS
      CALL TNMCDC( 'ENTIER' , NBSOM , MNOUIS )
      MNO = MNOUIS - 1
C
C     CALCUL DU CRITERE EN CHACUN DES SOMMETS DU MAILLAGE
      NBSOEX = 0
      DO 10 N=1,NBSOM
         IF( NU1STX .LE. N .AND. N .LE. NUDSTX ) THEN
C           UN SOMMET DE PLUS
            NBSOEX = NBSOEX + 1
            MCN( MNO+NBSOEX ) = N
         ENDIF
 10   CONTINUE
      IF( NBSOEX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'AUCUN SOMMET EXTRAIT'
         CALL LEREUR
         IERR = 4
         GOTO 9000
      ENDIF
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DU POINT FINAL
C     --------------------------------------------------
      CALL LXTNDC( NTLXPO , 'XYZSOMMET' , 'MOTS'   , WYZSOM+3*NBSOEX )
      CALL LXTSOU( NTLXPO , 'XYZSOMMET' ,  NTSOPO  , MNSOPO )
C
C     LE NOMBRE DE SOMMETS DU POINT
      MCN( MNSOPO + WNBSOM ) = NBSOEX
C
C     LES 3 COORDONNEES DES SOMMETS DU POINT
      MN0 = MNSOLG + WYZSOM - 3
      MN1 = MNSOPO + WYZSOM
      DO 500 I=1,NBSOEX
         N  = MCN(MNO+I)
         MN = MN0 + 3 * N
         RMCN(MN1  ) = RMCN(MN  )
         RMCN(MN1+1) = RMCN(MN+1)
         RMCN(MN1+2) = RMCN(MN+2)
         MN1 = MN1 + 3
 500  CONTINUE
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOPO) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOPO + WNBTGS ) = 0
      MCN( MNSOPO + WBCOOR ) = 3
      MCN( MNSOPO + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE
 9000 CALL TNMCDS( 'ENTIER' , NBSOM   , MNOUIS )
      END
