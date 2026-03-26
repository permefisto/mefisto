        SUBROUTINE SUEX06( NTLXSU , LADEFI ,
     %                     NTFASU , MNFASU , NTSOFA , MNSOFA , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN TRIANGLE TRANSFINI STRUCTURE
C -----    OU NON STRUCTURE AVEC OU SANS TANGENTES AUX ARETES
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU TRIANGLE
C LADEFI : TABLEAU DE DEFINITION DU TRIANGLE
C          CF $MEFISTO/td/d/a_surface__definition
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES SOMMETS DES EF
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES SOMMETS DES EF
C          CF $MEFISTO/td/d/a___nsef
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DU TRIANGLE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DU TRIANGLE
C          CF $MEFISTO/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS         OCTOBRE 1996
C234567..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           LADEFI(0:*)
      INTEGER           NBSOCT(3),MNSOCT(3),NUCOTE(3)
      INTEGER           NBARLI(3),MNARLI(3),
     %                  MNXYTG(3),MNNTGL(3),MNCGEF(3)
      REAL              XYZ(3,3,2)
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE'
      NOFOTI = NOFOTIEL()
C
C     VERIFICATION DES LIGNES DES 3 COTES DU TRIANGLE ALGEBRIQUE
      IERR  = 0
      NBTGS = 0
      DO 1 N=1,3
C
C        LE NUMERO DE LA LIGNE N
         NOLI = LADEFI(WU3COT-1+N)
C
C        RECHERCHE DE LA LIGNE DU COTE N ET DES TABLEAUX NSEF XYZSOMMET
C        ATTENTION: SUR CHAQUE COTE LA LIGNE DOIT ETRE STRUCTUREE
         CALL LISTRE( NOLI, NTLXLI, NTARLI, MNARLI(N),
     %                NTSOLI, MNSOCT(N), IERRP )
         IF( IERRP .NE. 0 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:2),'(I2)') N
            KERR(1) = 'ERREUR LA LIGNE DU COTE' // KERR(MXLGER)(1:2)
            KERR(2) = 'EST NON STRUCTUREE'
            CALL LEREUR
            IERR = N
            GOTO 1
         ENDIF
C
C        LE NOMBRE DE SOMMETS DE LA LIGNE
         NBSOCT(N) = MCN( MNSOCT(N) + WNBSOM )
C
C        CALCUL DES 2 NUMEROS DES 2 TANGENTES DES ARETES DE LA LIGNE
         CALL TGARLI( MNARLI(N), MNSOCT(N),
     %                NBARLI(N), MNXYTG(N), MNNTGL(N), MNCGEF(N) )
C
C        LE NOMBRE DE COTES AVEC DES TANGENTES
         IF( NBARLI(N) .GT. 0 ) NBTGS = NBTGS + 1
C
 1    CONTINUE
      IF( IERR .NE. 0 ) GOTO 9999
C
C     LES LIGNES SONT ORDONNEES POUR OBTENIR LE SENS SUIVANT
C     C1:S1S2 C2:S2S3 C3:S3S1
C
C                S3
C                |  \
C                |    \
C          C3   \/     /\ C2
C                |         \
C                |            \
C                S1---->-------S2
C                      C1
C
      DO 2 N=1,3
C
C        LE NUMERO DE LA LIGNE (1 A 3) DU COTE N
         NUCOTE(N) = N
C
C        LES COORDONNEES DU POINT INITIAL DU COTE N
         MN = MNSOCT(N) + WYZSOM
         XYZ(1,N,1)=RMCN(MN  )
         XYZ(2,N,1)=RMCN(MN+1)
         XYZ(3,N,1)=RMCN(MN+2)
C
C        LES COORDONNEES DU POINT FINAL DU COTE N
         MN = MNSOCT(N) + WYZSOM - 3 + 3 * NBSOCT(N)
         XYZ(1,N,2)=RMCN(MN  )
         XYZ(2,N,2)=RMCN(MN+1)
         XYZ(3,N,2)=RMCN(MN+2)
C
  2   CONTINUE
C
C     RECHERCHE DU COTE 2 S2->S3
      DO 3 N=2,3
         CALL XYZIDE( XYZ(1,1,2),XYZ(1,N,1),IDENT )
         IF (IDENT.EQ.1) THEN
C           LE COTE 2 EST LA LIGNE N (DE 1 A 3) DANS LES DONNEES UTILISATEUR
            NUCOTE(2)=N
C           INDICE DU POINT OPPOSE (INITIAL<->FINAL)
            I2 = 2
            GO TO 4
         ENDIF
         CALL XYZIDE( XYZ(1,1,2),XYZ(1,N,2),IDENT )
         IF (IDENT.EQ.1) THEN
            NUCOTE(2)=-N
C           INDICE DU POINT OPPOSE (INITIAL<->FINAL)
            I2 = 1
            GO TO 4
         ENDIF
 3    CONTINUE
      GOTO 9900
C
C     RECHERCHE DU COTE 3 S3->S1
 4    DO 6 N=2,3
         CALL XYZIDE( XYZ(1,1,1),XYZ(1,N,2),IDENT )
         IF (IDENT.EQ.1) THEN
            NUCOTE(3)=N
C           INDICE DU POINT OPPOSE (INITIAL<->FINAL)
            I3 = 1
            GO TO 8
         ENDIF
         CALL XYZIDE( XYZ(1,1,1),XYZ(1,N,1),IDENT )
         IF (IDENT.EQ.1) THEN
            NUCOTE(3)=-N
C           INDICE DU POINT OPPOSE (INITIAL<->FINAL)
            I3 = 2
            GO TO 8
         ENDIF
 6    CONTINUE
      GOTO 9900
C
C     VERIFICATION DE LA FERMETURE DU TRIANGLE COURBE
 8    N    = ABS( NUCOTE(2) )
      NOLI = ABS( NUCOTE(3) )
      CALL XYZIDE( XYZ(1,N,I2), XYZ(1,NOLI,I3), IDENT )
      IF (IDENT.NE.1) GOTO 9900
C
C     SI LA FONCTION TAILLE_IDEALE(x,y,z) EXISTE
C     ALORS TRIANGULATION NON STRUCTUREE FORCEE
      IF( NOFOTI .GT. 0 ) GOTO 100
C
      IF( NBSOCT(1) .EQ. NBSOCT(2) .AND. NBSOCT(1) .EQ. NBSOCT(3) ) THEN
C
C        TRIANGLE STRUCTURE ALGEBRIQUE AVEC OU SANS TG
C        =============================================
         CALL TRSTAT( NTLXSU, NUCOTE, NBTGS,
     %                NBSOCT(1), MNSOCT,
     %                NBARLI, MNXYTG, MNNTGL,
     %                NTFASU, MNFASU, NTSOFA, MNSOFA )
         GOTO 9000
      ENDIF
C
C     TRIANGLE NON STRUCTURE ALGEBRIQUE AVEC OU SANS TG
C     =================================================
 100  CALL TRNSAL( NTLXSU, NUCOTE, NBTGS,
     %             NBSOCT, MNSOCT,
     %             NBARLI, MNXYTG, MNNTGL,
     %             NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C
C     LE NOMBRE DE COORDONNEES D'UN SOMMET
 9000 MCN( MNSOFA + WBCOOR ) = 3
C
C     SUPPRESSION DES TANGENTES DOUBLES
C     =================================
      CALL MOINTG( NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C
C     IDENTIFICATION ET IMPOSITION DES COORDONNEES DES SOMMETS
C     DES LIGNES DU CONTOUR AU SOMMET LE PLUS PROCHE DE LA SURFACE
C     ============================================================
      IF( MNSOFA .GT. 0 ) THEN
         DO 9010 I=1,3
            IF( MNSOCT(I) .GT. 0 ) THEN
               CALL IDLISU( MNSOCT(I), MNSOFA )
            ENDIF
 9010    CONTINUE
      ENDIF
      GOTO 9999
C
C     ERREUR TRIANGLE NON FERME
 9900 NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR: LES 3 LIGNES NE FORMENT PAS UN TRIANGLE FERME'
      CALL LEREUR
      IERR = 6
C
 9999 RETURN
      END
