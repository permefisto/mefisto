       SUBROUTINE IDENT2ST( NDIM, MNXYZS1, MNNSEF1, NUMST1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER 2 SOMMETS EN UN SEUL POUR FAIRE UNE COUTURE DU
C -----    MAILLAGE, SA POSITION DEVIENT LE MILIEU DE CES 2 SOMMETS

C ENTREES:
C --------
C NDIM   : 2 MAILLAGE BIDIEMSIONNEL OU 3 MAILLAGE TRIDIMENSIONNEL
C NX1    : ABSCISSE DU 1-ER  POINT A IDENTIFIER
C NY1    : ORDONNEE DU 1-ER  POINT A IDENTIFIER
C NX2    : ABSCISSE DU 2-EME POINT A IDENTIFIER
C NY2    : ORDONNEE DU 2-EME POINT A IDENTIFIER
C MNXYZS1: ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF1: ADRESSE DU TABLEAU NSEF      DE LA SURFACE

C SORTIES:
C --------
C NUMST1 : >0 NUMERO DU SOMMET FINAL (PLUS PETIT DES 2 NUMEROS DE SOMMETS )
C             LE SOMMET NUMST1 DEVIENT LE MILIEU (NUMST1+NUMST2)/2
C          =0 SI ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & St PIERRE du PERRAY NOVEMBRE 2015
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NUMST1, NUMST2

C     NBTQ NOMBRE DES EF
      NBTQ   = MCN(MNNSEF1+WBEFOB)
      MNSOEL = MNNSEF1 + WUSOEF

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZS1+WNBSOM)
      NUMST1 = 0

      IF( NDIM .EQ. 2 ) THEN

C        2D CAS BIDIM
C        RECHERCHE DU PREMIER SOMMET CLIQUE
C        ----------------------------------
         CALL SAIPTC( NOTYEV, NX1, NY1, NOCHAR )
C        COORDONNEES OBJET DU POINT CLIQUE NX1,NY1
         X=XOB2PX(NX1)
         Y=YOB2PX(NY1)
C        RECHERCHE DU SOMMET LE PLUS PROCHE NUMST1
         MN1 = MNXYZS1+WYZSOM
         DISTMIN=(X-RMCN(MN1))**2 + (Y-RMCN(MN1+1))**2
         NUMST1=1
C        BOUCLE SUR LES POINTS DU MAILLAGE POUR TROUVER LE
C        SOMMET LE PLUS PROCHE DE (X,Y)
         MN1 = MNXYZS1+WYZSOM+3
         DO I=2,NBSOM
            DISTPS = (X-RMCN(MN1))**2 + (Y-RMCN(MN1+1))**2
            IF (DISTPS .LT. DISTMIN) THEN
               DISTMIN=DISTPS
               NUMST1=I
            ENDIF
            MN1 = MN1 + 3
         ENDDO

C        RECHERCHE DU SECOND SOMMET CLIQUE
C        ---------------------------------
         CALL SAIPTC( NOTYEV, NX2, NY2, NOCHAR )
C        COORDONNEES OBJET DU POINT CLIQUE NX2,NY2
         X=XOB2PX(NX2)
         Y=YOB2PX(NY2)
C        RECHERCHE DU SOMMET LE PLUS PROCHE NUMST2
         MN1 = MNXYZS1+WYZSOM
         DISTMIN = (X-RMCN(MN1))**2 + (Y-RMCN(MN1+1))**2
         NUMST2=1
C        BOUCLE SUR LES POINTS DU MAILLAGE POUR TROUVER LE
C        SOMMET LE PLUS PROCHE DE (X,Y)
         MN1 = MNXYZS1+WYZSOM+3
         DO I=2,NBSOM
            DISTPS = (X-RMCN(MN1))**2 + (Y-RMCN(MN1+1))**2
            IF (DISTPS .LT. DISTMIN) THEN
               DISTMIN=DISTPS
               NUMST2=I
            ENDIF
            MN1 = MN1 + 3
         ENDDO

      ELSE

C        3D CAS TRIDIM
C        RECHERCHE DU 1-ER SOMMET NUMST1 CLIQUE
C        --------------------------------------
         CALL SESTCLIC( RMCN(MNXYZS1+WYZSOM), ITST1, NUMST1 )
         IF( NUMST1 .LE. 0 ) THEN
C           NUMST1=0 SI LE POINT CLIQUE EST HORS MAILLAGE ou ABANDON DEMANDE
            IERR = 1
            GOTO 9999
         ENDIF
         MN = MNXYZS1 + WYZSOM + 3*NUMST1 -3
         PRINT*,'ident2st: SOMMET 1 CLIQUE: No=',NUMST1,
     %          ' X=',RMCN(MN),' Y=',RMCN(MN+1),' Z=',RMCN(MN+2)

C        RECHERCHE DU 2-EME SOMMET NUMST2 CLIQUE
C        ---------------------------------------
         CALL SESTCLIC( RMCN(MNXYZS1+WYZSOM), ITST2, NUMST2 )
         IF( NUMST2 .LE. 0 ) THEN
C           NUMST2=0 SI LE POINT CLIQUE EST HORS MAILLAGE ou ABANDON DEMANDE
            IERR = 1
            GOTO 9999
         ENDIF
         MN = MNXYZS1 + WYZSOM + 3*NUMST2 -3
         PRINT*,'s2p1tria: SOMMET 2 CLIQUE: No=',NUMST2,
     %          ' X=',RMCN(MN),' Y=',RMCN(MN+1),' Z=',RMCN(MN+2)

         IF( NUMST2 .EQ. NUMST1 ) THEN
            WRITE(IMPRIM,*)
     %     'ident2st: SOMMET 1 et 2 CLIQUES sont IDENTIQUES. A REDONNER'
            NUMST1 = 0
            GOTO 9999
         ENDIF

      ENDIF

C     TRAITEMENT AVEC NDIM=2 ou 3
C     ---------------------------
      IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,*)'ident2st: St',NUMST1,NUMST2,' sont IDENTIFIES'
      ELSE
        WRITE(IMPRIM,*)'ident2st: Vertices',NUMST1,NUMST2,
     %                 ' are IDENTIFIED'
      ENDIF

      IF( NUMST1 .EQ. NUMST2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,*)'ident2st: St',NUMST1,NUMST2,' sont IDENTIQUES'
         ELSE
         WRITE(IMPRIM,*)'ident2st: Vertices',NUMST1,NUMST2,' are SAME'
         ENDIF
      ENDIF

      IF( NUMST1 .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'ident2st: SOMMET 1 NON DETECTE'
         ELSE
            WRITE(IMPRIM,*) 'ident2st: VERTEX 1 NOT DETECTED'
         ENDIF
         GOTO 9999
      ENDIF

      IF( NUMST2 .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'ident2st: SOMMET 2 NON DETECTE'
         ELSE
            WRITE(IMPRIM,*) 'ident2st: VERTEX 2 NOT DETECTED'
         ENDIF
         GOTO 9999
      ENDIF

C     LE SOMMET NUMST1 DEVIENT LE MILIEU DE SOMMET1-SOMMET2 (NUMST1-NUMST2)
      CALL IDENT2SM( NUMST1, NUMST2, MNXYZS1, MNNSEF1 )

 9999 RETURN
      END
