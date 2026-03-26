      SUBROUTINE IDENT1STPP( NDIM, MNXYZ1, MNNSEF1,  NUMST1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER 1 SOMMET AU SOMMET LE PLUS PROCHE POUR FAIRE UNE
C -----    COUTURE DU MAILLAGE, SA POSITION DEVIENT LE MILIEU DE
C          CES 2 SOMMETS

C ENTREES:
C --------
C NDIM   : =2 MAILLAGE BIDIEMSIONNEL OU =3 MAILLAGE TRIDIMENSIONNEL
C NX1    : ABSCISSE DU POINT A IDENTIFIER A SON PLUS PROCHE VOISIN
C NY1    : ORDONNEE DU POINT A IDENTIFIER A SON PLUS PROCHE VOISIN
C MNXYZ1 : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF1: ADRESSE DU TABLEAU NSEF      DE LA SURFACE
C
C SORTIES:
C --------
C NUMST1 : >0 NUMERO DU SOMMET FINAL (PLUS PETIT DES 2 NUMEROS DE SOMMETS )
C          =0 SI ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC PARIS& St PIERRE du PERRAY Mai 2016
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

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZ1+WNBSOM)
      NUMST1 = 0

      IF( NDIM .EQ. 2 ) THEN

C        2D CAS BIDIM:  LE PREMIER SOMMET CLIQUE
C        ---------------------------------------
         CALL SAIPTC( NOTYEV, NX1, NY1, NOCHAR )
C        COORDONNEES OBJET DU POINT CLIQUE NX1,NY1
         X=XOB2PX(NX1)
         Y=YOB2PX(NY1)
C        RECHERCHE DU SOMMET LE PLUS PROCHE NUMST1
         MN1 = MNXYZ1 + WYZSOM
         DISTMIN=(X-RMCN(MN1))**2 + (Y-RMCN(MN1+1))**2
         NUMST1=1
C        ON BOUCLE SUR LES POINTS DU MAILLAGE
         MN1 = MNXYZ1 + WYZSOM + 3
         DO I=2,NBSOM
            DISTPS=(X-RMCN(MN1))**2 + (Y-RMCN(MN1+1))**2
            IF (DISTPS .LT. DISTMIN) THEN
               DISTMIN = DISTPS
               NUMST1=I
            ENDIF
            MN1 = MN1 + 3
         ENDDO

C        COORDONNEES DU SOMMET NUMST1
         MN1 = MNXYZ1+WYZSOM+3*NUMST1-3
         X = RMCN(MN1)
         Y = RMCN(MN1+1)

C        RECHERCHE DU SOMMET NUMST2 LE PLUS PROCHE DU SOMMET NUMST1
         NUMST2   = 0
         DISTMIN = 1E30
C        BOUCLE SUR LES POINTS DU MAILLAGE POUR TROUVER 
C        LE PLUS PROCHE SOMMET NUMST2 DU SOMMET NUMST1
         MN2 = MNXYZ1 + WYZSOM
         DO I=1,NBSOM
            IF( I .NE. NUMST1 ) THEN
               DISTPS = (X-RMCN(MN2))**2 + (Y-RMCN(MN2+1))**2
               IF (DISTPS .LE. DISTMIN) THEN
                  DISTMIN=DISTPS
                  NUMST2=I
               ENDIF
            ENDIF
            MN2 = MN2 + 3
         ENDDO

      ELSE

C        3D CAS TRIDIM:  RECHERCHE DU SOMMET NUMST1 CLIQUE
C        -------------------------------------------------
         CALL SESTCLIC( RMCN(MNXYZ1+WYZSOM), ITST1, NUMST1 )
         IF( NUMST1 .LE. 0 ) THEN
C           NUMST1=0 SI LE POINT CLIQUE EST HORS MAILLAGE ou ABANDON DEMANDE
            IERR = 1
            GOTO 9999
         ENDIF
         MN1 = MNXYZ1 + WYZSOM + 3*NUMST1 -3
         PRINT*,'ident1stpp: SOMMET 1 CLIQUE: No=',NUMST1,
     %          ' X=',RMCN(MN1),' Y=',RMCN(MN1+1),' Z=',RMCN(MN1+2)

C        COORDONNEES DU SOMMET NUMST1
         X=RMCN(MN1)
         Y=RMCN(MN1+1)
         Z=RMCN(MN1+2)

C        BOUCLE SUR LES POINTS DU MAILLAGE POUR TROUVER 
C        LE PLUS PROCHE SOMMET NUMST2 DU SOMMET NUMST1
         NUMST2  = 0
         DISTMIN = 1E30
         MN2 = MNXYZ1 + WYZSOM
         DO I=1,NBSOM
            IF( I .NE. NUMST1 ) THEN
               DISTPS = (X-RMCN(MN2)  ) **2
     %                + (Y-RMCN(MN2+1)) **2
     %                + (Z-RMCN(MN2+2)) **2
               IF (DISTPS .LE. DISTMIN) THEN
                  DISTMIN=DISTPS
                  NUMST2=I
               ENDIF
            ENDIF
            MN2 = MN2 + 3
         ENDDO

      ENDIF

C     EN DIMENSION 2 ou 3
      IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,*)'ident1stpp: St',NUMST1,NUMST2,' sont IDENTIFIES'
      ELSE
        WRITE(IMPRIM,*)'ident1stpp: Vertices',NUMST1,NUMST2,
     %                 ' are IDENTIFIED'
      ENDIF

      IF( NUMST1 .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'ident1stpp: SOMMET 1 NON DETECTE'
         ELSE
            WRITE(IMPRIM,*) 'ident1stpp: VERTEX 1 NOT DETECTED'
         ENDIF
         NUMST1 = 0
         GOTO 9999
      ENDIF

C     LE SOMMET NUMST1 DEVIENT LE MILIEU DE SOMMET1-SOMMET2 (NUMST1-NUMST2)
      CALL IDENT2SM( NUMST1, NUMST2, MNXYZ1, MNNSEF1 )

 9999 RETURN
      END
