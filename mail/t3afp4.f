      SUBROUTINE T3AFP4( NBARFR, NBFPLA, NOFAFR,
     %                   MOARFR, MXARFR, LAREFR,  NBCOOR, XYZPOI,
     %                   TRIANG, NOAXE,  XYZSFP,  TEMSFP,
     %                   NCAS0,  NCAS1,  NCAS,
     %                   TMIN,   TMAX,   NORMALE, NCOPLAN, NTLPLAN  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES FRONTALIERES ET EN COULEURS
C ----     LA TEMPERATURE DES FACES DANS DES PLANS DE SECTION SELON UN AXE
C          ET ORTHOGONALEMENT A CE PLAN
C
C ENTREES:
C --------
C NBARFR : NOMBRE DES ARETES FRONTALIERES
C NBFPLA : NOMBRE DE FACES INTERSECTIONS AVEC LES PLANS
C NOFAFR : NUMERO DES FACES SELON LEUR DISTANCE DECROISSANTE A L'OEIL
C MOARFR : NOMBRE DE MOTS PAR ARETE FRONTALIERE DU TABLEAU LAREFR
C MXARFR : NOMBRE DE FACES DU TABLEAU LAREFR
C LAREFR : TABLEAU NUMERO DES 2 SOMMETS ET LIEN
C          LAREFR(1,I)= NO DU 1-ER  SOMMET DE L'ARETE FRONTALIERE
C          LAREFR(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LAREFR(3,I)= 0 OU NUMERO DE L'ARETE FRONTALIERE SUIVANTE
C                       DANS LE HACHAGE
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS (3 ou 6)
C XYZPOI : COORDONNEES DES SOMMETS DES FACES FRONTALIERES (ET INTERNES)
C TRIANG : REEL VALEUR TEMOIN DE TRIANGLE DANS XYZSFP(4,4,*)
C NOAXE  : NUMERO COORDONNEE DU PLAN (1 a 3)
C XYZSFP : XYZ DES 4 SOMMETS DES FACES DES PLANS DE SECTION
C TEMSFP : TEMPERATURE AUX 4 SOMMETS DES FACES DES PLANS DE SECTION
C NCAS0  : NUMERO DU PREMIER CAS OU VECTEUR SOLUTION
C NCAS1  : NUMERO DU DERNIER CAS OU VECTEUR SOLUTION
C NCAS   : NUMERO DU CAS A TRACER
C TMIN   : TEMPERATURE MINIMALE POUR TOUS LES CAS TRAITES
C TMAX   : TEMPERATURE MAXIMALE POUR TOUS LES CAS TRAITES
C NORMALE: VECTEUR NORMAL AUX PLANS DE SECTION
C NCOPLAN: NUMERO DE LA COULEUR DES ARETES DES PLANS
C NTLPLAN: NUMERO DU TYPE DE LIGNE DE LA COULEUR DES ARETES DES PLANS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :PERRONNET ALAIN UPMC ANALYSE NUMERIQUE LJLL PARIS OCTOBRE 2003
C2345X7..............................................................012
      PARAMETER    ( LIGCON=0, LIGTIR=1 )
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      INTEGER        NOFAFR(1:NBFPLA*2+NBARFR)
      INTEGER        LAREFR(1:MOARFR, 1:MXARFR)
      REAL           XYZPOI(1:NBCOOR, 1:*)
      REAL           XYZSFP(1:4, 1:4, 1:NBFPLA)
      REAL           TEMSFP(1:4, NCAS0:NCAS1, 1:NBFPLA)
      REAL           NORMALE(3)
      REAL           COUL(4)
      REAL           XYZ1(3), XYZ2(3)
      REAL           TRIANG
      REAL           XYZ(3,4), XYZP(3)
C
cccC     LA PALETTE ARC EN CIEL
ccc      CALL PALCDE(11)
C
C     LE NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL
C
C     LES ARETES SONT TRACEES AVEC UNE EPAISSEUR
      CALL XVEPAISSEUR( 1 )
C
C     PLAN DE TRACE POSITIONNE AU MINIMUM OU ZERO DES TEMPERATURES
      IF( TMIN .LE. 0 .AND. 0 .LE. TMAX ) THEN
         TMILIEU = 0
      ELSE
ccc      TMILIEU = (TMIN+TMAX)/2
         TMILIEU = TMIN
      ENDIF
C
C     LE TRACE DES FACES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     =======================================================
      DO 50 NF = 1, NBFPLA * 2 + NBARFR
C
C        LE NUMERO DE LA FACE OU ARETE LA PLUS ELOIGNEE NON ENCORE TRACEE
         NAOF = NOFAFR( NF )
C
         IF( NAOF .GT. 0 ) THEN
C
C           TRACE DE LA FACE INTERSECTION AVEC UN PLAN
C           ------------------------------------------
C           LE NOMBRE DE SOMMETS DE LA FACE
            IF( NAOF .LE. NBFPLA ) THEN
                NUF = NAOF
            ELSE
                NUF = NAOF - NBFPLA
            ENDIF
            IF( XYZSFP(3,4,NUF) .EQ. TRIANG ) THEN
C              TRIANGLE
               NAF = 3
            ELSE
C              QUADRANGLE
               NAF = 4
            ENDIF
C
            IF( NAOF .LE. NBFPLA ) THEN
C
C              LA COULEUR DES NAF SOMMETS DE LA FACE DE SECTION OU PROFIL
               DO 10 I=1,NAF
C                 LA COULEUR
                  IF( TMIN .NE. TMAX ) THEN
                     COUL(I) = ( TEMSFP(I,NCAS,NAOF) - TMIN )
     %                       / (TMAX-TMIN) * NBCOUL + N1COUL
                  ELSE
                     COUL(I)=NDCOUL
                  ENDIF
                  IF( COUL(I) .GT. NDCOUL ) COUL(I)=NDCOUL
                  IF( COUL(I) .LT. N1COUL ) COUL(I)=N1COUL
 10            CONTINUE
C
               CALL XVTYPETRAIT( NTLAPL )
C
C              LES XYZ DES SOMMETS DE LA FACE SUR LE PROFIL
C              MODIFICATION XYZ -> XYZ + NORMALE * TEMPERATURE
               DO 30 I=1,NAF
                  DO 20 J=1,3
                     XYZ(J,I) = XYZSFP(J,I,NAOF)
     %                        + NORMALE(J)*(TEMSFP(I,NCAS,NAOF)-TMILIEU)
 20               CONTINUE
 30            CONTINUE
C
               IF( PREDUF .GT. 0 ) THEN
C
C                 REDUCTION DE LA FACE
                  REDUCF = PREDUF * 0.01
                  REDUC1 = 1.0 - REDUCF
C                 CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
                  CALL COBAPO( NAF, XYZ, XYZP )
C                 L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
                  DO 40 I=1,NAF
                     XYZ(1,I) = XYZ(1,I) * REDUC1 + XYZP(1) * REDUCF
                     XYZ(2,I) = XYZ(2,I) * REDUC1 + XYZP(2) * REDUCF
                     XYZ(3,I) = XYZ(3,I) * REDUC1 + XYZP(3) * REDUCF
 40               CONTINUE
C
               ENDIF
C
               IF( NAF .EQ. 3 ) THEN
C                 LE TRACE DU TRIANGLE
                  CALL TRIACOUL3DBORD( XYZ, COUL, NCOAPL, 1 )
               ELSE
C                 LE TRACE DU QUADRANGLE
                  CALL QUADCOUL3DBORD( XYZ, COUL, NCOAPL, 1 )
               ENDIF
C
            ELSE
C
               IF( NCOPLAN .GE. 0 ) THEN
C
C                 TRACE DES ARETES D'UNE FACE D'UN PLAN DE COORDONNEE NOAXE
                  CALL XVTYPETRAIT( NTLPLAN )
C                 LA PREMIERE EXTREMITE DE L'ARETE NAF
                  XYZ2(1)     = XYZSFP(1,NAF,NUF)
                  XYZ2(2)     = XYZSFP(2,NAF,NUF)
                  XYZ2(3)     = XYZSFP(3,NAF,NUF)
                  XYZ2(NOAXE) = XYZSFP(4,NAF,NUF)
C
                  DO 45 I=1,NAF
C                    LA SECONDE EXTREMITE DE L'ARETE NAF
                     XYZ1(1)     = XYZSFP(1,I,NUF)
                     XYZ1(2)     = XYZSFP(2,I,NUF)
                     XYZ1(3)     = XYZSFP(3,I,NUF)
C                    LA COORDONNEE NOAXE EST IMPOSEE
                     XYZ1(NOAXE) = XYZSFP(4,I,NUF)
                     CALL TRAIT3D( NCOPLAN, XYZ2, XYZ1 )
C
C                    LE SECOND DEVIENT LE PREMIER SOMMET
                     XYZ2(1) = XYZ1(1)
                     XYZ2(2) = XYZ1(2)
                     XYZ2(3) = XYZ1(3)
 45               CONTINUE
C
               ENDIF
C
            ENDIF
C
         ELSE
C
C           TRACE DE L'ARETE FRONTALIERE DU MAILLAGE
C           ----------------------------------------
C           LE NUMERO DE L'ARETE DANS LE TABLEAU LAREFR
            IF( NCOAFR .GE. 0 ) THEN
               CALL XVTYPETRAIT( NTLAFR )
               NAOF = -NAOF
               CALL TRAIT3D( NCOAFR, XYZPOI(1,LAREFR(1,NAOF)),
     %                               XYZPOI(1,LAREFR(2,NAOF)) )
            ENDIF
         ENDIF
 50   CONTINUE
C
C     RETOUR AU TRACE CONTINU DES LIGNES
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
