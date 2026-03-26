      SUBROUTINE AJOMIL19( COSMAX, MXFACO, LEFACO, NBFACO, N1FASC,
     %                     NBFAPE, MXFAPE, NOFAPE, N1TETS, NOTETR,
     %                     NBSOMM, MXSOMM, PTXYZD, NPSOFR, NOSTIS,
     %                     HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                     MXTRCF, NOTRCF, NOSTCF, MXARCF,N1ARCF,NOARCF,
     %                     MXETOI, NAETOI, LARETE,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER A TOUTE ARETE PERDUE SON MILIEU
C -----    DECOUPER EN 2 FACES TOUTES LES FACES DE LEFACO AYANT
C          CETTE ARETE PERDUE COMMUNE
C
C ENTREES :
C ---------
C COSMAX : SEUIL DU COSINUS DE L'ANGLE FORME PAR LES NORMALES AUX
C          2 FACES ET AU DESSUS DUQUEL LES FACES
C          SONT CONSIDEREES COPLANAIRES
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
C
C MODIFIES :
C ----------
C LEFACO : LES 3 SOMMETS, 2 MATERIAUX, 3 FACES VOISINES ET CHAINAGE
C          DES FACES TRIANGULAIRES DU CONTOUR ET INTERFACES
C NBFACO ; NOMBRE DE FACES ACTIVES DE LEFACO
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS
C NBFAPE : NOMBRE DE FACES PERDUES
C MXFAPE : NOMBRE MAXIMAL DE FACES PERDUES
C NOFAPE : NUMERO DES FACES PERDUES
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NBSOMM : NUMERO DU DERNIER SOMMET AJOUTE A LA TETRAEDRISATION
C          MIS A JOUR EN SORTIE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TETRAEDRISATION
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : NUMERO DES POINTS INITIAUX
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR UNE SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET RECONNU TROP PROCHE
C          = -3 SI LE POINT EST SOMMET NON SUPPRIMABLE REFUSE DANS TETRAEDRISATION
C          = -4 SI LE POINT EST SOMMET DE LA GRILLE REGULIERE SUPPRIMABLE
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE OU NO DE POINT INTERNE

C MXTRCF : NOMBRE MAXIMAL DECLARABLE D'ARETES OU TRIANGLES
C          OU SOMMETS DANS L'ETOILE
C
C NAETOI : LISTE DES ARETES DE L'ETOILE DU POINT COURANT
C          TABLEAU ENTIER(4,MXTRCF)
C LARETE : TABLEAU DE HACHAGE DES ARETES DES TRIANGLES PERDUS COPLANAIRES
C          LARETE(1,J) = SOMMET 1 DANS PTXYZD
C          LARETE(2,J) = SOMMET 2 DANS PTXYZD
C          LARETE(3,J) = NUMERO DANS LARETE DE L'ARETE SUIVANTE
C          LARETE(4,J) = NUMERO DU 1-ER  TRIANGLE
C          LARETE(5,J) = NUMERO DU 2-EME TRIANGLE
C          HACHAGE(NS1,NS2) = (NS1 + NS2) MODULO MXARET + 1
C ATTENTION : MEME MC QUE NAETOI ( TABLEAUX EN EQUIVALENCE A L'APPEL )

C N1ARCF : TABLEAU (0:MXTRCF) AUXILIAIRE
C NOARCF : NUMERO DES ARETES DE LA LIGNE DU CONTOUR FERME SELON UN SENS
C
C SORTIES:
C --------
C NOSTIS : NUMERO DES SOMMETS ISOLES N'APPARTENANT PAS AU CONTOUR
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DU CF
C NOSTCF : TABLEAU DU NUMERO DANS PTXYZD DES SOMMETS DU CF
C IERR   : 0 SI PAS D'ERREUR
C          1 SI SATURATION DES SOMMETS DE PTXYZD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du perray     AOUT 2012
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      PARAMETER       ( MXSMIN=1024 )
      INTEGER           NOSMIN(MXSMIN), NODSNO(MXSMIN)
      DOUBLE PRECISION  DISTMIN(MXSMIN)
      INTEGER           NBIPAV(3), N1SPAVE(0:*), NOPTSUIV(1:*)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3)
C
C     ANGLE AU DELA DUQUEL IL FAUT LE DECOUPER
      DOUBLE PRECISION  COSMIN
      PARAMETER        (COSMIN=-0.8D0)
C     SOIT  UN ANGLE ARCOS(-0.76)=139.46 DEGRES
C     SOIT  UN ANGLE ARCOS(-0.80)=143.13 DEGRES
C     SOIT  UN ANGLE ARCOS(-0.85)=148.21 DEGRES
      DOUBLE PRECISION  PTXYZD(1:4,1:*)
      INTEGER           NOTETR(8,*), N1TETS(MXSOMM)
      INTEGER           LEFACO(11,0:MXFACO),
     %                  N1FASC(1:*),
     %                  NOFAPE(1:MXFAPE),
     %                  NPSOFR(1:MXSOMM),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(1:3,1:MXARCF),
     %                  NAETOI(4,MXETOI),
     %                  LARETE(5,*),
     %                  NOTRCF(MXTRCF),
     %                  NOSTCF(MXTRCF),
     %                  NOSTIS(MXSOMM)
C
      INTEGER           NR(0:6),NS(2)
      DOUBLE PRECISION  D,DMAX,P(4),COS3PD,DIS2ST,DISTAR(4),C,
     %                  SURTRD,S,S1
      INTRINSIC         SQRT
C
C     DECLARATION DE LA QUEUE DES FACES POUR VD1F2F
      MNPIFA = 0
      MXPIFA = 1024
      CALL TNMCDC( 'ENTIER', MXPIFA, MNPIFA )
C
      NBFAPE0 = NBFAPE
      DO 1000 N=1,NBFAPE0
C
C        LE NUMERO LEFACO DE LA FACE PERDUE
         NF = NOFAPE( N )
         IF( NF .LE. 0 ) GOTO 1000

         if( nbsomm .eq. 1178 ) then
            print *,'ajomil19.f nbsomm=',nbsomm
         endif
C
C        FORMATION DU PREMIER CONTOUR FERME FORME DES
C        TRIANGLES PERDUS COPLANAIRES ADJACENTS A PARTIR DE NF
         CALL VDR1CF( COSMAX, NF    , NBSOMM, PTXYZD, N1TETS, NOTETR,
     %                NBFAPE, NOFAPE, MXFACO, LEFACO, N1FASC, NR,
     %                MXETOI, NAETOI,
     %                MXTRCF, NBTRCF, NOTRCF,
     %                MXARCF, NBCF,   N1ARCF, NOARCF,
     %                MXTRCF, NBSTCF, NOSTCF,
     %                MXSOMM, NBSTIS, NOSTIS,
     %                IERR )
         IF( IERR .NE. 0 ) GOTO 9900
C
C        ===============================================================
C        TRAITEMENT SELON LE CONTOUR FERME DES FACES PERDUES COPLANAIRES
C        ===============================================================

         IF( NBTRCF .NE. 1 ) GOTO 20

C        LE CF EST REDUIT A UN TRIANGLE PERDU
C        AJOUT DU BARYCENTRE DU TRIANGLE DE PLUS GRANDE SURFACE
C        ET SUBDIVISION EN 3 FACES DU PLUS GRAND TRIANGLE
C        ------------------------------------------------------
 5       DO I=1,4
            P(I) = ( PTXYZD(I,LEFACO(1,NF))
     %             + PTXYZD(I,LEFACO(2,NF))
     %             + PTXYZD(I,LEFACO(3,NF)) ) / 3D0
         ENDDO
C
C        RECHERCHE DES SOMMETS LES PLUS PROCHES A UNE DISTANCE D'UN PAVE
         CALL DMINSTPAV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE,
     %                   NOPTSUIV, P,      PTXYZD,
     %                   MXSMIN, NBSMIN, NOSMIN, NODSNO, DISTMIN )

C        UNE REFERENCE RELATIVE DE LONGUEUR D'ARETE
         D = SURTRD( PTXYZD(1,LEFACO(1,NF)),
     %               PTXYZD(1,LEFACO(2,NF)),
     %               PTXYZD(1,LEFACO(3,NF)) )
         IF( DISTMIN(1) .LT. D * 1D-4 ) THEN
C           LE BARYCENTRE DE LA FACE NF EST TROP PROCHE D'UN SOMMET
C           CETTE FACE EST OUBLIEE
            GOTO 1000
         ENDIF
C
         IF( NBSOMM .GE. MXSOMM ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:10),'(I10)') MXSOMM
            KERR(1) = 'MAXIMUM DE SOMMETS ' // KERR(MXLGER)(1:10)
     %                 // ' ATTEINT'
            CALL LEREUR
            IERR = 1
            GOTO 9900
         ENDIF

         NBSOMM = NBSOMM + 1
         DO I=1,4
            PTXYZD(I,NBSOMM) = P(I)
         ENDDO
         CALL VD1F3F( NBSOMM, NF, MXFACO, LEFACO, NBFACO,
     %                N1FASC, NBFAPE, MXFAPE, NOFAPE,
     %                NPSOFR, IERR )
         IF( IERR .NE. 0 ) GOTO 9900
C
C        LA FACE NF EST SUPPOSEE N'ETRE PLUS PERDUE
C        APRES L'AJOUT DU BARYCENTRE DE LA FACE NF
         NOFAPE(N) = -NF
         GOTO 1000


 20      IF( NBTRCF .EQ. 2 ) THEN
C           ---------------------------------------
C           DEUX FACES FORMANT UN CF QUADRANGULAIRE
C           ---------------------------------------
C           LES 2 FACES SONT ELLES DES DEMI-FACES C-A-D
C           UN DES 4 SOMMETS EST IL LE MILIEU DE 2 AUTRES ?
            NF  = NOTRCF(1)
            NF1 = NOTRCF(2)
            NR(1) = N1ARCF(1)
            DO 102 I=1,4
C              LE NUMERO NOARCF D'ARETE
               NR(I+1) = NOARCF(2,NR(I))
C              LE NUMERO PTXYZD DU SOMMET I
               NR( I ) = NOARCF(1,NR(I))
 102        CONTINUE
            NR(5) = NR(1)
            NR(6) = NR(2)
            NR(0) = NR(4)
C
            DO 104 I=1,4
               DO 103 JJ=1,3
                  D =( PTXYZD(JJ,NR(I-1)) + PTXYZD(JJ,NR(I+1)) ) * 0.5D0
                  IF( ABS( PTXYZD(JJ,NR(I)) - D ) .GT.
     %                ABS(D) * 1D-9 ) GOTO 104
C                     SOMMET I NON MILIEU DE I-1,I+1
 103           CONTINUE
C
C              ICI : LE SOMMET I EST MILIEU DE I-1 I+1
C              LE NUMERO DE L'ARETE OPPOSEE AU SOMMET I DANS NF ET NF1
               CALL NO1A1F( NR(I-1), NR(I+2), LEFACO(1,NF), NAR )
               IF( NAR .EQ. 0 ) THEN
                  CALL NO1A1F( NR(I+1), NR(I+2), LEFACO(1,NF ), NAR  )
                  CALL NO1A1F( NR(I-1), NR(I+2), LEFACO(1,NF1), NAR1 )
               ELSE
                  CALL NO1A1F( NR(I+1), NR(I+2), LEFACO(1,NF1), NAR1 )
               ENDIF
C
C              LAQUELLE DES 2 ARETES EST LA PLUS LONGUE?
               DISTAR(1) = DIS2ST( PTXYZD(1,LEFACO(NAR,NF)) ,
     %                             PTXYZD(1,LEFACO(MOD(NAR,3)+1,NF)) )
               DISTAR(2) = DIS2ST( PTXYZD(1,LEFACO(NAR1,NF1)) ,
     %                             PTXYZD(1,LEFACO(MOD(NAR1,3)+1,NF1)) )
C
C              AJOUT DU MILIEU DE LA PLUS GRANDE DES 2 ARETES
C              ----------------------------------------------
               IF( DISTAR(1) .LT. DISTAR(2) ) THEN
                  NAR = NAR1
                  NF  = NF1
               ENDIF
               CALL VD1F2F( 1, NAR, NF, MXFACO, LEFACO, NBFACO,
     %                      N1FASC, NBFAPE, NOFAPE,
     %                      NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                      MXPIFA, MNPIFA, IERR )
               IF( IERR .NE. 0 ) GOTO 9900
               GOTO 1000
 104        CONTINUE
C
C           ICI: AUCUN SOMMET N'EST LE MILIEU D'UNE ARETE
C           ---------------------------------------------
C           LE MILIEU A AJOUTER SUR L'ARETE COMMUNE CREERAIT IL
C           UN ANGLE OBTU TROP GRAND DANS L'UNE DES 2 FACES?
C           L'ARETE COMMUNE NE PEUT ETRE QUE NR(1-3) OU NR(2-4)
            N1 = 1
            N2 = 3
            CALL NO1A1F( NR(1), NR(3), LEFACO(1,NF) , NAR  )
            CALL NO1A1F( NR(1), NR(3), LEFACO(1,NF1), NAR1 )
            IF( NAR .LE. 0 .OR. NAR1 .LE. 0 ) THEN
               N1 = 2
               N2 = 4
               CALL NO1A1F( NR(2), NR(4), LEFACO(1,NF) , NAR  )
               CALL NO1A1F( NR(2), NR(4), LEFACO(1,NF1), NAR1 )
            ENDIF
C
C           ARETE COMMUNE : NAR DANS NF ET NAR1 DANS NF1
C           LES 3 COORDONNEES DU MILIEU DE L'ARETE
            DO 106 I=1,4
               P(I) = ( PTXYZD(I,NR(N1)) + PTXYZD(I,NR(N2)) ) * 0.5D0
 106        CONTINUE
C
C           N3 LE 3-EME SOMMET DE NF
            DO 107 I=1,3
               IF( LEFACO(I,NF) .NE. NR(N1) .AND.
     %             LEFACO(I,NF) .NE. NR(N2) ) GOTO 108
 107        CONTINUE
 108        N3 = LEFACO(I,NF)
C
C           N4 LE 3-EME SOMMET DE NF1
            DO 109 I=1,3
               IF( LEFACO(I,NF1) .NE. NR(N1) .AND.
     %             LEFACO(I,NF1) .NE. NR(N2) ) GOTO 110
 109        CONTINUE
 110        N4 = LEFACO(I,NF1)
C
C           CALCUL DU MINIMUM DES ANGLES AUTOUR DU POINT MILIEU P
            DISTAR(1) = COS3PD( P, PTXYZD(1,NR(N2)), PTXYZD(1,N3) )
            DISTAR(2) = COS3PD( P, PTXYZD(1,N3), PTXYZD(1,NR(N1)) )
            DISTAR(3) = COS3PD( P, PTXYZD(1,NR(N1)), PTXYZD(1,N4) )
            DISTAR(4) = COS3PD( P, PTXYZD(1,N4), PTXYZD(1,NR(N2)) )
            C = MIN( DISTAR(1), DISTAR(2), DISTAR(3), DISTAR(4) )
C
            IF( C .GE. COSMIN ) THEN
C
C              ANGLE PLUS PETIT QUE LE MAXIMUM PERMIS
C              LE MILIEU DE L'ARETE COMMUNE EST ACCEPTABLE
               DO NAR=1,3
                  IF( LEFACO(5+NAR,NF) .EQ. NF1 ) GOTO 170
               ENDDO
C
C              CE MILIEU D'ARETE EST IL TROP PROCHE D'UN SOMMET DE LA TETRA?
C              -------------------------------------------------------------
C              RECHERCHE DES SOMMETS LES PLUS PROCHES A UNE DISTANCE D'UN PAVE
 170           CALL DMINSTPAV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE,
     %                         NOPTSUIV, P,      PTXYZD,
     %                         MXSMIN, NBSMIN, NOSMIN, NODSNO, DISTMIN )
C
C              CARRE DE LA LONGUEUR DE L'ARETE COMMUNE
               D = (PTXYZD(1,NR(N2)) - PTXYZD(1,NR(N1))) ** 2
     %           + (PTXYZD(2,NR(N2)) - PTXYZD(2,NR(N1))) ** 2
     %           + (PTXYZD(3,NR(N2)) - PTXYZD(3,NR(N1))) ** 2
C
               IF( DISTMIN(1) .LT. D * 1D-4 ) GOTO 172
C
C              AJOUT DU MILIEU DE L'ARETE COMMUNE ET SUBDIVISION DES 2 FACES
C              -------------------------------------------------------------
               CALL VD1F2F( 1, NAR, NF, MXFACO, LEFACO, NBFACO,
     %                      N1FASC, NBFAPE, NOFAPE,
     %                      NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                      MXPIFA, MNPIFA, IERR )
               IF( IERR .NE. 0 ) GOTO 9900
               GOTO 1000
C
C              MILIEU TROP PROCHE D'UN SOMMET TETRAEDRISE
C              AJOUT DU BARYCENTRE DU TRIANGLE DE PLUS GRANDE SURFACE
C              ET SUBDIVISION EN 3 FACES DU PLUS GRAND TRIANGLE
C              ------------------------------------------------------
 172           S  = SURTRD( PTXYZD(1,LEFACO(1,NF)),
     %                      PTXYZD(1,LEFACO(2,NF)),
     %                      PTXYZD(1,LEFACO(3,NF)) )
               S1 = SURTRD( PTXYZD(1,LEFACO(1,NF1)),
     %                      PTXYZD(1,LEFACO(2,NF1)),
     %                      PTXYZD(1,LEFACO(3,NF1)) )
               NBPAS = 0
C
 174           IF( S1 .GT. S ) THEN
                  I   = NF
                  NF  = NF1
                  NF1 = I
               ENDIF
C
               IF( NBSOMM .GE. MXSOMM ) THEN
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:10),'(I10)') MXSOMM
                  KERR(1) = 'MAXIMUM DE SOMMETS ' // KERR(MXLGER)(1:10)
     %                    // ' ATTEINT'
                  CALL LEREUR
                  IERR = 1
                  GOTO 9900
               ENDIF
               NBSOMM = NBSOMM + 1
               DO I=1,4
                  PTXYZD(I,NBSOMM) = ( PTXYZD(I,LEFACO(1,NF))
     %                               + PTXYZD(I,LEFACO(2,NF))
     %                               + PTXYZD(I,LEFACO(3,NF)) ) / 3D0
               ENDDO
C
C              RECHERCHE DES SOMMETS LES PLUS PROCHES A UNE DISTANCE D'UN PAVE
               CALL DMINSTPAV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE,
     %                         NOPTSUIV, PTXYZD(1,NBSOMM),     PTXYZD,
     %                         MXSMIN, NBSMIN, NOSMIN, NODSNO, DISTMIN )
               IF( DISTMIN(1) .LT. D * 1D-4 ) THEN
                  IF( NBPAS .GE. 1 ) THEN
C                    LE BARYCENTRE DES 2 FACES EST TROP PROCHE D'UN SOMMET
C                    CAS A TRAITER
                     GOTO 178
                  ENDIF
                  NBSOMM = NBSOMM - 1
                  S1=S*2+1
                  NBPAS = NBPAS + 1
                  GOTO 174
               ENDIF
C
               CALL VD1F3F( NBSOMM, NF, MXFACO, LEFACO, NBFACO,
     %                      N1FASC, NBFAPE, MXFAPE, NOFAPE,
     %                      NPSOFR, IERR )
               IF( IERR .NE. 0 ) GOTO 9900
C
C              LA FACE NF1 EST SUPPOSEE N'ETRE PLUS PERDUE
C              APRES L'AJOUT DU BARYCENTRE DE LA FACE NF
               DO I=1,NBFAPE
                  IF( NOFAPE(I) .EQ. NF1 ) THEN
                     NOFAPE(I) = -NF1
                     GOTO 1000
                  ENDIF
               ENDDO
               GOTO 1000
C 
            ENDIF
C
C           ANGLE MAXIMUM AUTOUR DU MILIEU SUPERIEUR AU MAXIMUM PERMIS
C           AJOUT DU MILIEU DU COTE OPPOSE AU PLUS GRAND ANGLE
C           ----------------------------------------------------------
 178        DO 180 I=1,4
               IF( C .EQ. DISTAR(I) ) GOTO 190
 180        CONTINUE
 190        IF( I .EQ. 1 ) THEN
C              LES 2 SOMMETS DE L'ARETE OPPOSEE AU PLUS GRAND ANGLE
               N1 = NR(N2)
               N2 = N3
            ELSE IF( I .EQ. 2 ) THEN
               N2 = NR(N1)
               N1 = N3
            ELSE IF( I .EQ. 3 ) THEN
               N1 = NR(N1)
               N2 = N4
            ELSE
               N1 = N4
               N2 = NR(N2)
            ENDIF
            CALL NO1A1F( N1, N2, LEFACO(1,NF1) , NAR1 )
            IF( NAR1 .EQ. 0 ) THEN
               CALL NO1A1F( N1, N2, LEFACO(1,NF) , NAR1 )
               NF1 = NF
            ENDIF
            CALL VD1F2F( 1, NAR1, NF1, MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   MXPIFA, MNPIFA, IERR )
            IF( IERR .NE. 0 ) GOTO 9900
            GOTO 1000
         ENDIF
C
         IF( NBTRCF .EQ. 3 .AND. NBSTIS .EQ. 0 ) THEN
C           --------------------------------
C           TROIS FACES FORMANT UN PENTAGONE
C           --------------------------------
C           AJOUT DU BARYCENTRE DU PENTAGONE => 5 FACES
            CALL VDCFCX( PTXYZD, N1ARCF(1), NOARCF, SINMIN )
            IF( SINMIN .LT. 0.1 ) THEN
C              LE CONTOUR FERME EST JUGE NON CONVEXE OU
C              PRESENTANT DES ARETES DANS LE PROLONGEMENT L'UNE DE L'AUTRE
               GOTO 900
            ENDIF
            CALL VD3F5F( NOTRCF, N1ARCF, NOARCF,
     %                   MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR, IERR )
            IF( IERR .NE. 0 ) GOTO 9900
            GOTO 1000
         ENDIF
C
         IF( NBTRCF .GE. 4 .AND. NBTRCF .LE. 5 .AND. NBSTIS .EQ. 0) THEN
C           --------------------------------------------------------
C           4 A 5 FACES FORMENT UN POLYGONE SANS POINT INTERNE ISOLE
C           --------------------------------------------------------
            CALL VDCFCX( PTXYZD, N1ARCF(1), NOARCF, SINMIN )
            IF( SINMIN .LT. 0.1 ) THEN
C              LE CONTOUR FERME EST JUGE NON CONVEXE OU
C              PRESENTANT DES ARETES DANS LE PROLONGEMENT L'UNE DE L'AUTRE
               GOTO 900
            ENDIF
C           AJOUT DU BARYCENTRE DU POLYGONE =>
            CALL VDMFNF( NBTRCF, NOTRCF, N1ARCF, NOARCF,
     %                   MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR, IERR )
            IF( IERR .EQ. 10 ) THEN
C              LE BARYCENTRE N'ETOILE PAS LE POLYGONE
               IERR = 0
               GOTO 900
            ENDIF
            IF( IERR .NE. 0  ) GOTO 9900
            GOTO 1000
         ENDIF
C
         IF( NBTRCF .GE. 3 .AND. NBSTIS .EQ. 1 ) THEN
C           ----------------------------------------------------------------
C           AJOUT DU MILIEU DE LA PLUS GRANDE ARETE OPPOSEE AU POINT INTERNE
C           ----------------------------------------------------------------
            DMAX = 0D0
            DO 430 I=1,NBTRCF
               NF = NOTRCF( I )
               DO 410 J=1,3
C                 RECHERCHE DU SOMMET NOSTIS(1)
                  IF( LEFACO(J,NF) .EQ. NOSTIS(1) ) GOTO 420
 410           CONTINUE
               GOTO 430
C              L'ARETE OPPOSEE AU SOMMET INTERNE AU CF
 420           IF( J .EQ. 3 ) THEN
                  J = 1
               ELSE
                  J = J + 1
               ENDIF
               N1 = LEFACO(J,NF)
               N2 = LEFACO(MOD(J,3)+1,NF)
               D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %            + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %            + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
               IF( D .GT. DMAX ) THEN
                  DMAX  = D
                  NAR   = J
                  NFMAX = NF
               ENDIF
 430        CONTINUE
C           AJOUT DU MILIEU DE LA PLUS GRANDE ARETE DES TRIANGLES PERDUS
            CALL VD1F2F( 1, NAR, NFMAX,  MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   MXPIFA, MNPIFA, IERR )
            IF( IERR .NE. 0  ) GOTO 9900
            GOTO 1000
C
CCC
CCC            CALL VDCFCX( PTXYZD, N1ARCF(1), NOARCF, SINMIN )
CCC            IF( SINMIN .LT. 0.1 ) THEN
CCCC              LE CONTOUR FERME EST JUGE NON CONVEXE OU
CCCC              PRESENTANT DES ARETES DANS LE PROLONGEMENT L'UNE DE L'AUTRE
CCC               GOTO 900
CCC            ENDIF
CCCC           BARYCENTRE JOINT AUX 4 MILIEUX DES ARETES ET 2 SOMMETS OPPOSES
CCC            CALL VD4F8F( NOSTIS(1), NOTRCF, N1ARCF, NOARCF,
CCC     %                   MXFACO, LEFACO, NBFACO,
CCC     %                   N1FASC, NBFAPE, NOFAPE,
CCC     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
CCC     %                   MXPIFA, MNPIFA, IERR )
         ENDIF
C
C        --------------------------------
C        TRAITEMENT GENERAL DE LA FACE NF
C        --------------------------------
CCCC        AJOUT DU MILIEU DE LA PLUS GRANDE ARETE PERDUE DE LA PLUS
CCCC        GRANDE FACE PERDUE NF ET SI PAS D'ARETE PERDUE AJOUT DU
CCCC        BARYCENTRE DE LA PLUS GRANDE ARETE
CCCC
CCCC        TRI SELON LES SURFACES DECROISSANTES DES FACES PERDUES
CCC 900     IF( NBTRCF.GT.1 ) CALL TRIFAP( NBTRCF, NOTRCF, LEFACO, PTXYZD)
CCCC
CCCC        TRI CROISSANT DES ARETES DE LA FACE NF
CCC         NF = NOTRCF(1)
CCC         DO 910 I=1,3
CCC            DISTAR(I) = DIS2ST( PTXYZD(1,LEFACO(I,NF)) ,
CCC     %                          PTXYZD(1,LEFACO(MOD(I,3)+1,NF)) )
CCC 910     CONTINUE
CCCC
CCCC        RECHERCHE ET TRAITEMENT DE LA PLUS LONGUE ARETE NAR
CCC         D = 0
CCC         DO 920 I=1,3
CCC            IF( D .LT. DISTAR(I) ) THEN
CCC               D   = DISTAR(I)
CCC               NAR = I
CCC            ENDIF
CCC 920     CONTINUE
CCCC
CCCC        AJOUT DU MILIEU DE LA PLUS GRANDE ARETE DES TRIANGLES DU CF
CCC 900     DMAX = 0D0
CCC         DO 920 J=1,NBTRCF
CCC            NF = NOTRCF( J )
CCC            DO 910 I=1,3
CCC               N1 = LEFACO(I,NF)
CCC               N2 = LEFACO(MOD(I,3)+1,NF)
CCC               D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
CCC     %            + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
CCC     %            + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
CCC               IF( D .GT. DMAX ) THEN
CCC                  DMAX  = D
CCC                  NAR   = I
CCC                  NFMAX = NF
CCC               ENDIF
CCC 910        CONTINUE
CCC 920     CONTINUE
C
C        RECHERCHE DE LA PLUS GRANDE ARETE COMMUNE A 2 TRIANGLES
C        GENERATION DU TABLEAU DE HACHAGE
C        LARETE : TABLEAU DE HACHAGE DES ARETES DES OT
C        LARETE(1,J) = SOMMET 1 DANS PTXYZD
C        LARETE(2,J) = SOMMET 2 DANS PTXYZD
C        LARETE(3,J) = NUMERO DANS LARETE DE L'ARETE SUIVANTE
C        LARETE(4,J) = NUMERO DU 1-ER  TRIANGLE
C        LARETE(5,J) = NUMERO DU 2-EME TRIANGLE
C        HACHAGE(NS1,NS2) = (NS1 + NS2) MODULO MXARET + 1
C
C        NOMBRE MAXIMAL D'ARETES
 900     MXARET = 3 * NBTRCF
         DO 910 I=1,MXARET
            LARETE(1,I) = 0
            LARETE(3,I) = 0
            LARETE(5,I) = 0
 910     CONTINUE
C
C        LA 1-ERE ARETE LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
         LAVIDE = MXARET
C
         DO 930 J=1,NBTRCF
            NF = NOTRCF( J )
            DO 920 I=1,3
               NS(1) = LEFACO(I,NF)
               NS(2) = LEFACO(MOD(I,3)+1,NF)
               IF( NS(1) .GT. NS(2) ) THEN
                  K     = NS(1)
                  NS(1) = NS(2)
                  NS(2) = K
               ENDIF
C
C              HACHAGE DE L'ARETE POUR LA RETROUVER
               CALL HACHAG( 2, NS, 5, MXARET, LARETE, 3, LAVIDE, NOAR )
C
               IF( NOAR .EQ. 0 )  THEN
C                 SATURATION DES ARETES DE LARETE
                  NBLGRC(NRERR) = 1
                  KERR(1) ='SP AJOMIL19: TABLEAU LARETE SATURE'
                  CALL LEREUR
                  IERR = 3
                  GOTO 9900
               ENDIF
C
               IF( NOAR .LT. 0 ) THEN
C                 L'ARETE EST CREEE   1-ER TRIANGLE
                  LARETE( 4, -NOAR ) = NF
               ELSE
C                 L'ARETE EXISTE   2-EME TRIANGLE
                  LARETE( 5, NOAR ) = NF
               ENDIF
 920        CONTINUE
 930     CONTINUE
C
C        RECHERCHE DE LA PLUS GRANDE ARETE COMMUNE A 2 TRIANGLES PERDUS
C        ET COPLANAIRES
         DMAX = 0D0
         DO 950 I=1,MXARET
            IF( LARETE(1,I) .EQ. 0 ) GOTO 950
            IF( LARETE(5,I) .EQ. 0 ) GOTO 950
C           ARETE DOUBLE
            N1 = LARETE( 1, I )
            N2 = LARETE( 2, I )
            D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %         + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %         + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
            IF( D .GT. DMAX ) THEN
               DMAX = D
               NOAR = I
            ENDIF
 950     CONTINUE
C
         IF( DMAX .LE. 0D0 ) THEN
C           PAS D'ARETE DOUBLE
            DO 960 I=1,MXARET
               IF( LARETE(1,I) .EQ. 0 ) GOTO 960
               N1 = LARETE( 1, I )
               N2 = LARETE( 2, I )
               D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %            + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %            + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
               IF( D .GT. DMAX ) THEN
                  DMAX = D
                  NOAR = I
               ENDIF
 960        CONTINUE
         ENDIF
C
C        RECHERCHE DU NUMERO D'ARETE DANS LE TRIANGLE
         NF = LARETE( 4, NOAR )
         N1 = LARETE( 1, NOAR )
         N2 = LARETE( 2, NOAR )
         CALL NO1A1F( N1, N2, LEFACO(1,NF), NOAR )

C        LES 3 COORDONNEES DU MILIEU DE L'ARETE
         DO I=1,4
            P(I) = ( PTXYZD(I,N1) + PTXYZD(I,N2) ) * 0.5D0
         ENDDO

C        CE MILIEU D'ARETE EST IL TROP PROCHE D'UN SOMMET DE LA TETRA?
C        RECHERCHE DES SOMMETS LES PLUS PROCHES A UNE DISTANCE D'UN PAVE
         CALL DMINSTPAV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE,
     %                   NOPTSUIV, P,      PTXYZD,
     %                   MXSMIN, NBSMIN, NOSMIN, NODSNO, DISTMIN )
C
C        CARRE DE LA LONGUEUR DE L'ARETE COMMUNE
         D = (PTXYZD(1,N2) - PTXYZD(1,N1)) ** 2
     %     + (PTXYZD(2,N2) - PTXYZD(2,N1)) ** 2
     %     + (PTXYZD(3,N2) - PTXYZD(3,N1)) ** 2
C
         IF( DISTMIN(1) .LT. D * 1D-4 ) GOTO 5
C
C        AJOUT DU MILIEU DE LA PLUS GRANDE ARETE DES TRIANGLES PERDUS
         CALL VD1F2F( 1, NOAR, NF,    MXFACO, LEFACO, NBFACO,
     %                N1FASC, NBFAPE, NOFAPE,
     %                NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                MXPIFA, MNPIFA, IERR )
         IF( IERR .NE. 0 ) GOTO 9900

 1000 CONTINUE
C
C     DESTRUCTION DE LA QUEUE
 9900 CALL TNMCDS( 'ENTIER', MXPIFA, MNPIFA )
      RETURN
      END
