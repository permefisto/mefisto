      SUBROUTINE AJOMILFP( KNMVOL, COANPL, NFLPER,
     %                     MXFACO, LEFACO, NBFACO, N1FASC,
     %                     NBFAPE, MXFAPE, NOFAPE,
     %                     NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     MXTRCF, NBTRCF, NOTRCF, NBSTIS, NOSTIS,
     %                     MXARCF, N1ARCF, NOARCF,
     %                     MXETOI, LARETE, MXPIFA, MNPIFA, IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER A CERTAINES ARETES DES FACES PERDUES LEUR MILIEU et
C -----    DECOUPER EN 2 FACES TOUTES LES FACES DE LEFACO AYANT
C          CETTE ARETE PERDUE COMMUNE

C ENTREES :
C ---------
C KNMVOL : NOM DU VOLUME A TETRAEDRISER
C COANPL : SEUIL DU COSINUS DE L'ANGLE FORME PAR LES NORMALES AUX
C          2 FACES ET AU DESSUS DUQUEL LES FACES
C          SONT CONSIDEREES COPLANAIRES
C NFLPER : NUMERO LEFACO DE LA FACE PERDUE A TRAITER
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO

C MODIFIES :
C ----------
C LEFACO : LES 3 SOMMETS, 2 MATERIAUX, 3 FACES VOISINES ET CHAINAGE
C          DES FACES TRIANGULAIRES DU CONTOUR ET INTERFACES
C NBFACO ; NOMBRE DE FACES ACTIVES DE LEFACO
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS
C NBFAPE : NOMBRE DE FACES PERDUES
C NOFAPE : NUMERO DES FACES PERDUES
C NBSOMM : NUMERO DU DERNIER SOMMET AJOUTE A LA TETRAEDRISATION
C          MIS A JOUR EN SORTIE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TETRAEDRISATION
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : =  0  SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C LARETE : TABLEAU DE HACHAGE DES ARETES DES TRIANGLES PERDUS COPLANAIRES
C          LARETE(1,J) = SOMMET 1 DANS PTXYZD
C          LARETE(2,J) = SOMMET 2 DANS PTXYZD
C          LARETE(3,J) = NUMERO DANS LARETE DE L'ARETE SUIVANTE
C          LARETE(4,J) = NUMERO DU 1-ER  TRIANGLE
C          LARETE(5,J) = NUMERO DU 2-EME TRIANGLE
C          HACHAGE(NS1,NS2) = (NS1 + NS2) MODULO MXARET + 1
C ATTENTION: MEME MC QUE NAETOI ( TABLEAUX EN EQUIVALENCE A L'APPEL )
C MXTRCF : NOMBRE MAXIMAL DECLARABLE D'ARETES OU TRIANGLES
C          OU SOMMETS DANS L'ETOILE
C MXARCF : MAXIMUM D'ARETES DECLARABLES DANS N1ARCF et NOARCF
C N1ARCF : TABLEAU (0:MXARCF) AUXILIAIRE
C NOARCF : NUMERO DES ARETES DE LA LIGNE DU CONTOUR FERME SELON UN SENS
C NBTRCF : NOMBRE DE FACES TRIANGULAIRES PERDUES DU CF
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DU CF
C NBSTIS : NOMBRE DE SOMMETS ISOLES DANS LE CF
C NOSTIS : NUMERO DES SOMMETS ISOLES N'APPARTENANT PAS AU CONTOUR

C SORTIES:
C --------
C MXPIFA : NOMBRE D'ENTIERS DU TABLEAU PIFA
C MNPIFA : ADRESSE MCN DU TABLEAU PIFA
C IERR   : 0 SI PAS D'ERREUR
C          1 SI SATURATION DES SOMMETS DE PTXYZD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1993
C MODFS  : ALAIN PERRONNET  Saint PIERRE du PERRAY         Decembre 2019
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
C     ANGLE AU DELA DUQUEL IL FAUT LE DECOUPER
      DOUBLE PRECISION  COSMIN
      PARAMETER        (COSMIN=-0.8D0)
C     SOIT  UN ANGLE ARCOS(-0.76)=139.46 DEGRES
C     SOIT  UN ANGLE ARCOS(-0.80)=143.13 DEGRES
C     SOIT  UN ANGLE ARCOS(-0.85)=148.21 DEGRES

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      DOUBLE PRECISION  PTXYZD(1:4,1:MXSOMM)
      INTEGER           LEFACO(11,0:MXFACO),
     %                  N1FASC(1:MXSOMM),
     %                  NOFAPE(1:MXFAPE),
     %                  NPSOFR(1:MXSOMM),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(1:3,1:MXARCF),
     %                  LARETE(5,MXETOI),
     %                  NOTRCF(MXTRCF),
     %                  NOSTIS(MXSOMM)
      CHARACTER*24      KNMVOL
      CHARACTER*80      KTITRE
      INTEGER           NR(0:6), NS(2), NO0FAR(1)
      DOUBLE PRECISION  D, DMAX, P(3), COS3PD, DIS2ST, DISTAR(4), C

      TRACTE0 = TRACTE
      LORBITE = 1
      NBSOM0  = NBSOMM

ccc      NFLPER  = NOTRCF(1)
      PRINT*
      PRINT*,'ajomilfp: Debut Volume: ',KNMVOL,' Face perdue',NFLPER,
     %       ' NBTRCF=',NBTRCF,' COSMIN=',COSMIN,' COANPL=',COANPL,
     %       ' NBFAPE=',NBFAPE,' FACES PERDUES  NBSOM0=',NBSOM0


      if( nflper .eq. 1132298 ) then
         PRINT*,'ajomilfp: debog notetr(865330)'
         tracte=.true.
         KTITRE='ajomilfp: debog notetr(865330)'
         CALL TRDUCF( KTITRE, PTXYZD, NBTRCF, NOTRCF, NBSTIS, NOSTIS,
     %                LEFACO, NO0FAR )
         print*
      endif


cccC     TRACE DES FACES LEFACO AVANT AJOUT DE POINTS MILIEUX
ccc      tracte = .true.
ccc      CALL TRLEFACO( KNMVOL, MXFACO, LEFACO, PTXYZD )

C     ===============================================================
C     TRAITEMENT SELON LE CONTOUR FERME DES FACES PERDUES COPLANAIRES
C     ===============================================================

      IF( NBTRCF .EQ. 0 ) THEN
C           ---------------------
C           NON EXECUTE: BARYCENTRE EN TEST... AU DESSOUS
C           ---------------------
C           LAQUELLE DES 3 ARETES DE NFLPER EST LA PLUS LONGUE?
            DMAX = -1D0
            DO I=1,3
               IF( I .LT. 3 ) THEN
                  I1 = I + 1
               ELSE
                  I1 = 1
               ENDIF
               D = DIS2ST( PTXYZD( 1, LEFACO( I +1, NFLPER ) ),
     %                     PTXYZD( 1, LEFACO( I1+1, NFLPER ) ) )
               IF( D .GT. DMAX ) THEN
                  DMAX = D
                  NAR  = I
               ENDIF
            ENDDO

C           AJOUT DU POINT MILIEU DE LA PLUS GRANDE DES 3 ARETES DE NFLPER
C           --------------------------------------------------------------
            CALL VD1F2F( 0, NAR, NFLPER, MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   MXPIFA, MNPIFA, IERR )
ccc            IF( IERR .NE. 0 ) GOTO 9900
            GOTO 9900

      ENDIF


      IF( NBTRCF .EQ. 1 ) THEN
C        -------------------------------------------
C        AJOUT DU BARYCENTRE DE LA FACE TRIANGULAIRE
C        -------------------------------------------

C        CREATION DES COORDONNEES DU POINT MILIEU DE L'ARETE NAR DE NF
         IF( NBSOMM .GE. MXSOMM ) THEN
          PRINT*,'ajomilfp: SATURATION DU TABLEAU PTXYZD MXSOMM=',MXSOMM
            IERR = 1
            GOTO 9999
         ENDIF
         NBSOMM = NBSOMM + 1
         DO I=1,4
            PTXYZD(I,NBSOMM) = ( PTXYZD(I,LEFACO(1,NFLPER)) +
     %                           PTXYZD(I,LEFACO(2,NFLPER)) +
     %                           PTXYZD(I,LEFACO(3,NFLPER)) ) / 3D0
         ENDDO
C        BARYCENTRE D'UNE FACE LEFACO
         NPSOFR(NBSOMM) = 2

         CALL VD1F3F( NBSOMM, NFLPER, MXFACO, LEFACO, NBFACO,
     %                N1FASC, NBFAPE, MXFAPE, NOFAPE,
     %                NPSOFR, IERR )
         GOTO 9900

      ENDIF


      IF( NBTRCF .EQ. 2 ) THEN
C           ---------------------------------------
C           DEUX FACES FORMANT UN CF QUADRANGULAIRE
C           ---------------------------------------
C           LES 2 FACES SONT ELLES DES DEMI-FACES C-A-D
C           UN DES 4 SOMMETS EST IL LE MILIEU DE 2 AUTRES ?
            NF  = NOTRCF(1)
            NF1 = NOTRCF(2)
            NR(1) = N1ARCF(1)
            DO I=1,4
C              LE NUMERO NOARCF D'ARETE
               NR(I+1) = NOARCF(2,NR(I))
C              LE NUMERO PTXYZD DU SOMMET I
               NR( I ) = NOARCF(1,NR(I))
            ENDDO
            NR(5) = NR(1)
            NR(6) = NR(2)
            NR(0) = NR(4)

            DO 104 I=1,4
               DO JJ=1,3
                  D =( PTXYZD(JJ,NR(I-1)) + PTXYZD(JJ,NR(I+1)) ) * 0.5D0
                  IF( ABS( PTXYZD(JJ,NR(I)) - D ) .GT.
     %                ABS(D) * 1D-9 ) GOTO 104
C                     SOMMET I NON MILIEU DE I-1,I+1
               ENDDO
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
C              AJOUT DU POINT MILIEU DE LA PLUS GRANDE DES 2 ARETES
C              ----------------------------------------------------
               IF( DISTAR(1) .LT. DISTAR(2) ) THEN
                  NAR = NAR1
                  NF  = NF1
               ENDIF
               CALL VD1F2F( 0, NAR, NF, MXFACO, LEFACO, NBFACO,
     %                      N1FASC, NBFAPE, NOFAPE,
     %                      NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                      MXPIFA, MNPIFA, IERR )
ccc               IF( IERR .NE. 0 ) GOTO 9900
               GOTO 9900
 104        ENDDO

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

C           ARETE COMMUNE : NAR DANS NF ET NAR1 DANS NF1
C           LES 3 COORDONNEES DU MILIEU DE L'ARETE
            DO I=1,3
               P(I) = ( PTXYZD(I,NR(N1)) + PTXYZD(I,NR(N2)) ) * 0.5D0
            ENDDO

C           N3 LE 3-EME SOMMET DE NF
            DO I=1,3
               IF( LEFACO(I,NF) .NE. NR(N1) .AND.
     %             LEFACO(I,NF) .NE. NR(N2) ) GOTO 108
            ENDDO
 108        N3 = LEFACO(I,NF)

C           N4 LE 3-EME SOMMET DE NF1
            DO I=1,3
               IF( LEFACO(I,NF1) .NE. NR(N1) .AND.
     %             LEFACO(I,NF1) .NE. NR(N2) ) GOTO 110
            ENDDO
 110        N4 = LEFACO(I,NF1)

C           CALCUL DU MINIMUM DES ANGLES AUTOUR DU POINT MILIEU P
            DISTAR(1) = COS3PD( P, PTXYZD(1,NR(N2)), PTXYZD(1,N3) )
            DISTAR(2) = COS3PD( P, PTXYZD(1,N3), PTXYZD(1,NR(N1)) )
            DISTAR(3) = COS3PD( P, PTXYZD(1,NR(N1)), PTXYZD(1,N4) )
            DISTAR(4) = COS3PD( P, PTXYZD(1,N4), PTXYZD(1,NR(N2)) )
            C = MIN( DISTAR(1), DISTAR(2), DISTAR(3), DISTAR(4) )

            IF( C .GE. COSMIN ) THEN

C              ANGLE PLUS PETIT QUE LE MAXIMUM PERMIS
C              LE MILIEU DE L'ARETE COMMUNE EST ACCEPTABLE
               DO NAR=1,3
                  IF( LEFACO(5+NAR,NF) .EQ. NF1 ) GOTO 170
               ENDDO

C              AJOUT DU MILIEU DE L'ARETE COMMUNE ET SUBDIVISION DES 2 FACES
C              -------------------------------------------------------------
 170           CALL VD1F2F( 0, NAR, NF, MXFACO, LEFACO, NBFACO,
     %                      N1FASC, NBFAPE, NOFAPE,
     %                      NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                      MXPIFA, MNPIFA, IERR )
ccc               IF( IERR .NE. 0 ) GOTO 9900
               GOTO 9900
            ENDIF

C           ANGLE MAXIMUM AUTOUR DU MILIEU SUPERIEUR AU MAXIMUM PERMIS
C           AJOUT DU MILIEU DU COTE OPPOSE AU PLUS GRAND ANGLE
C           ----------------------------------------------------------
            DO I=1,4
               IF( C .EQ. DISTAR(I) ) GOTO 190
            ENDDO
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
            CALL VD1F2F( 0, NAR1, NF1, MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   MXPIFA, MNPIFA, IERR )
ccc            IF( IERR .NE. 0 ) GOTO 9900
            GOTO 9900
      ENDIF


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
ccc            IF( IERR .NE. 0 ) GOTO 9900
            GOTO 9900
      ENDIF


      IF( (NBTRCF .EQ. 4 .OR. NBTRCF .EQ. 5) .AND. NBSTIS .EQ. 0) THEN
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

            GOTO 9900
      ENDIF


ccc      IF( NBTRCF .GE. 3 .AND. NBSTIS .EQ. 1 ) THEN
      IF( NBTRCF .GE. 3 .AND. NBSTIS .LT. 0 ) THEN
C           inactif
C           -------------------------------------------------------------------
cccC           AJOUT DU MILIEU DE LA PLUS GRANDE ARETE OPPOSEE AU POINT INTERNE
C           AJOUT DU MILIEU DE LA PLUS GRANDE ARETE DES NBTRCF TRIANGLES DU CF
C           -------------------------------------------------------------------
            DMAX = 0D0
            DO 430 I=1,NBTRCF
               NF = NOTRCF( I )

ccc               DO J=1,3
cccC                 RECHERCHE DU SOMMET NOSTIS(1)
ccc                  IF( LEFACO(J,NF) .EQ. NOSTIS(1) ) GOTO 420
ccc               ENDDO
ccc               GOTO 430
cccC              L'ARETE OPPOSEE AU SOMMET INTERNE AU CF
ccc 420           IF( J .EQ. 3 ) THEN
ccc                  J = 1
ccc               ELSE
ccc                  J = J + 1
ccc               ENDIF

               DO J=1,3
                  N1 = LEFACO(J,NF)
                  N2 = LEFACO(MOD(J,3)+1,NF)
                  D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %               + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %               + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
                  IF( D .GT. DMAX ) THEN
                     DMAX  = D
                     NAR   = J
                     NFMAX = NF
                  ENDIF
               ENDDO

 430        ENDDO
C           AJOUT DU MILIEU DE LA PLUS GRANDE ARETE DES TRIANGLES PERDUS
            CALL VD1F2F( 0, NAR, NFMAX,  MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NBFAPE, NOFAPE,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   MXPIFA, MNPIFA, IERR )
ccc            IF( IERR .NE. 0  ) GOTO 9900
            GOTO 9900

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

C        ------------------------------------
C        TRAITEMENT GENERAL DE LA FACE NFLPER
C        ------------------------------------

CCCC        AJOUT DU MILIEU DE LA PLUS GRANDE ARETE PERDUE DE LA PLUS
CCCC        GRANDE FACE PERDUE NF ET SI PAS D'ARETE PERDUE AJOUT DU
CCCC        BARYCENTRE DE LA PLUS GRANDE ARETE
CCCC
CCCC        TRI SELON LES SURFACES DECROISSANTES DES FACES PERDUES
CCC 900     IF( NBTRCF.GT.1 ) CALL TRIFAP( NBTRCF, NOTRCF, LEFACO, PTXYZD)

CCCC        TRI CROISSANT DES ARETES DE LA FACE NF
CCC         NF = NOTRCF(1)
CCC         DO 910 I=1,3
CCC            DISTAR(I) = DIS2ST( PTXYZD(1,LEFACO(I,NF)) ,
CCC     %                          PTXYZD(1,LEFACO(MOD(I,3)+1,NF)) )
CCC 910     ENDDO

CCCC        RECHERCHE ET TRAITEMENT DE LA PLUS LONGUE ARETE NAR
CCC         D = 0
CCC         DO 920 I=1,3
CCC            IF( D .LT. DISTAR(I) ) THEN
CCC               D   = DISTAR(I)
CCC               NAR = I
CCC            ENDIF
CCC 920     ENDDO

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
CCC 910        ENDDO
CCC 920     ENDDO


 900  IF( NBSTIS .LT. 0 ) GOTO 970

C     RECHERCHE DE LA PLUS GRANDE ARETE COMMUNE A 2 TRIANGLES DU CF
C     -------------------------------------------------------------
C     GENERATION DU TABLEAU DE HACHAGE
C     LARETE : TABLEAU DE HACHAGE DES ARETES DES OT
C     LARETE(1,J) = SOMMET 1 DANS PTXYZD
C     LARETE(2,J) = SOMMET 2 DANS PTXYZD
C     LARETE(3,J) = NUMERO DANS LARETE DE L'ARETE SUIVANTE
C     LARETE(4,J) = NUMERO DU 1-ER  TRIANGLE
C     LARETE(5,J) = NUMERO DU 2-EME TRIANGLE
C     HACHAGE(NS1,NS2) = (NS1 + NS2) MODULO MXARET + 1
C     NOMBRE MAXIMAL D'ARETES
      MXARET = 3 * NBTRCF
      DO I=1,MXARET
         LARETE(1,I) = 0
         LARETE(3,I) = 0
         LARETE(5,I) = 0
      ENDDO

C     LA 1-ERE ARETE LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
      LAVIDE = MXARET

      DO J=1,NBTRCF
         NF = NOTRCF( J )
         DO I=1,3
            NS(1) = LEFACO(I,NF)
            IF( NS(1) .LE. 0 ) GOTO 9900
            NS(2) = LEFACO(MOD(I,3)+1,NF)
            IF( NS(1) .GT. NS(2) ) THEN
               K     = NS(1)
               NS(1) = NS(2)
               NS(2) = K
            ENDIF

C           HACHAGE DE L'ARETE POUR LA RETROUVER
            CALL HACHAG( 2, NS, 5, MXARET, LARETE, 3, LAVIDE, NOAR )

            IF( NOAR .EQ. 0 )  THEN
C              SATURATION DES ARETES DE LARETE
               NBLGRC(NRERR) = 1
               KERR(1) ='ajomilfp: TABLEAU LARETE SATURE a AUGMENTER'
               CALL LEREUR
               IERR = 3
               GOTO 9900
            ENDIF

            IF( NOAR .LT. 0 ) THEN
C              L'ARETE EST CREEE   1-ER TRIANGLE
               LARETE( 4, -NOAR ) = NF
            ELSE
C              L'ARETE EXISTE   2-EME TRIANGLE
               LARETE( 5, NOAR ) = NF
            ENDIF
         ENDDO
      ENDDO

C     RECHERCHE DE LA PLUS GRANDE ARETE COMMUNE A 2 TRIANGLES PERDUS
      NOAR = 0
      DMAX = 0D0
      DO 950 I=1,MXARET
         IF( LARETE(1,I) .EQ. 0 ) GOTO 950
         IF( LARETE(5,I) .EQ. 0 ) GOTO 950
C        ARETE DOUBLE
         N1 = LARETE( 1, I )
         N2 = LARETE( 2, I )
         D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %      + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %      + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
         IF( D .GT. DMAX ) THEN
            DMAX = D
            NOAR = I
         ENDIF
 950  ENDDO

      IF( DMAX .LE. 0D0 ) THEN
C        PAS D'ARETE DOUBLE
         DO 960 I=1,MXARET
            IF( LARETE(1,I) .EQ. 0 ) GOTO 960
            N1 = LARETE( 1, I )
            N2 = LARETE( 2, I )
            D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %         + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %         + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
            IF( D .GT. DMAX ) THEN
               DMAX = D
               NOAR = I
            ENDIF
 960     ENDDO
      ENDIF

C     RECHERCHE DU NUMERO D'ARETE DANS LE TRIANGLE
      NF = LARETE( 4, NOAR )
      N1 = LARETE( 1, NOAR )
      N2 = LARETE( 2, NOAR )
      CALL NO1A1F( N1, N2, LEFACO(1,NF), NOAR )
      GOTO 980


C     RECHERCHE DE LA PLUS GRANDE ARETE DES NBTRCF TRIANGLES DU CF
C     ------------------------------------------------------------
C     inactif
 970  DMAX = 0
      DO J=1,NBTRCF
         NF = NOTRCF( J )
         DO I=1,3
            N1 = LEFACO(I,NF)
            N2 = LEFACO(MOD(I,3)+1,NF)
            D  = ( PTXYZD(1,N2)-PTXYZD(1,N1) ) ** 2
     %         + ( PTXYZD(2,N2)-PTXYZD(2,N1) ) ** 2
     %         + ( PTXYZD(3,N2)-PTXYZD(3,N1) ) ** 2
            IF( D .GT. DMAX ) THEN
               DMAX  = D
               NOAR  = I
               NFMAX = NF
            ENDIF
         ENDDO
      ENDDO

C     CHOIX DU MILIEU DE LA PLUS GRANDE ARETE
      NF = NFMAX


C     AJOUT DU MILIEU DE LA PLUS GRANDE ARETE DES TRIANGLES PERDUS
C     ------------------------------------------------------------
 980  CALL VD1F2F( 0, NOAR,   NF,  MXFACO, LEFACO, NBFACO,
     %             N1FASC, NBFAPE, NOFAPE,
     %             NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %             MXPIFA, MNPIFA, IERR )
ccc      IF( IERR .NE. 0 ) GOTO 9900


 9900 PRINT*,'ajomilfp: Fin   Volume: ',KNMVOL,' Face perdue',NFLPER,
     %       ' AJOUT des MILIEUX',NBSOM0+1,' a',NBSOMM,' IERR=',IERR
      PRINT*

C     TRACE DES FACES LEFACO APRES AJOUT DES POINTS MILIEUX
ccc      tracte = .true.
      CALL TRFAPEPT( NBFAPE, NOFAPE, MXFACO, LEFACO,
     %               NBSOM0+1, NBSOMM, PTXYZD )

 9999 tracte = tracte0
      RETURN
      END
