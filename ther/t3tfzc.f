      SUBROUTINE T3TFZC( NBCOOR, NBPOIN, XYZPOI, NBTRFF, NSTRFF, XYZVQ2,
     %                   NCAS0,  NCAS1,  NCAS,   NTYP,   TEMPER, dptemp,
     %                   TMIN,   TMAX,   NUTRFF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES TRIANGLES FRONTALIERS PAR ZONES DE COULEURS
C ----     SELON LA TEMPERATURE ET L'ALGORITHME DU PEINTRE
C
C ENTREES:
C --------
C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT OU NOEUD (3 ou 6)
C NBPOIN : NOMBRE TOTAL DE NOEUDS=POINTS OU EST CONNUE LA TEMPERATURE
C XYZPOI : XYZ DES NBPOIN DU MAILLAGE
C NBTRFF : NOMBRE DE TRIANGLES A TRACER
C NSTRFF : NUMEROS DES 3 SOMMETS DES NBTRFF TRIANGLES
C          SI >NBPOIN ALORS A RECHERCHER DANS XYZVQ2
C          SINON  A CHERCHER DANS XYZPOI
C XYZVQ2 : XYZ ET TEMPERATURE DES FACES Q2
C NCAS0:NCAS1: CAS DE CALCULS DE LA TEMPERATURE
C NCAS   : NUMERO DU CAS A TRACER
C TEMPER : LES NBPOIN * NCAS0:NCAS1 TEMPERATURES
C TMIN   : TEMPERATURE MIN
C TMAX   : TEMPERATURE MAX
C NUTRFF : NUTRFF(1)=NUMERO DU TRIANGLE LE PLUS ELOIGNE DE L'OEIL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1994
C2345X7..............................................................012
      PARAMETER     ( LIGCON=0 )
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      INTEGER           NSTRFF(3,NBTRFF), NUTRFF(NBTRFF)
      REAL              XYZPOI(NBCOOR,NBPOIN), XYZVQ2(4,*)
      REAL              XYZ(3,3), XYZP(3), COUL(3)
      DOUBLE PRECISION  TEMPER(NBPOIN,NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp

C
C     LE NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL
C
C     TRAIT TRACES AVEC UNE EPAISSEUR DE 1
      CALL XVEPAISSEUR( 1 )
C
C     ARETES TIRETEES
      CALL XVTYPETRAIT( LIGCON )
C
C     LE TRACE DES FACES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     =======================================================
      DO 100 NF = 1, NBTRFF
C
C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NT = NUTRFF( NF )
C
C        LES COORDONNEES DES 3 SOMMETS DU TRIANGLE NT
         DO 30 I=1,3
C           LE NUMERO DU SOMMET DU TRIANGLE NT
            NS = NSTRFF(I,NT)
            IF( NS .LE. NBPOIN ) THEN
               DO 10 J=1,3
                  XYZ(J,I) = XYZPOI(J,NS)
 10            ENDDO
C              LA COULEUR
               IF( NTYP .EQ. 0 ) THEN
                  T = REAL( TEMPER(NS,NCAS) )
               ELSE
                  T = REAL( dptemp( NCAS )%dptab( NS ) )
               ENDIF
               IF( TMIN .NE. TMAX ) THEN
                  COULEUR = REAL( (T-TMIN) / (TMAX-TMIN) )
               ELSE
                  COULEUR = 1.0
               ENDIF
C              PROJECTION SI DEPASSEMENT AU DELA DU MIN OU MAX
               IF( COULEUR .GT. 1.0 ) COULEUR = 1.0
               IF( COULEUR .LT. 0.0 ) COULEUR = 0.0
               COUL(I) = N1COUL + NBCOUL * COULEUR
            ELSE
               NS = NS - NBPOIN
               DO 20 J=1,3
                  XYZ(J,I) = XYZVQ2(J,NS)
 20            ENDDO
C              LA COULEUR
               IF( TMIN .NE. TMAX ) THEN
                  COULEUR = (XYZVQ2(4,NS)-TMIN) / (TMAX-TMIN)
               ELSE
                  COULEUR = 1.0
               ENDIF
C              PROJECTION SI DEPASSEMENT AU DELA DU MIN OU MAX
               IF( COULEUR .GT. 1.0 ) COULEUR = 1.0
               IF( COULEUR .LT. 0.0 ) COULEUR = 0.0
               COUL(I) = N1COUL + NBCOUL * COULEUR
            ENDIF
 30      ENDDO
C
         IF( PREDUF .GT. 0 ) THEN

C           REDUCTION DE LA FACE
            REDUCF = PREDUF * 0.01
            REDUC1 = 1.0 - REDUCF

C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
            CALL COBAPO( 3, XYZ, XYZP )

C           L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
            DO 40 I=1,3
               XYZ(1,I) = XYZ(1,I) * REDUC1 + XYZP(1) * REDUCF
               XYZ(2,I) = XYZ(2,I) * REDUC1 + XYZP(2) * REDUCF
               XYZ(3,I) = XYZ(3,I) * REDUC1 + XYZP(3) * REDUCF
 40         ENDDO

         ENDIF

C        LE TRACE DU TRIANGLE
         CALL TRIACOUL3DBORD( XYZ, COUL, NCOUAF, 1 )

 100  ENDDO

      RETURN
      END
