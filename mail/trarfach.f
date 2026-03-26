      SUBROUTINE TRARFACH( NOCOUL, MOFACE, L1FACH, LFACES, XYZPOI )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DES FACES CHAINEES DEFINIES DANS LE TABLEAU
C -----    LFACES SELON LA COULEUR NOCOUL
C
C ENTREES:
C --------
C NOCOUL : NUMERO DE LA COULEUR DE TRACE DES FACES
C MOFACE : NOMBRE DE MOTS PAR FACE CHAINEE DU TABLEAU LFACES
C L1FACH : NUMERO DE LA PREMIERE FACE CHAINEE DANS LE TABLEAU LFACES
C LFACES : TABLEAU DES FACES DU MAILLAGE
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE OU 0 SI TRIANGLE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C          SI SOMMET 2 < DERNIER SOMMET  => FACE   DIRECTE DANS LE CUBE
C                      <                 => FACE INDIRECTE
C          UNE FACE DIRECTE EST VUE DE L EXTERIEUR AU CUBE
C          SOUS LA FORME DIRECTE
C          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CE CUBE
C                       <0 SI FACE INDIRECTE DANS CE CUBE
C          SI LA FACE APPARTIENT A 2 CUBES ALORS
C          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CE CUBE
C                       <0 SI FACE INDIRECTE DANS CE CUBE
C          SINON
C          LFACES(7,I)= NUMERO DE LA FACE FRONTALIERE ou INTERFACE SUIVANTE
C                       0 SI C'EST LA DERNIERE
C
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS TABLEAU NUTGFA
C XYZPOI : 3 COORDONNEES DES SOMMETS DES FACES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS St PIERRE du PERRAY JUILLET2009
C23456---------------------------------------------------------------012
      include"./incl/mecoit.inc"
      INTEGER  LFACES(1:MOFACE,1:*)
      REAL     XYZPOI(1:3,1:*)
      REAL     XYZ(3,4), XYZP(3)
C
C     BOUCLE SUR LES FACES CHAINEES
      N = L1FACH
C
      IF( PREDUF .LE. 0.0 ) THEN
C
C        PAS DE REDUCTION DEMANDEE DES FACES
C        ===================================
 10      IF( N .GT. 0 ) THEN
C
C           TRACE DE LA FACE CHAINEE N
            IF( LFACES(4,N) .EQ. 0 ) THEN
C              TRIANGLE
               NBSF = 3
            ELSE
C              QUADRANGLE
               NBSF = 4
            ENDIF
C
            J0 = NBSF
            DO 20 J=1,NBSF
               NS1 = LFACES(J0,N)
               NS2 = LFACES(J ,N)
               CALL TRAIT3D( NOCOUL, XYZPOI(1,NS1), XYZPOI(1,NS2) )
               J0 = J
 20         CONTINUE
C
C           PASSAGE A LA FACE SUIVANTE
            N = LFACES(7,N)
            GOTO 10
C
         ENDIF
C
      ELSE
C
C        REDUCTION DEMANDEE DES FACES
C        ============================
         REDUCF = PREDUF * 0.01
         REDUC1 = 1.0 - REDUCF
C
 100     IF( N .GT. 0 ) THEN
C
            IF( LFACES(4,N) .EQ. 0 ) THEN
C              TRIANGLE
               NBSF = 3
            ELSE
C              QUADRANGLE
               NBSF = 4
            ENDIF
C
            DO 110 J=1,NBSF
C              LES COORDONNEES DU SOMMET J
               NS1 = LFACES(J,N)
               XYZ(1,J) = XYZPOI(1,NS1)
               XYZ(2,J) = XYZPOI(2,NS1)
               XYZ(3,J) = XYZPOI(3,NS1)
 110        CONTINUE
C
C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
            CALL COBAPO( NBSF, XYZ, XYZP )
C
C           L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
            DO 120 J=1,NBSF
               XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
               XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
               XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
 120        CONTINUE
C
C           TRACE DE LA FACE CHAINEE N REDUITE
            J0 = NBSF
            DO 130 J=1,NBSF
               CALL TRAIT3D( NOCOUL, XYZ(1,J0), XYZ(1,J) )
               J0 = J
 130        CONTINUE
C
C           PASSAGE A LA FACE SUIVANTE
            N = LFACES(7,N)
            GOTO 100
C
         ENDIF
      ENDIF
C
      RETURN
      END
