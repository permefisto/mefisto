      SUBROUTINE TRARFAFR( MOFACE, MXFACE, LFACES, XYZPOI )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DES FACES FRONTALIERES DEFINIES
C -----    DANS LE TABLEAU LFACES SELON LA COULEUR NCOAFR

C ENTREES:
C --------
C NCOAFR : NUMERO DE LA COULEUR DE TRACE DES ARETES DES FACES

C ATTENTION: ICI LFACES N'EST PAS LE TABLEAU DU CAS GENERAL
C            LE TABLEAU a_objet_face a du etre modifie auparavant
C            cf les sp flui/parparv2.f parpartr.f
C MOFACE : NOMBRE DE MOTS PAR FACE CHAINEE DU TABLEAU LFACES
C MXFACE : NOMBRE DE LIGNES DU TABLEAU LFACES
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
C          LFACES(7,I)= 0   <== GENERALEMENT C'EST LE CHAINAGE NON NUL
C                               DES FACES FRONTALIERES!
C                               cf les sp flui/parparv2.f parpartr.f

C XYZPOI : 3 COORDONNEES DES SOMMETS DES FACES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS St PIERRE du PERRAY JUILLET2009
C23456---------------------------------------------------------------012
      include"./incl/mecoit.inc"
      INTEGER  LFACES(1:MOFACE,1:MXFACE)
      REAL     XYZPOI(1:3,1:*)
      REAL     XYZ(3,4), XYZP(3)

      IF( IAVFAC .LE. 0 .OR. NCOAFR .LT. 0 ) GOTO 9999

C     BOUCLE SUR LES FACES DU TABLEAU LFACES
C     ======================================
      CALL XVEPAISSEUR( NEPARF )
      DO 100 NF = 1, MXFACE

C        ELIMINATION DES FACES NON UTILISEES DANS MNFACE
         IF( LFACES(1,NF) .EQ. 0 ) GOTO 100

C        ELIMINATION DES FACES APPARTENANT A 2 CUBES
         IF( LFACES(7,NF) .NE. 0 ) GOTO 100

C        LA FACE NF EST FRONTALIERE
         IF( PREDUF .LE. 0.0 ) THEN

C           PAS DE REDUCTION DEMANDEE DES FACES
C           -----------------------------------
C           TRACE DE LA FACE FRONTALIERE NF
            IF( LFACES(4,NF) .EQ. 0 ) THEN
C              TRIANGLE
               NBSF = 3
            ELSE
C              QUADRANGLE
               NBSF = 4
            ENDIF

            J0 = NBSF
            DO J=1,NBSF
               NS1 = LFACES(J0,NF)
               NS2 = LFACES(J ,NF)
               CALL TRAIT3D( NCOAFR, XYZPOI(1,NS1), XYZPOI(1,NS2) )
               J0 = J
            ENDDO

         ELSE

C           REDUCTION DEMANDEE DES FACES FRONTALIERES
C           -----------------------------------------
            REDUCF = PREDUF * 0.01
            REDUC1 = 1.0 - REDUCF

            IF( LFACES(4,NF) .EQ. 0 ) THEN
C              TRIANGLE
               NBSF = 3
            ELSE
C              QUADRANGLE
               NBSF = 4
            ENDIF

            DO J=1,NBSF
C              LES COORDONNEES DU SOMMET J
               NS1 = LFACES(J,NF)
               XYZ(1,J) = XYZPOI(1,NS1)
               XYZ(2,J) = XYZPOI(2,NS1)
               XYZ(3,J) = XYZPOI(3,NS1)
            ENDDO

C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
            CALL COBAPO( NBSF, XYZ, XYZP )

C           L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
            DO J=1,NBSF
               XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
               XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
               XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
            ENDDO

C           TRACE DE LA FACE NF REDUITE
            J0 = NBSF
            DO J=1,NBSF
               CALL TRAIT3D( NCOAFR, XYZ(1,J0), XYZ(1,J) )
               J0 = J
            ENDDO

         ENDIF

 100  ENDDO

 9999 RETURN
      END
