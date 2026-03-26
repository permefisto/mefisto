      SUBROUTINE TRTETRA( COULARETE, NOSOTE, PTXYZD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      TRACE DES 6 ARETES D'UN TETRAEDRE DEFINI PAR LE NUMERO
C -----      DE SES 4 SOMMETS DANS PTXYZD
C
C ENTREES:
C --------
C COULARETE: NO COULEUR DE TRACE DES ARETES DU TRIANGLE
C NOSOTE   : NO DES 4 SOMMETS DANS LE TABLEAU PTXYZD
C PTXYZD   : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY     AVRIL 2008
C MODIFS: ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY SEPTEMBRE 2016
C2345X7..............................................................012
      include"./incl/mecoit.inc"

      INTEGER           COULARETE, NOSOTE(4)
      DOUBLE PRECISION  PTXYZD(4,*), XYZBAR(3)
      REAL              XYZ1(3), XYZ2(3)

C     PREDUF : pourcentage de reduction des faces
      IF( PREDUF .LE. 0.0 ) THEN

C        PAS DE REDUCTION DES FACES DU TETRAEDRE
C        ---------------------------------------
C        LES 3 ARETES DE LA FACE 1 DU TETRAEDRE
         NS1 = NOSOTE(3)
         DO I = 1, 3
            NS2 = NOSOTE(I)
            DO K=1,3
               XYZ1(K) = REAL( PTXYZD(K,NS1) )
               XYZ2(K) = REAL( PTXYZD(K,NS2) )
            ENDDO
            CALL TRAIT3D( COULARETE, XYZ1, XYZ2 )
            NS1 = NS2
         ENDDO

C        LES 3 ARETES NS4-NS1,NS2,NS3 DU TETRAEDRE
         NS1 = NOSOTE(4)
         DO K=1,3
            XYZ1(K) = REAL( PTXYZD(K,NS1) )
         ENDDO
         DO I = 1, 3
            NS2 = NOSOTE(I)
            DO K=1,3
               XYZ2(K) = REAL( PTXYZD(K,NS2) )
            ENDDO
            CALL TRAIT3D( COULARETE, XYZ1, XYZ2 )
         ENDDO

      ELSE

C        REDUCTION DES FACES DU TETRAEDRE
C        --------------------------------
C        CALCUL DU BARYCENTRE DU TETRAEDRE
         XYZBAR(1) = 0D0
         XYZBAR(2) = 0D0
         XYZBAR(3) = 0D0
         DO I=1,4
            NS1 = NOSOTE(I)
            DO K=1,3
               XYZBAR(K) = XYZBAR(K) +  PTXYZD(K,NS1)
            ENDDO
         ENDDO
         DO K=1,3
            XYZBAR(K) = XYZBAR(K) / 4
         ENDDO

C        LES 3 ARETES DE LA FACE 1 DU TETRAEDRE REDUIT DE PREDUF
         REDUC0 = PREDUF * 0.01
         REDUC1 = 1.0 - REDUC0

         NS1 = NOSOTE(3)
         DO I = 1, 3
            NS2 = NOSOTE(I)
            DO K=1,3
               XYZ1(K) = REAL( REDUC1*PTXYZD(K,NS1) + REDUC0*XYZBAR(K) )
               XYZ2(K) = REAL( REDUC1*PTXYZD(K,NS2) + REDUC0*XYZBAR(K) )
            ENDDO
            CALL TRAIT3D( COULARETE, XYZ1, XYZ2 )
            NS1 = NS2
         ENDDO

C        LES 3 ARETES NS4-NS1,NS2,NS3 DU TETRAEDRE
         NS1 = NOSOTE(4)
         DO K=1,3
            XYZ1(K) = REAL( REDUC1*PTXYZD(K,NS1) + REDUC0*XYZBAR(K) )
         ENDDO
         DO I = 1, 3
            NS2 = NOSOTE(I)
            DO K=1,3
               XYZ2(K) = REAL( REDUC1*PTXYZD(K,NS2) + REDUC0*XYZBAR(K) )
            ENDDO
            CALL TRAIT3D( COULARETE, XYZ1, XYZ2 )
         ENDDO

      ENDIF

      RETURN
      END
