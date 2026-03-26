      SUBROUTINE TRABOULE( RAYON, XYZCENTR, NCBOULE, NCARETE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     TRACER AVEC LA COULEUR NCBOULE LA TRIANGULATION D'UNE BOULE
C -----     DE RAYON ET DE CENTRE DE COORDONNEES XYZCENTR FACON COVID

C ENTREES:
C --------
C RAYON   : RAYON DE LA BOULE
C XYZCENTR: XYZ DU POINT CENTRE DE LA BOULE
C NCBOULE : NUMERO DE LA COULEUR DE TRACE DE LA SURFACE    DE LA BOULE
C NCARETE : NUMERO DE LA COULEUR DES ARETES DE L'ICOSAEDRE DE LA BOULE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY         Novembre 2020
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"

C     DECLARATION POUR CONSTRUIRE LA TRIANGULATION D'UNE SPHERE
C     NBDEAR NOMBRE DE SUBDIVISION D'UNE ARETE DE L'ICOSAEDRE
      PARAMETER        (NBDEAR=2)
      INTEGER           T (0:(NBDEAR+1)*(2*NBDEAR)),
     %                  T0(0:(NBDEAR+1)*(2*NBDEAR)),
     %                  TBG(0:3*NBDEAR+2),
     %                  TBD(0:3*NBDEAR+2),
     %                  NOSOSP(4,20*NBDEAR**2)
      REAL              XYZSSP(3,20*NBDEAR**2), XYZS(3,3),
     %                  XYZ1(3), XYZ2(3)

      REAL              RAYON, XYZCENTR(3)
      INTEGER           NCBOULE, NCARETE

C     GENERER LA TRIANGULATION REGULIERE DE LA SURFACE DE LA BOULE
      CALL SPH1S1( NBDEAR, RAYON, XYZCENTR(1),XYZCENTR(2),XYZCENTR(3),
     %             T,      T0,     TBG,    TBD,
     %             NBSTSP, NBTRSP, XYZSSP, NOSOSP )

C     BOUCLE SUR LES NBTRSP DE LA BOULE
      CALL XVEPAISSEUR( 0 )
      DO NT = 1, NBTRSP

C        LES 3 COORDONNEES DES 3 SOMMETS DU TRIANGLE NT
         DO J=1,3
            DO K=1,3
               XYZS( K, J ) = XYZSSP( K, NOSOSP(J,NT) )
            ENDDO
         ENDDO

C        LE TRACE EFFECTIF DU TRIANGLE NT SANS LES ARETES
         CALL FACE3D( NCBOULE, -1, 3, XYZS )

C        TRACE DES 3 ARETES SAUF LA MOITIE MEDIANE  FACON COVID
         J0 = 3
         DO J=1,3
C           ARETE J0-J
            DO K=1,3
               R = ( XYZS(K,J) - XYZS(K,J0) )/4
               XYZ1(K) = XYZS(K,J0) + R
               XYZ2(K) = XYZS(K,J)  - R
            ENDDO
            CALL TRAIT3D( NCARETE, XYZS(1,J0), XYZ1 )
            CALL TRAIT3D( NCARETE, XYZS(1,J),  XYZ2 )
            J0 = J
         ENDDO

      ENDDO

      RETURN
      END

