      SUBROUTINE EF4SEF( NBSOSI, XYZSTI, NBEFSI, NSEFSI,
     %                   MXSOSF, NBSOSF, XYZSTF, NBEFSF, NSEFSF,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 4 SOUS EF DES EF D'UNE SURFACE
C -----
C
C ENTREES:
C --------
C NBSOSI : NOMBRE DE SOMMETS DE LA SURFACE INITIALE
C XYZSTI : 3 COORDONNEES DES NBSOSI SOMMETS DE LA SURFACE INITIALE
C NBEFSI : NOMBRE D'EF DE LA SURFACE INITIALE
C NSEFSI : NO DES 4 SOMMETS DES EF DE LA SURFACE INITIALE
C MXSOSF : NOMBRE MAXIMAL DE SOMMETS DECLARABLES POUR LE MAILLAGE FINAL
C
C SORTIES:
C --------
C NBSOSF : NOMBRE DE SOMMETS DE LA SURFACE FINALE
C XYZSTF : 3 COORDONNEES DES NBSOSF SOMMETS DE LA SURFACE FINALE
C NBEFSF : NOMBRE D'EF DE LA SURFACE FINALE
C NSEFSF : NO DES 4 SOMMETS DES EF DE LA SURFACE FINALE
C IERR   : 0 SI PAS D'ERREUR, >0 SI MXSOSF TROP FAIBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC  et VEULETTES sur MER  Juillet 2008
C.......................................................................
      include"./incl/sotgar.inc"
      include"./incl/sotgfc.inc"
      REAL     XYZSTI(3,NBSOSI), XYZSTF(3,MXSOSF)
      INTEGER  NSEFSI(4,NBEFSI), NSEFSF(4,NBEFSF)

      INTEGER  NUSEF(9)

      INTEGER  NUSTRIA(3,4)
      INTEGER  NUSQUAD(4,4)
      DATA     NUSTRIA / 1,4,6,    4,2,5,   4,5,6,   6,5,3 /
      DATA     NUSQUAD / 1,5,9,8,  5,2,6,9, 8,9,7,4, 9,6,3,7 /

      IERR = 0

C     COPIE DES XYZ DES SOMMETS INITIAUX EN DEBUT DES FINAUX
      DO 20 J=1,NBSOSI
         DO 10 I=1,3
            XYZSTF(I,J) = XYZSTI(I,J)
 10      CONTINUE
 20   CONTINUE
      NBSOSF = NBSOSI

C     SUBDIVISION DES TRIANGLES et QUADRANGLES
      NBEFSF = 0
      DO 100 N=1,NBEFSI

         IF( NBSOSF .GE. MXSOSF ) THEN
            IERR = 1
            RETURN
         ENDIF

         IF( NSEFSI(4,N) .EQ. 0 ) THEN

C           TRIANGLE
C           ========
C           CONSTRUCTION DES NUMEROS DES 3 SOMMETS INITIAUX
            DO 22 I=1,3
               NUSEF(I) = NSEFSI(I,N)
 22         CONTINUE

C           CONSTRUCTION DES 3 MILIEUX DES ARETES ET DE LEUR NUMERO
            DO 24 J=1,3
               NBSOSF = NBSOSF + 1
               NUSEF(3+J) = NBSOSF
               NSA1 = NUSEF( J )
               IF( J .LT. 3 ) THEN
                  J1 = J+1
               ELSE
                  J1 = 1
               ENDIF
               NSA2 = NUSEF( J1 )
               XYZSTF(1,NBSOSF) = (XYZSTF(1,NSA1) + XYZSTF(1,NSA2))* 0.5
               XYZSTF(2,NBSOSF) = (XYZSTF(2,NSA1) + XYZSTF(2,NSA2))* 0.5
               XYZSTF(3,NBSOSF) = (XYZSTF(3,NSA1) + XYZSTF(3,NSA2))* 0.5
 24         CONTINUE

C           DECOUPAGE DU TRIANGLE EN 4 SOUS-TRIANGLES
            DO 29 J=1,4
               NBEFSF = NBEFSF + 1
               DO 25 I=1,3
C                 LE NUMERO DES 3 SOMMETS DU SOUS EF J
                  NSEFSF( I, NBEFSF ) = NUSEF( NUSTRIA(I,J) )
 25            CONTINUE
               NSEFSF( 4, NBEFSF ) = 0
 29         CONTINUE

         ELSE

C           QUADRANGLE
C           ==========
C           CONSTRUCTION DES NUMEROS DES 4 SOMMETS INITIAUX
            DO 62 I=1,4
               NUSEF(I) = NSEFSI(I,N)
 62         CONTINUE

C           CONSTRUCTION DES 4 MILIEUX DES ARETES ET DE LEUR NUMERO
            DO 64 J=1,4
               NBSOSF = NBSOSF + 1
               NUSEF(4+J) = NBSOSF
               NSA1 = NUSEF( J )
               IF( J .LT. 4 ) THEN
                  J1 = J+1
               ELSE
                  J1 = 1
               ENDIF
               NSA2 = NUSEF( J1 )
               XYZSTF(1,NBSOSF) = (XYZSTF(1,NSA1) + XYZSTF(1,NSA2))* 0.5
               XYZSTF(2,NBSOSF) = (XYZSTF(2,NSA1) + XYZSTF(2,NSA2))* 0.5
               XYZSTF(3,NBSOSF) = (XYZSTF(3,NSA1) + XYZSTF(3,NSA2))* 0.5
 64         CONTINUE

C           CONSTRUCTION DU BARYCENTRE DU QUADRANGLE
            NBSOSF = NBSOSF + 1
            NUSEF(9) = NBSOSF
            XYZSTF(1,NBSOSF) = ( XYZSTF(1,NUSEF(1)) + XYZSTF(1,NUSEF(2))
     %                         + XYZSTF(1,NUSEF(3)) + XYZSTF(1,NUSEF(4))
     %                         ) * 0.25
            XYZSTF(2,NBSOSF) = ( XYZSTF(2,NUSEF(1)) + XYZSTF(2,NUSEF(2))
     %                         + XYZSTF(2,NUSEF(3)) + XYZSTF(2,NUSEF(4))
     %                         ) * 0.25
            XYZSTF(3,NBSOSF) = ( XYZSTF(3,NUSEF(1)) + XYZSTF(3,NUSEF(2))
     %                         + XYZSTF(3,NUSEF(3)) + XYZSTF(3,NUSEF(4))
     %                         ) * 0.25

C           DECOUPAGE DU QUADRANGLE EN 4 SOUS-QUADRANGLES
            DO 69 J=1,4
               NBEFSF = NBEFSF + 1
               DO 68 I=1,4
C                 LE NUMERO DES 4 SOMMETS DU SOUS EF J
                  NSEFSF(I,NBEFSF) = NUSEF( NUSQUAD(I,J) )
 68            CONTINUE
 69         CONTINUE

         ENDIF

 100  CONTINUE

      RETURN
      END
