      SUBROUTINE VOEXH2( NBX,NBY,NBZ,NBXYZ,
     %                   XYZ1,XYZ2,XYZ3,XYZ4,XYZ5,XYZ6,
     %                   HEXYZ,CUXYZ,
     %                   NSENS,NPSOM,COSO)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECOPIER LES COORDONNEES DES SOMMETS DE CHAQUE FACE
C -----    DE L'HEXAEDRE A PARTIR DES DONNEES
C
C ENTREES:
C --------
C NBX,NBY,NBZ : NOMBRE DE SOMMETS DANS CHAQUE DIRECTION
C NBXYZ       : NOMBRE DE SOMMETS DANS CHAQUE FACE (MAJORATION)
C XYZ1 A XYZ6 : LES COORDONNEES DES SOMMETS DES 6 FACES
C NSENS(I)    : SENS DE PARCOURS DE LA FACE I PARMI LES DONNEES
C NPSOM(I)    : NUMERO DU PREMIER SOMMET LA FACE I
C
C SORTIES:
C --------
C HEXYZ       : LES COORDONNEES DES SOMMETS DES FACES DE L'HEXAEDRE
C CUXYZ       : LES COORDONNEES DES SOMMETS DES FACES DU CUBE UNITE
C COSO        : LES COORDONNEES DES POINTS DU BORD DE L'HEXAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS       DECEMBRE  1988
C.......................................................................
      INTEGER  NPSOM(6),NSENS(6)
      REAL     XYZ1(3,NBX*NBY),XYZ2(3,NBX*NBZ),XYZ3(3,NBY*NBZ),
     %         XYZ4(3,NBX*NBY),XYZ5(3,NBX*NBZ),XYZ6(3,NBY*NBZ),
     %         HEXYZ(3,6,NBXYZ),CUXYZ(3,6,NBXYZ),
     %         COSO(3,NBX*NBY*NBZ)
C
C     INITIALISATION
C     ==============
      NBS = NBX*NBY*NBZ
      DO NS=1,NBS
         DO J=1,3
            COSO(J,NS)=0.
         ENDDO
      ENDDO
      NBXNBY = NBX * NBY
C
C     CONSTRUCTION DE L'HEXAEDRE
C     ==========================
C
C     LA NUMEROTATION DANS L'HEXAEDRE EST NATURELLE :
C     NS = NX + (NY-1) * NBX + (NZ-1) * NBX * NBY
C
C     1) LA FACE 1
C
      NF = 1
      CALL VOEXH4(NBX,NBY,NPSOM(NF),NSENS(NF),N1,N2,N3)
      DO NX=1,NBX
         DO NY=1,NBY
            NSH = NX + (NY-1) * NBX
            NSF = N1 + N2 * NX + N3 * NY
            NSC = NSH
            DO J=1,3
               COSO(J,NSH)     = XYZ1(J,NSF)
               HEXYZ(J,NF,NSC) = XYZ1(J,NSF)
            ENDDO
         ENDDO
      ENDDO
      CALL VOEXH5(NF,HEXYZ,CUXYZ,NBX,NBY,NBXYZ)
C
C     2) LA  FACE 2
C
      NF = 2
      CALL VOEXH4(NBX,NBZ,NPSOM(NF),NSENS(NF),N1,N2,N3)
      DO NX=1,NBX
         DO NZ=1,NBZ
            NSH = NX + (NZ-1) * NBXNBY
            NSF = N1 + N2 * NX + N3 * NZ
            NSC = NX + (NZ-1) * NBX
            DO J=1,3
               COSO(J,NSH)     = XYZ2(J,NSF)
               HEXYZ(J,NF,NSC) = XYZ2(J,NSF)
            ENDDO
         ENDDO
      ENDDO
      CALL VOEXH5(NF,HEXYZ,CUXYZ,NBX,NBZ,NBXYZ)
C
C     3) LA  FACE 3
C
      NF = 3
      CALL VOEXH4(NBY,NBZ,NPSOM(NF),NSENS(NF),N1,N2,N3)
      DO NY=1,NBY
         DO NZ=1,NBZ
            NSH = 1  + (NY-1) * NBX + (NZ-1) * NBXNBY
            NSF = N1 + N2 * NY + N3 * NZ
            NSC = NY + (NZ-1) * NBY
            DO J=1,3
               COSO(J,NSH)     = XYZ3(J,NSF)
               HEXYZ(J,NF,NSC) = XYZ3(J,NSF)
            ENDDO
         ENDDO
      ENDDO
      CALL VOEXH5(NF,HEXYZ,CUXYZ,NBY,NBZ,NBXYZ)
C
C     4) LA  FACE 4
C
      NF = 4
      CALL VOEXH4(NBX,NBY,NPSOM(NF),NSENS(NF),N1,N2,N3)
      NDECAL = (NBZ-1) * NBXNBY
      DO NX=1,NBX
         DO NY=1,NBY
C           NSH = NX + (NY-1) * NBX + (NBZ-1) * NBXNBY
            NSH = NX + (NY-1) * NBX + NDECAL
            NSF = N1 + N2 * NX + N3 * NY
            NSC = NX + (NY-1) * NBX
            DO J=1,3
               COSO(J,NSH)     = XYZ4(J,NSF)
               HEXYZ(J,NF,NSC) = XYZ4(J,NSF)
            ENDDO
         ENDDO
      ENDDO
      CALL VOEXH5(NF,HEXYZ,CUXYZ,NBX,NBY,NBXYZ)
C
C     5) LA  FACE 5
C
      NF = 5
      CALL VOEXH4(NBX,NBZ,NPSOM(NF),NSENS(NF),N1,N2,N3)
      NDECAL = (NBY-1) * NBX
      DO NX=1,NBX
         DO NZ=1,NBZ
C           NSH = NX + (NBY-1) * NBX + (NZ-1) * NBXNBY
            NSH = NX + NDECAL + (NZ-1) * NBXNBY
            NSF = N1 + N2 * NX + N3 * NZ
            NSC = NX + (NZ-1) * NBX
            DO J=1,3
               COSO(J,NSH)     = XYZ5(J,NSF)
               HEXYZ(J,NF,NSC) = XYZ5(J,NSF)
            ENDDO
         ENDDO
      ENDDO
      CALL VOEXH5(NF,HEXYZ,CUXYZ,NBX,NBZ,NBXYZ)
C
C     6) LA  FACE 6
C
      NF = 6
      CALL VOEXH4(NBY,NBZ,NPSOM(NF),NSENS(NF),N1,N2,N3)
      DO NY=1,NBY
         DO NZ=1,NBZ
            NSH =  NBX + (NY-1) * NBX + (NZ-1) * NBXNBY
            NSF = N1 + N2 * NY + N3 * NZ
            NSC = NY + (NZ-1) * NBY
            DO J=1,3
               COSO(J,NSH)     = XYZ6(J,NSF)
               HEXYZ(J,NF,NSC) = XYZ6(J,NSF)
            ENDDO
         ENDDO
      ENDDO
      CALL VOEXH5(NF,HEXYZ,CUXYZ,NBY,NBZ,NBXYZ)
C
      RETURN
      END
