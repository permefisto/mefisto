      SUBROUTINE XYZBF10( PTXYZD, N1FEOC, NPFEOC, NFETOI,
     %                    NBFETO, XYZMIX, XYZBAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :   CALCULER LES COORDONNEES XYZBAR DU BARYCENTRE DES FACES
C  -----   N1FEOC A NPFEOC DE L'ETOILE DEFINIE DANS NFETOI

C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NPFEOC : NUMERO NFETOI DE LA PREMIERE FACE DE L'ETOILE A NE PAS TRAITER
C NFETOI : LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST DIRIGE VERS L'INTERIEUR DE L'ETOILE
C          5: CHAINAGE SUR LA FACE SUIVANTE OCCUPEE OU VIDE

C SORTIES:
C --------
C NBFETO : NOMBRE DE FACES DE L'ETOILE RECALCULE
C XYZMIX : COORDOONNEES EXTREMES DES SOMMETS DES FACES DE L'ETOILE
C XYZBAR : 3 COORDONNEES XYZ DU BARYCENTRE DES FACES N1FEOC A NPFEOC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER   NOVEMBRE 2014
C2345X7..............................................................012
      DOUBLE PRECISION  PTXYZD(4,*), XYZMIX(3,2), XYZBAR(3), S
      INTEGER           NFETOI(5,*)

C     COORDONNEES DES POINTS EXTREMES ET DU BARYCENTRE DES SOMMETS
C     DES FACES DE L'ETOILE
      DO K=1,3
         XYZBAR(K) = 0D0
         XYZMIX(K,1) = 1D100
         XYZMIX(K,2) =-1D100
      ENDDO

      NBFETO = 0
      NBS    = 0
      NF1    = N1FEOC

 20   IF( NF1 .NE. NPFEOC ) THEN

C        UNE FACE DE PLUS
         NBFETO = NBFETO + 1

C        BOUCLE SUR LES 3 SOMMETS DE LA FACE NF1
         DO K=1,3

C           NUMERO DU SOMMET K DE LA FACE NF1
            NS = ABS( NFETOI(1+K,NF1) )

            DO L=1,3

C              COORDONNEE L DU SOMMET NS
               S = PTXYZD(L,NS)
               IF( S .LT. XYZMIX(L,1) )  XYZMIX(L,1) = S
               IF( S .GT. XYZMIX(L,2) )  XYZMIX(L,2) = S

C              BARYCENTRE
               XYZBAR(L) = XYZBAR(L) + S

            ENDDO

         ENDDO
         NBS = NBS + 3

C        LA FACE SUIVANTE
         NF1 = NFETOI(5,NF1)
         GOTO 20

      ENDIF

      DO L=1,3
         XYZBAR(L) = XYZBAR(L) / NBS
      ENDDO

ccc      print *,'xyzbf10: Xmin=',XYZMIX(1,1),' X Barycentre=',XYZBAR(1),
ccc     %                ' XMax=',XYZMIX(1,2)
ccc      print *,'xyzbf10: Ymin=',XYZMIX(2,1),' Y Barycentre=',XYZBAR(2),
ccc     %                ' YMax=',XYZMIX(2,2)
ccc      print *,'xyzbf10: Zmin=',XYZMIX(3,1),' Z Barycentre=',XYZBAR(3),
ccc     %                ' ZMax=',XYZMIX(3,2)

      RETURN
      END
