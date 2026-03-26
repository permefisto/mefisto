      SUBROUTINE TGSPHE( RAYON,  XC, YC, ZC, NBTRIA, NOSOTR,
     %                   XYZSOM, XYZTGS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 6*NBTRIA TANGENTES
C -----    (ON POURRAIT DIVISER PAR 2 LE STOCKAGE)
C
C ENTREES:
C -------
C RAYON  : RAYON DE LA SPHERE
C XC     : COORDONNEE X DU CENTRE DE LA SPHERE
C YC     : COORDONNEE Y DU CENTRE DE LA SPHERE
C ZC     : COORDONNEE Z DU CENTRE DE LA SPHERE
C NBTRIA : NOMBRE DE TRIANGLES
C NOSOTR : NUMERO DES 3 SOMMETS PUIS 0 DES NBTRIA TRIANGLES
C XYZSOM : TABLEAU DES COORDONNES DES SOMMETS DES TRIANGLES
C
C SORTIES:
C -------
C XYZTGS : LES 3 COORDONNEES DES TANGENTES DES EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET  ANALYSE NUMERIQUE UPMC  PARIS         JUIN 1996
C234567--------------------------------------------------------------012
      INTEGER  NOSOTR(4,NBTRIA)
      REAL     XYZSOM(3,*), XYZTGS(3,*)
      REAL     A(3,3), ANGLE(3)
C
C     BOUCLE SUR LES TRIANGLES
      R2   = RAYON ** 2
      NUTG = 0
      DO 100 N=1,NBTRIA
C
C        CALCUL DES 3 ANGLES AU CENTRE (C SI, C SI+1)
         DO 10 I=1,3
C           NUMERO DU SOMMET I DU TRIANGLE N
            NS = NOSOTR(I,N)
            IF( I .LT. 3 ) THEN
               NS1 = I+1
            ELSE
               NS1 = 1
            ENDIF
C           NUMERO DU SOMMET I+1 DU TRIANGLE N
            NS1 = NOSOTR(NS1,N)
C
C           ANGLE AU CENTRE (C SI, C SI+1)
            A(1,1) = XYZSOM(1,NS) - XC
            A(2,1) = XYZSOM(2,NS) - YC
            A(3,1) = XYZSOM(3,NS) - ZC
C
            A(1,2) = XYZSOM(1,NS1) - XC
            A(2,2) = XYZSOM(2,NS1) - YC
            A(3,2) = XYZSOM(3,NS1) - ZC
C
C           LE VECTEUR NORMAL
            CALL PROVER( A(1,1), A(1,2), A(1,3) )
C           SA NORME
            S = SQRT( A(1,3)**2 + A(2,3)**2 + A(3,3)**2 )
C           L'ANGLE POSITIF AU CENTRE (C SI, C SI+1)
            ANGLE(I) = ASIN( S / R2 )
 10      CONTINUE
C
C        LES 2 VECTEURS TANGENTS AUX 3 SOMMETS DU TRIANGLE N
         I0  = 3
         NS0 = NOSOTR(3,N)
         DO 20 I=1,3
C           NUMERO DU SOMMET I DU TRIANGLE N
            NS = NOSOTR(I,N)
            IF( I .LT. 3 ) THEN
               NS1 = I+1
            ELSE
               NS1 = 1
            ENDIF
C           NUMERO DU SOMMET I+1 DU TRIANGLE N
            NS1 = NOSOTR(NS1,N)
C           NUMERO DES 2 TANGENTES DU SOMMET SI
            NTG1 = NUTG + 1
            NTG2 = NUTG + 2
            NUTG = NTG2
C
C           LE VECTEUR RAYON C SI
            A(1,1) = XYZSOM(1,NS) - XC
            A(2,1) = XYZSOM(2,NS) - YC
            A(3,1) = XYZSOM(3,NS) - ZC
C
C           LE VECTEUR ORTHOGONAL AU RAYON ET SI SI+1
            A(1,2) = XYZSOM(1,NS1) - XYZSOM(1,NS)
            A(2,2) = XYZSOM(2,NS1) - XYZSOM(2,NS)
            A(3,2) = XYZSOM(3,NS1) - XYZSOM(3,NS)
            CALL PROVER( A(1,1), A(1,2), A(1,3) )
C           LE VECTEUR ORTHOGONAL AU RAYON DANS LE PLAN SI SI+1
            CALL PROVER( A(1,3), A(1,1), XYZTGS(1,NTG1) )
C           NORMALISATION A 1 DU VECTEUR TANGENT SI SI+1 EN SI
            CALL NORMER( 3, XYZTGS(1,NTG1), J )
C           NORME MISE A  RAYON * ANGLE(CSI, CSI+1)
            S = RAYON * ANGLE(I)
            XYZTGS(1,NTG1) = XYZTGS(1,NTG1) * S
            XYZTGS(2,NTG1) = XYZTGS(2,NTG1) * S
            XYZTGS(3,NTG1) = XYZTGS(3,NTG1) * S
C
C           LE VECTEUR ORTHOGONAL AU RAYON ET SI SI+1
            A(1,2) = XYZSOM(1,NS0) - XYZSOM(1,NS)
            A(2,2) = XYZSOM(2,NS0) - XYZSOM(2,NS)
            A(3,2) = XYZSOM(3,NS0) - XYZSOM(3,NS)
            CALL PROVER( A(1,1), A(1,2), A(1,3) )
C           LE VECTEUR ORTHOGONAL AU RAYON DANS LE PLAN SI SI-1
            CALL PROVER( A(1,3), A(1,1), XYZTGS(1,NTG2) )
C           NORMALISATION A 1 DU VECTEUR TANGENT SI SI-1 EN SI
            CALL NORMER( 3, XYZTGS(1,NTG2), J )
C           NORME MISE A  RAYON * ANGLE(CSI, CSI-1)
            S = RAYON * ANGLE(I0)
            XYZTGS(1,NTG2) = XYZTGS(1,NTG2) * S
            XYZTGS(2,NTG2) = XYZTGS(2,NTG2) * S
            XYZTGS(3,NTG2) = XYZTGS(3,NTG2) * S
C
C           MISE A JOUR
            I0  = I
            NS0 = NS
 20      CONTINUE
 100  CONTINUE
      END
