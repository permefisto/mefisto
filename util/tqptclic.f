      SUBROUTINE TQPTCLIC( NX, NY, NDIM, XYZSOM, NBTQ, NOSOTQ,  NUMTQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO DU TRIANGLE OU QUADRANGLE 3D VISIBLE
C -----    CONTENANT LE CLIC SOURIS ET LE PLUS PROCHE DE L'OEIL

C ENTREES:
C --------
C NX     : ABSCISSE ECRAN DU POINT CLIQUE
C NY     : ORDONNEE ECRAN DU POINT CLIQUE
C NDIM   : DIMENSION (2 ou 3) DE L'ESPACE DE LA SURFACE
C XYZSOM : XYZ 3 COORDONNEES CARTESIENNES DES SOMMETS DU MAILLAGE
C NBTQ   : NOMBRE D'EF (TRIANGLES OU QUADRANGLES) DU MAILLAGE
C NOSOTQ : NO DES 4 SOMMETS DES TRIANGLES ET/OU QUADRANGLES DU MAILLAGE

C SORTIES:
C --------
C NUMTQ  : >0 NUMERO DU TRIANGLE OU QUADRANGLE CONTENANT LE POINT CLIQUE
C          =0 SI AUCUN  TRIANGLE OU QUADRANGLE CONTIENT  LE POINT CLIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      REAL         XYZSOM(1:3,*)
      INTEGER      NOSOTQ(1:4,1:NBTQ)
C
      REAL         XYZAX(3,4), XYZAPT(3)
      EQUIVALENCE (XYZAX(1,1),X1),(XYZAX(2,1),Y1),(XYZAX(3,1),Z1),
     %            (XYZAX(1,2),X2),(XYZAX(2,2),Y2),(XYZAX(3,2),Z2),
     %            (XYZAX(1,3),X3),(XYZAX(2,3),Y3),(XYZAX(3,3),Z3),
     %            (XYZAX(1,4),X4),(XYZAX(2,4),Y4),(XYZAX(3,4),Z4)
      EQUIVALENCE (XYZAPT(1),XAPT), (XYZAPT(2),YAPT), (XYZAPT(3),ZAPT)
C
      NUMTQ = 0

C     COORDONNEES AXONOMETRIQUES DU POINT CLIQUE EN NX NY PIXELS
C     PAS D'INFORMATION DE LA 3-EME COORDONNEE DANS LA DIRECTION DE VISEE
C     ===================================================================
      XAPT = XOB2PX( NX )
      YAPT = YOB2PX( NY )

C     BOUCLE SUR LES TRIANGLES OU QUADRANGLES DU MAILLAGE
C     ===================================================
      ZAMAX = -1E28
      DO I=1,NBTQ

C       LES COORDONNEES DES SOMMETS DU TRIANGLE ou QUADRANGLE NUMERO I
         IF( NDIM .EQ. 2 ) THEN
            DO J=1,4
               N = NOSOTQ(J,I)
               XYZAX(1,J) = XYZSOM(1,N)
               XYZAX(2,J) = XYZSOM(2,N)
               XYZAX(3,J) = XYZSOM(3,N)
            ENDDO
         ELSE
            DO J=1,4
               N = NOSOTQ(J,I)
               IF( N .GT. 0 ) THEN
C                 CALCUL DES 3 COORDONNEES AXONOMETRIQUES XYZAX
C                 A PARTIR DES 3 COORDONNEES XYZSOM DU SOMMET N
                  CALL XYZAXO( XYZSOM(1,N), XYZAX(1,J) )
               ENDIF
            ENDDO
         ENDIF
C
C        LE POINT (NX,NY) EST DANS LE TRIANGLE ou QUADRANGLE NUMTQ
C        CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU POINT (XAPT,YAPT)
C        DANS LE TRIANGLE NUMTQ ou LE TRIANGLE 123 DU QUADRANGLE NUMTQ
         D   =    (X2-X1 )*(Y3-Y1 )  -  (X3-X1 )*(Y2-Y1)
         CB1 = ( (X2-XAPT)*(Y3-YAPT) - (X3-XAPT)*(Y2-YAPT) ) / D
         CB2 = ( (X3-XAPT)*(Y1-YAPT) - (X1-XAPT)*(Y3-YAPT) ) / D
         CB3 = ( (X1-XAPT)*(Y2-YAPT) - (X2-XAPT)*(Y1-YAPT) ) / D
         IF( CB1 .GE. 0.0 .AND. CB1 .LE. 1.0 .AND.
     %       CB2 .GE. 0.0 .AND. CB2 .LE. 1.0 .AND.
     %       CB3 .GE. 0.0 .AND. CB3 .LE. 1.0 ) THEN
C
C           POINT INTERIEUR AU TRIANGLE ou QUADRANGLE I
            IF( NDIM .EQ. 2 ) THEN
C              LE PREMIER EF ATTEINT EST LE BON
               NUMTQ = I
               GOTO 9000
            ELSE
C              PARCOURS DE TOUS LES EF POUR TROUVER L'EF VISIBLE
               ZAXOBA = (Z1+Z2+Z3) / 3
               IF( ZAXOBA .GT. ZAMAX ) THEN
                  NUMTQ = I
                  ZAMAX = ZAXOBA
               ENDIF
            ENDIF
C
         ELSE
C
C           POINT EXTERIEUR AU TRIANGLE I 123
            IF( NOSOTQ(4,I) .GT. 0 ) THEN
C
C              I EST UN QUADRANGLE 1234
C              CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU POINT (XAPT,YAPT)
C              DANS LE TRIANGLE 1 3 4 DU QUADRANGLE NUMTQ
               D   =    (X3-X1 )*(Y4-Y1)   -  (X4-X1 )*(Y3-Y1)
               CB1 = ( (X3-XAPT)*(Y4-YAPT) - (X4-XAPT)*(Y3-YAPT) ) / D
               CB2 = ( (X4-XAPT)*(Y1-YAPT) - (X1-XAPT)*(Y4-YAPT) ) / D
               CB3 = ( (X1-XAPT)*(Y3-YAPT) - (X3-XAPT)*(Y1-YAPT) ) / D
               IF( CB1 .GE. 0.0 .AND. CB1 .LE. 1.0 .AND.
     %             CB2 .GE. 0.0 .AND. CB2 .LE. 1.0 .AND.
     %             CB3 .GE. 0.0 .AND. CB3 .LE. 1.0 ) THEN
C
C                 POINT INTERIEUR AU TRIANGLE 134 DU QUADRANGLE I
                  IF( NDIM .EQ. 2 ) THEN
C                    LE PREMIER EF ATTEINT EST LE BON
                     NUMTQ = I
                     GOTO 9000
                  ELSE
C                    PARCOURS DE TOUS LES EF POUR TROUVER L'EF VISIBLE
                     ZAXOBA = (Z1+Z3+Z4) / 3
                     IF( ZAXOBA .GT. ZAMAX ) THEN
                        NUMTQ = I
                        ZAMAX = ZAXOBA
                     ENDIF
                  ENDIF
C
               ENDIF
            ENDIF
         ENDIF
C
      ENDDO
C
 9000 IF( NUMTQ .GT. 0 ) THEN
         PRINT*,'tqptclic: TQANGLE CLIQUE: No',NUMTQ,' St:',
     %          (NOSOTQ(K,NUMTQ),K=1,4)
      ENDIF

      RETURN
      END
