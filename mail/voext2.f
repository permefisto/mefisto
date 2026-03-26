        SUBROUTINE VOEXT2( NBSA,   NBSF,
     %                     XYZ1,   XYZ2,  XYZ3,  XYZ4, XYZ,
     %                     NSENS,  NPSOM,
     %                     XYSFTR, COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :         RECOPIER LES COORDONNEES DES SOMMETS DE CHACUNE DES
C -----         4 FACES DU TETRAEDRE COURBE A PARTIR DES DONNEES
C               POSITIONNER TOUS LES SOMMETS DES 4 FACES
C               TRIANGULAIRES RECTANGLES UNITE DU TETRAEDRE UNITE
C
C ENTREES:
C --------
C NBSA        : NOMBRE DE SOMMETS SUR CHAQUE ARETE
C NBSF        : NOMBRE DE SOMMETS DANS CHAQUE FACE
C XYZ1 A XYZ4 : LES COORDONNEES DES SOMMETS DES 4 FACES
C XYZ         : TABLEAU DE TRAVAIL
C NSENS(I)    : SENS DE PARCOURS DE LA FACE I PARMI LES DONNEES
C NPSOM(I)    : NUMERO DU PREMIER SOMMET LA FACE I
C
C SORTIES:
C --------
C XYSFTR      : LES 3 COORDONNEES DES SOMMETS DES 4 FACES DU TETRAEDRE UNITE
C COSO        : LES 3 COORDONNEES DES SOMMETS DES 4 FACES DU TETRAEDRE COURBE
C               DANS LA NUMEROTATION GLOBALE DES SOMMETS DU TETRAEDRE COURBE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY       ANALYSE NUMERIQUE UPMC PARIS     MARS  1989
C MODIFS : ALAIN  PERRONNET  ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1997
C23456---------------------------------------------------------------012
      INTEGER  NPSOM(4), NSENS(4)
      REAL     XYZ1(3,NBSF), XYZ2(3,NBSF), XYZ3(3,NBSF), XYZ4(3,NBSF),
     %         XYZ(3,NBSF)
      REAL     XYSFTR(2,NBSF,1:4), COSO(3,*)
C
C     INITIALISATION
C     ==============
      NBS = NBSA*(NBSA+1)*(NBSA+2)/6
C
C     CONSTRUCTION DU TETRAEDRE
C     =========================
C
C     1) LA FACE 1
C
      NF = 1
      CALL VOEXT4( NBSA, NPSOM(NF), NSENS(NF), XYZ1, XYZ )
      DO N1=1,NBSA
         NT = (N1-1)*N1*(N1+1)/6
         DO N2=1,N1
            NSTE = NT + N2*(N2-1)/2 + 1
            NSFA = N1*(N1-1)/2 + N2
            DO J=1,3
               COSO(J,NSTE) = XYZ(J,NSFA)
            ENDDO
         ENDDO
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS DE LA FACE 1 DU TETRAEDRE UNITE
      CALL TRCTR1( NBSA, XYZ, XYSFTR(1,1,1) )
C
C     2) LA  FACE 2
C
      NF = 2
      CALL VOEXT4( NBSA, NPSOM(NF), NSENS(NF), XYZ2, XYZ )
      DO N1=1,NBSA
         NT = (N1-1)*N1*(N1+1)/6 + N1*(N1-1)/2
         DO N2=1,N1
            NSTE = NT + N2
            NSFA = N1*(N1-1)/2 + N2
            DO J=1,3
               COSO(J,NSTE) = XYZ(J,NSFA)
            ENDDO
         ENDDO
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS DE LA FACE 2 DU TETRAEDRE UNITE
      CALL TRCTR1( NBSA, XYZ, XYSFTR(1,1,2) )
C
C     3) LA  FACE 3
C
      NF = 3
      CALL VOEXT4( NBSA, NPSOM(NF), NSENS(NF), XYZ3, XYZ )
      DO N1=1,NBSA
         NT = (N1-1)*N1*(N1+1)/6
         DO N2=1,N1
            NSTE = NT + N2*(N2+1)/2
            NSFA = N1*(N1-1)/2 + N1 + 1 - N2
            DO J=1,3
               COSO(J,NSTE) = XYZ(J,NSFA)
            ENDDO
         ENDDO
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS DE LA FACE 3 DU TETRAEDRE UNITE
      CALL TRCTR1( NBSA, XYZ, XYSFTR(1,1,3) )
C
C     4) LA  FACE 4
C
      NF = 4
      CALL VOEXT4( NBSA, NPSOM(NF), NSENS(NF), XYZ4, XYZ )
      NT = (NBSA-1)*NBSA*(NBSA+1)/6
      DO N1=1,NBSA
         DO N2=1,N1
            NSFA = (N1-1)*N1/2 + N2
            NSTE = NT + NSFA
            DO J=1,3
               COSO(J,NSTE) = XYZ(J,NSFA)
            ENDDO
         ENDDO
      ENDDO

C     CALCUL DES COORDONNEES DES SOMMETS DE LA FACE 4 DU TETRAEDRE UNITE
      CALL TRCTR1( NBSA, XYZ, XYSFTR(1,1,4) )
C
      RETURN
      END
