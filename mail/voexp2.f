      SUBROUTINE VOEXP2( NBSA,   NBSH,
     %                   XYZ1,   XYZ2,   XYZ3,   XYZ4,  XYZ5,
     %                   NSENS,  NPSOM,
     %                   XYZFAC, XYSTR1, XYSTR5, XYSCA1,
     %                   COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :         RECOPIER LES COORDONNEES DES SOMMETS DE CHAQUE FACE
C -----         DU PENTAEDRE COURBE A PARTIR DES DONNEES
C               POSITIONNER TOUS LES POINTS DU BORD DU PENTAEDRE UNITE
C
C ENTREES:
C --------
C NBSA,NBSH   : NOMBRE DE SOMMETS DANS CHAQUE DIRECTION
C XYZ1 A XYZ5 : LES COORDONNEES DES SOMMETS DES 5 FACES AVANT RENUMEROTATION
C NSENS(I)    : SENS DE PARCOURS DE LA FACE I PARMI LES DONNEES
C NPSOM(I)    : NUMERO DU PREMIER SOMMET LA FACE I
C
C TABLEAU AUXILIAIRE:
C -------------------
C XYZFAC : LES 3 COORDONNEES DES SOMMETS D'UNE FACE DU PENTAEDRE
C
C SORTIES:
C --------
C XYSTR1 : LES 2 COORDONNEES DES POINTS SUR LE TRIANGLE UNITE
C          DE LA FACE 1 DU PENTAEDRE
C XYSTR5 : LES 2 COORDONNEES DES POINTS SUR LE TRIANGLE UNITE
C          DE LA FACE 5 DU PENTAEDRE
C XYSCA1 : LES 2 COORDONNEES DES POINTS SUR LE CARRE UNITE
C          DES FACES 2 3 4 DU PENTAEDRE
C COSO   : LES COORDONNEES DES POINTS DES 5 FACES DU PENTAEDRE COURBE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1988
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1997
C23456---------------------------------------------------------------012
      INTEGER  NPSOM(5),  NSENS(5)
      REAL     XYZ1(3,*), XYZ2(3,*), XYZ3(3,*), XYZ4(3,*), XYZ5(3,*)
      REAL     XYZFAC(3,*)
      REAL     XYSTR1(2,*), XYSTR5(2,*)
      REAL     XYSCA1(2,NBSA,NBSH,2:4)
      REAL     COSO(3,*)
C
C     INITIALISATION
C     ==============
      NBST = NBSA * (NBSA+1) / 2
      NBS  = NBST * NBSH
C
C     CONSTRUCTION DU PENTAEDRE COURBE ET DES FACES RECTANGLE UNITE
C     =============================================================
C
C     1) LA FACE 1 TRIANGULAIRE DE SOMMETS S1 S2 S3
C        ------------------------------------------
      NF = 1
C     LE SOMMETS DE LA FACE SONT REORDONNES DANS LA CONFIGURATION DU PENTAEDRE
      CALL VOEXT4( NBSA, NPSOM(NF), NSENS(NF), XYZ1, XYZFAC )
C     LES SOMMETS DE LA FACE SONT PLONGES DANS LE PENTAEDRE
      DO NB1=1,NBSA
         DO NB2=1,NB1
            NSPE = (NB1-1)*NB1/2+NB2
            NSFA = NSPE
            DO J=1,3
               COSO(J,NSPE) = XYZFAC(J,NSFA)
            ENDDO
         ENDDO
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS INTERNES A LA FACE 1 DU PENTAEDRE UNITE
      CALL TRCTR1( NBSA, XYZFAC, XYSTR1 )
C
C     2) LA  FACE 2 QUADRANGULAIRE DE SOMMETS S1 S2 S5 S4
C        ------------------------------------------------
      NF = 2
      CALL VOEXP4(NBSA,NBSH,NPSOM(NF),NSENS(NF),N1,N2,N3)
      NSPLAN = 0
      NSF    = 0
      DO NB1=1,NBSH
         DO NB2=1,NBSA
            NSPE = NSPLAN + NB2*(NB2-1)/2+1
            NSFA = N1 + N2 * NB1 + N3 * NB2
            NSF  = NSF + 1
            DO J=1,3
               COSO(J,NSPE) = XYZ2(J,NSFA)
               XYZFAC(J,NSF) = XYZ2(J,NSFA)
            ENDDO
         ENDDO
         NSPLAN = NSPLAN + NBST
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS INTERNES A LA FACE 2 DU PENTAEDRE UNITE
      CALL QUCCA1( NBSA, NBSH, XYZFAC, XYSCA1(1,1,1,2) )
C
C     3) LA  FACE 3 QUADRANGULAIRE DE SOMMETS S2 S3 S6 S5
C        ------------------------------------------------
      NF = 3
      CALL VOEXP4(NBSA,NBSH,NPSOM(NF),NSENS(NF),N1,N2,N3)
      NSPLAN = 0
      NSF    = 0
      NDECAL = (NBSA-1)*NBSA/2
      DO NB1=1,NBSH
         DO NB2=1,NBSA
            NSF  = NSF + 1
            NSPE = NSPLAN + NDECAL + NB2
            NSFA = N1 + N2 * NB1 + N3 * NB2
            DO J=1,3
               COSO(J,NSPE) = XYZ3(J,NSFA)
               XYZFAC(J,NSF) = XYZ3(J,NSFA)
            ENDDO
         ENDDO
         NSPLAN = NSPLAN + NBST
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS INTERNES A LA FACE 3 DU PENTAEDRE UNITE
      CALL QUCCA1( NBSA, NBSH, XYZFAC, XYSCA1(1,1,1,3) )
C
C     4) LA  FACE 4 QUADRANGULAIRE DE SOMMETS S3 S1 S4 S6
C        ------------------------------------------------
      NF = 4
      CALL VOEXP4(NBSA,NBSH,NPSOM(NF),NSENS(NF),N1,N2,N3)
      NSF    = 0
      NSPLAN = 0
      DO NB1=1,NBSH
         DO NB2=1,NBSA
            NSF  = NSF + 1
            NB3  = NBSA + 1 - NB2
            NSPE = NSPLAN + NB3*(NB3+1)/2
            NSFA = N1 + N2 * NB1 + N3 * NB2
            DO J=1,3
               COSO(J,NSPE) = XYZ4(J,NSFA)
               XYZFAC(J,NSF) = XYZ4(J,NSFA)
            ENDDO
         ENDDO
         NSPLAN = NSPLAN + NBST
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS INTERNES A LA FACE 4 DU PENTAEDRE UNITE
      CALL QUCCA1( NBSA, NBSH, XYZFAC, XYSCA1(1,1,1,4) )
C
C     5) LA  FACE 5 TRIANGULAIRE DE SOMMETS S4 S5 S6
C        -------------------------------------------
      NF = 5
      CALL VOEXT4( NBSA, NPSOM(NF), NSENS(NF), XYZ5, XYZFAC )
      NSTP = NBST * (NBSH-1)
      DO NB1=1,NBSA
         DO NB2=1,NB1
            NSFA = (NB1-1)*NB1/2 + NB2
            NSPE = NSTP + NSFA
            DO J=1,3
               COSO(J,NSPE) = XYZFAC(J,NSFA)
            ENDDO
         ENDDO
      ENDDO
C     CALCUL DES COORDONNEES DES SOMMETS INTERNES A LA FACE 5 DU PENTAEDRE UNITE
      CALL TRCTR1( NBSA, XYZFAC, XYSTR5 )
C
      RETURN
      END
