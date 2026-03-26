      SUBROUTINE XYZS6C( AR6CUB, NA6CUB, CUB6DI, CUB6RG, NC6CUB, NBSTCH,
     %                   NUSTDF, XYZDF,  XYZCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 6 COORDONNEES DES SOMMETS D'UN MAILLAGE D'UN 6CUBE
C -----    CENTRE EN SOUS 6-CUBES FORMES PAR NA6CUB DIFFERENCES FINIES ET
C          NC6CUB COUCHES HOMOTHETIQUES EN PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C AR6CUB : LONGUEUR D'UNE ARETE DU NOYAU DANS UNE DIRECTION
C NA6CUB : NOMBRE D'ARETES DANS UNE DIRECTION DES DIFFERENCES FINIES
C CUB6DI : LARGEUR TOTALE DU 6-CUBE
C CUB6RG : PROGRESSION GEOMETRIQUE DES ARETES DES COUCHES
C NC6CUB : NOMBRE DE COUCHES HOMOTHETIQUES
C NBSTCH : NOMBRE DE SOMMETS D'UNE COUCHE
C
C SORTIES:
C --------
C NUSTDF : NO GLOBAL DES SOMMETS DE LA 1ERE COUCHE AU DELA DU NOYAU
C XYZDF  : 6 COORDONNEES DES SOMMETS DIFFERENCES FINIES
C XYZCH  : 6 COORDONNEES DES SOMMETS DES COUCHES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire JLL UMPC PARIS      SEPTEMBRE 2006
C MODIF : ALAIN PERRONNET LJLL UMPC & Saint Pierre du Perray   MARS 2009
C.......................................................................
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C     CONSTRUCTION DU TABLEAU XYZUVW DU MAILLAGE EN 6-CUBES
      REAL    XYZDF(6,0:NA6CUB,0:NA6CUB,0:NA6CUB,
     %                0:NA6CUB,0:NA6CUB,0:NA6CUB)
      INTEGER NUSTDF(0:NA6CUB,0:NA6CUB,0:NA6CUB,
     %               0:NA6CUB,0:NA6CUB,0:NA6CUB)
      REAL    XYZCH(6,1:NBSTCH,1:NC6CUB), X0
C
C     CALCUL DES 6 COORDONNEES DES SOMMETS DIFFERENCES FINIES
C     LE CUBE6 GENERAL EST CENTRE D'OU LA TRANSLATION DE LA
C     DEMI LARGEUR DU NOYAU
      X0 = - ( NA6CUB * AR6CUB ) / 2
C
C     NOYAU CENTRE: LES SOMMETS DIFFERENCES FINIES ET DE LA FRONTIERE
      NBSTDF = (NA6CUB+1)**6
      DO 60 N=0,NA6CUB
         DO 50 M=0,NA6CUB
            DO 40 L=0,NA6CUB
               DO 30 K=0,NA6CUB
                  DO 20 J=0,NA6CUB
                     DO 10 I=0,NA6CUB
C                       LES 6 COORDONNEES DU SOMMET DIFFERENCES FINIES
                        XYZDF(1,I,J,K,L,M,N) = X0 + I * AR6CUB
                        XYZDF(2,I,J,K,L,M,N) = X0 + J * AR6CUB
                        XYZDF(3,I,J,K,L,M,N) = X0 + K * AR6CUB
                        XYZDF(4,I,J,K,L,M,N) = X0 + L * AR6CUB
                        XYZDF(5,I,J,K,L,M,N) = X0 + M * AR6CUB
                        XYZDF(6,I,J,K,L,M,N) = X0 + N * AR6CUB
C
C                       LE SOMMET EST IL SUR LA FRONTIERE DES DIFFERENCES FINIES
                        IF( I.EQ.0 .OR. I.EQ.NA6CUB .OR.
     %                      J.EQ.0 .OR. J.EQ.NA6CUB .OR.
     %                      K.EQ.0 .OR. K.EQ.NA6CUB .OR.
     %                      L.EQ.0 .OR. L.EQ.NA6CUB .OR.
     %                      M.EQ.0 .OR. M.EQ.NA6CUB .OR.
     %                      N.EQ.0 .OR. N.EQ.NA6CUB ) THEN
C                           SOMMET SUR LA FRONTIERE
                            NBSTDF = NBSTDF + 1
C                           LE NUMERO GLOBAL DU SOMMET SUR LA PREMIERE COUCHE
                            NUSTDF(I,J,K,L,M,N) = NBSTDF
                        ENDIF
 10                  CONTINUE
 20               CONTINUE
 30            CONTINUE
 40         CONTINUE
 50      CONTINUE
 60   CONTINUE
C
C     NOMBRE DE SOMMETS FRONTALIERS OU D'UNE COUCHE
      WRITE(IMPRIM,*)'NB SOMMETS FRONTALIERS DES DIFFERENCES FINIES=',
     %                NBSTCH, (NA6CUB+1)**6-(NA6CUB-1)**6
      IF( NC6CUB .LE. 0 ) GOTO 9999
C
C     CALCUL DES 6 COORDONNEES DES SOMMETS DES NC6CUB COUCHES HOMOTHETIQUES
C     => NBSDCH SOMMETS NOUVEAUX PAR COUCHE
      CUBES2 = CUB6DI / 2
      DEPS   = 0.001*CUB6DI
      NUSDF  = 0
      DO 160 N=0,NA6CUB
         DO 150 M=0,NA6CUB
            DO 140 L=0,NA6CUB
               DO 130 K=0,NA6CUB
                  DO 120 J=0,NA6CUB
                     DO 110 I=0,NA6CUB
C
C                       LE SOMMET EST IL SUR LA FRONTIERE DES DIFFERENCES FINIES
                        IF( I.EQ.0 .OR. I.EQ.NA6CUB .OR.
     %                      J.EQ.0 .OR. J.EQ.NA6CUB .OR.
     %                      K.EQ.0 .OR. K.EQ.NA6CUB .OR.
     %                      L.EQ.0 .OR. L.EQ.NA6CUB .OR.
     %                      M.EQ.0 .OR. M.EQ.NA6CUB .OR.
     %                      N.EQ.0 .OR. N.EQ.NA6CUB ) THEN
C
C                           SOMMET SUR LA FRONTIERE
                            NUSDF = NUSDF + 1
C
C                           LES NC6CUB-1 PREMIERES COUCHES
                            DO 108 NBC=1,NC6CUB-1
                               AMPLI = CUB6RG**NBC
                               DO 106 NC=1,6
                                  R = XYZDF(NC,I,J,K,L,M,N) * AMPLI
                                  XYZCH(NC,NUSDF,NBC) = R
 106                           CONTINUE
 108                        CONTINUE
C
C                           LA DERNIERE COUCHE POUR EVITER LES ERREURS D'ARRONDI
                            AMPLI = CUB6RG**NC6CUB
                            DO 109 NC=1,6
                               R = XYZDF(NC,I,J,K,L,M,N) * AMPLI
                               IF( ABS(ABS(R)-CUBES2) .LE. DEPS ) THEN
                                  IF( R .GT. 0 ) THEN
                                     R =  CUBES2
                                  ELSE
                                     R = -CUBES2
                                  ENDIF
                               ENDIF
                               XYZCH(NC,NUSDF,NC6CUB) = R
 109                        CONTINUE
C
                        ENDIF
C
 110                 CONTINUE
 120              CONTINUE
 130           CONTINUE
 140        CONTINUE
 150     CONTINUE
 160  CONTINUE
C
 9999 WRITE(IMPRIM,*)'NB TOTAL de SOMMETS du 6-CUBE =',
     %               (NA6CUB+1)**6+NBSTCH*NC6CUB
      RETURN
      END
