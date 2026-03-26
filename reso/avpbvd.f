C/UPDATE ADD NAME=AVPBVD,SSI=82121622
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                    MODULE AVPBVD
C                    -------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C BUT : FORMER  BS COMBINAISON LINEAIRE DES TABLEAUX B1 ET B2
C ----- BS(J) <= SOMME ALFA1D(I,J) * B1(I) + SOMME ALFA2D(I,J) * B2(I)
C                I=1,...,NDSM1              I=1,...,NDSM2
C       POUR J=1,...,NDSMS
C       VERSION REELLE DOUBLE PRECISION
C
C PARAMETRES D ENTREE :
C --------------------
C KK1   : NO DE LA PREMIERE COLONNE DE BS A TRAITER
C KK2   : NO DE LA DERNIERE COLONNE DE BS A TRAITER
C NTDL  : VALEUR MAXIMALE DU 2-EME INDICE DES TABLEAUX B1,B2,BS
C ALFA1D: TABLEAU REEL DOUBLE PRECISION
C NDSM1 : VALEUR MAXIMALE DU 1-ER INDICE DU TABLEAU B1
C K1    : NO DE LA COLONNE QUI PRECEDE LA 1-ERE DE LA PAGE DE B1
C B1    : TABLEAU B1(NDSM1,NTDL)
C ALFA2D: TABLEAU REEL DOUBLE PRECISION
C NDSM2 : VALEUR MAXIMALE DU 1-ER INDICE DU TABLEAU B2
C K2    : NO DE LA COLONNE QUI PRECEDE LA 1-ERE DE LA PAGE DE B2
C B2    : TABLEAU B2(NDSM2,NTDL)
C NDSMS : VALEUR MAXIMALE DU 1-ER INDICE DU TABLEAU BS
C KS    : NO DE LA COLONNE QUI PRECEDE LA 1-ERE DE LA PAGE DE BS
C BS    : TABLEAU BS(NDSMS,NTDL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS  ET IRIA NOVEMBRE 1979
C.......................................................................
       SUBROUTINE AVPBVD(KK1,KK2,NTDL,ALFA1D,NDSM1,K1,B1,
     &                                ALFA2D,NDSM2,K2,B2,
     &                                       NDSMS,KS,BS)
C.......................................................................
      DOUBLE PRECISION S,B1,B2,BS,ALFA1D,ALFA2D
      DIMENSION ALFA1D(NDSM1,NDSMS),ALFA2D(NDSM2,NDSMS),
     &          B1(NDSM1,NTDL),B2(NDSM2,NTDL),BS(NDSMS,NTDL)
C
            DO 1 K=KK1,KK2
            KKK1 = K - K1
            KKK2 = K - K2
            KKKS = K - KS
                 DO 2 J=1,NDSMS
                 S = 0.D0
                      DO 3 I=1,NDSM1
                      S = S + ALFA1D(I,J) * B1(I,KKK1)
    3                 CONTINUE
C
                      DO 4 I=1,NDSM2
                      S = S + ALFA2D(I,J) * B2(I,KKK2)
    4                 CONTINUE
                 BS(J,KKKS) = S
    2            CONTINUE
    1       CONTINUE
C
      RETURN
      END
