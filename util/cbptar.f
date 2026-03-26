      SUBROUTINE CBPTAR( S1, S2, PT, CBPTA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA COORDONNEE BARYCENTRIQUE DU POINT PT
C -----    SUR L'ARETE S1-S2 (PT DOIT ETRE SUR LA DROITE S1-S2)
C
C ENTREES:
C --------
C S1,S2  : LES 2 POINTS QUI DEFINISSENT L'ARETE
C PT     : LE POINT DE LA DROITE S1-S2 COORDONNEE BARYCENTRIQUE A CALCULER
C
C SORTIE :
C --------
C CBPTA  : LA COORDONNEE BARYCENTRIQUE DU POINT PT DANS L'ARETE S1-S2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC& St Pierre du Perray SEPTEMBRE 2011
C2345X7..............................................................012
      DOUBLE PRECISION  S1(3), S2(3), PT(3), CBPTA, X1
C
C     CBPTA = PT-S1 / S2-S1
      CBPTA = 0D0
      NBC   = 0
C
      DO K = 1, 3
         X1 = S1(K)
         IF( X1 .NE. S2(K) ) THEN
            NBC   = NBC + 1
            CBPTA = CBPTA + ( PT(K)-X1 ) / ( S2(K)-X1 )
         ENDIF
      ENDDO
C
C     VALEUR MOYENNE
      CBPTA = CBPTA / NBC
C
      RETURN
      END
