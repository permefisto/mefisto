      SUBROUTINE NUTDTS( KNOMTS , NUTD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETROUVER LE NUMERO DU TD A PARTIR DU NOM DE TMS ET
C ----- DU DICTIONNAIRE DES NOMS DE TABLEAUX DESCRIPTEURS
C
C   EXEMPLE: TS: ~>POINT>P1>XYZSOMMET"MT10 => TD: ~>POINT>>XYZSOMMET
C
C ENTREE :
C --------
C KNOMTS : NOM DU TABLEAU TMS
C
C SORTIE :
C --------
C NUTD   : NUMERO DU TABLEAU DESCRIPTEUR DANS LE DICTIONNAIRE DICOTD
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      CHARACTER*(*)     KNOMTS
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     TRANSFORMATION EN MINUSCULES DU NOM DU TMS
      CALL MINUSC( KNOMTS )
C
C     LE 1-ER CARACTERE DIFFERENT DE ~>
      N = 1
      IF( KNOMTS(1:1) .EQ. '~' ) N = 2
      IF( KNOMTS(N:N) .EQ. '>' ) N = N + 1
C
C     LE NOMBRE DE CARACTERES DE KNOMTS  ET DES TD
      LTS = LEN( KNOMTS )
      LTD = 48
C
C     BOUCLE SUR LES TABLEAUX DESCRIPTEURS
      DO 100 NUTD = 1 , NBTD
C        LE 1-ER CARACTERE DE NOM DE TMS A ANALYSER
         N1 = N
C        LE 1-ER CARACTERE DE NOM DE TD A ANALYSER
         M1 = 3
C        ANALYSE CARACTERE APRES CARACTERE
C        RECHERCHE DU PROCHAIN > DANS LE KNOMTS
 10      N2 = INDEX( KNOMTS(N1:LTS) , '>' )
C        RECHERCHE DU PROCHAIN > DANS LE KNOMTD
         M2 = INDEX( DICOTD(NUTD)(M1:LTD) , '>' )
         IF( N2 .GT. 0 ) THEN
            N2 = N1 + N2 - 2
C           S'IL EXISTE UN SUFFIXE, ALORS IL EST SUPPRIME
            NSUF = INDEX( KNOMTS(N1:N1+N2-1) , KSUFIX )
            IF( NSUF .GT. 0 ) THEN
               N3 = N1 + NSUF - 2
            ELSE
               N3 = N2
            ENDIF
C           KNOMTS(N1:N3) EST UN CHAMP DU NOM DE TMS
            IF( M2 .GT. 0 ) THEN
               M2 = M1 + M2 - 2
C              NOMTD(M1:M2) EST UN CHAMP DU NOM DE TD
               IF( M2 .LT. M1 ) THEN
C                 >> DANS LE TD
C                 PASSAGE AU CHAMP DU NOMTS SUIVANT
                  N1 = N2 + 2
                  M1 = M1 + 1
                  GOTO 10
               ELSE
C                 2 CHAMPS QUI DOIVENT ETRE IDENTIQUES
                  IF( KNOMTS(N1:N3) .EQ. DICOTD(NUTD)(M1:M2) ) THEN
C                    CHAMPS EGAUX. PASSAGE AU SUIVANT
                     N1 = N2 + 1
                     M1 = M2 + 1
                     GOTO 10
                  ELSE
C                    CHAMPS DIFFERENTS
                     GOTO 100
                  ENDIF
               ENDIF
            ELSE
C              IL EXISTE UN CHAMP TMS ET PAS DANS TD
               GOTO 100
            ENDIF
         ELSE
C           C'EST LE DERNIER CHAMP DU NOMTS
            IF( M2 .GT. 0 ) THEN
C              CE N'EST PAS LE DERNIER CHAMP DE NOMTD
               GOTO 100
            ELSE
C              C'EST AUSSI LE DERNIER CHAMP DU NOMTD
C              LES 2 DOIVENT ETRE IDENTIQUES
               N2 = INDEX( KNOMTS(N1:LTS) , ' ' )
               IF( N2 .EQ. 0 ) THEN
                  N2 = LTS
               ELSE
                  N2 = N1 + N2 - 2
               ENDIF
C              S'IL EXISTE UN SUFFIXE, ALORS IL EST SUPPRIME
               NSUF = INDEX( KNOMTS(N1:N2) , KSUFIX )
               IF( NSUF .GT. 0 ) N2 = N1 + NSUF - 2
               M2 = INDEX( DICOTD(NUTD)(M1:LTD) , ' ' )
               IF( M2 .EQ. 0 ) THEN
                  M2 = LTD
               ELSE
                  M2 = M1 + M2 - 2
               ENDIF
               IF( KNOMTS(N1:N2) .EQ. DICOTD(NUTD)(M1:M2) ) THEN
C                 LE NOMTS CORRESPOND AU NOMTD
C                 TRANSFORMATION EN MAJUSCULES DU NOM DU TMS
                  CALL MAJUSC( KNOMTS )
                  RETURN
               ELSE
                  GOTO 100
               ENDIF
            ENDIF
         ENDIF
 100  CONTINUE
C
C     AUCUN NOM DE TD N'EST CORRECT POUR CE TMS
      NUTD = 0
C
      RETURN
      END
