      SUBROUTINE CALPN2(NTDL,LPLIGC,LPCOLC,LPDILU)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULE LE POINTEUR LPDILU : ADRESSE DANS AGC DU DERNIER
C ---   COEFFICIENT DE LA PARTIE TRIANGULAIRE INFERIEURE
C
C PARAMETRES D ENTREE :
C ---------------------
C LPLIGC,LPCOLC : LES POINTEURS ASSOCIES
C
C PARAMETRE DE SORTIE :
C --------------------
C LPDILU        : LE POINTEUR SUR LE DERNIER COEFFICIENT DE L
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1989
C23456---------------------------------------------------------------012
      DIMENSION LPLIGC(NTDL+1),LPCOLC(1),LPDILU(NTDL)
C
      K1=1
      DO 1 I=1,NTDL
         K2=LPLIGC(I+1)
         LPDILU(I)=K1-1
         DO 2 K=K1,K2-1
            IF(LPCOLC(K).LT.I) THEN
               LPDILU(I)=K
            ELSE
               GO TO 3
            END IF
2        CONTINUE
3        K1=K2+1
1     CONTINUE
C
      RETURN
      END
