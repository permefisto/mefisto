      SUBROUTINE CALPN4(NTDL,LPLIGC,LPCOLC,LPDILU)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULE LES TABLEAUX LPLIGC ET LPCOLC
C ---   DANS LE CAS D'UNE MATRICE AGC SYMETRIQUE
C
C PARAMETRES D ENTREE :
C ---------------------
C NTDL                 : LE NOMBRE TOTAL D'INCONNUES
C LPLIGC,LPCOLC,LPDILU : LES POINTEURS ASSOCIES
C
C PARAMETRE DE SORTIE  :
C --------------------
C LPLIGC,LPCOLC        : LES POINTEURS ASSOCIES
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1989
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DIMENSION LPLIGC(NTDL+1),LPCOLC(1),LPDILU(NTDL)
C
      LLPCOL=0
      K1=1
      DO 1 I=1,NTDL
         K2=LPDILU(I)
         DO 2 K=K1,K2
            LLPCOL=LLPCOL+1
            LPCOLC(LLPCOL)=LPCOLC(K)
2        CONTINUE
         LLPCOL=LLPCOL+1
         LPCOLC(LLPCOL)=I
         K1=LPLIGC(I+1)+1
         LPLIGC(I+1)=LLPCOL
1     CONTINUE
C
      RETURN
      END
