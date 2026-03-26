      SUBROUTINE DRREGC( NTDL, LPLIGC, LPDILU, LPCOLC, AGC, R, Z )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CE SP DESCEND LE SYSTEME TRIANGULAIRE INFERIEUR:
C -----                     (( L AGC )) (Z) = (R)
C       REMONTE LE SYSTEME TRIANGULAIRE SUPERIEUR:
C                           (( U AGC )) (Z) = (Z)
C
C   PARAMETRES D'ENTREE:
C   -------------------
C   AGC,LPLIGC,LPDILU,LPCOLC : MATRICE DE PRECONDITIONNEMENT
C                              ET POINTEURS ASSOCIES
C   R                        : SECOND MEMBRE
C   NTDL                     : NOMBRE D'INCONNUES
C
C   PARAMETRE DE SORTIE:
C   -------------------
C   Z                        : VECTEUR SOLUTION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Pascal JOLY LABORATOIRE D'ANALYSE NUMERIQUE PARIS 6  MAI 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AGC, Z, R, SOM
      DIMENSION         AGC(*), Z(NTDL), R(NTDL)
      DIMENSION         LPLIGC(NTDL+1), LPDILU(NTDL), LPCOLC(*)
C
C     LA DESCENTE  ( L AGC )  (Z) = (R)
C     ---------------------------------
      DO 1 I = 1, NTDL
         SOM = R(I)
         DO 2 K = LPLIGC(I)+1, LPDILU(I)
            SOM = SOM - AGC(K) * Z( LPCOLC(K) )
2        CONTINUE
         Z(I) = SOM
1     CONTINUE
C
C     LA REMONTEE   ( U AGC ) (Z) = (R)
C     ---------------------------------
      DO 3 I = NTDL, 1, -1
         SOM = Z(I)
         DO 4 K = LPDILU(I)+1, LPLIGC(I+1)-1
            SOM = SOM - AGC(K) * Z( LPCOLC(K) )
4        CONTINUE
         Z(I) = SOM * AGC( LPLIGC(I+1) )
3     CONTINUE
C
      RETURN
      END
