      SUBROUTINE ORTVVP( NBDL, NCODSB, MUB, B, COMXB,
     %                   NCVALP, NBVALP, VALP, VECP, V, W )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    B-ORTHONORMALISER LES VECTEURS DE VECP PAR RAPPORT A EUX-MEMES
C -----
C          SI UN VECTEUR EST NUL APRES B-ORTHONORMALISATION TOUTES LES
C          VALEURS ET VECTEURS PROPRES SONT DECALES
C
C ENTREES:
C --------
C NBDL   : NOMBRE DE LIGNES ET COLONNES DE LA MATRICE B
C NCODSB : CODE 0 => B DIAGONALE , 1=> B PROFIL SYMETRIQUE
C MUB    : POINTEUR DU PROFIL DE B
C B      : LA MATRICE (DE MASSE) B
C COMXB  : MAXIMUM ACTUEL DES COEFFICIENTS DE LA MATRICE M
C
C NCVALP : NOMBRE DE VECTEURS DECLARES
C
C MODIFIES:
C ---------
C NBVALP : NOMBRE DE VALEURS PROPRES
C VALP   : LES VALEURS PROPRES
C VECP   : VECTEURS A B-ORTHONORMALISES
C V,W    : TABLEAUX AUXILIAIRES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET ANALYSE NUMERIQUE PARIS ET TEXAS A&M    AOUT 2003
C ......................................................................
      include"./incl/langue.inc"
      include"./incl/epsvvp.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  B(*), VALP(*), VECP(NCVALP,NBDL), V(*), W(*)
      DOUBLE PRECISION  COMXB, S, SMAX
      INTEGER           MUB(*)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'ENTREE DE ORTVVP: NOMBRE VALEURS PROPRES=',
     %                   NBVALP
      ELSE
         WRITE(IMPRIM,*)'START of ORTVVP: NUMBER of EIGENVALUES=',
     %                   NBVALP
      ENDIF
      NBV  = 0
      SMAX = 0D0
C
      DO 100 J=1,NBVALP
         DO 20 L=1,NBDL
            V(L) = VECP(J,L)
 20      CONTINUE
C
C        W = B * V  ET V=VECP(J)
         CALL APBEBD( NCODSB,MUB,B,V,1,NBDL, W )
C
C        B-ORTHOGONALISATION DE V PAR RAPPORT AUX VECP(1),...,VECP(NBV)
C        --------------------------------------------------------------
         DO 50 I=1,NBV
C           S=( VECT(I), B VECT(J) )
            S=0D0
            DO 30 L=1,NBDL
               S = S + VECP(I,L) * W(L)
 30         CONTINUE
C           VECT(J) = VECT(J) - S * VECP(I)
            DO 40 L=1,NBDL
               V(L) = V(L) - S * VECP(I,L)
 40         CONTINUE
 50      CONTINUE
C
C        B-ORTHONORMALISATION DE V ET MISE DANS VECP(J)
         CALL APBEBD( NCODSB,MUB,B,V,1,NBDL, W )
         S=0D0
         DO 60 L=1,NBDL
            S = S + V(L) * W(L)
 60      CONTINUE
C
         IF( S .GT. SMAX ) SMAX = S
CCC         WRITE(IMPRIM,10060) J, S, SMAX
CCC10060    FORMAT('ORTVVP: NORME tV M V',I5,' = ',G15.7,' max=',G15.7)
C
CCC         IF( S .GT. PRECIS*SMAX ) THEN
         IF( S .GT. PRECIS*COMXB ) THEN
C
C           B-ORTHONORMALISATION DE V ET MISE DANS VECP(NBV)
C           ------------------------------------------------
            NBV = NBV + 1
            S   = 1D0 / SQRT( ABS(S) )
            DO 90 L=1,NBDL
               VECP(NBV,L) = V(L) * S
 90         CONTINUE
C           DECALAGE DE LA VALEUR PROPRE
            VALP(NBV) = VALP(J)
C
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
CCC               WRITE(IMPRIM,10010) J,S,SMAX
               WRITE(IMPRIM,*) 'VALEUR PROPRE PERDUE=', VALP(J)
               WRITE(IMPRIM,10010) J,S,COMXB
            ELSE
CCC               WRITE(IMPRIM,20010) J,S,SMAX
               WRITE(IMPRIM,*) 'LOST EIGENVALUE=', VALP(J)
               WRITE(IMPRIM,20010) J,S,COMXB
            ENDIF
         ENDIF
 100  CONTINUE
C
10010 FORMAT('ORTVVP: APRES M-ORTHONORMALISATION NORME VECTEUR',I5,
     %' =',G15.7,' <',G15.7,' SUPPOSE NUL => IL EST SUPPRIME')
20010 FORMAT('ORTHOD: After M-ORTHONORMALISATION, VECTOR NORM',I5,
     %' =',G15.7,' <',G15.7,' SUPPOSED NULL => It is DELETED')
C
C     LE NOMBRE REEL DE VALEURS ET VECTEURS PROPRES
      NBVALP = NBV
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'SORTIE DE ORTVVP: NOMBRE VALEURS PROPRES=',
     %                   NBVALP
      ELSE
         WRITE(IMPRIM,*)'EXIT  of ORTVVP: NUMBER of EIGENVALUES=',
     %                   NBVALP
      ENDIF
      END
