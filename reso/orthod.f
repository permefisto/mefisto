      SUBROUTINE ORTHOD( NBDL, MUB, B, NCODSB, COMXB,
     %                   VECPI, NCI, NVI,
     %                   VECP,  NC,
     %                   V, W,
     %                   VALP, R )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    B-ORTHONORMALISER LES VECTEURS DE VECP PAR RAPPORT A EUX-MEMES
C -----    ET AUX NCI VECTEURS DE VECPI AVEC LE RESULTAT DANS R
C
C          SI UN VECTEUR EST NUL APRES B-ORTHONORMALISATION TOUTES LES
C          VALEURS ET VECTEURS PROPRES SONT DECALES
C          ET LA DERNIERE VALEUR PROPRE EST FIXEE A 100 * PRECEDENTE
C          ET LE VECTEUR ASSOCIE EST GENERE DE FACON ALEATOIRE
C
C ENTREES:
C --------
C NBDL   : NOMBRE DE LIGNES ET COLONNES DE LA MATRICE B
C MUB    : POINTEUR DU PROFIL DE B
C B      : LA MATRICE (DE MASSE) B
C NCODSB : CODE 0 => B DIAGONALE , 1=> B PROFIL SYMETRIQUE
C COMXB  : COEFFICIENT MAX DE LA MATRICE B
C
C VECPI  : VECTEURS B-ORTHONORMALISES DEJA CALCULES
C NCI    : NOMBRE DECLARE DE VECTEURS VECPI, 0 EST POSSIBLE
C NVI    : NOMBRE DE VECTEURS DEJA CALCULES A PRENDRE EN COMPTE
C VECP   : VECTEURS A B-ORTHONORMALISER
C NC     : NOMBRE DE TELS VECTEURS VECP A B-ORTHONORMALISER
C
C V,W    : TABLEAUX AUXILIAIRES
C
C SORTIES:
C --------
C VALP   : LES VALEURS PROPRES
C R      : LES NC VECTEURS B-ORTHONORMALISES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET ANALYSE NUMERIQUE PARIS ET INRIA    NOVEMBRE 1989
C AUTEUR : ALAIN PERRONNET     TEXAS A & M UNIVERSITY       JUILLET 2003
C ......................................................................
      include"./incl/langue.inc"
      include"./incl/epsvvp.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           MUB(1+NBDL)
      DOUBLE PRECISION  S, COMXB, B(*), VECPI(NCI,NBDL), VECP(NC,NBDL),
     %                  R(NC,NBDL), V(NBDL), W(NBDL), VALP(NC)
      DOUBLE PRECISION  VV0, VV1
ccc      DOUBLE PRECISION  SMAX
C
ccc      SMAX = 0D0
      N1   = 1
C     NOMBRE DE RETRAIT DE VALEURS PROPRES
      NBRETR = 0
C
      DO 100 J=1,NC
C
C        LE VECTEUR J DE PRODUIT SCALAIRE VV0
 1       VV0 = 0D0
         DO 2 L=1,NBDL
            V(L) = VECP(J,L)
            VV0  = VV0 + V(L) ** 2
 2       CONTINUE
C
C        W = B * V  ET V=VECP(J)
         CALL APBEBD( NCODSB, MUB, B, V, 1, NBDL,  W )
C
C        B-ORTHOGONALISATION DE V=VECP(J) AUX VECPI(1),...,VECPI(NVI)
C        ------------------------------------------------------------
         DO 5 I=1,NVI
            S=0D0
            DO 3 L=1,NBDL
               S = S + VECPI(I,L) * W(L)
    3       CONTINUE
            DO 4 L=1,NBDL
               V(L) = V(L) - S * VECPI(I,L)
    4       CONTINUE
    5    CONTINUE
C
C        B-ORTHOGONALISATION DE V PAR RAPPORT AUX VECP(1),...,VECP(J-1)
C        --------------------------------------------------------------
         DO 8 I=1,J-1
            S=0D0
            DO 6 L=1,NBDL
               S = S + R(I,L) * W(L)
 6          CONTINUE
            DO 7 L=1,NBDL
               V(L) = V(L) - S * R(I,L)
 7          CONTINUE
 8       CONTINUE
C
C        PRODUIT SCALAIRE DU V B-ORTHONORMALISE
         VV1 = 0D0
         DO 9 L=1,NBDL
            VV1 = VV1 + V(L) ** 2
 9       CONTINUE
C
C        B-ORTHONORMALISATION DE V ET MISE DANS R(J)
C        -------------------------------------------
         CALL APBEBD( NCODSB, MUB, B, V, 1, NBDL,  W )
         S=0D0
         DO 10 L=1,NBDL
            S = S + V(L) * W(L)
 10      CONTINUE
ccc
ccc         IF( S .GT. SMAX ) SMAX = S
C
CCC         IF( LANGAG .EQ. 0 ) THEN
CCC            WRITE(IMPRIM,10009) J, S, SMAX, COMXB
CCC         ELSE
CCC            WRITE(IMPRIM,20009) J, S, SMAX, COMXB
CCC         ENDIF
CCC10009 FORMAT('ORTHOD: NORME tVP M VP',I5,'=',G15.7,' smax=',G15.7,
CCC     %' comxb=',G15.7)
CCC20009 FORMAT('ORTHOD: NORM tVP M VP',I5,'=',G15.7,' smax=',G15.7,
CCC     %' comxb=',G15.7)
cccC
CCC         IF( S .LE. PRECIS*SMAX ) THEN
ccc
         IF( VV1 .LE. VV0*PRECIS .AND. S .LE. PRECIS*COMXB ) THEN
C
C           SI PLUS DE VALEURS PROPRES AU DELA DE J => RETOUR
            IF( J+NBRETR .GT. NC ) RETURN
            NBRETR = NBRETR + 1
C
C           VECTEUR NUL APRES B-ORTHONORMALISATION => REPRISE
C           -------------------------------------------------
            IF( LANGAG .EQ. 0 ) THEN
CCC               WRITE (IMPRIM,10020) J, S, PRECIS*SMAX
              WRITE (IMPRIM,10020) VV0, VV1, VALP(J), J, S, PRECIS*COMXB
            ELSE
CCC               WRITE (IMPRIM,20020) J, S, PRECIS*SMAX
              WRITE (IMPRIM,20020) VV0, VV1, VALP(J), J, S, PRECIS*COMXB
            ENDIF
10020 FORMAT('ORTHOD: VV0=',G12.4,' VV1=',G12.4,' tVP [M] VP=',G12.4,
     %' VECTEUR',I5,' NORME=',G12.4,
     %'<',G12.4,' NUL => NOUVELLE INITIALISATION ALEATOIRE')
20020 FORMAT('ORTHOD: VV0=',G12.4,' VV1=',G12.4,' tVP [M] VP=',G12.4,
     %' VECTOR',I5,' NORM=',G12.4,
     %'<',G12.4,' NULL => NEW RANDOM INITIALISATION')
C
C           TRANSLATION DES VALEURS ET VECTEURS PROPRES
C           POUR SUPPRIMER LA VALEUR ET VECTEUR PROPRE J
            DO 40 I=J+1,NC
               I1 = I - 1
               VALP(I1) = VALP(I)
               DO 30 L=1,NBDL
                  VECP(I1,L) = VECP(I,L)
 30            CONTINUE
 40         CONTINUE
C
C           GENERATION ALEATOIRE DU DERNIER VECTEUR PROPRE
            VALP(NC) = 100 * VALP(NC-1)
            N1=N1+2
            IF( N1 .GT. NBDL ) N1 = 3
            CALL ALEAD( N1, NBDL, V )
            DO 50 L=1,NBDL
               VECP(NC,L) = V(L)
 50         CONTINUE
            VECP(NC,N1) = 0.1D0 * VECP(NC,N1)
C
            GOTO 1
         ENDIF
C
C        B-ORTHONORMALISATION DE V ET MISE DANS R
C        ----------------------------------------
ccc         print *,'orthod: vecteur',j,' tV [M] V=',s
         S = 1D0 / SQRT( ABS(S) )
         DO 60 L=1,NBDL
            R(J,L) = V(L) * S
 60      CONTINUE
 100  CONTINUE
C
      RETURN
      END
