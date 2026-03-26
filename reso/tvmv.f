      SUBROUTINE TVMV( NBC,  NCODSB, MUB, B,
     %                 NBVD, NBV, VECTEUR, VAUX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER tVECT [B] VECT DES NBV VECTEURS
C -----    LES RESULTATS DOIVENT ETRE EGAUX A 1
C          CONDITION DE B-ORTHONORMALISATION
C
C ENTREES:
C --------
C NBC    : NOMBRE DE LIGNES ET COLONNES DE LA MATRICE B
C NCODSB : CODE 0 => B DIAGONALE , 1=> B PROFIL SYMETRIQUE
C MUB    : POINTEUR DU PROFIL DE B
C B      : LA MATRICE B
C NBVD   : NOMBRE DE VECTEURS DECLARES (NBVD,NBC)
C NBV    : NOMBRE DES PREMIERS VECTEURS A TRAITER
C VECTEUR: VECTEURS(NBVD,NBC)
C VAUX   : 2 VECTEURS AUXILIAIRES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2009
C ......................................................................
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           MUB(0:NBC)
      DOUBLE PRECISION  B(*), VECTEUR(NBVD,NBC), VAUX(NBC,2), PROSCD
C
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*)
     %'VERIFICATION: tV [M] V = 1? ------------------------------'
      DO 10 NV=1,NBV
C
C        COPIE DE VECTEUR(NV) DANS VAUX(1)
         DO 5 I=1,NBC
            VAUX(I,1) = VECTEUR(NV,I)
 5       CONTINUE
C
C        VAUX(2) = B * VECTEUR(NV)
         CALL MAPRVE( 0, 1D0, NBC, NCODSB, MUB, B, VAUX(1,1),
     %                VAUX(1,2) )
C
C       (VECTEUR(NV), B VECTEUR(NV))
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'VECTEUR PROPRE',NV,' tV [M] V=',
     %                       PROSCD( VAUX(1,1), VAUX(1,2), NBC )
         ELSE
            WRITE(IMPRIM,*) 'EIGENVECTOR',NV,' tV [M] V=',
     %                       PROSCD( VAUX(1,1), VAUX(1,2), NBC )
         ENDIF
 10   CONTINUE
C
      RETURN
      END
