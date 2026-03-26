      SUBROUTINE NORMTQ( NBSF, XYZSTQ, COORNO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES COMPOSANTES DU VECTEUR NORMAL A UNE FACE
C -----    TRIANGULAIRE OU QUADRANGULAIRE
C
C ENTREES:
C --------
C NBSF   : NOMBRE D ARETES DE LA FACE ( 3 OU 4 )
C XYZSTQ : COORDONNEES DES 3 OU 4 SOMMETS
C
C SORTIES:
C --------
C COORNO : LES 3 COMPOSANTES DU VECTEUR NORMALE A LA FACE DE NORME 1
C IERR   : =0 SI PAS D'ERREUR
C          >0 EN CAS D'ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS          OCTOBRE 1994
C ......................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      REAL              XYZSTQ(1:3,1:NBSF),
     %                  COORNO(3),
     %                  V(3,2)
C
      IERR   = 0
      NUMPAS = 0
 5    D      = 0
C
      IF( NBSF .EQ. 3 ) THEN
C
C        SI LA FACE EST TRIANGULAIRE LA NORMALE EST LE PRODUIT VECTORIEL
C        DU COTE 1 PAR LE COTE 3
C        ---------------------------------------------------------------
         DO 20 J=1,2
            DO 10 I=1,3
               V(I,J) = XYZSTQ(I,J+1) - XYZSTQ(I,1)
               D      = D + V(I,J) ** 2
 10         CONTINUE
 20      CONTINUE
C
      ELSE IF( NBSF .EQ. 4 ) THEN
C
C        SI LA FACE EST QUADRANGULAIRE LA NORMALE EST LE PRODUIT VECTORIEL
C        DES 2 DIAGONALES 13 24
C        -----------------------------------------------------------------
         DO 50 J=1,2
            DO 40 I=1,3
               V(I,J) = XYZSTQ(I,J+2) - XYZSTQ(I,J)
               D      = D + V(I,J) ** 2
 40         CONTINUE
 50      CONTINUE
C
      ELSE
C
C        ERREUR : LA FACE N EST NI TRIANGULAIRE NI QUADRANGULAIRE
C        --------------------------------------------------------
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBSF
         KERR(1) = 'NORMTQ:FACE DE '//KERR(MXLGER)(1:4)
     %          //' COTES AU LIEU DE 3 OU 4'
         CALL LEREUR
C
C        NORMALE IMPOSEE AU VECTEUR X
         COORNO(1)=1
         COORNO(2)=0
         COORNO(3)=0
         IERR = 1
         GOTO 9999
      ENDIF
C
C     LE PRODUIT VECTORIEL ET LA NORMALISATION A 1.
C     ---------------------------------------------
      CALL PROVER( V(1,1) , V(1,2) , COORNO )
      R = PROSCR( COORNO , COORNO , 3 )
C
      IF( R .LT. 1E-10*D ) THEN
         IF( NBSF .EQ. 4 ) THEN
            IF( NUMPAS .EQ. 0 ) THEN
C              QUADRANGLE A 2 ``DIAGONALES'' PARALLELES
C              ESSAI AVEC LA PREMIERE AUTRE POSSIBILITE 12 34
               DO 111 I=1,3
                  R = XYZSTQ(I,2)
                  XYZSTQ(I,2) = XYZSTQ(I,3)
                  XYZSTQ(I,3) = R
 111           CONTINUE
               NUMPAS = 1
               GOTO 5
            ELSE IF( NUMPAS .EQ. 1 ) THEN
C              QUADRANGLE A 2 ``DIAGONALES'' PARALLELES
C              ESSAI AVEC LA SECONDE AUTRE POSSIBILITE 14 32
               DO 112 I=1,3
                  R = XYZSTQ(I,3)
                  XYZSTQ(I,3) = XYZSTQ(I,4)
                  XYZSTQ(I,4) = R
 112           CONTINUE
               NUMPAS = 2
               GOTO 5
            ENDIF
         ENDIF
C
C        NORMALE NULLE IMPOSEE AU VECTEUR X
         COORNO(1)=1
         COORNO(2)=0
         COORNO(3)=0
         IERR = 2
CCC
CCC         NBLGRC(NRERR) = 1
CCC         KERR(1) = 'PAS DE NORMALE A UNE FACE REDUITE A UNE ARETE'
CCC         CALL LEREUR
CCC         WRITE(IMPRIM,10100) (J,(XYZSTQ(I,NOSOEL(J)),I=1,3),J=1,NBSF)
CCC10100 FORMAT(' PB NORMTQ:FACE REDUITE A UNE ARETE'/
CCC     %   ( ' SOMMET',I1,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7/))
CCC
         GOTO 9999
      ENDIF
C
C     NORMALISATION A 1
      R = 1.0/ SQRT( R )
      DO 120 I=1,3
         COORNO( I ) = COORNO( I ) * R
 120  CONTINUE
C
 9999 RETURN
      END
