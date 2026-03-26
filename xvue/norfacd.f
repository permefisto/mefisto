      SUBROUTINE NORFACD( NAR, COORSO, COORNO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES COMPOSANTES DU VECTEUR NORMAL UNITE A UNE FACE
C -----   TRIANGULAIRE OU QUADRANGULAIRE
C         VERSION AVEC DES REELS EN DOUBLE PRECISION

C ENTREES:
C --------
C NAR    : NOMBRE D ARETES DE LA FACE ( 3 OU 4 )
C COORSO : COORDONNEES DES SOMMETS DE LA FACE ( COORSO(3,3 OU 4) )

C SORTIES:
C --------
C COORNO : LES 3 COMPOSANTES DE LA NORMALE A LA FACE
C          VECTEUR NUL EN CAS D'ERREUR ( FACE REDUITE A UNE ARETE )
C IERR   : =0 SI NORMALE CORRECTEMENT CALCULEE
C          =1 SI FACE NON TRIANGULAIRE et NON QUADRANGULAIRE
C          =2 SI FACE REDUITE A UNE ARETE => COORNO=0,0,0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS             MARS 1987
C ......................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  COORSO(3,NAR), COORNO(3), V(3,2), D, PROSCD
      INTRINSIC         SQRT

C     SI LA FACE EST TRIANGULAIRE LA NORMALE EST LE PRODUIT VECTORIEL
C     DU COTE 1 PAR LE COTE 3
C     ---------------------------------------------------------------
      IERR = 0
      IF( NAR .EQ. 3 ) THEN
         DO J=1,2
            DO I=1,3
               V(I,J) = COORSO(I,J+1) - COORSO(I,1)
            ENDDO
         ENDDO
         GOTO 100
      ENDIF

C     SI LA FACE EST QUADRANGULAIRE LA NORMALE EST LE PRODUIT VECTORIEL
C     DES 2 DIAGONALES
C     -----------------------------------------------------------------
      IF( NAR .EQ. 4 ) THEN
         DO J=1,2
            DO I=1,3
               V(I,J) = COORSO(I,J+2) - COORSO(I,J)
            ENDDO
         ENDDO
         GOTO 100
      ENDIF

C     ERREUR : LA FACE N EST NI TRIANGULAIRE NI QUADRANGULAIRE
C     --------------------------------------------------------
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NAR
      KERR(1) = 'norfacd: FACE de '//KERR(MXLGER)(1:4)
     %           //' COTES. AU LIEU de 3 ou 4'
      CALL LEREUR
      IERR = 1
      GOTO 9000

C     LE PRODUIT VECTORIEL ET LA NORMALISATION A 1.
C     ---------------------------------------------
 100  CALL PROVEC( V(1,1), V(1,2), COORNO )
      D = PROSCD( COORNO, COORNO, 3 )
C
      IF( D .EQ. 0D0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FACE REDUITE A UNE ARETE => PAS DE NORMALE'
         CALL LEREUR
         WRITE(IMPRIM,10100) (J,(COORSO(I,J),I=1,3),J=1,NAR)
10100 FORMAT(' norfacd: Probleme une FACE est REDUITE a UNE ARETE'/
     %   ( ' SOMMET',I1,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7))
         IERR = 2
         GOTO 9000
      ENDIF

C     LES 3 COMPOSANTES DE LA NORMALE UNITAIRE
      D = 1D0 / SQRT( D )
      DO I=1,3
         COORNO( I ) = COORNO( I ) * D
      ENDDO
      GOTO 9999


C     LES 3 COMPOSANTES NULLES DE LA NORMALE INCORRECTE
 9000 COORNO(1) = 0D0
      COORNO(2) = 0D0
      COORNO(3) = 0D0

 9999 RETURN
      END
