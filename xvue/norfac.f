      SUBROUTINE NORFAC( NAR, COORSO, COORNO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES COMPOSANTES DU VECTEUR NORMAL UNITE A UNE FACE
C -----   TRIANGULAIRE OU QUADRANGULAIRE
C         VERSION AVEC DES REELS EN SIMPLE PRECISION

C ENTREES:
C --------
C NAR    : NOMBRE D ARETES DE LA FACE ( 3 OU 4 )
C COORSO : COORDONNEES DES SOMMETS DE LA FACE ( COORSO(3,3 OU 4) )

C SORTIES:
C --------
C COORNO : LES 3 COMPOSANTES DE LA NORMALE A LA FACE DE NORME 1
C          VECTEUR NUL EN CAS D'ERREUR ( FACE REDUITE A UNE ARETE )
C IERR   : =0 SI NORMALE CORRECTEMENT CALCULEE
C          =1 SI FACE NON TRIANGULAIRE et NON QUADRANGULAIRE
C          =2 SI FACE REDUITE A UNE ARETE => COORNO=0,0,0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS             MARS 1987
C ......................................................................
      REAL              COORSO(3,NAR), COORNO(3)
      DOUBLE PRECISION  V(3,2), VECNOR(3), D, PROSCD

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
      PRINT*,'norfac: FACE de ',NAR,' ARETES. Au LIEU de 3 ou 4'
      IERR = 1
      GOTO 9900

C     LE PRODUIT VECTORIEL ET LA NORMALISATION A 1.
C     ---------------------------------------------
 100  CALL PROVEC( V(1,1), V(1,2), VECNOR )
      D = PROSCD( VECNOR, VECNOR, 3 )
C
      IF( D .EQ. 0D0 ) THEN
         PRINT 10100, (J,(COORSO(I,J),I=1,3),J=1,NAR)
10100 FORMAT(' norfac: Probleme: une FACE REDUITE A UNE ARETE => PAS de 
     %VECTEUR NORMAL'/
     %   ( ' norfac: Sommet ',I1,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7))
         IERR = 2
         GOTO 9900
      ENDIF

C     LES 3 COMPOSANTES DE LA NORMALE UNITAIRE
      D = SQRT( D )
      DO I=1,3
         COORNO( I ) = REAL( VECNOR( I ) / D )
      ENDDO
      GOTO 9999


C     LES 3 COMPOSANTES NULLES DE LA NORMALE INCORRECTE
 9900 COORNO(1) = 0
      COORNO(2) = 0
      COORNO(3) = 0

 9999 RETURN
      END
