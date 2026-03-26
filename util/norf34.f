      SUBROUTINE NORF34( NBSF, NOSOEL, XYZSOM, COORNO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES COMPOSANTES DU VECTEUR NORMAL UNITE A UNE FACE
C -----    TRIANGULAIRE OU QUADRANGULAIRE
C ENTREES:
C --------
C NBSF   : NOMBRE D ARETES DE LA FACE ( 3 OU 4 )
C NOSOEL : NUMERO DES NBSF SOMMETS DE LA FACE
C XYZSOM : 3 COORDONNEES DES NBSF SOMMETS

C SORTIES:
C --------
C COORNO : LES 3 COMPOSANTES DE LA NORMALE UNITE A LA FACE
C IERR   : =0 SI PAS D'ERREUR
C          =1 VECTEUR NORMAL NON CALCULABLE -> VECTEUR UNITE SELON X
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS             MARS 1987
C ......................................................................
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NOSOEL(1:NBSF)
      REAL              XYZSOM(1:3,1:*),COORNO(3),V(3,2)

      IERR = 0

C     SI LA FACE EST TRIANGULAIRE LA NORMALE EST LE PRODUIT VECTORIEL
C     DU COTE 1 PAR LE COTE 3
C     ---------------------------------------------------------------
      IF( NBSF .NE. 3 ) GOTO 10
      DO J=1,2
         DO I=1,3
            V(I,J) = XYZSOM(I,NOSOEL(J+1)) - XYZSOM(I,NOSOEL(1))
         ENDDO
      ENDDO
      GOTO 20

C     SI LA FACE EST QUADRANGULAIRE LA NORMALE EST LE PRODUIT VECTORIEL
C     DES 2 DIAGONALES
C     -----------------------------------------------------------------
 10   IF( NBSF .NE. 4 ) GOTO 9900
      DO J=1,2
         DO I=1,3
            V(I,J) = XYZSOM(I,NOSOEL(J+2)) - XYZSOM(I,NOSOEL(J))
         ENDDO
      ENDDO

C     LE PRODUIT VECTORIEL DES 2 VECTEURS ET LA NORMALISATION A 1.
C     ------------------------------------------------------------
 20   CALL PROVER( V(1,1), V(1,2), COORNO )
      R = PROSCR( COORNO, COORNO, 3 )
      IF( R .LE. 0. ) GOTO 9990

C     NORMALISATION A 1. DU VECTEUR NORMAL
      R = 1./ SQRT( R )
      DO I=1,3
         COORNO( I ) = COORNO( I ) * R
      ENDDO
      GOTO 9999


C     ERREUR : LA FACE N EST NI TRIANGULAIRE NI QUADRANGULAIRE
C     --------------------------------------------------------
 9900 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NBSF
      KERR(1) = 'norf34: FACE de '//KERR(MXLGER)(1:4)
     %           //' SOMMETS AU LIEU de 3 OU 4'
      CALL LEREUR
C     AFFICHAGE DES NBSF SOMMETS
      GOTO 9997


C     ERREUR: VECTEUR NORMAL NON CALCULABLE -> VECTEUR UNITE SELON X
C     -------------------------------------
 9990 NBLGRC(NRERR) = 1
      KERR(1) = 'norf34: PAS de VECTEUR NORMAL a UNE FACE DEGENEREE'
      CALL LEREUR
      PRINT 19000
19000 FORMAT(' norf34: Probleme la FACE est DEGENEREE')


 9997 DO J=1,NBSF
         NS = NOSOEL(J)
         IF( NS .GT. 0 ) THEN
            PRINT 19001, NS, (XYZSOM(I,NS),I=1,3)
         ELSE
            PRINT 19001, NS
         ENDIF
      ENDDO
19001 FORMAT(' norf34: SOMMET',I9,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)

C     LE VECTEUR NORMAL est IMPOSE -> VECTEUR UNITE SELON X
      COORNO(1)=1
      COORNO(2)=0
      COORNO(3)=0
      IERR = 1
      GOTO 9999

 9999 RETURN
      END
