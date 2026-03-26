      SUBROUTINE TRIENT( N , X )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRIER DE FACON CROISSANTE LES N ENTIERS DE LA FILE X
C ----- TRI VALABLE POUR N FAIBLE ET FILE PRESQUE TRIEE

C PARAMETRE D ENTREE :
C --------------------
C N      : NOMBRE DE ENTIERS A TRIER

C PARAMETRE MODIFIE :
C -------------------
C X      : ENTIERS A TRIER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS            AVRIL 1987
C ......................................................................
      INTEGER  X( N ) , Y

C     TRI DIT DES FICHES
C     ==================
      DO I=2,N
         DO J1=I-1,1,-1
            J = J1 + 1
            IF( X(J1) .GT. X(J) ) THEN
               Y       = X( J  )
               X( J  ) = X( J1 )
               X( J1 ) = Y
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
