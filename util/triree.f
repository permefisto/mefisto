      SUBROUTINE TRIREE( N , X )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRIER DE FACON CROISSANTE LES N REELS DE LA FILE X
C ----- TRI VALABLE POUR N FAIBLE ET FILE PRESQUE TRIEE
C
C PARAMETRE D ENTREE :
C --------------------
C N      : NOMBRE DE REELS A TRIER
C
C PARAMETRE MODIFIE :
C -------------------
C X      : REELS A TRIER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  AVRIL 1987
C ......................................................................
      REAL  X( N ) , Y
C
C     TRI DIT DES FICHES
C     ==================
      DO 20 I=2,N
         DO 10 J1=I-1,1,-1
            J = J1 + 1
            IF( X(J1) .GT. X(J) ) THEN
               Y       = X( J  )
               X( J  ) = X( J1 )
               X( J1 ) = Y
            ENDIF
   10    CONTINUE
   20 CONTINUE
      END
