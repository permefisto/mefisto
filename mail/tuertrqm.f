      SUBROUTINE TUERTRQM( XYZSOM, QUATRQM, NBTRIA, NSTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER de la TRIANGULATION TOUT TRIANGLE de
C -----    QUALITE < QUATRQM

C ENTREES:
C --------
C XYZSOM : XYZ DES SOMMETS DE LA TRIANGULATION
C QUATRQM: QUALITE AU DESSOUS DE LAQUELLE LES TRIANGLES SONT SUPPRIMES

C MODIFIES:
C ---------
C NBTRIA : NOMBRE DE TRIANGLES DU TABLEAU NSTRIA AVANT PUIS APRES
C NSTRIA : NUMERO DES 4 SOMMETS DE CHACUN DES NBTRIA FACES
C          UNE FACE DONT TOUTES LES ARETES N'ONT PAS DE FACE OPPOSEE
C          A TOUS SES NUMEROS DE SOMMETS REMIS A ZERO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY               Mai 2020
C2345X7..............................................................012
      REAL       XYZSOM(3,*)
      INTEGER    NSTRIA(4,NBTRIA)

      NBTR = 0
      DO 10 NTR = 1, NBTRIA

         IF( NSTRIA(1,NTR) .EQ. 0 ) GOTO 10

C        CALCUL DE LA QUALITE DU TRIANGLE NTR
         CALL QUATRI( NSTRIA(1,NTR), XYZSOM, QUALTR )

         IF( QUALTR .GE. QUATRQM ) THEN

C           UN TRIANGLE CONSERVE
            NBTR = NBTR + 1
            DO K=1,4
               NSTRIA( K, NBTR ) = NSTRIA( K, NTR )
            ENDDO

         ENDIF

 10   ENDDO

      PRINT*,'tuertrqm:',NBTRIA,' TRIANGLES INITIAUX'
      PRINT*,'tuertrqm:',NBTR,  ' TRIANGLES FINAUX, soit',
     %        NBTRSUP,' TRIANGLES SUPPRIMES de QUALITE <',QUATRQM

      NBTRSUP = NBTRIA - NBTR
      NBTRIA = NBTR

      RETURN
      END
