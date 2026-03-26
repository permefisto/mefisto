      SUBROUTINE AFNOARCF( NBCF, N1ARCF, NOARCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES ARETES DES CONTOURS FERMES DU TABLEAU NOARCF
C -----
C ENTREES:
C --------
C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C N1ARCF : NUMERO DU 1-ER SOMMET OU 1-ERE ARETE DU CONTOUR FERME
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY          Fevrier 2018
C2345X7..............................................................012
      INTEGER  N1ARCF(0:NBCF), NOARCF(3,*)

      DO 20 NCF = 1, NBCF

         PRINT*,'afnoarcf: NO du CF=',NCF

C        LA PREMIERE ARETE DU CF NCF
         NOA  = 0
         NA00 = N1ARCF( NCF )
         NA0  = NA00
         NA1  = NA0

C        L'ARETE NA0-NA1
 10      NA1 = NOARCF( 2, NA0 )

C        SES 2 SOMMETS
         NS0 = NOARCF( 1, NA0 )
         NS1 = NOARCF( 1, NA1 )

C        NUMERO DE L'ARETE DU CF NCF
         NOA = NOA + 1
         PRINT*,'afnoarcf: ARETE(',NOA,')=',NA0,' Sommets:',NS0, NS1,
     %          ' Suivante:',NA1,' Triangle OPPOSE', NOARCF(3,NA0)
         IF( NS0 .EQ. NS1 ) THEN
            PRINT*,'afnoarcf: ATTENTION SOMMETS INCORRECTS. A CORRIGER'
            GOTO 9999
         ENDIF

C        PASSAGE A L'ARETE SUIVANTE DU CF NCF
         IF( NA1 .NE. NA00 ) THEN
            NA0 = NA1
            GOTO 10
         ENDIF

 20   ENDDO

 9999 RETURN
      END
