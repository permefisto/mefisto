      SUBROUTINE NO1F1T( NOSOTR, NOSOTE, NUFACTET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NO LOCAL NUFACTET DE LA FACE DE SOMMETS NOSOTR
C -----    PARMI LES 4 SOMMETS NOSOTE D'UN TETRAEDRE
C
C ENTREES:
C --------
C NOSOTR : NUMERO DES 3 SOMMETS SUPPOSES CROISSANTS DE LA FACE(TRIANGLE)
C                               ===================
C NOSOTE : NUMERO DES 4 SOMMETS DU TETRAEDRE
C
C SORTIE :
C --------
C NUFACTET: >0 NUMERO DE 1 A 4 DE LA FACE DANS LE TETRAEDRE SI RETROUVEE
C              FACE 1: 1,3,2   2: 2,3,4   3: 3,1,4   4: 4,1,2
C           =0 SI FACE NON RETROUVEE DANS LE TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1991
C....................................................................012
      INTEGER      NOSOTR(1:3), NOSOTE(1:4), NS(1:3),
     %             NOSOFATE(1:3,1:4)
C                  NUMERO LOCAL DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      DATA         NOSOFATE/ 1,3,2, 2,3,4, 3,1,4, 4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

C     PARCOURS DES 4 FACES DU TETRAEDRE
      DO NUFACTET=1,4
C
C        LES 3 SOMMETS DE LA FACE NUFACTET DU TETRAEDRE
         NS(1) = NOSOTE( NOSOFATE(1,NUFACTET) )
         NS(2) = NOSOTE( NOSOFATE(2,NUFACTET) )
         NS(3) = NOSOTE( NOSOFATE(3,NUFACTET) )
C
C        TRI CROISSANT DES 3 NUMEROS DE SOMMETS NS
C        SI   NS(1) > NS(2) ILS SONT PERMUTES
         N = NS(1)
         IF( N .GT. NS(2) )  THEN
            NS(1) = NS(2)
            NS(2) = N
         ENDIF
C        ICI NS(1) < NS(2)
C
C        SI  NS(2) > NS(3) ILS SONT PERMUTES
         N = NS(2)
         IF( N .GT. NS(3) ) THEN
            NS(2) = NS(3)
            NS(3) = N
C           ICI NS(2) < NS(3)
C
C           SI  NS(1) > NS(2) ILS SONT PERMUTES
            N = NS(1)
            IF( N .GT. NS(2) )  THEN
               NS(1) = NS(2)
               NS(2) = N
            ENDIF
         ENDIF
C
C        EST ELLE LA BONNE FACE ?
         IF( NS(1) .EQ. NOSOTR(1) ) THEN
            IF( NS(2) .EQ. NOSOTR(2) ) THEN
               IF( NS(3) .EQ. NOSOTR(3) ) GOTO 9999
            ENDIF
         ENDIF
      ENDDO
C
C     FACE NON RETROUVEE
      NUFACTET = 0

 9999 RETURN
      END
