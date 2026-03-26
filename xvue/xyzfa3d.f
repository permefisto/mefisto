      SUBROUTINE XYZFA3D( NBS, NOSOEL, XYZSOM, XYZST )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE TABLEAU XYZST DES XYZ DES NBS SOMMETS
C -----    POUR APPELER ENSUITE FACE3D

C ENTREES:
C --------
C NBS    : NOMBRE DE SOMMETS DE L'EF
C NOSOEL : NUMERO XYZSOM DES NBS SOMMETS DE L'EF
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE

C SORTIE :
C --------
C XYZST  : XYZ DES NBS SOMMETS DE L'EF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN  Saint PIERRE du PERRAY               Mai 2020
C2345X...............................................................012
      INTEGER   NOSOEL(NBS)
      REAL      XYZSOM(3,*), XYZST(3,NBS)

      DO N = 1, NBS
         NS = NOSOEL( N )
         DO K = 1, 3
            XYZST( K, N ) = XYZSOM( K, NS )
         ENDDO
      ENDDO

      RETURN
      END
