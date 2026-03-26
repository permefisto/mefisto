      SUBROUTINE MOBMAX(NCODSA,NDSM,NTDL,LPLIGN,LPCOLO,AG,X,Y)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                          SP MOBMAX
C                          ---------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CALCULER LE RESIDU Y = Y - A * X DU SYSTEME LINEAIRE
C ----
C      Y , X TABLEAUX (NDSM,NTDL)  A MATRICE MORSE SYMETRIQUE
C                                  OU DIAGONALE
C                                  OU NON SYMETRIQUE
C
C PARAMETRES D ENTREE:
C --------------------
C NCODSA : 0 SI LA MATRICE A EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE PROFIL
C         -1 SI LA MATRICE EST NON-SYMETRIQUE PROFIL
C NDSM   : NOMBRE DE CAS OU SECONDS MEMBRES
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE OU ORDRE DE LA MATRICE
C LPLIGN :
C LPCOLO : LES POINTEURS DE LA MATRICE MORSE AG
C AG     :
C X      : LE VECTEUR DONNE
C
C PARAMETRE MODIFIE  :
C --------------------
C Y      : LE VECTEUR RESULTAT
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR: P. JOLY  ANALYSE NUMERIQUE PARIS      MAI 1990
C.......................................................................
      DIMENSION         LPLIGN(NTDL+1),LPCOLO(1:*)
      DOUBLE PRECISION  AG(1:*),X(NTDL,NDSM),Y(NTDL,NDSM)
C
C     AIGUILLAGE SELON NCODSA
C     -----------------------
C
      IF ( NCODSA .GT. 0 ) THEN
C        ------   MATRICE SYMETRIQUE   ------
         J1 = 1
         DO I=1,NTDL
            J2=LPLIGN(I+1) - 1
            IF(J1.GT.J2) GOTO 3
            DO J=J1,J2
               NOCOL = LPCOLO(J)
               DO ND=1,NDSM
                  Y(I,ND) = Y(I,ND) - AG(J) * X(NOCOL,ND)
                  Y(NOCOL,ND) = Y(NOCOL,ND) - AG(J) * X(I,ND)
               ENDDO
            ENDDO
C           ---   LE COEFFICIENT DIAGONAL   ---
 3          J2 = J2 + 1
            DO ND=1,NDSM
               Y(I,ND) = Y(I,ND) - AG(J2) *  X(I,ND)
            ENDDO
            J1 = J2 + 1
         ENDDO
      ELSE
C        ------   MATRICE NON SYMETRIQUE   ------
         J1 = 1
         DO I=1,NTDL
            J2 = LPLIGN(I+1)
            DO J=J1,J2
               NOCOL = LPCOLO(J)
               DO ND=1,NDSM
                  Y(I,ND) = Y(I,ND) - AG(J) * X(NOCOL,ND)
               ENDDO
            ENDDO
            J1 = J2 + 1
         ENDDO
      END IF

      RETURN
      END
