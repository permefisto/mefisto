      SUBROUTINE XYZSTG( MNS, MNTG, NBS, NOSOEL, NUEFTG,
     %                   X, Y, Z )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER OU CALCULER LES COORDONNEES DES NBS SOMMETS
C -----    ET DES 2 NBS TANGENTES DE L'EF AVEC TG
C
C ENTREES:
C --------
C MNS    : ADRESSE MCN -3 DE LA PREMIERE COORDONNEE DU PREMIER SOMMET
C MNTG   : ADRESSE MCN -3 DE LA PREMIERE COORDONNEE DE LA PREMIERE TG
C NBS    : NOMBRE DE SOMMETS ( 3 ou 4 ) DE LA FACE A TRACER
C NOSOEL : NUMERO DES NBSOEF SOMMETS
C          SUIVI EVENTUELLEMENT DES NBTGEF TANGENTES DE L'EF NUELEM
C          (=> 12 AU PLUS)
C          SI PAS DE TANGENTES (NBTGEF=0):
C             TRIANGLE  : NO SOMMET1 , NS2 , NS3, 0
C             QUADRANGLE: NO SOMMET1 , NS2 , NS3, NS4
C
C          S'IL EXISTE DES TANGENTES (NBTGEF>0) ELLES SONT RANGEES PAR SOMMET:
C             TRIANGLE  : NO  SOMMET1, NS2 , NS3, 0,
C                         NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                         NO TANGENTE5(S3S1), NT6(S3S2),   0,         0
C             QUADRANGLE: NO  SOMMET1, NS2 , NS3, NS4,
C                         NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                         NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C          CE CHOIX PERMET UNE BOUCLE SUR LES TANGENTES PAR LES SOMMETS
C
C NUEFTG : NUMERO DE L'EF A TG ET SINON 0
C
C SORTIES:
C --------
C X      : ABSCISSE OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C Y      : ORDONNEE OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C Z      : COTE     OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1996
C2345X...............................................................012
      INTEGER        NOSOEL(1:12)
      REAL           X(12), Y(12), Z(12)
      include"./incl/pp.inc"
      COMMON         RMCN(MOTMCN)
C
C     LES 3 COORDONNEES DES NBS SOMMETS DE LA FACE
      DO 10 J=1,NBS
         MN   = MNS + 3 * NOSOEL(J)
         X(J) = RMCN( MN     )
         Y(J) = RMCN( MN + 1 )
         Z(J) = RMCN( MN + 2 )
 10   CONTINUE
C
C     LA FACE EST ELLE UN EF AVEC DES TANGENTES?
      IF( NUEFTG .GT. 0 ) THEN
C
C        OUI : LES NUMEROS DES 8 TGS DU QUADRANGLE
C              SONT RANGEES DANS NOSOEL(5:12)
         NT = 4
         NJ = NBS
         DO 50 J=1,NBS
C
C           LE NUMERO DE LA PREMIERE TANGENTE DU SOMMET J
            NT  = NT + 1
            NTG = NOSOEL( NT )
            NJ  = NJ + 1
            IF( NTG .NE. 0 ) THEN
C              IL EXISTE UNE TANGENTE
               MN = MNTG + 3 * ABS(NTG)
               IF( NTG .GT. 0 ) THEN
                  X(NJ) = RMCN( MN     )
                  Y(NJ) = RMCN( MN + 1 )
                  Z(NJ) = RMCN( MN + 2 )
               ELSE
                  X(NJ) = -RMCN( MN     )
                  Y(NJ) = -RMCN( MN + 1 )
                  Z(NJ) = -RMCN( MN + 2 )
               ENDIF
            ELSE
C              PAS DE TANGENTE: COTE DROIT
               IF( J .LT. NBS ) THEN
                  J1 = J + 1
               ELSE
                  J1 = 1
               ENDIF
               X(NJ) = X(J1) - X(J)
               Y(NJ) = Y(J1) - Y(J)
               Z(NJ) = Z(J1) - Z(J)
            ENDIF
C
C           LE NUMERO DE LA SECONDE TANGENTE DU SOMMET J
            NT  = NT + 1
            NTG = NOSOEL( NT )
            NJ  = NJ + 1
            IF( NTG .NE. 0 ) THEN
C              IL EXISTE UNE TANGENTE
               MN = MNTG + 3 * ABS(NTG)
               IF( NTG .GT. 0 ) THEN
                  X(NJ) = RMCN( MN     )
                  Y(NJ) = RMCN( MN + 1 )
                  Z(NJ) = RMCN( MN + 2 )
               ELSE
                  X(NJ) = -RMCN( MN     )
                  Y(NJ) = -RMCN( MN + 1 )
                  Z(NJ) = -RMCN( MN + 2 )
               ENDIF
            ELSE
C              PAS DE TANGENTE: COTE DROIT
               IF( J .NE. 1 ) THEN
                  J1 = J - 1
               ELSE
                  J1 = NBS
               ENDIF
               X(NJ) = X(J1) - X(J)
               Y(NJ) = Y(J1) - Y(J)
               Z(NJ) = Z(J1) - Z(J)
            ENDIF
 50      CONTINUE
      ENDIF
      END
