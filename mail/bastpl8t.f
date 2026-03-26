      SUBROUTINE BASTPL8T( COSMAXPL, MXSOMM, NBSOMM, XYZSOM,
     %                     MXTRIA,   NBTRIA, NOTRIA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUT DU BARYCENTRE DU TRIANGLE DE PLUS PETIT RAPPORT des
C -----    SURFACES DES PLUS DE 8 TRIANGLES D'UN SOMMET

C ENTREES:
C --------
C COSMAXPL:COSINUS DE L'ANGLE DIEDRE ENTRE LES 2 PLANS DES TRIANGLES
C          AU DESSOUS DUQUEL LES 2 TRIANGLES SONT CONSIDERES NON-COPLANAIRES
C MXSOMM : NOMBRE MAXIMUM DE SOMMETS DANS LE TABLEAU XYZSOM
C NBSOMM : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : XYZ DES NBSOMM SOMMETS DE LA TRIANGULATION
C MXTRIA : NOMBRE MAXIMUM DE TRIANGLES DU TABLEAU NOTRIA
C NOTRIA : NUMERO DE 1 A NBSOMM DES SOMMETS DES TRIANGLES DE LA TRIANGULATION
C          1 2 3 LE NO XYZSOM DES 3 SOMMETS DU TRIANGLE
C          4 5 6 LE NO DU TRIANGLE ADJACENT PAR L'ARETE 1 2 3

C MODIFIES:
C ---------
C NBSOMM : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : XYZ DES NBSOMM SOMMETS + BARYCENTRES AJOUTES DE LA TRIANGULATION
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint Pierre du Perray          Fevrier 2020
C23456...............................................................012
      PARAMETER       ( MXTRST=256 )
C     MXTRST=MAXIMUM DE TRIANGLES D'UN MEME SOMMET
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      INTEGER  NOTRIA(6,MXTRIA), NOTRST(MXTRST)
      REAL     XYZSOM(3,NBSOMM), XYZBAR(3)

      TRACTE0 = TRACTE

C     NOMBRE D'AJOUTS DE BARYCENTRES
      NBAJBA  = 0

C     1-ERE ITERATION EN NUMEROS DE SOMMET DECROISSANTS
      NSMIN = NBSOMM
      NSMAX = 1
      LEPAS =-1

      DO 200 ITER = 1, 1
      NB2DIA  = 0
      NBSOMM0 = NBSOMM

      DO 100 NSt = NSMIN, NSMAX, LEPAS

C        NOMBRE DE TRIANGLES DE SOMMET NST
         NBTRST = 0

C        RECHERCHE DES NBTRST TRIANGLES DE SOMMET NST
         DO 10 NTR = 1, NBTRIA

            IF( NOTRIA( 1, NTR ) .LE. 0 ) GOTO 10

            DO I = 1, 3
C              NUMERO DU SOMMET I DU TRIANGLE NTR
               NS = NOTRIA( I, NTR )
               IF( NS .EQ. NST ) THEN
C                 NTR A POUR SOMMET NST
                  IF( NBTRST .LT. MXTRST ) THEN
                     NBTRST = NBTRST + 1
                     NOTRST( NBTRST ) = NTR
                  ENDIF
                  GOTO 10
               ENDIF
            ENDDO

C           PASSAGE AU TRIANGLE DE SOMMET NSt SUIVANT
 10      ENDDO

         IF( NBTRST .LE. 8 ) GOTO 100

C        IL Y A PLUS DE 8 TRIANGLES DE SOMMET NSt
C        ----------------------------------------
C        TRACE DES TRIANGLES NOTRST(1:NBTNEW) DE NOTRIA, LEURS VOISINS
C        et LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
         CALL TRTRIAN( 'bastpl8t 0', XYZSOM, 6, MXTRST, NBTRST, NOTRST,
     %                  NOTRIA )

         IF( NBSOMM   .GE. MXSOMM ) GOTO 9999
         IF( NBTRIA+2 .GT. MXTRIA ) GOTO 9999

         IF( ITER .EQ. 1 ) THEN

C           RECHERCHE DU TRIANGLE DE SOMMET NSt DE PLUS PETIT
C           RAPPORT des SURFACES
C           -------------------------------------------------
            NTRMIN = 0
            SMIN   = 1E27
            DO 20 N = 1, NBTRST

               NTR = NOTRST( N )
C              COMPARAISON DES SURFACES du TRIANGLE NTR et
C              du TRIANGLE ADJACENT d'ARETE SANS SOMMET NSt
               STR = SURTRR( XYZSOM( 1, NOTRIA(1,NTR) ),
     %                       XYZSOM( 1, NOTRIA(2,NTR) ),
     %                       XYZSOM( 1, NOTRIA(3,NTR) ) )
               IF( STR .LE. 0. ) GOTO 20

C              RECHERCHE DE L'ARETE DE NTR SANS SOMMET NSt
               DO I=1,3
                  NS1 = NOTRIA( I, NTR )
                  IF( NS1 .NE. NSt ) THEN
                     IF( I .EQ. 3 ) THEN
                        I1 = 1
                     ELSE
                        I1 = I+1
                     ENDIF
                     NS2 = NOTRIA( I1, NTR )
                     IF( NS2 .NE. NSt ) THEN

C                       LE TRIANGLE ADJACENT DE NTR PAR L'ARETE
C                       PERIPHERIQUE I
                        NTRADJ = NOTRIA( 3+I, NTR )
                        IF( NTRADJ .LE. 0 ) GOTO 20

C                       COMPARAISON DES SURFACES DES TRIANGLES NTR NTRADJ
                        S = SURTRR( XYZSOM( 1, NOTRIA(1,NTRADJ) ),
     %                              XYZSOM( 1, NOTRIA(2,NTRADJ) ),
     %                              XYZSOM( 1, NOTRIA(3,NTRADJ) ) )
                        IF( S .LE. 0. ) GOTO 20

C                       LE RAPPORT DES SURFACES ENTRE 0 et 1
                        IF( STR .LE. S ) THEN
                           RAP2TR = STR / S
                        ELSE
                           RAP2TR = S / STR
                        ENDIF

                        IF( RAP2TR .LT. SMIN ) THEN
                           SMIN   = RAP2TR
                           NTRMIN = NTR
                        ENDIF
                        GOTO 20
                     ENDIF
                  ENDIF
               ENDDO

 20         ENDDO

         ELSE

C           ITER>1: RECHERCHE DU TRIANGLE DE SOMMET NSt DE PLUS
C                   GRANDE SURFACE
C           ---------------------------------------------------
            NTRMIN = 0
            SMIN = 1E27
            DO N = 1, NBTRST
               NTR = NOTRST( N )
               S = SURTRR( XYZSOM( 1, NOTRIA(1,NTR) ),
     %                     XYZSOM( 1, NOTRIA(2,NTR) ),
     %                     XYZSOM( 1, NOTRIA(3,NTR) ) )
               IF( S .LT. SMIN ) THEN
                  SMIN   = S
                  NTRMIN = NTR
               ENDIF
            ENDDO

         ENDIF


C        NTRMIN LE TRIANGLE DE PLUS PETIT RAPPORT d'AIRES ou GRANDE AIRE
C        EST DECOUPE EN 3 SOUS-TRIANGLES A PARTIR DE SON BARYCENTRE
C        ---------------------------------------------------------------
         DO K=1,3
            XYZBAR( K ) = 0.
         ENDDO
         DO K=1,3
            NS = NOTRIA( K, NTRMIN )
            DO I=1,3
               XYZBAR(I) = XYZBAR(I) + XYZSOM(I,NS)
            ENDDO
         ENDDO

C        AJOUT DU BARYCENTRE DU TRIANGLE NTRMIN
         NBAJBA = NBAJBA + 1
         NBSOMM = NBSOMM + 1

         DO K=1,3
            XYZSOM( K, NBSOMM ) = XYZBAR( K ) / 3
         ENDDO

C        DECOUPAGE EN 3 SOUS-TRIANGLES DU TRIANGLE NTRMIN
         NTR1 = NTRMIN

C        LES 3 SOMMETS DE NTRMIN
         NS1 = NOTRIA( 1, NTRMIN )
         NS2 = NOTRIA( 2, NTRMIN )
         NS3 = NOTRIA( 3, NTRMIN )

C        LES 3 TRIANGLES ADJACENTS DE NTRMIN PAR LES ARETES
         NTR1AD = NOTRIA( 4, NTRMIN )
         NTR2AD = NOTRIA( 5, NTRMIN )
         NTR3AD = NOTRIA( 6, NTRMIN )

C        LES 2 TRIANGLES CREES DE SOMMET NBSOMM
         NTR2 = NBTRIA + 1
         NTR3 = NBTRIA + 2
         NBTRIA = NTR3

C        CONSTRUCTION DU TRIANGLE NTR1: NS1-NS2-NBSOMM
         NOTRIA( 1, NTR1 ) = NS1
         NOTRIA( 2, NTR1 ) = NS2
         NOTRIA( 3, NTR1 ) = NBSOMM

C        LES TRIANGLES ADJACENTS DE NTR1
         NOTRIA( 4, NTR1 ) = NTR1AD
         IF( NTR1AD .GT. 0 .AND. NOTRIA(1,NTR1AD) .GT. 0 ) THEN
            CALL NUAR2STR( NS1, NS2, NOTRIA(1,NTR1AD), NATA )
            NOTRIA( 3+NATA, NTR1AD ) = NTR1
         ENDIF
         NOTRIA( 5, NTR1 ) = NTR2
         NOTRIA( 6, NTR1 ) = NTR3

C        CONSTRUCTION DU TRIANGLE NTR2: NS2-NS3-NBSOMM
         NOTRIA( 1, NTR2 ) = NS2
         NOTRIA( 2, NTR2 ) = NS3
         NOTRIA( 3, NTR2 ) = NBSOMM

C        LES TRIANGLES ADJACENTS DE NTR2
         NOTRIA( 4, NTR2 ) = NTR2AD
         IF( NTR2AD .GT. 0 .AND. NOTRIA(1,NTR2AD) .GT. 0 ) THEN
            CALL NUAR2STR( NS2, NS3, NOTRIA(1,NTR2AD), NATA )
            NOTRIA( 3+NATA, NTR2AD ) = NTR2
         ENDIF
         NOTRIA( 5, NTR2 ) = NTR3
         NOTRIA( 6, NTR2 ) = NTR1

C        CONSTRUCTION DU TRIANGLE NTR3: NS3-NS1-NBSOMM
         NOTRIA( 1, NTR3 ) = NS3
         NOTRIA( 2, NTR3 ) = NS1
         NOTRIA( 3, NTR3 ) = NBSOMM

C        LES TRIANGLES ADJACENTS DE NTR3
         NOTRIA( 4, NTR3 ) = NTR3AD
         IF( NTR3AD .GT. 0 .AND. NOTRIA(1,NTR3AD) .GT. 0 ) THEN
            CALL NUAR2STR( NS3, NS1, NOTRIA(1,NTR3AD), NATA )
            NOTRIA( 3+NATA, NTR3AD ) = NTR3
         ENDIF
         NOTRIA( 5, NTR3 ) = NTR1
         NOTRIA( 6, NTR3 ) = NTR2

C        TENTATIVE D'ECHANGE DES DIAGONALES DES ARETES D'ADJACENCE
C        DES TRIANGLES DE SOMMET NSt + NBSOMM
         NBTRST = NBTRST + 1
         NOTRST( NBTRST ) = NTR2
         NBTRST = NBTRST + 1
         NOTRST( NBTRST ) = NTR3

         DO NT = NBTRST, 1, -1
            NTR = NOTRST( NT )
            DO K = 1, 3
               CALL EC2DIA(COSMAXPL, NTR, K, NOTRIA, NBSOM, XYZSOM, NT0)
               IF( NT0 .GT. 0 ) THEN
C                 L'ECHANGE A ETE FAIT
                  NB2DIA = NB2DIA + 1
               ENDIF
            ENDDO
         ENDDO

C        TRACE DES TRIANGLES NOTRST(1:NBTNEW) DE NOTRIA, LEURS VOISINS et
C        LEURS VOISINS DE VOISINS POUR VOIR LEUR ENVIRONNEMENT
ccc         tracte = .true.
         CALL TRTRIAN( 'bastpl8t 1', XYZSOM, 6, MXTRST, NBTRST, NOTRST,
     %                  NOTRIA )
         tracte = tracte0

 100  ENDDO

      PRINT*,'bastpl8t: Iteration',ITER,
     %       ': Nb de SOMMETS avec PLUS de 8 TRIANGLES',
     %        NBAJBA,' BARYCENTRES AJOUTES et',NB2DIA,
     %       ' ECHANGES de DIAGONALES 2Tr<->2Tr'

      IF( NBAJBA .EQ. 0 .AND. NB2DIA .EQ. 0 ) GOTO 9999

C     SECONDE ITERATION EN NUMEROS DE SOMMET CROISSANT ...
      I     = NSMIN
      NSMIN = NSMAX
      NSMAX = I
      LEPAS = -LEPAS

 200  ENDDO


 9999 TRACTE = TRACTE0
      RETURN
      END
