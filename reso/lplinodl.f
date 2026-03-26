      SUBROUTINE LPLINODL( NBNOEU,  NTDL,   NDDLNO, NCODSA, LPLIGNO,
     %                     LPLIGDL, NBCOAG, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     A PARTIR DU POINTEUR SUR LA DIAGONALE DU PROFIL AUX NOEUDS
C -----     GENERER LE POINTEUR SUR LA DIAGONALE DES DEGRES DE LIBERTE
C           DE CHAQUE NOEUD DEFINI PAR NDDLNO POINTEUR SUR LE NO DU
C           DERNIER DL DE CHAQUE NOEUD

C ENTREES :
C ---------
C NBNOEU  : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C NTDL    : NOMBRE TOTAL DE DL DU MAILLAGE
C NDDLNO  : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD
C           CE TABLEAU EST DIMENSIONNE (0:NBNOEU)
C NCODSA  : CODE DE STOCKAGE D'UNE MATRICE GLOBALE
C           -1 MATRICE NON SYMETRIQUE MAIS POSITION SYMETRIQUE
C              DES COEFFICIENTS NON NULS
C            0 MATRICE DIAGONALE
C            1 MATRICE SYMETRIQUE
C LPLIGNO : TABLEAU DES POINTEURS SUR LES COEFFICIENTS DE LA DIAGONALE
C           DE LA MATRICE PROFIL EN SUPPOSANT QUE LES NOEUDS ONT 1 SEUL DL

C SORTIES :
C ---------
C LPLIGDL : TABLEAU DES POINTEURS SUR LA DIAGONALE POUR TOUS LES DL
C           A PARTIR DU TABLEAU NDDLNO POINTEUR SUR LE DERNIER DL
C           DE CHAQUE NOEUD. CE TABLEAU EST DIMENSIONNE (0:NTDL)
C NBCOAG  : NOMBRE DE COEFFICIENTS DE LA MATRICE = LPLIGDL(NTDL)
C IERR    : =1 DEBORDEMENT DE 2**31 VARIABLES POUR LE PROFIL
C           =0 SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: SOFIANE BENHAMADOUCHE ANALYSE NUMERIQUE UPMC PARIS Janvier2000
C MODIFS: ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Avril 2010
C MODIFS: ALAIN PERRONNET Saint Pierre du Perray            Fevrier 2022
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      INTEGER  NDDLNO(0:NBNOEU), LPLIGNO(0:NBNOEU), LPLIGDL(0:NTDL)

      IERR   = 0
      NODLIB = 0
      LPLIGDL(0) = 0

      IF( NCODSA .GT. 0 ) THEN

C        MATRICE SYMETRIQUE
         DO I = 1, NBNOEU

            NBTENN = 0
            DO J = I-(LPLIGNO(I)-LPLIGNO(I-1)-1), I-1
C              NOMBRE DE DL DU NOEUD J
               NDLJ = NDDLNO(J) - NDDLNO(J-1)
C              NOMBRE DE TERMES NON NULS DU PROFIL DE LA MATRICE
               NBTENN = NBTENN + NDLJ
            ENDDO

C           NOMBRE DE DL DU NOEUD I
            NDLI = NDDLNO(I) - NDDLNO(I-1)
            DO K = 1,NDLI
C              LIGNE DU DEGRE DE LIBERTE NODLIB DE LA MATRICE PAR DL
               NODLIB = NODLIB + 1
               LPLIGDL(NODLIB) = LPLIGDL(NODLIB-1) + NBTENN + K
            ENDDO

         ENDDO

      ELSE

C        MATRICE NON SYMETRIQUE
         DO I = 1, NBNOEU

            NBTENN = 0
            IH     = ( LPLIGNO(I) - LPLIGNO(I-1) ) / 2
            DO J = I-IH, I-1
C              NOMBRE DE DL DU NOEUD J NON DIAGONAL
               NDLJ = NDDLNO(J) - NDDLNO(J-1)
C              NOMBRE DE TERMES NON NULS DU PROFIL DE LA MATRICE
               NBTENN = NBTENN + NDLJ * 2
            ENDDO

C           NOMBRE DE DL DU NOEUD I
            NDLI = NDDLNO(I) - NDDLNO(I-1)
            DO K = 1,NDLI

C              LIGNE DU DEGRE DE LIBERTE NODLIB DE LA MATRICE PAR DL
               NODLIB = NODLIB + 1
               LPLIGDL(NODLIB) = LPLIGDL(NODLIB-1) + NBTENN + 2*K-1

C              DETECTION DE PASSAGE DES 2**31 DEPASSEMENT MEMOIRE
               IF( LPLIGDL(NODLIB) .LT. LPLIGDL(NODLIB-1) ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT *,'lplinodl: PROFIL de la MATRICE > 2**31 COE
     %FFICIENTS => DEPASSEMENT MEMOIRE'
                  ELSE
                     PRINT *,'lplinodl: SKYLINE MATRIX > 2**31 COEFFICIE
     %NTS => SEGMENTATION FAULT'
                  ENDIF
                  IERR = 1
                  GOTO 9999
               ENDIF

            ENDDO
         ENDDO

      ENDIF

      NBCOAG = LPLIGDL( NTDL )
      IF( LANGAG .EQ. 0 ) THEN
         print *,'lplinodl: MATRICE PROFIL de', NBCOAG,
     %           ' Reels Double Precision'
      ELSE
         print *,'lplinodl: SKYLINE MATRIX of', NBCOAG,
     %           ' Real Double Precision'
      ENDIF

 9999 RETURN
      END
