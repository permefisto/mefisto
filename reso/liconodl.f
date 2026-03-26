      SUBROUTINE LICONODL( NBNOEU, NTDL,   NDDLNO, NCODSA,
     %                     LPLINO, LPCONO,
     %                     LPLIDL, LPCODL, NBCOAG, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     POUR TRAITER une MATRICE MORSE CONDENSEE
C -----     PASSAGE des POINTEURS sur les NOEUDS 
C                   aux POINTEURS SUR LES DEGRES DE LIBERTE:
C           A PARTIR DES NOEUDS le POINTEUR SUR LE COEFF DIAGONAL LPLINO
C           le NUMERO des COLONNES LPCONO de la MATRICE
C           GENERER LE POINTEUR SUR LES DEGRES DE LIBERTE 
C           du POINTEUR sur le COEFF DIAGONAL LPLIDL et
C           le NO des COLONNES LPCODL et CELA
C           a PARTIR DE NDDLNO POINTEUR SUR LE NO DU DERNIER DL
C           DE CHAQUE NOEUD

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
C LPLINO  : TABLEAU DES POINTEURS SUR LES COEFFICIENTS DE LA DIAGONALE
C           DE LA MATRICE MORSE EN SUPPOSANT QUE LES NOEUDS ONT 1 SEUL DL
C LPCONO  : TABLEAU DES NUMEROS DE COLONNES DES COEFFICIENTS DE LA MATRICE
C           PAR NOEUD

C SORTIES :
C ---------
C LPLIDL  : TABLEAU DES POINTEURS SUR LA DIAGONALE POUR TOUS LES DL
C           A PARTIR DU TABLEAU NDDLNO POINTEUR SUR LE DERNIER DL
C           DE CHAQUE NOEUD. CE TABLEAU EST DIMENSIONNE (0:NTDL)
C LPCODL  : TABLEAU DES NUMEROS DE COLONNES DES COEFFICIENTS DE LA MATRICE
C           PAR DL
C NBCOAG  : NOMBRE DE COEFFICIENTS COLONNE DE LA MATRICE = LPLIDL(NTDL)
C IERR    : =1 DEBORDEMENT DE 2**31 VARIABLES POUR LE PROFIL
C           =0 SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Saint Pierre du Perray            Fevrier 2022
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      INTEGER  NDDLNO(0:NBNOEU), LPLINO(0:NBNOEU), LPLIDL(0:NTDL)
      INTEGER  LPCONO(1:*), LPCODL(1:*)

      IERR  = 0

      IF( NCODSA .GT. 0 ) THEN

C        MATRICE SYMETRIQUE
         LPLIDL(0) = 0
         NOCODL = 0

         DO I = 1, NBNOEU

            DO NODLI = NDDLNO(I-1)+1, NDDLNO(I)

               DO NCONO = LPLINO(I-1)+1, LPLINO(I)

C                 NO DE LA COLONNE PAR NOEUD
                  NOCONO = LPCONO(NCONO)

C                 PARCOURS DES DL DE LA COLONNE NLCONO
                  DO NDL = NDDLNO(NOCONO-1)+1, MIN(NDDLNO(NOCONO),NODLI)
                     NOCODL = NOCODL + 1
                     LPCODL( NOCODL ) = NDL
                  ENDDO

               ENDDO

C              DERNIERE COLONNE DL DU DL NODLI    
               LPLIDL(NODLI) = NOCODL

C              DETECTION DE PASSAGE DES 2**31
               IF( LPLIDL(NODLI) .LT. LPLIDL(NODLI-1) ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT *,'liconodl: MATRICE CONDENSEE > 2**31 COEFFI
     %CIENTS'
                  ELSE
                     PRINT *,'liconodl: CONDENSED MATRIX > 2**31 COEFFIC
     %IENTS'
                  ENDIF
                  IERR = 1
                  GOTO 9999
               ENDIF

            ENDDO

         ENDDO

      ELSE

         IF( LANGAG .EQ. 0 ) THEN
            print *,'liconodl: MATRICE DIAGONALE ou NON SYMETRIQUE NON P
     %ROGRAMMEE'
         ELSE
            print *,'liconodl: NON SYMETRIQUE or DIAGONAL MATRIX NOT COM
     %PUTED'
         ENDIF
      ENDIF

      NBCOAG = LPLIDL( NTDL )
      IF( LANGAG .EQ. 0 ) THEN
         print *,'liconodl: MATRICE CONDENSEE de', NBCOAG,
     %           ' Reels Double Precision'
      ELSE
         print *,'liconodl: CONDENSED MATRIX of', NBCOAG,
     %           ' Real Double Precision'
      ENDIF

 9999 RETURN
      END
