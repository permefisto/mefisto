      SUBROUTINE PRGCLC( NBNOE,  NBDLNO,  NCODSA, LPTDVOI, LISTVOI,
     %                   LPLIGN, MOLPCOI, LPCOLO, MOLPCOF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE STOCKAGE CONDENSE DE LA MATRICE DU MAILLAGE POUR
C -----    NBDLNO DEGRES DE LIBERTE ( CONSTANT EN TOUT NOEUD ) PAR NOEUD

C ENTREES:
C --------
C NBNOE  : NOMBRE DE NOEUDS DU MAILLAGE
C NBDLNO : NOMBRE CONSTANT DE DEGRES DE LIBERTE PAR NOEUD
C NCODSA : CODE DE STOCKAGE DE LA MATRICE PROFIL
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C LPTDVOI: NUMERO DANS LISTVOI DU DERNIER VOISIN DE CHAQUE NOEUD
C LISTVOI: NO DES VOISINS DE CHAQUE NOEUD TRIES CROISSANTS
C MOLPCOI: NOMBRE INITIAL D'ENTIERS DU TABLEAU LPCOLO

C SORTIES:
C --------
C LPLIGN : POINTEUR SUR LE DERNIER COEFFICIENT STOCKE DE CHAQUE LIGNE
C          DE LA MATRICE MORSE GLOBALE EF
C LPCOLO : LISTE DES NUMEROS DE COLONNES DES COEFFICIENTS NON NULS
C          DE LA MATRICE
C MOLPCOF: NOMBRE FINAL D'ENTIERS INITIALISES DU TABLEAU LPCOLO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY       ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1989
C MODIFS : ALAIN  PERRONNET  Saint PIERRE du PERRAY         FEVRIER 2021
C MODIFS : ALAIN  PERRONNET  Saint PIERRE du PERRAY         AVRIL   2023
C23456---------------------------------------------------------------012
      INTEGER  LPTDVOI(0:*), LISTVOI(*), LPLIGN(0:*),  LPCOLO(MOLPCOI)

C     NOMBRE TOTAL DE DEGRES DE LIBERTE OU DE LIGNES DE LA MATRICE
      NTDL = NBNOE * NBDLNO

C     FORMATION DES TABLEAUX LPLIGN ET LPCOLO
C     =======================================
      LPLIGN(0) = 0
      LPC0 = 0
      LPC  = 0
      IF( NCODSA .GT. 0 ) THEN

C        MATRICE SYMETRIQUE
C        ------------------
         LPV1 = LPTDVOI(0) + 1
         DO NOE=1,NBNOE

            LPCVOI = 0
            LPV2   = LPTDVOI(NOE)

C           CALCUL DU NOMBRE DE VOISINS DU NOEUD NOE DE NUMERO < NOE
            NBV = 0
            DO LPV=LPV1,LPV2
C              NUMERO DU NOEUD VOISIN DE NOE
               NUMV = LISTVOI( LPV )
               IF( NUMV .LT. NOE ) THEN
C                 LE NO EST INFERIEUR A NOE DONC RETENU
                  NBV = NBV + 1
               ENDIF
            ENDDO

C           PARCOURS DES NBV NOEUDS VOISINS DE NO < NOE
            IF( NBDLNO .EQ. 1 ) THEN

C              1 DL PAR NOEUD
               DO LPV=LPV1,LPV1+NBV-1
C                 NUMERO DU NOEUD VOISIN
                  NUMV = LISTVOI( LPV )
C                 LE DEGRE DE LIBERTE DU NOEUD NUMV
                  LPC = LPC + 1
                  LPCOLO(LPC) = NUMV
               ENDDO

C              L'UNIQUE DEGRE DE LIBERTE DU NOEUD NOE
               LPC = LPC + 1
               LPCOLO(LPC) = NOE
               LPLIGN(NOE) = LPC

            ELSE

C              PLUSIEURS DL PAR NOEUD
               DO LPV=LPV1,LPV1+NBV-1
C                 NUMERO DU NOEUD VOISIN
                  NUMV = LISTVOI(LPV)
C                 LES NBDLNO DEGRES DE LIBERTE DU NOEUD NUMV
                  LPLIV = ( NUMV - 1 ) * NBDLNO
                  DO IN=1,NBDLNO
                     DO JN=1,NBDLNO
                        LPCV = LPC0 + JN + ( ( IN - 1 ) * IN ) / 2
     %                       + NBV * NBDLNO * ( IN - 1 )
     %                       + LPCVOI
                        LPCOLO(LPCV) = LPLIV + JN
                     ENDDO
                  ENDDO
                  LPCVOI = LPCVOI + NBDLNO
               ENDDO

C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD NOE
               LPC  = LPC0
               LPLI = ( NOE - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,IN
                     LPC = LPC0 + JN + ( ( IN - 1 ) * IN ) / 2
     %                   + NBV * NBDLNO *  IN
                     LPCOLO(LPC) = LPLI + JN
                  ENDDO
                  LPLIGN(LPLI+IN) = LPC
               ENDDO

            ENDIF

            LPC0 = LPC
            LPV1 = LPV2 + 1

         ENDDO

      ELSE

C        MATRICE NON SYMETRIQUE
C        ----------------------
         LPV1 = 1
         DO NOE=1,NBNOE

            LPCVOI = 0
            LPV2 = LPTDVOI(NOE)
            NBV  = LPV2 - LPV1 + 1
            DO LPV=LPV1,LPV2
               NUMV = LISTVOI(LPV)
C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD NUMV
               LPLIV = ( NUMV - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,NBDLNO
                     LPCV = LPC0 + JN + LPCVOI
     %                    + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     LPCOLO(LPCV) = LPLIV + JN
                  ENDDO
               ENDDO
               LPCVOI = LPCVOI + NBDLNO
            ENDDO

C           LES NBDLNO DEGRES DE LIBERTE DU NOEUD NOE
            LPLI = ( NOE - 1 ) * NBDLNO
            DO IN=1,NBDLNO
               KN = 0
               DO JN=1,NBDLNO
                  IF (JN.NE.IN) THEN
                     KN = KN + 1
C                    PEUT ETRE UN PROBLEME SUR LPLIGN(NOE) QUI SUIT ???
                     LPC = LPC0 + KN + LPLIGN(NOE) * NBDLNO
     %                   + ( NBV + 1 ) * NBDLNO * ( IN - 1 )
                     LPCOLO(LPC) = LPLI + JN
                  ENDIF
               ENDDO
               LPC = LPC0 + ( NBV + 1 ) * NBDLNO *  IN
               LPCOLO(LPC) = LPLI + IN
               LPLIGN(LPLI+IN) = LPC
            ENDDO
            LPC0 = LPC
            LPV1 = LPV2 + 1

         ENDDO

      ENDIF

C     MISE A JOUR DU NUMERO DU DERNIER NUMERO DE COLONNE
      MOLPCOF = LPLIGN( NTDL )

      RETURN
      END
