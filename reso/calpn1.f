      SUBROUTINE CALPN1(NTDL,NCODSA,LPLIGN,LPCOLO,LPLIGC,LPCOLC,LPDILU)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISE LES POINTEURS LPLIGC ET LPCOLC
C ---
C
C PARAMETRES D ENTREE :
C ---------------------
C NTDL          : LE NOMBRE TOTAL D'INCONNUES
C NCODSA        : LE MODE DE STOCKAGE DE LA MATRICE MORSE
C LPLIGC,LPCOLC : LES POINTEURS ASSOCIES A LA MATRICE MORSE AGC
C LPDILU        : TABLEAU AUXILIAIRE
C
C PARAMETRE DE SORTIE :
C --------------------
C LPLIGC,LPCOLC : LES POINTEURS ASSOCIES A LA MATRICE MORSE AGC
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1989
C23456---------------------------------------------------------------012
      DIMENSION LPLIGN(NTDL+1),LPCOLO(*)
      DIMENSION LPLIGC(NTDL+1),LPCOLC(*),LPDILU(NTDL)

C     INITIALISATION
C     --------------
      DO I=1,NTDL+1
         LPLIGC(I)=LPLIGN(I)
      ENDDO

      IF(NCODSA.LT.0) THEN

C        STOCKAGE NON SYMETRIQUE
C        -----------------------
         DO I=1,LPLIGN(NTDL+1)
            LPCOLC(I)=LPCOLO(I)
         ENDDO
         RETURN

      ELSE

C        STOCKAGE SYMETRIQUE
C        -------------------
         DO I=1,NTDL
            LPDILU(I)=0
         ENDDO
         DO I=1,NTDL
            K1=LPLIGN(I)+1
            K2=LPLIGN(I+1)-1
            DO K=K1,K2
               J=LPCOLO(K)
               LPDILU(J)=LPDILU(J)+1
            ENDDO
         ENDDO

C        LE POINTEUR LPLIGC
         LLP=0
         DO I=1,NTDL
            LLP=LLP+LPDILU(I)
            LPLIGC(I+1)=LPLIGC(I+1)+LLP
            LPDILU(I)=LPLIGC(I)
         ENDDO

C        LE TABLEAU LPCOLC
         DO I=1,NTDL
            DO K=LPLIGN(I)+1,LPLIGN(I+1)-1
               J=LPCOLO(K)
               LI=LPDILU(I)+1
               LPCOLC(LI)=J
               LPDILU(I)=LI
               LJ=LPDILU(J)+1
               LPCOLC(LJ)=I
               LPDILU(J)=LJ
            ENDDO
            LPCOLC(LPLIGC(I+1))=I
         ENDDO

C        RECLASSEMENT DE LPCOLC
         DO I=1,NTDL
           K1=LPLIGC(I)+1
           K2=LPLIGC(I+1)-1
           DO K=K1,K2
              DO L=K,K2

ccc              IF(LPCOLC(K)-LPCOLC(L)) 10,10,11

                 IF(LPCOLC(K) .GT. LPCOLC(L)) THEN
                    LPCOLCK  =LPCOLC(K)
                    LPCOLC(K)=LPCOLC(L)
                    LPCOLC(L)=LPCOLCK
                 ENDIF

              ENDDO
           ENDDO
        ENDDO

      ENDIF

      RETURN
      END
