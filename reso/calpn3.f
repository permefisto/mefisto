      SUBROUTINE CALPN3( NTDL,   LPLIGC, LPCOLC, LPDILU,
     &                   LPDIRE, LPCORE, INDIC,  LONG,   IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CE SP CALCULE LES POINTEURS LPLIGC,LPCOLC
C ----  CORRESPONDANTS AU 'REMPLISSAGE' DE LA
C       MATRICE EN COURS DE FACTORISATION
C
C PARAMETRES D'ENTREE:
C -------------------
C LPLIGC,LPCOLC,LPDILU : LES POINTEURS DE LA MATRICE INITIALE
C INDIC                : TABLEAU UTILITAIRE
C LONG                 : TAILLE MAXIMALE DU TABLEAU LPCOLC
C
C PARAMETRE DE SORTIE:
C -------------------
C LPDIRE,LPCORE       : LES POINTEURS DE REMPLISSAGE
C IERR                : CODE D'ERREUR >0 SI ERREUR
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : P.JOLY  LABORATOIRE D'ANALYSE NUMERIQUE  PARIS 6
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      DIMENSION LPLIGC(NTDL+1),LPCOLC(1:*),LPDILU(NTDL)
      DIMENSION LPDIRE(NTDL+1),LPCORE(1:*),INDIC(NTDL)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)

C     CALCUL DU REMPLISSAGE : TABLEAUX LPDIRE ET LPCORE
      DO I=1,NTDL
         INDIC(I)=0
      ENDDO

      MAXCOE = 0
      LPDIRE(1)=0
      IC1=0
      K1 =1

      DO I=1,NTDL
C        TENTATIVE DE PROTECTION CONTRE LE DEBORDEMENT (COMPROMIS)
         IF( IC1+MAXCOE .GT. LONG ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'TAILLE DE LA MEMOIRE INSUFFISANTE'
               KERR(2) = 'AUGMENTER MOTMCN OU REDUIRE LE MAILLAGE'
            ELSE
               KERR(1) = 'NOT ENOUGH MEMORY'
               KERR(2) = 'AUGMENT MOTMCN or REDUCE the MESH'
            ENDIF
            CALL LEREUR
            IERR = 1
            RETURN
         ENDIF
         K2=LPLIGC(I+1)
         DO K=K1,K2
            J=LPCOLC(K)
            INDIC(J)=I
         ENDDO
         DO K=K1,LPDILU(I)
            J=LPCOLC(K)
            DO L=LPDILU(J)+1,LPLIGC(J+1)
               JJ=LPCOLC(L)
               IF(INDIC(JJ).NE.I) THEN
                  IC1=IC1+1
                  LPCORE(IC1)=JJ
                  INDIC(JJ)=I
               END IF
            ENDDO
         ENDDO
         LPDIRE(I+1)=IC1
C        LE NOMBRE MAXIMAL DE COEFFICIENTS D'UNE LIGNE PRECEDENTE
         K = IC1-LPDIRE(I)
         IF( K .GT. MAXCOE ) MAXCOE = K
         K1=K2+1
      ENDDO

C     REUNION DE LPCOLC ET LPCORE DANS LPCOLC
      IF(LPDIRE(NTDL+1).GT.0) THEN

C     TEST SUR LA PLACE MEMOIRE
      KC2=LPLIGC(NTDL+1)+LPDIRE(NTDL+1)
      IF( KC2.GT.LONG ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,1000) KC2,LONG
         ELSE
            WRITE(IMPRIM,2000) KC2,LONG
         ENDIF
         IERR = 1
         RETURN
      END IF

      DO I=NTDL,1,-1
         KC2A=LPLIGC(I+1)
         KD2=LPDIRE(I+1)+1
         LC=LPLIGC(I+1)-LPLIGC(I)-1
         LD=LPDIRE(I+1)-LPDIRE(I)
         LPLIGC(I+1)=KC2
         LPCOLC(KC2)=I
         DO J=1,LD
            LPCOLC(KC2-J)=LPCORE(KD2-J)
         ENDDO
         KC2=KC2-LD
         DO J=1,LC
            LPCOLC(KC2-J)=LPCOLC(KC2A-J)
         ENDDO
         KC2=KC2-LC

C        ON ORDONNE LES D.L. DE LA LIGNE I PAR ORDRE CROISSANT
         KC=LPLIGC(I+1)-1
         DO K=KC2,KC
            DO L=K,KC

ccc            IF(LPCOLC(K)-LPCOLC(L)) 9,9,10
               IF( LPCOLC(K) .GT. LPCOLC(L) ) THEN
                  MACK     =LPCOLC(K)
                  LPCOLC(K)=LPCOLC(L)
                  LPCOLC(L)=MACK
               ENDIF

            ENDDO
         ENDDO
         KC2=KC2-1
      ENDDO

      LPLIGC(1)=0
      ENDIF
C
      RETURN
1000  FORMAT(' NOMBRE de COEFFICIENTS NON NULS de AGC',T41,I10/,
     %       ' SUPERIEUR a la LIMITE AUTORISEE ',T41,I10)
2000  FORMAT(' NUMBER of NO NULL  COEFFICIENTS of AGC',T41,I10/,
     %       ' GREATEUR than the AUTORIZED LIMIT ',T41,I10)
      END
