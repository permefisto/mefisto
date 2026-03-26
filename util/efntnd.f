      SUBROUTINE EFNTND( NOOBVC , NOOBSF , NOOBLA , NOOBPS ,
     %                   NOTYOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C -----     EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, LES SURFACES ET LES VOLUMES
C
C ENTREES :
C ---------
C NOOBVC : NUMERO DU VOLUME  DE L'EF
C NOOBSF : NUMERO DE SURFACE DES NFACE  FACES   DE L'EF
C NOOBLA : NUMERO DE LIGNE   DES NARET  ARETES  DE L'EF
C NOOBPS : NUMERO DE POINT   DES NBNSOM SOMMETS DE L'EF
C
C SORTIES :
C ---------
C NOTYOB : NUMERO DU TYPE ET OBJET DE CHAQUE NOEUD DE L'EF
C          NOTYOB(1,*)=NO DE 1 A 4 ( 1:POINT , 2:LIGNE ,  ... )
C          NOTYOB(2,*)=NO DE L'OBJET DE CE TYPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/ponoel.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      INTEGER           NOOBSF(1:*),NOOBLA(1:*),NOOBPS(1:*)
      INTEGER           NOTYOB(2,*)
C
      IF( NOOBVC .LE. 0 ) THEN
C
C        LA VALEUR PAR DEFAUT EST ZERO
C        -----------------------------
         DO 10 I=1,NBNOE
            NOTYOB(1,I) = 0
            NOTYOB(2,I) = 0
 10      CONTINUE
C
      ELSE
C
C        LE VOLUME
C        ---------
         DO 20 I=1,NBNOE
C           LE NOEUD A LE NUMERO DE SON VOLUME
            NOTYOB(1,I) = 4
            NOTYOB(2,I) = NOOBVC
 20      CONTINUE
C
      ENDIF
C
C     LES SURFACES DES FACES
C     ----------------------
      DO 39 I=1,NFACE
         IF( NOOBSF(I) .GT. 0 ) THEN
C
C          LA FACE EST UNE PARTIE D'UNE SURFACE
C
C          LES NOEUDS INTERNES A LA FACE I
           DO 32 J=1,NBNOFA(I)
              L = NONOFA(J,I)
              NOTYOB(1,L) = 3
              NOTYOB(2,L) = NOOBSF(I)
 32        CONTINUE
C
C          LES NOEUDS INTERNES AUX ARETES DE LA FACE I
           DO 36 J=1,NBARFA(I)
C             LE NUMERO DE L'ARETE J DE LA FACE I
              NUAR = NOARFA(J,I)
              DO 34 K=1,NBNOAR(NUAR)
C                LES NOEUDS DE L'ARETE I
                 L = NONOAR(K,NUAR)
                 NOTYOB(1,L) = 3
                 NOTYOB(2,L) = NOOBSF(I)
 34           CONTINUE
 36        CONTINUE
C
C          LES SOMMETS DE LA FACE I
           DO 38 J=1,NBSOFA(I)
              L = NOSOFA(J,I)
              NOTYOB(1,L) = 3
              NOTYOB(2,L) = NOOBSF(I)
 38        CONTINUE
         ENDIF
 39   CONTINUE
C
C     LES LIGNES DES ARETES
C     ---------------------
      DO 80 I=1,NARET
         IF( NOOBLA(I) .GT. 0 ) THEN
C
C           L'ARETE EST SUR UNE LIGNE
C
            DO 75 J=1,NBNOAR(I)
C              LES NOEUDS DE L'ARETE I
               L = NONOAR(J,I)
               NOTYOB(1,L) = 2
               NOTYOB(2,L) = NOOBLA(I)
 75         CONTINUE
C
C           LES SOMMETS DE L'ARETE I
            DO 76 J=1,2
               L = NOSOAR(J,I)
               NOTYOB(1,L) = 2
               NOTYOB(2,L) = NOOBLA(I)
 76         CONTINUE
C
         ENDIF
 80   CONTINUE
C
C     LES POINTS DES SOMMETS
C     ----------------------
      DO 90 I=1,NBNSOM
         IF( NOOBPS(I) .GT. 0 ) THEN
            NOTYOB(1,I) = 1
            NOTYOB(2,I) = NOOBPS(I)
         ENDIF
 90   CONTINUE
C
      RETURN
      END
