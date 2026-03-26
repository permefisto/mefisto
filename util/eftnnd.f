      SUBROUTINE EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                   NOMTMS, MXDOTH, LPNTMS,
     %                   NOTYOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     LE CALCUL DU TYPE PLSV DE CHAQUE NOEUD DE L'ELEMENT FINI
C -----     QUI SUPPORTE UN TABLEAU "TMS" NOMTMS
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, LES SURFACES ET LES VOLUMES
C
C ENTREES:
C --------
C NOOBVC : NUMERO DU VOLUME  DE L'EF
C NOOBSF : NUMERO DE SURFACE DES NFACE  FACES   DE L'EF
C NOOBLA : NUMERO DE LIGNE   DES NARET  ARETES  DE L'EF
C NOOBPS : NUMERO DE POINT   DES NBNSOM SOMMETS DE L'EF
C
C NUMIOB : NUMERO MINIMAL DES PLSV
C MNDOEL : ADRESSE MCN DES DONNEES (THERMIQUE ou ELASTICITE ...) DE L'OBJET
C
C NOMTMS : NOM DU TABLEAU A SELECTIONNER FILTRE POUR LES TYPES DES NOEUDS
C MXDOTH : NOMBRE MAXIMAL DE TABLEAUX TMS DES DONNEES THERMIQUE ou ELASTICITE ...
C LPNTMS : NUMERO DE "NOMTMS" PARMI LES TABLEAUX TMS DES DONNEES
C
C SORTIES:
C --------
C NOTYOB : NUMERO DU TYPE ET OBJET DE CHAQUE NOEUD DE L'EF SUPPORTANT "NOMTMS"
C          NOTYOB(1,*)=NO DE 1 A 4 ( 1:POINT , 2:LIGNE ,  ... )
C          NOTYOB(2,*)=NO DE L'OBJET DE CE TYPE
C          NOTYOB(3,*)=ADRESSE MCN DU TABLEAU "NOMTMS"
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C23456---------------------------------------------------------------012
      include"./incl/ponoel.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      CHARACTER*(*)     NOMTMS
      INTEGER           NOOBSF(1:*),NOOBLA(1:*),NOOBPS(1:*)
      INTEGER           NUMIOB(4),MNDOEL(4)
      INTEGER           NOTYOB(3,*)
C
      IF( NOOBVC .LT. 0 )
     %    WRITE(IMPRIM,*) 'EFTNND: NOOBVC=',NOOBVC,'  NOMTMS=',NOMTMS
C
      IF( NOOBVC .EQ. 0 ) THEN
C
C        LA VALEUR PAR DEFAUT EST ZERO
C        -----------------------------
         DO 10 I=1,NBNOE
            NOTYOB(1,I) = 0
            NOTYOB(2,I) = 0
            NOTYOB(3,I) = 0
 10      CONTINUE
C
      ELSE
C
C        LE VOLUME
C        ---------
         MN = MNDOEL(4)
     %      + MXDOTH * ( NOOBVC - NUMIOB(4) ) - 1
     %      + LPNTMS
         MN = MCN( MN )
C
C        LE TMS "NOMTMS" EST SUPPORTE PAR CE VOLUME
C        => TOUS LES NOEUDS DE L'EF SONT DU TYPE DU VOLUME
         DO 20 I=1,NBNOE
C           LE NOEUD A LE NUMERO DE SON VOLUME
            NOTYOB(1,I) = 4
            NOTYOB(2,I) = NOOBVC
            NOTYOB(3,I) = MN
 20      CONTINUE
C
      ENDIF
C
C     LES SURFACES DES FACES SUPPORTANT "NOMTMS"
C     ------------------------------------------
      DO 39 I=1,NFACE
         IF( NOOBSF(I) .GT. 0 ) THEN
C
C          LA FACE I EST UNE PARTIE D'UNE SURFACE
           MN = MNDOEL(3)
     %        + MXDOTH * ( NOOBSF(I) - NUMIOB(3) ) - 1
     %        + LPNTMS
           MN = MCN( MN )
           IF( MN .LE. 0 ) GOTO 39
C
C          LES NOEUDS INTERNES A LA FACE I
           DO 32 J=1,NBNOFA(I)
              L = NONOFA(J,I)
              NOTYOB(1,L) = 3
              NOTYOB(2,L) = NOOBSF(I)
              NOTYOB(3,L) = MN
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
                 NOTYOB(3,L) = MN
 34           CONTINUE
 36        CONTINUE
C
C          LES SOMMETS DE LA FACE I
           DO 38 J=1,NBSOFA(I)
              L = NOSOFA(J,I)
              NOTYOB(1,L) = 3
              NOTYOB(2,L) = NOOBSF(I)
              NOTYOB(3,L) = MN
 38        CONTINUE
         ENDIF
 39   CONTINUE
C
C     LES LIGNES DES ARETES SUPPORTANT "NOMTMS"
C     -----------------------------------------
      DO 80 I=1,NARET
         IF( NOOBLA(I) .GT. 0 ) THEN
C
C           L'ARETE EST SUR UNE LIGNE
            MN = MNDOEL(2)
     %         + MXDOTH * ( NOOBLA(I) - NUMIOB(2) ) - 1
     %         + LPNTMS
            MN = MCN( MN )
            IF( MN .LE. 0 ) GOTO 80
C
            DO 75 J=1,NBNOAR(I)
C              LES NOEUDS DE L'ARETE I
               L = NONOAR(J,I)
               NOTYOB(1,L) = 2
               NOTYOB(2,L) = NOOBLA(I)
               NOTYOB(3,L) = MN
 75         CONTINUE
C
C           LES SOMMETS DE L'ARETE I
            DO 76 J=1,2
               L = NOSOAR(J,I)
               NOTYOB(1,L) = 2
               NOTYOB(2,L) = NOOBLA(I)
               NOTYOB(3,L) = MN
 76         CONTINUE
C
         ENDIF
 80   CONTINUE
C
C     LES POINTS DES SOMMETS SUPPORTANT "NOMTMS"
C     ------------------------------------------
      DO 90 I=1,NBNSOM
         IF( NOOBPS(I) .GT. 0 ) THEN
            MN = MNDOEL(1)
     %         + MXDOTH * ( NOOBPS(I) - NUMIOB(1) ) - 1
     %         + LPNTMS
            MN = MCN( MN )
            IF( MN .LE. 0 ) GOTO 90
            NOTYOB(1,I) = 1
            NOTYOB(2,I) = NOOBPS(I)
            NOTYOB(3,I) = MN
         ENDIF
 90   CONTINUE
C
      RETURN
      END
