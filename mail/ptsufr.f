      SUBROUTINE PTSUFR( NMVOLU, NBSOM, MNSOFR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  A PARTIR DES ELEMENTS FINIS 3D FORMATION POUR CHAQUE SOMMET
C -----  DE SON CODE INTERNE (0), FRONTALIER (1) ou INTERFACE (2)
C        SOMMET FRONTIERE SI SOMMET D'UNE FACE VUE UNE SEULE FOIS
C        SOMMET INTERFACE SI SOMMET D'UNE FACE VUE 2 FOIS DANS 2 EF
C        DE MATERIAUX DIFFERENTS
C        SOMMET INTERNE SI AUCUN DES 2 PRECEDENTS CAS

C ENTREES:
C --------
C NMVOLU : NOM DE L'OBJET A TRACER

C SORTIES :
C ---------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE DU VOLUME NMVOLU
C MNSOFR : TABLEAU SOFR : = 0 SI SOMMET INTERNE
C                         = 1 SI SOMMET FRONTALIER
C                         = 2 SI SOMMET SUR INTERFACE 2 MATERIAUX
C          SI SOMMET FRONTALIER ET INTERFACE => IL EST FORCE INTERFACE!
C          ADRESSE MCN DU TABLEAU SOFR OU 0 EN PRESENCE D'UNE ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    Decembre 1991
C MODIF  : PERRONNET ALAIN LJLL UPMC ET ST PIERRE DU PERRAY     Mai 2008
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___face.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMVOLU

C     GENERATION EVENTUELLE PAR HACHAGE DES FACES DES EF VOLUMIQUES
C     CHAINAGE DES FACES FRONTALIERES EN POSITION 7
C     AVEC UN LIEN NEGATIF POUR LES FACES FRONTALIERES
C     =============================================================
      MNSOFR = 0
      CALL HAFAVO( NMVOLU, 3, NTFAVO, MNFAVO, NBSOM, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ptsufr: VOLUME: ' // NMVOLU
         KERR(2) = 'IMPOSSIBLE DE CREER SES FACES'
         CALL LEREUR
         RETURN
      ENDIF

C     LE NOMBRE D'ENTIERS PAR FACE
      MOFACE = MCN( MNFAVO + WOFACE )
C     LE NOMBRE DE FACES FRONTALIERES
      NBFAFR = MCN( MNFAVO + WBFAFR )
C     LE NUMERO DE LA PREMIERE FACE FRONTALIERE
      L1FAFR = MCN( MNFAVO + W1FAFR )
C     LE NOMBRE DE FACES INTERFACES
      NBFA2M = MCN( MNFAVO + WBFA2M )
C     LE NUMERO DE LA PREMIERE FACE INTERFACE
      L1FA2M = MCN( MNFAVO + W1FA2M )

C     RESERVATION D'UN TABLEAU POUR LE NUMERO DE CODE FRONTIERE DES SOMMETS
C     =====================================================================
      CALL TNMCDC( 'ENTIER', NBSOM, MNSOFR )
      CALL AZEROI( NBSOM, MCN(MNSOFR) )

C     PARCOURS DES FACES FRONTALIERES (VUES 1 SEULE FOIS )
C     ====================================================
      MNFACE = MNFAVO + WFACES - MOFACE - 1
      DO 20 NF = NBFAFR, 1, -1

C        L'ADRESSE DU PREMIER SOMMET DE LA FACE FRONTALIERE
         MN = MNFACE + L1FAFR * MOFACE
C
C        LE NOMBRE DE SOMMETS DE LA FACE
         IF( MCN( MN+4 ) .GT. 0 ) THEN
            NAF = 4
         ELSE
            NAF = 3
         ENDIF

         DO 10 J=1,NAF
C           LE NUMERO DU J-EME SOMMET DE LA FACE FRONTALIERE
            NS = MCN( MN + J )
C           1 COMME FRONTALIER
            MCN( MNSOFR - 1 + NS ) = 1
 10      CONTINUE

C        LA FACE FRONTIERE SUIVANTE
         L1FAFR = MCN( MN + 7 )
 20   CONTINUE

C     PARCOURS DES FACES INTERFACES ENTRE 2 MATERIAUX
C     ===============================================
      DO 50 NF = NBFA2M, 1, -1

C        L'ADRESSE DU PREMIER SOMMET DE LA FACE INTERFACE
         MN = MNFACE + L1FA2M * MOFACE

C        LE NOMBRE DE SOMMETS DE LA FACE INTERFACE
         IF( MCN( MN+4 ) .GT. 0 ) THEN
            NAF = 4
         ELSE
            NAF = 3
         ENDIF

         DO 40 J=1,NAF
C           LE NUMERO DU J-EME SOMMET DE LA FACE INTERFACE
            NS = MCN( MN + J )
C           2 COMME SUR UN INTERFACE ENTRE 2 MATERIAUX
            MCN( MNSOFR - 1 + NS ) = 2
 40      CONTINUE

C        LA FACE INTERFACE SUIVANTE
         L1FA2M = MCN( MN + 8 )
 50   CONTINUE

      RETURN
      END
