      SUBROUTINE TS1LAG( X,      PENALI, NBJEUX, JEU,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NBPOLY, NPI,    POLY,
     %                   F1,     POIDEL, DP,
     %                   NOPART, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES SECONDS MEMBRES DES EF 1D LAGRANGE ISOPARAMETRIQUES
C -----
C
C ENTREES:
C --------
C X      : COORDONNEE DES NBPOLY POINTS DE L'EF
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NOOBLA : NUMERO DE LIGNE DE L'ELEMENT FINI
C SOURCE : FLUX NORMAL DE TEMPERATURE A LA PAROI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS LIGNES
C
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE SUR LE RECTANGLE
C POLY   : POLY(I, L) VALEUR DU POLYNOME I AU L-EME POINT D'INTEGRATION
C
C F1     : ABSCISSE DES NPI POINTS D'INTEGRATION
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C DP     : DP(NBPOLY, NPI) GRADIENT AUX POINTS D 'INTEGRATION DES
C          FONCTIONS  DE BASE ISOPARAMETRIQUES
C X      : COORDONNEES RAYON ET COTE DES NBPOLY POINTS DE L'EF
C NOPART : POUR NLSE SEULEMENT AU NIVEAU DE SOURCE=FORCE et CONTACT=FIXATION
C          1 SI PARTIE REELLE TRAITEE ou 2 SI PARTIE IMAGINAIRE TRAITEE
C          0 SI INACTIF (CAS THERMIQUE STANDARD D'UNE SOURCE)
C
C SORTIE :
C --------
C BE     : BE(NBPOLY) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     JUIN 2009
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  SOURCE(2),
     %                  POLY(NBPOLY, NPI),
     %                  POIDEL(NPI),
     %                  DP(NBPOLY, NPI), F1(NPI), XYZ(3),
     %                  BE(NBPOLY),
     %                  VITEFL(1),
     %                  PENALI
      REAL              X(NBPOLY)
      INTEGER           NOOBPS( 1:2 )
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      DOUBLE PRECISION  D, PROSCD
C
C     MISE A ZERO DE BE LE VECTEUR ELEMENTAIRE
C     ----------------------------------------
      CALL AZEROD( NBPOLY, BE )
C
C     ====================================
C     CONTRIBUTION DES SOURCES SUR L'ARETE
C     ====================================
      IF( LTDELI(LPSOUR,JEU,NOOBLA) .GT. 0 ) THEN
C
         IF( TESTNL .GE. 1 ) THEN
C           PB NON LINEAIRE:
C           LES SOURCES DEPENDENT DE LA TEMPERATURE
C           RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
            MN  = (MNTHET-1)/2
            MNT = (MNTHDL-1)/2
            DO 2 I=1,NBPOLY
               DMCN( MNT+I ) = DMCN( MN + MCN(MNNODL+I-1) )
 2          CONTINUE
         ENDIF
C
         DO 9 L=1,NPI
C
C           LA CONTRIBUTION DES SOURCES DE CHALEUR
            IF( TESTNL .GE. 1 ) THEN
C               PB NON LINEAIRE:
C               CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL = PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
            ENDIF
C
C           LA VALEUR DES SOURCES LINEIQUES EN CE POINT D'INTEGRATION
            XYZ(1) = F1(L)
            XYZ(2) = 0D0
            XYZ(3) = 0D0
            IF( TESTNL .LE. 5 ) THEN
               CALL RESOUR( 2, NOOBLA, 3, XYZ,
     %                      LTDELI(LPSOUR,JEU,NOOBLA), SOURCE )
            ELSE
               CALL REFORC( 2, NOOBLA, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                    0D0,   0D0,   0D0,
     %                      LTDELI(LPSOUR,JEU,NOOBLA), SOURCE )
               SOURCE(1) = SOURCE(NOPART)
            ENDIF
C
            DO 7 I=1,NBPOLY
               BE(I) = BE(I) + POIDEL(L) * POLY(I,L) * SOURCE(1)
 7          CONTINUE
 9       CONTINUE
C
      ENDIF
C
C     LA CONTRIBUTION DU TRANSPORT: - VITESSE * GRADIENT TEMPERATURE
C     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
C     -----------------------------------------------------------------
      IF( LTDELI(LPVIFL,JEU,NOOBLA) .GT. 0 ) THEN
C
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN  = (MNTHET-1)/2
         MNT = (MNTHDL-1)/2
         DO 12 I=1,NBPOLY
            DMCN( MNT+I ) = DMCN( MN + MCN(MNNODL+I-1) )
 12      CONTINUE
C
         DO 18 L=1,NPI
C
            IF( TESTNL .GE. 1 ) THEN
C              PB NON LINEAIRE:
C              CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL = PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
            ENDIF
C
C           LA VALEUR DE lA VITESSE DU FLUIDE AU POINT D'INTEGRATION L
            XYZ(1) = F1(L)
            XYZ(2) = 0D0
            XYZ(3) = 0D0
            CALL REVIFL( 2, NOOBLA, 1, 3, XYZ,
     %                   LTDELI(LPVIFL,JEU,NOOBLA), VITEFL )
C
C           - t[P] [V] [DP] [TEMPERATURE EF]
            DO 17 I=1,NBPOLY
C
               D = 0D0
               DO 15 K=1,NBPOLY
                  D = D + VITEFL(1) * DP(K,L) * DMCN(MNT+K)
 15            CONTINUE
C
               BE(I) = BE(I) - POLY(I,L) * D * POIDEL(L)
C
 17         CONTINUE
 18      CONTINUE
C
      ENDIF
C
C
C     ==============================================================
C     CONTRIBUTION D'UNE SOURCE OU UN CONTACT PENALISE AUX 2 SOMMETS
C     ==============================================================
      DO 50 K=1,2
C
C        LE NUMERO DE POINT DU SOMMET K
         NOOB = NOOBPS(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT
            IESOUR = 0
C
C           POINT SUPPORT D'UNE SOURCE DE CHALEUR?
            IF( LTDEPO(LPSOUR,JEU,NOOB) .GT. 0  ) IESOUR = 1
C
C           POINT SUPPORT D'UN CONTACT PENALISE ?
            IF( LTDEPO(LPCONT,JEU,NOOB) .GT. 0 .AND.
     %          PENALI .NE. 0D0                 ) IESOUR = 2
C
            IF( IESOUR .EQ. 0 ) GOTO 50
C
            IF( TESTNL .GE. 1 ) THEN
C              PB NON LINEAIRE: RECUPERATION DE LA TEMPERATURE AU SOMMET K
               TEMPEL = DMCN( (MNTHET-1)/2 + MCN(MNNODL+K-1) )
            ENDIF
C
            XYZ(1) = X(K)
            XYZ(2) = 0D0
            XYZ(3) = 0D0
            IF( IESOUR .EQ. 1 ) THEN
C
C              FLUX NORMAL DE CHALEUR: CONDITION NEUMANN OU FOURIER
C              UN TABLEAU SOURCE EXISTE POUR CE POINT
               IF( TESTNL .LE. 5 ) THEN
                  CALL RESOUR( 1, NOOB, 3, XYZ,
     %                         LTDEPO(LPSOUR,JEU,NOOB), SOURCE )
               ELSE
                  CALL REFORC( 1, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                     0D0,   0D0,   0D0,
     %                         LTDEPO(LPSOUR,JEU,NOOB), SOURCE )
                  SOURCE(1) = SOURCE(NOPART)
               ENDIF
C
            ELSE
C
               IF( TESTNL .LE. 5 ) THEN
C                 CONTACT PENALISE = PENALI x TEMPERATURE
                  CALL RECONT( 1, NOOB, 3, XYZ,
     %                         LTDEPO(LPCONT,JEU,NOOB), SOURCE )
                  SOURCE(1) = SOURCE(1) * PENALI
C
               ELSE
C                 FIXATION(2) PENALISEE
                  MN = LTDEPO(LPCONT,JEU,NOOB)
                  CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, SOURCE )
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     IF( NU .EQ. NOPART ) THEN
                        SOURCE(1) = SOURCE(I) * PENALI
                     ENDIF
                  ENDDO
               ENDIF
C
            ENDIF
C
C           CONTRIBUTION DU SOMMET K AU SECOND MEMBRE
            BE(K) = BE(K) + SOURCE(1)
C
         ENDIF
 50   CONTINUE
C
      RETURN
      END
