      SUBROUTINE ADARFR( TARMAX, MXSOMM, NBSOMM, NUPLIS, NOSUTR, NUISOP,
     %                   MXTRIA, N1TRVI, NOTRIA, NOTRSO,
     %                   NUMILF, NUMXLF, L1LGFR,
     %                   MXARCH, NBARCH, N1ARCH, NSARCH,
     %                   PXYD  , NBTRIA, NBARDE, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUT DU MILIEU DES ARETES FRONTALIERES DE TAILLE
C -----    PRESQUE EGALE AU DEMI PERIMETRE DU TRIANGLE OU
C          SUPERIEURE A LA TAILLE MAXIMALE PERMISE
C
C ENTREES:
C --------
C TARMAX : TAILLE DE L'ARETE MAXIMALE PERMISE SUR UNE FRONTIERE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C NUPLIS : NUMERO DE POINT OU LIGNE DE CHACUN DES SOMMETS
C          -2 000 000 - NUPLIS() ANCIEN SI LE SOMMET A DEJA ETE SUPPRIME
C          -NP SI NP EST LE NUMERO DU POINT UTILISATEUR DE CE SOMMET
C          -1 000 000 - NL1 - 1000 * NL2 SI LE SOMMET APPARTIENT
C                     AUX 2 LIGNES NL1 > NL2
C          NU LIGNE SI LE SOMMET EST SUR UNE LIGNE UTILISATEUR
C          0        SINON
C NOSUTR : NUMERO DE SURFACE DE CHAQUE TRIANGLE
C NUISOP : NUMERO DE L'ISO DU POINT ET 0 SINON
C N1TRVI : POINTE DANS NOTRIA VERS LE PREMIER TRIANGLE VIDE
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NOTRIA
C NOTRIA : LISTE CHAINEE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                           ADJACENT PAR L'ARETE I
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C
C NUMILF : NUMERO MINIMAL DES LIGNES FRONTIERES DE L'OBJET
C NUMXLF : NUMERO MAXIMAL DES LIGNES FRONTIERES DE L'OBJET
C L1LGFR : POUR CHAQUE LIGNE DE NUMILF A NUMXLF POINTE SUR LA PREMIERE
C          ARETE CHAINEE DE LA LIGNE
C
C MXARCH : NOMBRE MAXIMAL D'ARETES DECLARABLES DANS NSARCH
C NBARCH : NOMBRE D'ARETES CHAINEES DANS NSARCH
C N1ARCH : NUMERO DANS NSARCH DE LA PREMIERE ARETE VIDE
C NSARCH : NUMERO DES 2 SOMMETS DE L'ARETE, ARETE PRECEDENTE ET SUIVANTE
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE ACTUEL DE SOMMETS ( INTERNES ET EXTERNES )
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C NBTRIA : NOMBRE DE TRIANGLES DANS NOTRIA
C
C SORTIE :
C --------
C NBARDE : NOMBRE D'ARETES FRONTALIERES DECOUPEES
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1995
C....................................................................012
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION(CF DVTR2D)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NOTRSO(MXSOMM),
     %                  NUPLIS(MXSOMM),
     %                  NOSUTR(MXTRIA),
     %                  NUISOP(*),
     %                  L1LGFR(NUMILF:NUMXLF),
     %                  NSARCH(4,MXARCH)
      DOUBLE PRECISION  PXYD(3,MXSOMM),RLONG(3),P
C
      IERR   = 0
      NBARDE = 0
      DO 100 NT=1,MXTRIA
         IF( NOTRIA(1,NT) .GT. 0 ) THEN
C        LE TRIANGLE EXISTE
            DO 90 KK=1,3
               NTOP = NOTRIA(3+KK,NT)
               IF( NTOP .EQ. 0 ) THEN
C                 ARETE FRONTALIERE CALCUL DU PERIMETRE DU TRIANGLE NT
                  DO 10 K=1,3
                     IF( K .LT. 3 ) THEN
                        K1 = K + 1
                     ELSE
                        K1 = 1
                     ENDIF
                     NS1 = NOTRIA(K, NT)
                     NS2 = NOTRIA(K1,NT)
                     RLONG(K) = SQRT( (PXYD(1,NS2)-PXYD(1,NS1))**2
     %                              + (PXYD(2,NS2)-PXYD(2,NS1))**2 )
 10               CONTINUE
                  P = RLONG(1)+RLONG(2)+RLONG(3)
                  IF( P-RLONG(KK) .LT. RLONG(KK) * 1.2 .OR.
     %                RLONG(KK) .GT. TARMAX ) THEN
C
C                    AJOUT DU MILIEU DE L'ARETE FRONTALIERE (KK DANS NT)
                     IF( KK .LT. 3 ) THEN
                        K1 = KK + 1
                     ELSE
                        K1 = 1
                     ENDIF
                     IF( K1 .LT. 3 ) THEN
                        K2 = K1 + 1
                     ELSE
                        K2 = 1
                     ENDIF
C                    LES 3 SOMMETS DE NT
                     NS1 = NOTRIA(KK, NT)
                     NS2 = NOTRIA(K1, NT)
                     NS3 = NOTRIA(K2, NT)
                     IF( NBSOMM .GE. MXSOMM ) THEN
                        NBLGRC(NRERR) = 1
                        KERR(1) = 'TABLEAU PXYD SATURE'
                        CALL LEREUR
                        IERR = 1
                        RETURN
                     ENDIF
                     NBSOMM = NBSOMM + 1
                     PXYD(1,NBSOMM) = ( PXYD(1,NS1)+PXYD(1,NS2) )*0.5D0
                     PXYD(2,NBSOMM) = ( PXYD(2,NS1)+PXYD(2,NS2) )*0.5D0
                     PXYD(3,NBSOMM) = ( PXYD(3,NS1)+PXYD(3,NS2) )*0.5D0
C                    CE NOUVEAU POINT N'EST PAS SUR UNE ISO
                     NUISOP(NBSOMM) = 0
C
C                    CALCUL DU NUMERO DE LIGNE FRONTIERE DU POINT NBSOMM
                     IF( NUPLIS(NS1) .GT. - 2 000 000 .AND.
     %                   NUPLIS(NS1) .LT. - 1 000 000  ) THEN
C                       NS1 APPARTIENT A 2 LIGNES NL1 > NL2
C                       NUPLIS(NS1) = - 1 000 000 - NL1 - 1000 * NL2
                        NL1 = - NUPLIS(NS1) - 1 000 000
                        NL2 = NL1 / 1000
                        NL1 = NL1 - NL2 * 1000
                     ELSE
                        NL1 = NUPLIS(NS1)
                        NL2 = 0
                     ENDIF
C
                     IF( NUPLIS(NS2) .GT. - 2 000 000 .AND.
     %                   NUPLIS(NS2) .LT. - 1 000 000  ) THEN
C                       NS2 APPARTIENT A 2 LIGNES ML1 > ML2
C                       NUPLIS(NS2) = - 1 000 000 - ML1 - 1000 * ML2
                        ML1 = - NUPLIS(NS2) - 1 000 000
                        ML2 = ML1 / 1000
                        ML1 = ML1 - ML2 * 1000
                     ELSE
                        ML1 = NUPLIS(NS2)
                        ML2 = 0
                     ENDIF
C
C                    LE NUMERO DE LIGNE FRONTIERE DU POINT NBSOMM
                     NL = 0
                     IF( NL1 .EQ. ML1 .OR. NL1 .EQ. ML2 ) NL = NL1
                     IF( NL  .EQ. 0  .AND.
     %                  (NL2 .EQ. ML1 .OR. NL2 .EQ. ML2) ) NL = NL2
                     IF( NL .EQ. 0 .AND. NL1 .GT. 0 ) NL = NL1
                     IF( NL .EQ. 0 .AND. ML1 .GT. 0 ) NL = ML1
                     IF( NL .LE. 0 ) THEN
                        WRITE(IMPRIM,*)
     %                 'PT MILIEU FRONTALIER NON SUR UNE LIGNE'
                        NBSOMM = NBSOMM - 1
                        GOTO 90
                     ENDIF
                     NUPLIS(NBSOMM) = NL
C
C                    RECHERCHE DE L'ARETE NS1-NS2 DANS LE CHAINAGE DE LA LIGNE N
                     L = L1LGFR(NL)
 30                  IF( L .GT. 0 ) THEN
                        IF( NS1 .EQ. NSARCH(1,L) .AND.
     %                      NS2 .EQ. NSARCH(2,L) ) THEN
                            GOTO 50
                        ELSE IF( NS2 .EQ. NSARCH(1,L) .AND.
     %                           NS1 .EQ. NSARCH(2,L) ) THEN
                            GOTO 50
                        ENDIF
C                       PASSAGE A L'ARETE SUIVANTE
                        L = NSARCH(4,L)
                        GOTO 30
                     ENDIF
C                    ARETE NON RETROUVEE
                     GOTO 90
C
C                    MODIFICATION DE L'ARETE EN ARETES NS1-NBSOMM + NBSOMM-NS2
 50                  NSARCH(1,L) = NS1
                     NSARCH(2,L) = NBSOMM
C                    AJOUT DE LA NOUVELLE
                     CALL AJARCH( NBSOMM, NS2, L1LGFR(NL),
     %                            MXARCH, NBARCH, N1ARCH, NSARCH, IERR )
                     IF( IERR .NE. 0 ) GOTO 90
C
C                    DECOUPAGE DU TRIANGLE NT EN 2 TRIANGLES
                     NT1 = NOTRIA(3+K1,NT)
                     NT2 = NOTRIA(3+K2,NT)
                     IF( N1TRVI .LE. 0 ) THEN
C                       SATURATION DES TRIANGLES
                        IERR = 3
                        NBLGRC(NRERR) = 1
                        KERR(1) = 'SATURATION DES TRIANGLES'
                        CALL LEREUR
                        GOTO 90
                     ENDIF
C                    LE TRIANGLE A AJOUTER
                     NTOP = N1TRVI
C                    LE NOUVEAU TRIANGLE VIDE
                     N1TRVI = NOTRIA(4,N1TRVI)
C
C                    MISE A JOUR DU TRIANGLE NT
                     NOTRIA(KK,  NT) = NBSOMM
                     NOTRIA(3+K2,NT) = NTOP
C
C                    MISE A JOUR DU TRIANGLE NTOP
                     NOTRIA(1, NTOP) = NS1
                     NOTRIA(2, NTOP) = NBSOMM
                     NOTRIA(3, NTOP) = NS3
                     NOTRIA(4, NTOP) = 0
                     NOTRIA(5, NTOP) = NT
                     NOTRIA(6, NTOP) = NT2
                     NOSUTR(NTOP) = NOSUTR(NT)
C
C                    MISE A JOUR DU TRIANGLE NT2
                     IF( NT2 .GT. 0 ) THEN
                        DO 80 L=4,6
                           IF( NOTRIA(L,NT2) .EQ. NT ) THEN
                              NOTRIA(L,NT2) = NTOP
                              GOTO 85
                           ENDIF
 80                     CONTINUE
                     ENDIF
C
C                    MISE A JOUR DE NOTRSO
 85                  NOTRSO(NBSOMM) = NT
                     NOTRSO(NS1   ) = NTOP
C
C                    UNE ARETE FRONTALIERE DECOUPEE DE PLUS
                     NBARDE = NBARDE + 1
                     NBTRIA = NBTRIA + 1
C
C                    TRACE DES 2 TRIANGLES
                     CALL DVTRTR(PXYD, NOTRIA, NT,   NCJAUN, NCNOIR )
                     CALL DVTRTR(PXYD, NOTRIA, NTOP, NCJAUN, NCNOIR )
                  ENDIF
               ENDIF
 90         CONTINUE
         ENDIF
 100  CONTINUE
C
C     RECHERCHE D'UNE LIGNE FRONTALIERE REDUITE A UNE ARETE CHAINEE
C     POUR DECOUPAGE EN 2 ARETES FRONTALIERES
C     =============================================================
      DO 200 NL = NUMILF, NUMXLF
         NBAR = 0
         L    = L1LGFR(NL)
 105     IF( L .GT. 0 ) THEN
C           PASSAGE A L'ARETE SUIVANTE
            L    = NSARCH(4,L)
            NBAR = NBAR + 1
            GOTO 105
         ENDIF
         IF( NBAR .EQ. 1 ) THEN
C
C           UNE ARETE SUR LA LIGNE FRONTALIERE NL
            L   = L1LGFR(NL)
            NS1 = NSARCH(1,L)
            NS2 = NSARCH(2,L)
C
C           RECHERCHE DE CETTE ARETE DANS LE MAILLAGE
            DO 150 NT=1,MXTRIA
               IF( NOTRIA(1,NT) .GT. 0 ) THEN
C                 LE TRIANGLE EXISTE
                  DO 140 KK=1,3
C                    LE TRIANGLE OPPOSE A L'ARETE KK
                     NTOP = NOTRIA(3+KK,NT)
                     IF( NTOP .EQ. 0 ) THEN
C                       ARETE FRONTALIERE CALCUL DU PERIMETRE DU TRIANGLE NT
                        IF( KK .LT. 3 ) THEN
                           K1 = KK + 1
                        ELSE
                           K1 = 1
                        ENDIF
                        IF( NOTRIA(KK,NT) .EQ. NS1 ) THEN
                           IF( NOTRIA(K1,NT) .EQ. NS2 ) THEN
C                             ARETE RETROUVEE
                              GOTO 160
                           ENDIF
                        ENDIF
                        IF( NOTRIA(KK,NT) .EQ. NS2 ) THEN
                           IF( NOTRIA(K1,NT) .EQ. NS1 ) THEN
C                             ARETE RETROUVEE
                              GOTO 160
                           ENDIF
                        ENDIF
                     ENDIF
 140              CONTINUE
               ENDIF
 150        CONTINUE
C           ARETE NON RETROUVEE
            GOTO 200
C
C           ARETE RETROUVEE
C           AJOUT DU MILIEU DE L'ARETE FRONTALIERE (KK DANS NT)
 160        IF( K1 .LT. 3 ) THEN
               K2 = K1 + 1
            ELSE
               K2 = 1
            ENDIF
C           LES 3 SOMMETS DE NT
            NS1 = NOTRIA(KK, NT)
            NS2 = NOTRIA(K1, NT)
            NS3 = NOTRIA(K2, NT)
            IF( NBSOMM .GE. MXSOMM ) THEN
               NBLGRC(NRERR) = 1
                KERR(1) = 'TABLEAU PXYD SATURE'
               CALL LEREUR
               IERR = 1
               RETURN
            ENDIF
            NBSOMM = NBSOMM + 1
            PXYD(1,NBSOMM) = ( PXYD(1,NS1)+PXYD(1,NS2) )*0.5D0
            PXYD(2,NBSOMM) = ( PXYD(2,NS1)+PXYD(2,NS2) )*0.5D0
            PXYD(3,NBSOMM) = ( PXYD(3,NS1)+PXYD(3,NS2) )*0.5D0
C           CE NOUVEAU POINT N'EST PAS SUR UNE ISO
            NUISOP(NBSOMM) = 0
C
C           GESTION DE L'ARETE FRONTALIERE
            NUPLIS(NBSOMM) = NL
            L = L1LGFR(NL)
C           MODIFICATION DE L'ARETE EN ARETES NS1-NBSOMM + NBSOMM-NS2
            NSARCH(1,L) = NS1
            NSARCH(2,L) = NBSOMM
C           AJOUT DE LA NOUVELLE
            CALL AJARCH( NBSOMM, NS2, L1LGFR(NL),
     %                   MXARCH, NBARCH, N1ARCH, NSARCH, IERR )
            IF( IERR .NE. 0 ) GOTO 200
C
C           DECOUPAGE DU TRIANGLE NT EN 2 TRIANGLES
            NT1 = NOTRIA(3+K1,NT)
            NT2 = NOTRIA(3+K2,NT)
            IF( N1TRVI .LE. 0 ) THEN
C              SATURATION DES TRIANGLES
               IERR = 3
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES TRIANGLES'
               CALL LEREUR
               IERR = 2
               GOTO 200
            ENDIF
C           LE TRIANGLE A AJOUTER
            NTOP = N1TRVI
C           LE NOUVEAU TRIANGLE VIDE
            N1TRVI = NOTRIA(4,N1TRVI)
C
C           MISE A JOUR DU TRIANGLE NT
            NOTRIA(KK,  NT) = NBSOMM
            NOTRIA(3+K2,NT) = NTOP
C
C           MISE A JOUR DU TRIANGLE NTOP
            NOTRIA(1, NTOP) = NS1
            NOTRIA(2, NTOP) = NBSOMM
            NOTRIA(3, NTOP) = NS3
            NOTRIA(4, NTOP) = 0
            NOTRIA(5, NTOP) = NT
            NOTRIA(6, NTOP) = NT2
            NOSUTR(NTOP) = NOSUTR(NT)
C
C           MISE A JOUR DU TRIANGLE NT2
            IF( NT2 .GT. 0 ) THEN
               DO 180 L=4,6
                  IF( NOTRIA(L,NT2) .EQ. NT ) THEN
                     NOTRIA(L,NT2) = NTOP
                     GOTO 185
                  ENDIF
 180           CONTINUE
            ENDIF
C
C           MISE A JOUR DE NOTRSO
 185        NOTRSO(NBSOMM) = NT
            NOTRSO(NS1   ) = NTOP
C
C           UNE ARETE FRONTALIERE DECOUPEE DE PLUS
            NBARDE = NBARDE + 1
            NBTRIA = NBTRIA + 1
C
C           TRACE DES 2 TRIANGLES
            CALL DVTRTR(PXYD, NOTRIA, NT,   NCMAGE, NCBLAN )
            CALL DVTRTR(PXYD, NOTRIA, NTOP, NCMAGE, NCBLAN )
         ENDIF
 200  CONTINUE
      END
