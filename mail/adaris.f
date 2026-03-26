      SUBROUTINE ADARIS( TARMAX, MXSOMM, NBSOMM, NUPLIS, NOSUTR, NUISOP,
     %                   MXTRIA, N1TRVI, NOTRIA, NOTRSO,
     %                   PXYD  , NBTRIA, NBARDE, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUT DU MILIEU DES ARETES ISOVALEURS DE TAILLE
C -----    SUPERIEURE A LA TAILLE MAXIMALE PERMISE
C
C ENTREES:
C --------
C TARMAX : TAILLE DE L'ARETE MAXIMALE PERMISE SUR UNE FRONTIERE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C NUPLIS : NUMERO DE POINT OU LIGNE OU ISO DE CHACUN DES SOMMETS
C          -2 000 000 SI LE SOMMET A DEJA ETE SUPPRIME
C          -NP SI NP EST LE NUMERO DU POINT UTILISATEUR DE CE SOMMET
C          -1 234 567 SI LE SOMMET APPARTIENT A 2 LIGNES (SOMMET INITIAL)
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
     %                  NUISOP(*)
      DOUBLE PRECISION  PXYD(3,MXSOMM), RLONG
C
      IERR   = 0
      NBARDE = 0
      DO 100 NT=1,MXTRIA
         IF( NOTRIA(1,NT) .GT. 0 ) THEN
C        LE TRIANGLE EXISTE
            DO 90 KK=1,3
               NTOP = NOTRIA(3+KK,NT)
               IF( NTOP .GT. 0 ) THEN
C                 ARETE NON FRONTALIERE
                  IF( KK .LT. 3 ) THEN
                     K1 = KK + 1
                  ELSE
                     K1 = 1
                  ENDIF
                  NS1 = NOTRIA(KK,NT)
                  NS2 = NOTRIA(K1,NT)
C
                  IF( NUISOP(NS1) .EQ. NUISOP(NS2) .AND.
     %                NUISOP(NS1) .GT. 0 ) THEN
C
C                    ARETE SUR L'ISO NUPLIS(NS1)
                     RLONG = SQRT( (PXYD(1,NS2)-PXYD(1,NS1))**2
     %                           + (PXYD(2,NS2)-PXYD(2,NS1))**2 )
                     IF( RLONG .GT. TARMAX ) THEN
C
C                       AJOUT DU MILIEU DE L'ARETE FRONTALIERE (KK DANS NT)
                        IF( K1 .LT. 3 ) THEN
                           K2 = K1 + 1
                        ELSE
                           K2 = 1
                        ENDIF
C                       LE 3-EME SOMMET DU TRIANGLE NT
                        NS3 = NOTRIA(K2, NT)
C
C                       LES TRIANGLES OPPOSES AUX ARETES DE NTOP
                        DO 20 KL=1,3
                           IF( NOTRIA(KL,NTOP) .EQ. NS2 ) THEN
                              IF( KL .LT. 3 ) THEN
                                 K3 = KL + 1
                              ELSE
                                 K3 = 1
                              ENDIF
                              IF( NOTRIA(K3,NTOP) .EQ. NS1 ) THEN
C                                ARETE RETROUVEE
                                 GOTO 30
                              ENDIF
                           ENDIF
 20                     CONTINUE
                        GOTO 90
C
 30                     IF( K3 .LT. 3 ) THEN
                           K4 = K3 + 1
                        ELSE
                           K4 = 1
                        ENDIF
                        NS4 = NOTRIA(K4,NTOP)
C
                        IF( NBSOMM .GE. MXSOMM ) THEN
                           NBLGRC(NRERR) = 1
                           KERR(1) = 'TABLEAU PXYD SATURE'
                           CALL LEREUR
                           IERR = 1
                           RETURN
                        ENDIF
                        NBSOMM = NBSOMM + 1
                        PXYD(1,NBSOMM) = (PXYD(1,NS1)+PXYD(1,NS2))*0.5D0
                        PXYD(2,NBSOMM) = (PXYD(2,NS1)+PXYD(2,NS2))*0.5D0
                        PXYD(3,NBSOMM) = (PXYD(3,NS1)+PXYD(3,NS2))*0.5D0
C
C                       CE NOUVEAU POINT EST SUR L'ISO NUISOP(NS1)
                        NUISOP(NBSOMM) = NUISOP(NS1)
C
C                       ARETE INTERNE
                        NUPLIS(NBSOMM) = 0
C
C                       DECOUPAGE DU TRIANGLE NT   EN 2 TRIANGLES
C                       DECOUPAGE DU TRIANGLE NTOP EN 2 TRIANGLES
C                       LES TRIANGLES OPPOSES
                        NT1 = NOTRIA(3+K1,NT)
                        NT2 = NOTRIA(3+K2,NT)
                        NT3 = NOTRIA(3+K3,NTOP)
                        NT4 = NOTRIA(3+K4,NTOP)
                        IF( N1TRVI .LE. 0 ) THEN
C                          SATURATION DES TRIANGLES
                           IERR = 3
                           NBLGRC(NRERR) = 1
                           KERR(1) = 'SATURATION DES TRIANGLES'
                           CALL LEREUR
                           GOTO 90
                        ENDIF
C                       LE 1-ER TRIANGLE A AJOUTER
                        NTPL = N1TRVI
C                       LE NOUVEAU TRIANGLE VIDE
                        N1TRVI = NOTRIA(4,N1TRVI)
                        IF( N1TRVI .LE. 0 ) THEN
C                          SATURATION DES TRIANGLES
                           IERR = 3
                           NBLGRC(NRERR) = 1
                           KERR(1) = 'SATURATION DES TRIANGLES'
                           CALL LEREUR
                           GOTO 90
                        ENDIF
C                       LE SECOND TRIANGLE A AJOUTER
                        NTPP = N1TRVI
C                       LE NOUVEAU TRIANGLE VIDE
                        N1TRVI = NOTRIA(4,N1TRVI)
C
C                       MISE A JOUR DU TRIANGLE NT
                        NOTRIA(K1,  NT) = NBSOMM
                        NOTRIA(3+K1,NT) = NTPL
C
C                       MISE A JOUR DU TRIANGLE NTOP
                        NOTRIA(KL,  NTOP) = NBSOMM
                        NOTRIA(3+K4,NTOP) = NTPP
C
C                       MISE A JOUR DU TRIANGLE NTPL
                        NOTRIA(1, NTPL) = NBSOMM
                        NOTRIA(2, NTPL) = NS2
                        NOTRIA(3, NTPL) = NS3
                        NOTRIA(4, NTPL) = NTPP
                        NOTRIA(5, NTPL) = NT1
                        NOTRIA(6, NTPL) = NT
                        NOSUTR(NTPL) = NOSUTR(NT)
C
C                       MISE A JOUR DU TRIANGLE NTPP
                        NOTRIA(1, NTPP) = NBSOMM
                        NOTRIA(2, NTPP) = NS4
                        NOTRIA(3, NTPP) = NS2
                        NOTRIA(4, NTPP) = NTOP
                        NOTRIA(5, NTPP) = NT4
                        NOTRIA(6, NTPP) = NTPL
                        NOSUTR(NTPP) = NOSUTR(NTOP)
C
C                       MISE A JOUR DU TRIANGLE NT1
                        IF( NT1 .GT. 0 ) THEN
                           DO 40 L=4,6
                              IF( NOTRIA(L,NT1) .EQ. NT ) THEN
                                 NOTRIA(L,NT1) = NTPL
                                 GOTO 45
                              ENDIF
 40                        CONTINUE
                        ENDIF
C
C                       MISE A JOUR DU TRIANGLE NT4
 45                     IF( NT4 .GT. 0 ) THEN
                           DO 50 L=4,6
                              IF( NOTRIA(L,NT4) .EQ. NTOP ) THEN
                                 NOTRIA(L,NT4) = NTPP
                                 GOTO 55
                              ENDIF
 50                        CONTINUE
                        ENDIF
C
C                       MISE A JOUR DE NOTRSO
 55                     NOTRSO(NBSOMM) = NT
                        NOTRSO(NS2   ) = NTPL
C
C                       UNE ARETE ISOVALEUR DECOUPEE DE PLUS
                        NBARDE = NBARDE + 1
                        NBTRIA = NBTRIA + 2
C
C                       TRACE DES 4 TRIANGLES
                        CALL DVTRTR(PXYD, NOTRIA, NT,   NCCYAN, NCNOIR )
                        CALL DVTRTR(PXYD, NOTRIA, NTOP, NCCYAN, NCNOIR )
                        CALL DVTRTR(PXYD, NOTRIA, NTPL, NCCYAN, NCNOIR )
                        CALL DVTRTR(PXYD, NOTRIA, NTPP, NCCYAN, NCNOIR )
                     ENDIF
                  ENDIF
               ENDIF
 90         CONTINUE
         ENDIF
 100  CONTINUE
      END
