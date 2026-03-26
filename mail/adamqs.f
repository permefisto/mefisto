      SUBROUTINE ADAMQS( MXSOMM, NBSOMM, NUPLIS, NOSUTR, NUISOP,
     %                   MXTRIA, NOTRIA, NOTRSO,
     %                   MXTRNS, NOTRNS, NOARNS,
     %                   PXYD   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MODIFICATION DE LA TOPOLOGIE POUR AVOIR AU MIEUX 6 VOISINS
C -----    POUR CHAQUE SOMMET DE LA TRIANGULATION
C
C ENTREES:
C --------
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
C MXTRNS : >21 NOMBRE DE VARIABLES DU TABLEAU NOTRNS ( ET NOARNS )
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE ACTUEL DE SOMMETS ( INTERNES ET EXTERNES )
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C
C AUXILIAIRES:
C ------------
C NOTRNS, NOARNS : TABLEAUX DE MXTRNS ENTIERS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1995
C....................................................................012
      include"./incl/trvari.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION(CF DVTR2D)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NOTRSO(MXSOMM),
     %                  NUPLIS(MXSOMM),
     %                  NOSUTR(MXTRIA),
     %                  NUISOP(*)
      DOUBLE PRECISION  PXYD(3,MXSOMM)
      INTEGER           NOTRNS(MXTRNS),NOARNS(MXTRNS)
C
      IERR   = 0
C
      DO 1000 NS=1,NBSOMM
C
C        NT UN TRIANGLE CONTENANT LE SOMMET NS
         NT = NOTRSO(NS)
         IF( NT .GT. 0 )  THEN
C
            IF( NOTRIA(1,NT) .GT. 0 )  THEN
C
C              LE POINT NS EST INTERNE AU MAILLAGE
               IF( NUPLIS(NS) .EQ. 0 ) THEN
C
C                 PARCOURS DES TRIANGLES DE SOMMET NS
C                 -----------------------------------
C                 LE NUMERO DES SOMMETS DE NT
                  DO 10 I=1,3
                     IF( NOTRIA(I,NT) .EQ. NS ) GOTO 20
   10             CONTINUE
   20             J1 = I + 1
                  IF( J1 .EQ. 4 ) J1 = 1
C
C                 LE TRIANGLE DE L'AUTRE COTE
                  NT1 = NOTRIA(I+3,NT)
                  NS1 = NOTRIA(J1 ,NT)
C
C                 LA CONTRIBUTION DE NS1
                  NBS  = 1
                  NOTRNS( 1 ) = NT
                  NOARNS( 1 ) = J1
C
C                 PARCOURS PAR L'ARETE NS NS1
   25             IF( NT1 .NE. NT ) THEN
C                    RECHERCHE DE L'ARETE NS-NS1
                     DO 30 J=1,3
                        IF( NOTRIA(J,NT1) .EQ. NS1 ) GOTO 40
   30                CONTINUE
C
C                    L'ARETE SUIVANTE DE SOMMET NS
   40                J1 = J + 1
                     IF( J1 .EQ. 4 ) J1=1
                     J2 = J1 + 1
                     IF( J2 .EQ. 4 ) J2=1
                     NS2 = NOTRIA(J2,NT1)
C
C                    LA CONTRIBUTION AU BARYCENTRE DE NS2
                     NBS  = NBS + 1
                     IF( NBS .GT. MXTRNS ) GOTO 1000
                     NOTRNS( NBS ) = NT1
                     NOARNS( NBS ) = J2
C
C                    PASSAGE AU TRIANGLE SUIVANT
                     NT1 = NOTRIA(J1+3,NT1)
                     IF( NT1 .LE. 0 ) THEN
                        WRITE(IMPRIM,*) 'ADAMQS:ANOMALIE A CORRIGER'
                        GOTO 1000
                     ENDIF
                     NS1 = NS2
                     GOTO 25
                  ENDIF
C
C                 VERIFICATION QUE TOUS LES TRIANGLES APPARTIENNENT A LA MEME SU
                  NOSURF = NOSUTR( NOTRNS(1) )
                  DO 16 I=2,NBS
                     IF( NOSUTR( NOTRNS(I) ) .NE. NOSURF ) GOTO 1000
 16               CONTINUE
C
C                 BARYCENTRAGE EVENTUEL DU POINT NS SUR L'ISOVALEUR
                  IF( NUISOP(NS) .GT. 0 ) THEN
                     X  = REAL( PXYD(1,NS) )
                     Y  = REAL( PXYD(2,NS) )
                     NA = 0
                     DO 18 I=1,NBS
                        NS1 = NOTRIA( NOARNS(I), NOTRNS(I) )
                        IF( NUISOP(NS1) .EQ. NUISOP(NS) ) THEN
                           NA = NA + 1
                           X  = REAL( X + PXYD(1,NS1) )
                           Y  = REAL( Y + PXYD(2,NS1) )
                        ENDIF
 18                  CONTINUE
                     IF( NA .EQ. 2 ) THEN
                        PXYD(1,NS) = X / 3.0
                        PXYD(2,NS) = Y / 3.0
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 1000 CONTINUE
      END
