      SUBROUTINE DVAMQS( ITER,   NBITAQ, MXSOMM, NBSOMM, NLSOFR,
     %                   MXTRIA, NDTRIA, N1TRVI, NOTRIA, CETRIA, NOTRSO,
     %                   MXTRNS, NOTRNS, NOARNS,
     %                   PXYD   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    UNE ITERATION DE BARYCENTRAGE DES POINTS INTERNES
C -----    MODIFICATION DE LA TOPOLOGIE POUR AVOIR AU MIEUX 6 VOISINS
C          POUR CHAQUE SOMMET DE LA TRIANGULATION
C
C ENTREES:
C --------
C ITER   : NUMERO DE L'ITERATION D'AMELIORATION
C NBITAQ : NOMBRE D'ITERATIONS D'AMELIORATION DE LA QUALITE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C NLSOFR : NUMERO DE LA LIGNE FERMEE (1 A NBLFTR) DE CHAQUE SOMMET
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE
C          0 SI LE POINT EST INTERNE OU EXTERNE
C N1TRVI : POINTE DANS NOTRIA VERS LE PREMIER TRIANGLE VIDE
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NOTRIA
C NDTRIA : NUMERO DU PLUS GRAND TRIANGLE UTILISE DANS NOTRIA
C NOTRIA : LISTE CHAINEE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                           ADJACENT PAR L'ARETE I
C CETRIA : COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT ET
C          CARRE DU RAYON
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
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1994
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
     %                  NLSOFR(MXSOMM)
      DOUBLE PRECISION  PXYD(3,MXSOMM),
     %                  CETRIA(1:3,1:MXTRIA)
      INTEGER           NOSOQU(5),NOTRNS(MXTRNS),NOARNS(MXTRNS),NTR(3)
      DOUBLE PRECISION  PONDER, PONDE1, XBAR, YBAR, SURTD2
      REAL              COMXMI(2,2)
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES POUR LE TRACE
      COMXMI(1,1) =  1E28
      COMXMI(2,1) =  1E28
      COMXMI(1,2) = -1E28
      COMXMI(2,2) = -1E28
C
C     COEFFICIENT DE PONDERATION CROISSANT AVEC LES ITERATIONS
      PONDER = MIN( 1D0, ( 50 + (50*ITER)/NBITAQ ) * 0.01D0 )
      PONDE1 = 1D0 - PONDER
      IERR   = 0
C
      NBSOM0 = NBSOMM
      DO 1000 NS=1,NBSOM0
C        NT UN TRIANGLE CONTENANT LE SOMMET NS
         NT = NOTRSO(NS)
         IF( NT .GT. 0 )  THEN
            IF( NOTRIA(1,NT) .GT. 0 )  THEN
C           LE POINT NS EST INTERNE AU MAILLAGE
            IF( NLSOFR(NS) .EQ. 0 ) THEN
C
C              LE POINT N'EST PAS FRONTALIER NI INTERNE IMPOSE PAR L'UTILISATEUR
               COMXMI(1,1) =  1E28
               COMXMI(2,1) =  1E28
               COMXMI(1,2) = -1E28
               COMXMI(2,2) = -1E28
C
C              PARCOURS DES TRIANGLES DE SOMMET NS
C              -----------------------------------
C              LE NUMERO DES SOMMETS DE NT
               DO 10 I=1,3
                  IF( NOTRIA(I,NT) .EQ. NS ) GOTO 20
   10          CONTINUE
   20          J1 = I + 1
               IF( J1 .EQ. 4 ) J1 = 1
C
C              LE TRIANGLE DE L'AUTRE COTE
               NT1 = NOTRIA(I+3,NT)
               NS1 = NOTRIA(J1 ,NT)
C
C              LA CONTRIBUTION DE NS1
               NBS  = 1
               XBAR = PXYD(1,NS1)
               YBAR = PXYD(2,NS1)
               NOTRNS( 1 ) = NT
               NOARNS( 1 ) = J1
C
C              PARCOURS PAR L'ARETE NS NS1
   25          IF( NT1 .NE. NT ) THEN
C
C        RECHERCHE DU MIN ET MAX DES COORDONNEES POUR LE TRACE
         DO 8 J=1,3
            NS3 = NOTRIA(J,NT1)
            DO 6 L=1,2
               D = REAL( PXYD(L,NS3) )
               IF( D .GT. COMXMI(L,2) ) COMXMI(L,2) = D
               IF( D .LT. COMXMI(L,1) ) COMXMI(L,1) = D
 6          CONTINUE
 8       CONTINUE
C
C                 RECHERCHE DE L'ARETE NS-NS1
                  DO 30 J=1,3
                     IF( NOTRIA(J,NT1) .EQ. NS1 ) GOTO 40
   30             CONTINUE
C
C                 L'ARETE SUIVANTE DE SOMMET NS
   40             J1 = J + 1
                  IF( J1 .EQ. 4 ) J1=1
                  J2 = J1 + 1
                  IF( J2 .EQ. 4 ) J2=1
                  NS2 = NOTRIA(J2,NT1)
C
C                 LA CONTRIBUTION AU BARYCENTRE DE NS2
                  NBS  = NBS + 1
                  IF( NBS .GT. MXTRNS ) GOTO 1000
                  NOTRNS( NBS ) = NT1
                  NOARNS( NBS ) = J2
                  XBAR = XBAR + PXYD(1,NS2)
                  YBAR = YBAR + PXYD(2,NS2)
C
C                 PASSAGE AU TRIANGLE SUIVANT
                  NT1 = NOTRIA(J1+3,NT1)
                  IF( NT1 .LE. 0 ) THEN
                     WRITE(IMPRIM,*) 'DVAMQS:ANOMALIE A CORRIGER'
                     GOTO 1000
                  ENDIF
                  NS1 = NS2
                  GOTO 25
               ENDIF
C
C              LE BARYCENTRE
C              =============
               XBAR = XBAR / NBS
               YBAR = YBAR / NBS
C
C        LE CADRE OBJET GLOBAL EN UNITES UTILISATEUR
      IF( TRATRI ) THEN
         X = COMXMI(1,2) - COMXMI(1,1)
         Y = COMXMI(2,2) - COMXMI(2,1)
         CALL ISOFENETRE( COMXMI(1,1)-X/20, COMXMI(1,2)+X/20,
     %                    COMXMI(2,1)-Y/20, COMXMI(2,2)+Y/20 )
         DO 15 I=1,NBS
            CALL DVTRTR( PXYD, NOTRIA, NOTRNS(I), NCBLAN, NCNOIR )
 15      CONTINUE
      ENDIF
C
C              PONDERATION POUR EVITER LES DEGENERESCENSES
               PXYD(1,NS) = PONDE1 * PXYD(1,NS) + PONDER * XBAR
               PXYD(2,NS) = PONDE1 * PXYD(2,NS) + PONDER * YBAR
C
      IF( TRATRI ) THEN
         DO 17 I=1,NBS
            CALL DVTRTR( PXYD, NOTRIA, NOTRNS(I), NCNOIR, NCBLAN )
 17      CONTINUE
      ENDIF
C
C              REMISE A JOUR DU RAYON CIRCONSCRIT DES TRIANGLES
               DO 45 J=1,NBS
                  IER = 1
                  J1  = NOTRNS( J )
                  CALL CENCED( PXYD(1,NOTRIA(1,J1)),
     %                         PXYD(1,NOTRIA(2,J1)),
     %                         PXYD(1,NOTRIA(3,J1)), CETRIA(1,J1), IER )
                  IF( IER .NE. 0 ) WRITE(IMPRIM,*) '45: TRIANGLE ',J1
 45            CONTINUE
C
               IF( NBS .EQ. 3 ) THEN
C
C                 3 TRIANGLES DEVIENNENT 1 TRIANGLE
C                 =================================
C                 NOTRNS(1:3) NUMEROS DES 3 TRIANGLES DE SOMMET COMMUN NS
C                 NOARNS(1:3) NUMERO DANS LE TRIANGLE DE L'ARETE OPPOSEE AU SOMM
C                 NOSOQU(1:3) NUMERO DES 3 SOMMETS DU TRIANGLE A FORMER
C                 LES 3 TRIANGLES SONT RANGES SELON LE SENS DIRECT
                  J         = NOTRNS(2)
                  NOTRNS(2) = NOTRNS(3)
                  NOTRNS(3) = J
                  J         = NOARNS(2)
                  NOARNS(2) = NOARNS(3)
                  NOARNS(3) = J
C
                  DO 70 J=1,3
                     NOSOQU(J) = NOTRIA( NOARNS(J), NOTRNS(J) )
 70               CONTINUE
C
C                 LE NOUVEAU TRIANGLE ENSEMBLE DES 3 AUTRES TRIANGLES
                  NT1 = NOTRNS(1)
                  DO 73 J1 = 1 , 3
C                    LE TRIANGLE NT AU DELA DE NOTRNS(J1)
                     NT = NOTRIA( 3+NOARNS(J1), NOTRNS(J1) )
                     IF( NT .GT. 0 ) THEN
                        IF( NOTRIA(1,NT) .GT. 0 ) THEN
C                          LE SOMMET J1 + 1 DU TRIANGLE NT1
                           IF( J1 .LT. 3 ) THEN
                              J3 = J1 + 1
                           ELSE
                              J3 = 1
                           ENDIF
                           NS1 = NOSOQU(J3)
                           DO 71 J2 = 1 , 3
                              IF( NOTRIA(J2,NT) .EQ. NS1 ) GOTO 72
 71                        CONTINUE
 72                        NOTRIA(3+J2,NT) = NT1
                        ELSE
                           NT = 0
                        ENDIF
                     ELSE
                        NT = 0
                     ENDIF
                     NTR(J1) = NT
 73               CONTINUE
C
                  DO 75 J=1,3
C                    LE NUMERO DU SOMMET J DU TRIANGLE NT1
                     NOTRIA(   J, NT1 ) = NOSOQU( J )
C                    LE NUMERO DU TRIANGLE AU DELA DU TRIANGLE NT1
                     NOTRIA( 3+J, NT1 ) = NTR( J )
C                    LE NUMERO DU TRIANGLE DE SOMMET NOSOQU(J)
                     NOTRSO( NOSOQU(J) ) = NT1
 75               CONTINUE
C
C                 LE CENTRE ET RAYON DU CERCLE CIRCONSCRIT A NT1
                  IER = 1
                  CALL CENCED( PXYD(1,NOTRIA(1,NT1)),
     %                         PXYD(1,NOTRIA(2,NT1)),
     %                         PXYD(1,NOTRIA(3,NT1)),
     %                         CETRIA(1,NT1), IER )
                  IF( IER .NE. 0 ) WRITE(IMPRIM,*) '75: TRIANGLE ',NT1
C
C                 SUPPRESSION DANS NOTRIA DES 2 TRIANGLES DERNIERS TRIANGLES
                  NOTRIA( 1, NOTRNS(2) ) = 0
                  NOTRIA( 1, NOTRNS(3) ) = 0
C
C                 LES 2 TRIANGLES NOTRNS(2:3) SONT LIBRES
                  NOTRIA(4,NOTRNS(2)) = N1TRVI
                  N1TRVI = NOTRNS(2)
                  NOTRIA(4,NOTRNS(3)) = N1TRVI
                  N1TRVI = NOTRNS(3)
C
C                 LE SOMMET NS EST PERDU
                  NOTRSO( NS ) = 0
C
      IF( TRATRI ) THEN
            CALL DVTRTR( PXYD, NOTRIA, NT1, NCVERT, NCJAUN )
      ENDIF
C
               ELSE IF( NBS .EQ. 4 ) THEN
C
C                 DECOUPAGE DES 4 TRIANGLES DE SOMMET COMMUN NS
C                 EN 2 TRIANGLES ENGLOBANTS
C                 =============================================
C                 LES 4 SOMMETS DE NOSOQU SONT RANGES SELON LE
C                 MOUVEMENT DES AIGUILLES D'UNE MONTRE
                  DO 50 J=1,4
                     NOSOQU(J) = NOTRIA( NOARNS(J), NOTRNS(J) )
 50               CONTINUE
C
C                 TEST POUR SAVOIR SI LE QUADRANGLE FORME A PARTIR
C                 DES 4 TRIANGLES EST CONVEXE OU NON
                  SURT0 = REAL(
     %                    ABS( SURTD2( PXYD(1,NOSOQU(1)),
     %                                 PXYD(1,NOSOQU(3)),
     %                                 PXYD(1,NOSOQU(2)) ) )
     %                  + ABS( SURTD2( PXYD(1,NOSOQU(1)),
     %                                 PXYD(1,NOSOQU(4)),
     %                                 PXYD(1,NOSOQU(3)) ) ) )
C
                  SURT1 = REAL(
     %                    ABS( SURTD2( PXYD(1,NOSOQU(1)),
     %                                 PXYD(1,NOSOQU(4)),
     %                                 PXYD(1,NOSOQU(2)) ) )
     %                  + ABS( SURTD2( PXYD(1,NOSOQU(2)),
     %                                 PXYD(1,NOSOQU(4)),
     %                                 PXYD(1,NOSOQU(3)) ) ) )
C
                   IF( ABS(SURT0-SURT1) .GT. 1E-4*SURT0 ) THEN
C
C                     QUADRANGLE 1432 EN FORME DE CHEVRON
C                     IL EST PREFERABLE DE TRAITER 3T=>1T
                      GOTO 1000
CCCC
CCCC                     PERMUTATION CIRCULAIRE POUR QUADRANGLE 2143
CCC                      J1 = NOSOQU(4)
CCC                      J2 = NOARNS(4)
CCC                      J3 = NOTRNS(4)
CCC                      DO 55 J=3,1,-1
CCC                         NOSOQU(J+1) = NOSOQU(J)
CCC                         NOARNS(J+1) = NOARNS(J)
CCC                         NOTRNS(J+1) = NOTRNS(J)
CCC 55                   CONTINUE
CCC                      NOSOQU(1) = J1
CCC                      NOARNS(1) = J2
CCC                      NOTRNS(1) = J3
                  ENDIF
C
C                 LES 2 TRIANGLES MODIFIES = 4 TRIANGLES ANCIENS
                  J1 = NOTRNS(1)
                  J2 = NOTRNS(2)
C
C                 LE TRIANGLE AU DELA DE NOTRNS(1)
                  NT = NOTRIA(NOARNS(1)+3,NOTRNS(1))
                  IF( NT .GT. 0 ) THEN
                     IF( NOTRIA(1,NT) .GT. 0 ) THEN
                        DO 61 J=1,3
                           IF( NOTRIA(J,NT) .EQ. NOSOQU(4) ) GOTO 62
 61                     CONTINUE
                        WRITE(IMPRIM,*) 'PB DVAMQS 61: ST ',NOSOQU(4)
                        GOTO 1000
C
 62                     NOTRIA(J+3,NT) = J2
                     ENDIF
                  ELSE
                     NT = 0
                  ENDIF
C
C                 LE TRIANGLE AU DELA DE NOTRNS(2)
                  NT1 = NOTRIA(NOARNS(2)+3,NOTRNS(2))
                  IF( NT1 .GT. 0 ) THEN
                     IF( NOTRIA(1,NT1) .GT. 0 ) THEN
                        DO 63 J=1,3
                           IF( NOTRIA(J,NT1) .EQ. NOSOQU(1) ) GOTO 64
 63                     CONTINUE
                        WRITE(IMPRIM,*) 'PB DVAMQS 63: ST ',NOSOQU(1)
                        GOTO 1000
C
 64                     NOTRIA(J+3,NT1) = J1
                     ENDIF
                  ELSE
                     NT1 = 0
                  ENDIF
C
C                 LE TRIANGLE AU DELA DE NOTRNS(3)
                  NT2 = NOTRIA(NOARNS(3)+3,NOTRNS(3))
                  IF( NT2 .GT. 0 ) THEN
                     IF( NOTRIA(1,NT2) .GT. 0 ) THEN
                        DO 65 J=1,3
                           IF( NOTRIA(J,NT2) .EQ. NOSOQU(2) ) GOTO 66
 65                     CONTINUE
                        WRITE(IMPRIM,*) 'PB DVAMQS 65: ST ',NOSOQU(2)
                        GOTO 1000
C
 66                     NOTRIA(J+3,NT2) = J1
                     ENDIF
                  ELSE
                     NT2 = 0
                  ENDIF
C
C                 LE TRIANGLE AU DELA DE NOTRNS(4)
                  NT3 = NOTRIA(NOARNS(4)+3,NOTRNS(4))
                  IF( NT3 .GT. 0 ) THEN
                     IF( NOTRIA(1,NT3) .GT. 0 ) THEN
                        DO 67 J=1,3
                           IF( NOTRIA(J,NT3) .EQ. NOSOQU(3) ) GOTO 68
 67                     CONTINUE
                        WRITE(IMPRIM,*) 'PB DVAMQS 67: ST ',NOSOQU(3)
                        GOTO 1000
C
 68                     NOTRIA(J+3,NT3) = J2
                     ENDIF
                  ELSE
                     NT3 = 0
                  ENDIF
C
C                 LES TRIANGLES AU DELA DES TRIANGLES J1 et J2
                  NOTRIA(4,J1) = J2
                  NOTRIA(5,J1) = NT2
                  NOTRIA(6,J1) = NT1
                  NOTRIA(4,J2) = NT
                  NOTRIA(5,J2) = NT3
                  NOTRIA(6,J2) = J1
C
C                 LES SOMMETS DES 2 TRIANGLES
                  NOTRIA(1,J1) = NOSOQU(1)
                  NOTRIA(2,J1) = NOSOQU(3)
                  NOTRIA(3,J1) = NOSOQU(2)
C
                  NOTRIA(1,J2) = NOSOQU(1)
                  NOTRIA(2,J2) = NOSOQU(4)
                  NOTRIA(3,J2) = NOSOQU(3)
C
C                 LE CENTRE ET RAYON DU CERCLE CIRCONSCRIT A J1
                  IER = 1
                  CALL CENCED( PXYD(1,NOTRIA(1,J1)),
     %                         PXYD(1,NOTRIA(2,J1)),
     %                         PXYD(1,NOTRIA(3,J1)), CETRIA(1,J1), IER )
                  IF( IER .NE. 0 ) WRITE(IMPRIM,*) '68: TRIANGLE ',J1
C
C                 LE CENTRE ET RAYON DU CERCLE CIRCONSCRIT A J2
                  IER = 1
                  CALL CENCED( PXYD(1,NOTRIA(1,J2)),
     %                         PXYD(1,NOTRIA(2,J2)),
     %                         PXYD(1,NOTRIA(3,J2)), CETRIA(1,J2), IER )
                  IF( IER .NE. 0 ) WRITE(IMPRIM,*) '69: TRIANGLE ',J2
C
C                 LE SOMMET NS EST PERDU
                  NOTRSO( NS ) = 0
C                 LE NUMERO DE TRIANGLE DES 4 SOMMETS DE J1 J2
                  NOTRSO( NOSOQU(1) ) = J1
                  NOTRSO( NOSOQU(2) ) = J1
                  NOTRSO( NOSOQU(3) ) = J2
                  NOTRSO( NOSOQU(4) ) = J2
C
C                 SUPPRESSION DES TRIANGLES NOTRNS(3:4)
                  NOTRIA(1,NOTRNS(3)) = 0
                  NOTRIA(1,NOTRNS(4)) = 0
C
C                 LES 2 TRIANGLES NOTRNS(3:4) SONT LIBRES
                  NOTRIA(4,NOTRNS(3)) = N1TRVI
                  N1TRVI = NOTRNS(3)
                  NOTRIA(4,NOTRNS(4)) = N1TRVI
                  N1TRVI = NOTRNS(4)
C
      IF( TRATRI ) THEN
            CALL DVTRTR( PXYD, NOTRIA, J1, NCBLEU, NCJAUN )
            CALL DVTRTR( PXYD, NOTRIA, J2, NCBLEU, NCJAUN )
      ENDIF
C
               ELSE IF( NBS .EQ. 5 ) THEN
C
C                 RECHERCHE DU CAS DE 2 SOMMETS D'UNE MEME ARETE
C                 CHACUN AYANT 5 SOMMETS VOISINS
C                 DANS CE CAS SUPPRESSION DE 2 TRIANGLES, 2 SOMMETS
C                 ET CREATION DU MILIEU DE L'ARETE
C                 =================================================
                  DO 110 J=1,5
                     NS1 = NOTRIA( NOARNS(J), NOTRNS(J) )
                     IF( NLSOFR(NS1) .EQ. 0 ) THEN
C                       NS1 EST INTERNE ET NON IMPOSE PAR L'UTILISATEUR
C                       AU DELA DE 16 TRIANGLES DE SOMMET NS1 RETOUR
                        CALL DVRETS( NS1,  NOTRIA, NOTRSO,
     %                               NBS1, 16, NOTRNS(6), NOARNS(6) )
                        IF( NBS1 .EQ. 5 ) GOTO 115
                     ENDIF
 110              CONTINUE
                  GOTO 1000
C
C                 CONFIGURATION RETROUVEE AVEC 2 SOMMETS DE 5 SOMMETS VOISINS
 115              NT = NOTRNS(J)
                  J  = NOARNS(J)
C                 RECHERCHE DE L'ARETE DE SOMMETS NS - NS1
                  IF( J .GT. 1 ) THEN
                     J2 = J - 1
                  ELSE
                     J2 = 3
                  ENDIF
C                 RECHERCHE DE L'ARETE DE SOMMETS NS1 - NS2
                  IF( J2 .GT. 1 ) THEN
                     J1 = J2 - 1
                  ELSE
                     J1 = 3
                  ENDIF
C
C                 LE TRIANGLE OPPOSE D'ARETE NS-NS1
                  NT1 = NOTRIA(3+J2,NT)
                  IF( NT1 .LE. 0 ) GOTO 1000
                  IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 1000
                  DO 120 I=1,3
                     IF( NOTRIA(I,NT1) .EQ. NS ) GOTO 125
 120              CONTINUE
                  WRITE(IMPRIM,*) 'PB DVAMQS 120: TR ',NT,' ',NT1
                  GOTO 1000
C
C                 LES 2 ARETES NS-NS3 ET NS3-NS1
 125              IF( I .GT. 1 ) THEN
                     I2 = I - 1
                  ELSE
                     I2 = 3
                  ENDIF
C                 RECHERCHE DE L'ARETE DE SOMMETS NS3-NS1
                  IF( I2 .GT. 1 ) THEN
                     I1 = I2 - 1
                  ELSE
                     I1 = 3
                  ENDIF
C
C                 RECHERCHE DES 2 ARETES DU TRIANGLE NT NON NS-NS1
                  NT2 = NOTRIA(3+J,NT)
                  IF( NT2 .LE. 0 ) GOTO 1000
                  IF( NOTRIA(1,NT2) .LE. 0 ) GOTO 1000
                  DO 130 K1=4,6
                      IF( NOTRIA(K1,NT2) .EQ. NT ) GOTO 135
 130              CONTINUE
                  WRITE(IMPRIM,*) 'PB DVAMQS 130: TR ',NT,' ',NT2
                  GOTO 1000
C
C
 135              NT3 = NOTRIA(3+J1,NT)
                  IF( NT3 .LE. 0 ) GOTO 1000
                  IF( NOTRIA(1,NT3) .LE. 0 ) GOTO 1000
                  DO 140 K2=4,6
                      IF( NOTRIA(K2,NT3) .EQ. NT ) GOTO 145
 140              CONTINUE
                  WRITE(IMPRIM,*) 'PB DVAMQS 140: TR ',NT,NT3
                  GOTO 1000
C
C                 RECHERCHE DES 2 ARETES DU TRIANGLE NT1 NON NS-NS1
 145              NT12 = NOTRIA(3+I,NT1)
                  IF( NT12 .LE. 0 ) GOTO 1000
                  IF( NOTRIA(1,NT12) .LE. 0 ) GOTO 1000
                  DO 150 K12=4,6
                      IF( NOTRIA(K12,NT12) .EQ. NT1 ) GOTO 155
 150              CONTINUE
                  WRITE(IMPRIM,*) 'PB DVAMQS 150: TR ',NT1,' ',NT12
                  GOTO 1000
C
 155              NT13 = NOTRIA(3+I1,NT1)
                  IF( NT13 .LE. 0 ) GOTO 1000
                  IF( NOTRIA(1,NT13) .LE. 0 ) GOTO 1000
                  DO 160 K13=4,6
                      IF( NOTRIA(K13,NT13) .EQ. NT1 ) GOTO 165
 160              CONTINUE
                  WRITE(IMPRIM,*) 'PB DVAMQS 160: TR ',NT1,' ',NT13
                  GOTO 1000
C
C                 ICI IL EXISTE BIEN 2 COUPLES (NT2,NT3) ET (NT12,NT13)
C                 LES 2 TRIANGLES NT ET NT1 PEUVENT DISPARAITRE
 165              NOTRIA(K1,NT2) = NT3
                  NOTRIA(K2,NT3) = NT2
C                 MISE A JOUR DE NOTRSO DU 3-EME SOMMET NS2 CAR NT DISPARAIT
                  NOTRSO( NOTRIA(J1,NT) ) = NT3
C
                  NOTRIA(K12,NT12) = NT13
                  NOTRIA(K13,NT13) = NT12
C                 MISE A JOUR DE NOTRSO DU 3-EME SOMMET NS3 CAR NT1 DISPARAIT
                  NOTRSO( NOTRIA(I1,NT1) ) = NT13
C
C                 MISE A JOUR DE NOTRSO DU SOMMET NS
                  NOTRSO( NS ) = NT13
C
C                 REMPLACEMENT DE NS1 PAR NS DANS LES TRIANGLES DE SOMMET NS1
                  DO 175 K1=6,10
                     NT2 = NOTRNS( K1 )
                     DO 170 K2=1,3
                        IF( NOTRIA(K2,NT2) .EQ. NS1 ) THEN
                           NOTRIA(K2,NT2) = NS
                           GOTO 175
                        ENDIF
 170                 CONTINUE
 175              CONTINUE
C
C                 MISE A JOUR DES COORDONNEES DE NS = MILIEU DE NS-NS1
                  DO 180 K1=1,3
                     PXYD(K1,NS) = ( PXYD(K1,NS) + PXYD(K1,NS1) ) * 0.5
 180              CONTINUE
C
C                 SUPPRESSION DANS NOTRIA DES 2 TRIANGLES NT ET NT1
                  NOTRIA(1,NT ) = 0
                  NOTRIA(1,NT1) = 0
C
C                 LES 2 TRIANGLES NT ET NT1 SONT LIBRES
                  NOTRIA(4,NT ) = N1TRVI
                  N1TRVI = NT
                  NOTRIA(4,NT1) = N1TRVI
                  N1TRVI = NT1
C
C                 MISE A JOUR DE CETRIA (ICI NOTRIA(1,NT et NT1)=0)
                  DO 190 K1=1,10
                     NT2 = NOTRNS(K1)
                     IF( NOTRIA(1,NT2) .LE. 0 ) GOTO 190
C
       IF( TRATRI ) THEN
            CALL DVTRTR( PXYD, NOTRIA, NT2, NCCYAN, NCJAUN )
       ENDIF
C                    LE CENTRE ET RAYON DU CERCLE CIRCONSCRIT A J2
                     IER = 1
                     CALL CENCED( PXYD(1,NOTRIA(1,NT2)),
     %                            PXYD(1,NOTRIA(2,NT2)),
     %                            PXYD(1,NOTRIA(3,NT2)),
     %                            CETRIA(1,NT2), IER )
                  IF( IER .NE. 0 ) WRITE(IMPRIM,*) '190: TRIANGLE ',NT2
 190              CONTINUE
C
C                 MISE A JOUR NOTRSO POUR NS1
                  NOTRSO( NS1 ) = 0
C
               ELSE IF( NBS .GE. 8 ) THEN
C
C                 AJOUTER UN TRIANGLE CENTRAL ABC + 3 ADJACENTS A LA
C                 PLACE DE L'ETOILE DES TRIANGLES DE SOMMET NS
C                 ==================================================
                  CALL DVAJ4T( NS,     NBS,    NOTRNS, NOARNS,
     %                         MXSOMM, NBSOMM, PXYD,   NLSOFR,
     %                         MXTRIA, NDTRIA, N1TRVI, NOTRIA,
     %                         CETRIA, NOTRSO, IERR )
                  IF( IERR .NE. 0 ) THEN
                     IERR = 0
                     GOTO 1000
                  ENDIF
C
               ENDIF
            ENDIF
            ENDIF
         ENDIF
 1000  CONTINUE
      END
