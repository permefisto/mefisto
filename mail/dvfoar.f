      SUBROUTINE DVFOAR( NS1ARE, NS2ARE,
     %                   PXYD  , NLSOFR,
     %                   N1TRVI, NOTRIA, NOTRSO,
     %                   MXETRI, NAETOI, NARMIN,
     %                   MXARCF, N1ARCF, NOARCF, NOTRCF,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   FORCER L'ARETE NS1-NS2 DANS LA TRIANGULATION ACTUELLE
C -----   TRIANGULATION FRONTALE POUR LA REOBTENIR
C
C ENTREES:
C --------
C NS1ARE,NS2ARE: NUMERO PXYD DES 2 SOMMETS EXTREMITES DE L'ARETE A FORCER
C NBLFTR : NOMBRE DE LIGNES FERMEES LIMITANT LA SURFACE A TRIANGULER
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C NLSOFR : NUMERO(1 A NBLFTR) DE LA LIGNE FERMEE DU POINT
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE
C          0 SI LE POINT EST INTERNE OU EXTERNE NON IMPOSE
C MXARCF : NOMBRE DE VARIABLES DES TABLEAUX N1ARCF, NOARCF,NOTRCF
C
C ENTREES ET SORTIES :
C --------------------
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                          ADJACENT PAR L'ARETE i
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C
C AUXILIAIRES :
C -------------
C NAETOI : TABLEAU (4,MXETRI) AUXILIAIRE
C NARMIN : TABLEAU (MXETRI)   AUXILIAIRE
C N1ARCF : TABLEAU (0:MXARCF) AUXILIAIRE
C NOARCF : TABLEAU (3,MXARCF) AUXILIAIRE
C NOTRCF : TABLEAU (1:MXARCF) AUXILIAIRE
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          2 NS1 DANS AUCUN TRIANGLE
C          9 TABLEAU NOSOAR DE TAILLE INSUFFISANTE CAR TROP D'ARETES
C            A PROBLEME
C          10 UN DES TABLEAUX N1ARCF, NOARCF NOTRCF EST SATURE
C             AUGMENTER A L'APPEL MXARCF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C....................................................................012
      PARAMETER        (MXPITR=32)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
      INTEGER           NOTRIA(6,*),
     %                  NLSOFR(*),
     %                  NOTRSO(*),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(3,MXARCF)
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NAETOI(4,MXETRI),
     %                  NARMIN(MXETRI),
     %                  NOTRCF(MXARCF)
      INTEGER           LAPITR(MXPITR)
      DOUBLE PRECISION  X1,Y1,X2,Y2,D12,D3,D4,X,Y,D,DMIN
      INTEGER           NS(2),NSS(1:2),NSS1,NSS2
      EQUIVALENCE      (NSS(1),NSS1),(NSS(2),NSS2)
C
C     TRAITEMENT DE CETTE ARETE PERDUE
      NS1  = NS1ARE
      NS2  = NS2ARE
C
      IF( TRATRI ) THEN
C        LES TRACES SONT DEMANDES
         CALL EFFACE
C        LE CADRE OBJET GLOBAL EN UNITES UTILISATEUR
         XX1 = REAL( MIN( PXYD(1,NS1), PXYD(1,NS2) ) )
         XX2 = REAL( MAX( PXYD(1,NS1), PXYD(1,NS2) ) )
         YY1 = REAL( MIN( PXYD(2,NS1), PXYD(2,NS2) ) )
         YY2 = REAL( MAX( PXYD(2,NS1), PXYD(2,NS2) ) )
         IF( XX1 .GE. XX2 ) XX2 = XX1 + (YY2-YY1)
         IF( YY1 .GE. YY2 ) YY2 = YY1 + (XX2-XX1)*0.5
         CALL ISOFENETRE( XX1-(XX2-XX1), XX2+(XX2-XX1),
     %                    YY1-(YY2-YY1), YY2+(YY2-YY1) )
      ENDIF
CCC
CCC      NSA1 = NS1
CCC      NSA2 = NS2
C
C     RECHERCHE DU TRIANGLE VOISIN DANS LE SENS INDIRECT
      NSENS = 0
C
C     TRACE DE L'ARETE PERDUE
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
C
C     RECHERCHE DES TRIANGLES INTERSECTES PAR LE SEGMENT NS1-NS2
      IPAS = 0
C
 3    X1  = PXYD(1,NS1)
      Y1  = PXYD(2,NS1)
      X2  = PXYD(1,NS2)
      Y2  = PXYD(2,NS2)
      D12 = (X2-X1)**2 + (Y2-Y1)**2
C
C     RECHERCHE DU NO LOCAL DU SOMMET NS1 DANS L'UN DE SES TRIANGLES
      NT0  = NOTRSO( NS1 )
      IF( NT0 .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'DVFOAR: SOMMET ',NS1,' SANS AUCUN TRIANGLE'
         IERR = 2
         RETURN
      ENDIF
      IF( NOTRSO(NS2) .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'DVFOAR: SOMMET ',NS2,' SANS AUCUN TRIANGLE'
         IERR = 2
         RETURN
      ENDIF
C
 20   DO 22 NA00=1,3
         IF( NOTRIA(NA00,NT0) .EQ. NS1 ) GOTO 26
 22   CONTINUE
C
 25   IF( IPAS .EQ. 0 ) THEN
C        TENTATIVE D'INVERSION DES 2 SOMMETS EXTREMITES DE L'ARETE A FORCER
         NA00 = NS1
         NS1  = NS2
         NS2  = NA00
         IPAS = 1
         GOTO 3
      ELSE
         WRITE(IMPRIM,*)'DVFOAR:ARETE ',NS1ARE,' - ',NS2ARE,' A IMPOSER'
         WRITE(IMPRIM,*)'DVFOAR:ANOMALIE SOMMET ',NS1,
     %   'NON DANS LE TRIANGLE DE SOMMETS ',(NOTRIA(I,NT0),I=1,3)
         NA00 = 1
      ENDIF
C
C     INTERSECTION AVEC L'ARETE OPPOSEE A NS1
 26   NA0 = NOSUI3( NA00 )
      NS3 = NOTRIA( NA0, NT0 )
      NA1 = NOPRE3( NA00 )
      NS4 = NOTRIA( NA1, NT0 )
C
C     TRACE DU TRIANGLE NT0 ET DE L'ARETE PERDUE
      CALL DVTRTR( PXYD, NOTRIA, NT0, NCBLAN, NCBLEU )
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
C
      CALL INT1SD( NS1, NS2, NS3, NS4, PXYD, LINTER, X1, Y1 )
      IF( LINTER .LE. 0 ) THEN
C        PAS D'INTERSECTION : ROTATION AUTOUR DU POINT NS1
         IF( NSENS .EQ. 0 ) THEN
            J = NA00
         ELSE
            IF( NA00 .NE. 1 ) THEN
               J = NA00 - 1
            ELSE
               J = 3
            ENDIF
         ENDIF
         NT0 = NOTRIA( 3+J, NT0 )
         IF( NT0 .GT. 0 ) GOTO 20
C        LE PARCOURS SORT DU DOMAINE
C        IL FAUT TOURNER DANS L'AUTRE SENS AUTOUR DE NS1
         IF( NSENS .EQ. 0 ) THEN
            NSENS = 1
            NT0   = NOTRSO( NS1 )
            GOTO 20
         ENDIF
C
C        DANS LES 2 SENS, PAS D'INTERSECTION => IMPOSSIBLE
C        ESSAI AVEC L'ARETE INVERSEE NS1 <-> NS2
         IF( IPAS .EQ. 0 ) GOTO 25
         WRITE(IMPRIM,*) 'DVFOAR: ARETE ',NS1,' ',NS2,
     %  ' SANS INTERSECTION AVEC LES TRIANGLES ACTUELS'
         WRITE(IMPRIM,*)'REVOYEZ LES LIGNES DU CONTOUR'
         IERR = 11
         RETURN
      ENDIF
C
C     IL EXISTE UNE INTERSECTION AVEC L'ARETE OPPOSEE
C     NBTRCF : NOMBRE DE TRIANGLES DU CF
      NBTRCF = 1
      NOTRCF( 1 ) = NT0
C
C     LE TRIANGLE OPPOSE A L'ARETE NA0 DE NT0
 28   NT1 = NOTRIA( 3+NA0, NT0 )
C
C     TRACE DU TRIANGLE NT1 ET DE L'ARETE PERDUE
      CALL DVTRTR( PXYD, NOTRIA, NT1, NCJAUN, NCMAGE )
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
C
C     LE TRIANGLE NT1 CONTIENT IL NS2 ?
      DO 32 J=1,3
         IF( NOTRIA(J, NT1 ) .EQ. NS2 ) GOTO 70
 32   CONTINUE
C
C     RECHERCHE DE L'ARETE NA1 DANS NT1 DE L'ARETE NA0 DE NT0
      IF( NOTRIA(4,NT1) .EQ. NT0 ) THEN
         NA1 = 1
      ELSE IF( NOTRIA(5,NT1) .EQ. NT0 ) THEN
         NA1 = 2
      ELSE
         NA1 = 3
      ENDIF
C
C     RECHERCHE DE L'INTERSECTION DE NS1-NS2 AVEC LES 2 AUTRES ARETES
C     TRACE DU TRIANGLE NT1 ET DE L'ARETE PERDUE
      CALL DVTRTR( PXYD, NOTRIA, NT1, NCJAUN, NCMAGE )
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
      NA2 = NA1
      DO 50 I1 = 1,2
C        L'ARETE SUIVANTE
         NA2 = NOSUI3(NA2)
C        LES 2 SOMMETS DE L'ARETE NA2 DE NT1
         NS3 = NOTRIA( NA2 , NT1 )
         NS4 = NOTRIA( NOSUI3(NA2), NT1 )
         CALL INT1SD( NS1, NS2, NS3, NS4, PXYD, LINTER, X , Y )
         IF( LINTER .GT. 0 ) THEN
C
C           LES 2 ARETES S'INTERSECTENT EN (X,Y)
C           DISTANCE DE (X,Y) A NS3 ET NS4
            D3=(PXYD(1,NS3)-X)**2+(PXYD(2,NS3)-Y)**2
            D4=(PXYD(1,NS4)-X)**2+(PXYD(2,NS4)-Y)**2
C           NSP EST LE POINT LE PLUS PROCHE DE (X,Y)
            IF( D3 .LT. D4 ) THEN
               NSP = NS3
               D   = D3
            ELSE
               NSP = NS4
               D   = D4
            ENDIF
            IF( D .GT. 1D-5*D12 ) GOTO 60
C
C           ICI LE SOMMET NSP EST TROP PROCHE DE NS1-NS2
            IF( NLSOFR(NSP) .NE. 0 ) THEN
C              POINT UTILISATEUR OU FRONTALIER NON SUPPRIMABLE
               GOTO 60
            ENDIF
C
C           LE POINT INTERNE NSP EST SUPPRIME EN METTANT TOUS LES TRIANGLES
C           L'AYANT COMME SOMMET DANS LA PILE DES TRIANGLES
            WRITE(IMPRIM,*) 'DVFOAR:SOMMET ',NSP,' SUPPRIME'
            CALL TRC1ST( NSP, NOTRSO(NSP), NOTRIA,
     %                   MXPITR, NBT, LAPITR )
C           AJOUT DES TRIANGLES A NOTRCF
            NBTRC0 = NBTRCF
            DO 38 J=1,ABS(NBT)
               NT = LAPITR(J)
               DO 35 K=NBTRCF,1,-1
                  IF( NT .EQ. NOTRCF(K) ) GOTO 38
 35            CONTINUE
C              TRIANGLE AJOUTE
               NBTRCF = NBTRCF + 1
               NOTRCF( NBTRCF ) = NT
               CALL DVTRTR( PXYD, NOTRIA, NT, NCJAUN, NCMAGE )
               CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
 38         CONTINUE
C           SOMMET SUPPRIME APPARTENANT A AUCUN TRIANGLE
            NOTRSO( NSP ) = 0
C
C           NS2 EST-IL UN SOMMET DES TRIANGLES EMPILES?
            DO 40 NT=NBTRC0+1,NBTRCF
               NT1 = NOTRCF( NT )
               DO 39 K=1,3
C                 LE SOMMET K DE NT1
                  IF( NOTRIA( K, NT1 ) .EQ. NS2 ) THEN
C                    BUT ATTEINT
                     GOTO 80
                  ENDIF
 39            CONTINUE
 40         CONTINUE
C
C           RECHERCHE DU PLUS PROCHE POINT D'INTERSECTION DE NS1-NS2
C           PAR RAPPORT A NS2 AVEC LES ARETES DES TRIANGLES AJOUTES
            NT0  = 0
            DMIN = D12 * 10000
            DO 48 NT=NBTRC0+1,NBTRCF
               NT1 = NOTRCF( NT )
               DO 42 K=1,3
C                 LES 2 SOMMETS DE L'ARETE K DE NT
                  NS3 = NOTRIA( K, NT1 )
                  NS4 = NOTRIA( NOSUI3(K), NT1 )
                  CALL INT1SD( NS1, NS2, NS3, NS4, PXYD,
     %                         LINTER, X , Y )
                  IF( LINTER .GT. 0 ) THEN
C                    LES 2 ARETES S'INTERSECTENT EN (X,Y)
                     D = (X-X2)**2+(Y-Y2)**2
                     IF( D .LT. DMIN ) THEN
                        NT0  = NT1
                        NA0  = K
                        DMIN = D
                     ENDIF
                  ENDIF
 42            CONTINUE
 48         CONTINUE
C
C           REDEMARRAGE AVEC LE TRIANGLE NT0 ET L'ARETE NA0
            IF( NT0 .GT. 0 ) GOTO 28
            WRITE(IMPRIM,*) 'DVFOAR: ALGORITHME DEFAILLANT'
            IERR = 12
            RETURN
         ENDIF
 50   CONTINUE
C
C     PAS D'INTERSECTION DIFFERENTE DE L'INITIALE => SOMMET SUR NS1-NS2
C     ROTATION AUTOUR DU SOMMET PAR L'ARETE SUIVANT NA1
CCC      WRITE(IMPRIM,*) 'DVFOAR 50: POURQUOI? A ECLAIRCIR'
      WRITE(IMPRIM,*) 'DVFOAR 50: REVOYEZ VOS DONNEES'
      WRITE(IMPRIM,*) 'LES LIGNES FERMEES DOIVENT ETRE DISJOINTES'
      WRITE(IMPRIM,*) 'VERIFIEZ SI C''EST LE CAS'
      IERR = 13
      NBTRCF = NBTRCF + 1
      NOTRCF( NBTRCF ) = NT1
      NT0 = NT1
      NA0 = NOSUI3( NA1 )
      GOTO 28
C
C     CAS SANS PROBLEME : INTERSECTION DIFFERENTE DE CELLE INITIALE
 60   NBTRCF = NBTRCF + 1
      NOTRCF( NBTRCF ) = NT1
C     PASSAGE AU TRIANGLE SUIVANT
      NA0 = NA2
      NT0 = NT1
      GOTO 28
C
C     ----------------------------------------------------------
C     ICI TOUTES LES INTERSECTIONS DE NS1-NS2 ONT ETE PARCOURUES
C     ----------------------------------------------------------
 70   NBTRCF = NBTRCF + 1
      NOTRCF( NBTRCF ) = NT1
C
C     FORMATION DES ARETES UNIQUES DU CF AUTOUR DE L'ARETE NS1-NS2
C     ------------------------------------------------------------
C     REINITIALISATION A VIDE DES ARETES DE L'ETOILE
C     FORMEE DES ARETES VUES UNE FOIS DANS LES TRIANGLES DE L'ETOILE
 80   N1AEVI = 1
      N1AEOC = 0
      MMETRI = MIN(8*NBTRCF,MXETRI)
      DO 90 I=1,MMETRI
C        NUMERO DANS NAETOI DE L'ARETE SUIVANTE
         NAETOI(4,I) = I+1
 90   CONTINUE
      NAETOI(4,MMETRI) = 0
C
      DO 98 I=1,NBTRCF
C        AJOUT OU RETRAIT DES 3 ARETES DU TRIANGLE NOTRCF(I) A L'ETOILE
         CALL AJTRET( NOTRCF(I), NOTRIA, N1AEVI, N1AEOC, NAETOI )
 98   CONTINUE
C
C     MODIFICATION DU CONTENU DU TABLEAU NAETOI
C     LE NO DE TRIANGLE => NO TRIANGLE AU DELA DE L'ARETE
C     LE NO LOCAL DANS LE TRIANGLE => NO 1-ER SOMMET DE L'ARETE
C     LE NO INUTILISE              => NO 2-ME SOMMET DE L'ARETE
C     ---------------------------------------------------------
      NA1 = N1AEOC
C     BOUCLE SUR LES ARETES
 110  IF( NA1 .GT. 0 ) THEN
C        LE NO DU TRIANGLE ET LOCAL DE L'ARETE
         NT   = NAETOI(1,NA1)
         I    = ABS(NAETOI(2,NA1))
C        LE NUMERO DU TRIANGLE AU DELA DE L'ARETE
         NTOP = NOTRIA(3+I,NT)
         NAETOI(1,NA1) = NTOP
C        LE NUMERO DES 2 SOMMETS DE L'ARETE I DU TRIANGLE NT
         IF( I .EQ. 3 ) THEN
            NS4 = 1
         ELSE
            NS4 = I + 1
         ENDIF
C        NUMERO DU SOMMET 1 DE L'ARETE
         NS3 = NOTRIA(I,NT)
         NAETOI(2,NA1) = NS3
C        NUMERO DU TRIANGLE CONTENANT ENCORE NS3
         NOTRSO( NS3 ) = NTOP
C        NUMERO DU SOMMET 2 DE L'ARETE
         NS4 = NOTRIA(NS4,NT)
         NAETOI(3,NA1) = NS4
C        NUMERO DU TRIANGLE CONTENANT ENCORE NS4
         NOTRSO( NS4 ) = NTOP
C        PASSAGE A L'ARETE SUIVANTE
         NA1 = NAETOI(4,NA1)
         GOTO 110
      ENDIF
C
C     LES ARETES SONT REORDONNEES POUR FORMER UNE LIGNE FERMEE
C     ========================================================
      NA1 = N1AEOC
C     LA PREMIERE ARETE
      NS0 = NAETOI(2,NA1)
      NS1 = NAETOI(3,NA1)
C
C     LE 1-ER SOMMET OU ARETE DU CONTOUR FERME
      N1ARCF( 1 ) = 1
C     LE NOMBRE DE SOMMETS DU CONTOUR FERME DE L'ETOILE
      NBARCF = 1
C     LE PREMIER SOMMET DE L'ETOILE
      NOARCF( 1, NBARCF ) = NS0
C     LE SOMMET SUIVANT
      NOARCF( 2, NBARCF ) = NBARCF + 1
C     LE NUMERO DU TRIANGLE DE L'AUTRE COTE DE CETTE ARETE
      NOARCF( 3, NBARCF ) = NAETOI(1,NA1)
C
C     TRACE DE L'ARETE
      CALL DVTRAR( PXYD, NS0, NS1, NCVERT, NCBLAN )
C
C     L'ARETE SUIVANTE
      NA1    = NAETOI(4,NA1)
      N1AEOC = NA1
C
 120  IF( N1AEOC .GT. 0 ) THEN
C
C        RECHERCHE DE L'ARETE DE 1-ER SOMMET NS1
         NA0 = N1AEOC
         NA1 = NA0
 160     IF( NA1 .GT. 0 ) THEN
C
C           LE NUMERO DU PREMIER SOMMET DE L'ARETE
            IF ( NS1 .NE. NAETOI(2,NA1) .AND.
     %           NS1 .NE. NAETOI(3,NA1) ) THEN
C              PASSAGE A L'ARETE SUIVANTE
               NA0 = NA1
               NA1 = NAETOI(4,NA1)
               GOTO 160
            ENDIF
C
C           ARETE PERIPHERIQUE RETROUVEE
            NBARCF = NBARCF + 1
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            IF( NS1 .EQ. NAETOI(2,NA1) ) THEN
               NS1 = NAETOI(3,NA1)
               NS2 = 2
            ELSE
               NS1 = NAETOI(2,NA1)
               NS2 = 3
            ENDIF
            NOARCF( 1, NBARCF ) = NAETOI(NS2,NA1)
C           L'ARETE SUIVANTE
            NOARCF( 2, NBARCF ) = NBARCF + 1
C           LE NUMERO DU TRIANGLE DE L'AUTRE COTE
            NOARCF( 3, NBARCF ) = NAETOI(1,NA1)
C
C           TRACE DE L'ARETE
            CALL DVTRAR( PXYD, NAETOI(2,NA1), NAETOI(3,NA1),
     %                   NCVERT, NCBLAN )
C
C           SUPPRESSION DE L'ARETE
            IF( N1AEOC .EQ. NA1 ) THEN
                N1AEOC = NAETOI(4,NA1)
            ELSE
                NAETOI(4,NA0) = NAETOI(4,NA1)
            ENDIF
            GOTO 120
         ENDIF
      ENDIF
C
      IF( NS1 .NE. NS0 ) THEN
C        ARETE NON RETROUVEE : L'ETOILE NE SE REFERME PAS
CCC         WRITE(IMPRIM,*) 'DVFOAR: ANOMALIE DVFOAR A CORRIGER'
CCC         WRITE(IMPRIM,*) 'DVFOAR: CF NON REFERME'
         WRITE(IMPRIM,*) 'DVFOAR 125: REVOYEZ VOS DONNEES'
         WRITE(IMPRIM,*) 'LES LIGNES FERMEES DOIVENT ETRE DISJOINTES'
         WRITE(IMPRIM,*) 'VERIFIEZ SI C''EST LE CAS'
         IERR = 14
         RETURN
      ENDIF
C
C     LE SOMMET SUIVANT DU DERNIER SOMMET EST LE PREMIER
C     CHAINAGE CIRCULAIRE
      NOARCF( 2, NBARCF ) = 1
C
C     CHAINAGE DES SOMMETS DE ARCF VIDES ( 2 SONT NECESSAIRES ENSUITE )
      N1ARCF(0) = NBARCF+3
      MMARCF = MIN(8*NBARCF,MXARCF)
      DO 125 I=NBARCF+3,MMARCF
         NOARCF(2,I) = I+1
 125  CONTINUE
      NOARCF(2,MMARCF) = 0
CCCC
CCCC     EXISTE-T-IL UN SOMMET TROP PROCHE DE L'ARETE ?
CCCC     SI OUI: AJOUT AU CF DE TOUS LES TRIANGLES LE CONTENANT
CCCC             ET SUPPRESSION DU POINT
CCCC     ------------------------------------------------------
CCC      NA0 = NBTRCF
CCCC     LONGUEUR DE L'ARETE
CCC      DD = SQRT( (PXYD(1,NSA2)-PXYD(1,NSA1))**2 +
CCC     %                 (PXYD(2,NSA2)-PXYD(2,NSA1))**2 ) * 0.1
CCC      DO 150 I=1,NBARCF
CCC         NS2 = NOARCF(1,I)
CCC         IF( NLSOFR(NS2) .EQ. 0 ) THEN
CCCC           POINT NON SUR UNE LIGNE FRONTALIERE NI IMPOSE
CCC            IF( NS2 .NE. NSA1 .AND. NS2 .NE. NSA2 ) THEN
CCCC
CCCC              DISTANCE DU POINT I DU CF A L'ARETE NSA1 NSA2
CCC               D = DIPTDR( PXYD(1,NS2),PXYD(1,NSA1),PXYD(1,NSA2) )
CCC               IF( D .LT. DD ) THEN
CCCC
CCCC                 POINT TROP PROCHE DE L'ARETE
CCC                  WRITE(IMPRIM,*)'DVFOAR: SOMMET ',NS2,' SUPPRIME'
CCCC                 RECHERCHE DES TRIANGLES DE CE SOMMET
CCC                  NT00 = NOTRSO( NS2 )
CCC                  NT0  = NT00
CCCC
CCCC                 NT0 EST IL DEJA DANS NOTRCF ?
CCC 128                    DO 130 J=1,NBTRCF
CCC                     IF( NT0 .EQ. NOTRCF(J) ) GOTO 135
CCC 130                    CONTINUE
CCCC                 AJOUT DE NT0
CCC                  NBTRCF = NBTRCF + 1
CCC                  NOTRCF( NBTRCF ) = NT0
CCCCCC                  CALL DVTRTR( PXYD, NOTRIA, NT0, NCROUG, NCJAUN )
CCCC
CCCC                 RECHERCHE DU TRIANGLE SUIVANT
CCC 135                    DO 140 J=1,3
CCC                     IF( NS2 .EQ. NOTRIA(J,NT0) ) THEN
CCCC                       NS2 EST LE J-EME SOMMET
CCCC                       LE TRIANGLE SUIVANT
CCC                        NT1 = NOTRIA( J+3, NT0 )
CCC                        IF( NT1 .NE. NT00 ) THEN
CCC                           NT0 = NT1
CCC                           GOTO 128
CCC                        ENDIF
CCC                     ENDIF
CCC 140                    CONTINUE
CCC               ENDIF
CCC            ENDIF
CCC         ENDIF
CCC 150        CONTINUE
CCC      IF( NBTRCF .GT. NA0 ) GOTO 80
C
C     DESTRUCTION DES TRIANGLES DU CF
C     -------------------------------
      DO 166 I=1,NBTRCF
         NT0 = NOTRCF( I )
C        MISE DE NT0 DANS LES TRIANGLES VIDES
         NOTRIA( 1, NT0 ) = 0
         NOTRIA( 4, NT0 ) = N1TRVI
         N1TRVI = NT0
 166  CONTINUE
C
C     REPERAGE DES SOMMETS NS1 NS2 DANS LE CF  ( CF EQUIVALENCE )
C     ---------------------------------------
      NS1 = NS1ARE
      NS2 = NS2ARE
      NS(1) = NS1
      NS(2) = NS2
      DO 170 I=1,2
         NA0 = N1ARCF(1)
  167    IF( NOARCF(1,NA0) .EQ. NS(I) ) GOTO 168
         NA0 = NOARCF( 2, NA0 )
         GOTO 167
 168     NSS(I) = NA0
 170  CONTINUE
C
C     FORMATION DES 2 CF A PARTIR DE L'ARETE NS1-NS2
C     ----------------------------------------------
C     SAUVEGARDE DU SUIVANT DE NS1
      NA0 = NOARCF( 2, NSS1 )
      NT1 = NOARCF( 3, NSS1 )
C
C     LE PREMIER CF
      N1ARCF( 1 ) = NSS1
      NOARCF( 2, NSS1 ) = NSS2
      NOARCF( 3, NSS1 ) = -1
C
C     LE SECOND CF
C     L'ARETE DOUBLEE
      N1 = NBARCF + 1
      N2 = NBARCF + 2
      NOARCF( 1, N1 ) = NS2
      NOARCF( 2, N1 ) = N2
      NOARCF( 3, N1 ) = -1
      NOARCF( 1, N2 ) = NS1
      NOARCF( 2, N2 ) = NA0
      NOARCF( 3, N2 ) = NT1
      N1ARCF( 2 ) = N1
C
C     RECHERCHE DU PRECEDENT DE NSS2
 175  NA1 = NOARCF( 2, NA0 )
      IF( NA1 .NE. NSS2 ) THEN
C        PASSAGE A L'ARETE SUIVANTE
         NA0 = NA1
         GOTO 175
      ENDIF
C     NA0 PRECEDE NSS2 => IL PRECEDE N1
      NOARCF( 2, NA0 ) = N1
C
C     DEPART AVEC 2 CF
      NBET   = 2
      NBTRC1 = 0
C
C     TANT QUE LE NOMBRE DE CF CREES EST NON NUL FAIRE
C     TRIANGULATION FRONTALE DU CF
C     ================================================
 180  IF( NBET .GT. 0 ) THEN
C
C        L'ETOILE EN HAUT DE PILE A POUR PREMIERE ARETE
         NA01 = N1ARCF( NBET )
         NA1  = NOARCF( 2, NA01 )
C
C        LE CF A T IL TOUTES CES ARETES SUR LA FRONTIERE ?
 190     IF( NOARCF(3,NA1) .GT. 0 ) THEN
C           ARETE NON FRONTALIERE : RETOUR A L'ARETE INITIALE
            NA1  = NOARCF( 2, NA01 )
            GOTO 200
         ENDIF
         IF( NA1 .NE. NA01 ) THEN
C           PASSAGE A l'ARETE SUIVANTE
            NA1 = NOARCF(2,NA1)
            GOTO 190
         ENDIF
C        TOUTES LES ARETES ONT ETE VUES SUR LA FRONTIERE
         NBET = NBET - 1
         GOTO 180
C
C        CHOIX DU SOMMET DE L'ETOILE A RELIER A L'ARETE NA1
 200     CALL TRCHTD( PXYD, NA01, NA1, NOARCF,
     %                NA03, NA3, NARMIN )
         IF( NA3 .EQ. 0 ) THEN
            IERR = 2
            RETURN
         ENDIF
C
C        L'ARETE SUIVANTE DE NA1
         NA02 = NA1
         NA2  = NOARCF( 2, NA1 )
C
C        FORMATION DU TRIANGLE ARETE NA1 - SOMMET NOARCF(1,NA3)
         CALL TRTRCF( NBET,   NA01, NA1, NA02, NA2, NA03, NA3,
     %                N1TRVI, NOTRIA, NOTRSO,
     %                MMARCF, N1ARCF, NOARCF, NF )
         IF( NF .LE. 0 ) THEN
C           SATURATION DU TABLEAU NOTRIA OU NOARCF OU N1ARCF
            IERR = 1
            RETURN
         ENDIF
C
C        AJOUT DU TRIANGLE CREE A SA LISTE
         IF( NBTRC1 .GE. MXARCF ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) ='SATURATION DES TRIANGLES DU CF (NOTRCF)'
            CALL LERESU
            IERR = 1
            RETURN
         ENDIF
         NBTRC1 = NBTRC1 + 1
         NOTRCF( NBTRC1 ) = NF
         GOTO 180
      ENDIF
C
C     MISE A JOUR DU CHAINAGE DES FACES ADJACENTES PAR LES ARETES
C     -----------------------------------------------------------
      DO 700 NFP0 = 1, NBTRC1
C        LE NUMERO DE LA FACE DANS NOTRIA
         NF0 = NOTRCF( NFP0 )
C        BOUCLE SUR SES 3 ARETES
         DO 600 I=1,3
C           SEULE UNE ARETE SANS TRIANGLE OPPOSE EST TRAITEE
            IF( NOTRIA( 3+I, NF0 ) .GT. 0 ) GOTO 600
C           LES 2 SOMMETS
            NS1 = NOTRIA(I,NF0)
            IF( I .EQ. 3 ) THEN
               NS2 = 1
            ELSE
               NS2 = I + 1
            ENDIF
            NS2 = NOTRIA(NS2,NF0)
C
C           RECHERCHE DE L'ARETE NS1-NS2 DANS LES FACES NOTRIA
            DO 530 NFP1=NFP0+1,NBTRC1
               NF = NOTRCF( NFP1 )
               DO 510 J=1,3
C                 LES SOMMETS DE L'ARETE
                  NS3 = NOTRIA(J,NF)
                  IF( J .EQ. 3 ) THEN
                     NS4 = 1
                  ELSE
                     NS4 = J + 1
                  ENDIF
                  NS4 = NOTRIA(NS4,NF)
                  IF( NS3 .EQ. NS2 ) THEN
                     IF( NS4 .EQ. NS1 ) GOTO 550
                  ENDIF
C                 ARETE NON RETROUVEE DANS LE TRIANGLE NF
 510           CONTINUE
C              PASSAGE A LA FACE SUIVANTE DU CF
 530        CONTINUE
C           ARETE NON RETROUVEE => ELLE EST SUPPOSEE FRONTALIERE
            NOTRIA( 3+I, NF0 ) = 0
            GOTO 600
C
C           ARETE RETROUVEE ( NF0 , I ) <-> ( NF , J )
 550        NOTRIA(3+I,NF0) = NF
            NOTRIA(3+J,NF ) = NF0
 600     CONTINUE
 700  CONTINUE
C
      RETURN
      END
