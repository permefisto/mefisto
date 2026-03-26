      SUBROUTINE LIERCF( COSMIN, COSMAX, COANPL,
     %                   MXFAPE, NBFAPE, NOFAPE,
     %                   MXTETR, NOTETR, N1TETS,
     %                   MXFACO, NBFACO, LEFACO, N1FASC,
     %                   MXSOMM, NBSOMM, PTXYZD, NOSTIS, 
     %                   MXTRCF, NOTRCF, NOSTCF, MXARCF, N1ARCF, NOARCF,
     %                   MXETOI, NAETOI, MXTEFA, NOTEFA,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MODIFICATION DE LA TRIANGULATION LEFACO DES SURFACES CONTOUR
C -----    DE L'OBJET POUR INTEGRER LE MAXIMUM D'ARETES DES FACES DE LA
C          TETRAEDRISATION ACTUELLE
C
C ENTREES:
C --------
C COSMIN : COSINUS MINIMAL DES ANGLES DES TRIANGLES DE LA TRIANGULATION
C COSMAX : COSINUS MAXIMAL DES ANGLES DES TRIANGLES DE LA TRIANGULATION
C COANPL : SEUIL DU COSINUS DE L'ANGLE FORME PAR LES NORMALES AUX
C          2 FACES ET AU DESSUS DUQUEL LES FACES
C          SONT CONSIDEREES COPLANAIRES
C MXFAPE : NOMBRE MAXIMAL DE FACES PERDUES DE LEFACO DANS NOFAPE
C MXFACO : NOMBRE MAXIMAL DE TRIANGLES PERMIS POUR LA TRIANGULATION
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DE LA TETRAEDRISATION
C NBSOMM : NOMBRE ACTUEL  DE SOMMETS DE LA TETRAEDRISATION
C MXTRCF : NOMBRE MAXIMAL DECLARABLE D'ARETES OU TRIANGLES
C          OU SOMMETS DANS L'ETOILE OU NOMBRE DE CF
C MXARMI : NOMBRE MAXIMAL DE SOMMETS A MEME DISTANCE D'UNE ARETE
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C
C ENTREES ET SORTIES:
C -------------------
C NBFAPE : NOMBRE DE  FACES PERDUES DE   LEFACO
C NOFAPE : NUMERO DES FACES PERDUES DANS LEFACO
C NBFACO ; NOMBRE DE FACES ACTIVES DE LEFACO
C LEFACO : LES 3 SOMMETS, 2 MATERIAUX, 3 FACES VOISINES ET CHAINAGE
C          DES FACES TRIANGULAIRES DU CONTOUR ET INTERFACES
C N1FASC : N1FASC(I) NUMERO D'UN TRIANGLE DE LEFACO DE SOMMET I
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS ET DISTANCE SOUHAITEE
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C
C AUXILIAIRES :
C -------------
C NAETOI : LISTE DES ARETES DE L'ETOILE DU POINT COURANT
C          TABLEAU ENTIER(4,MXTRCF)
C N1ARCF : TABLEAU (0:MXTRCF) AUXILIAIRE
C NOARCF : NUMERO DES ARETES DE LA LIGNE DU CONTOUR FERME SELON UN SENS
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DE L'ETOILE
C
C SORTIES:
C --------
C NOSTIS : NUMERO DES SOMMETS ISOLES N'APPARTENANT PAS AU CONTOUR
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DU CF
C NOSTCF : TABLEAU DU NUMERO DANS PTXYZD DES SOMMETS DU CF
C IERR   : =0 SI PAS D'ERREUR
C          >0 SATURATION DU TABLEAU LEFACO OU NOARCF OU N1ARCF
C             OU IMPOSSIBILITE DE CONTINUER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1993
C2345X7..............................................................012
      DOUBLE PRECISION   QTRIMI
      PARAMETER        ( QTRIMI=0.25D0 )
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  PTXYZD(1:4,1:*),PT(3)
      DOUBLE PRECISION  SURTRD, D, DMIN, DINFO, Q, QQ, DMIMA
      DOUBLE PRECISION  Q1, Q2, Q3, Q4, S0, S1
      INTRINSIC         SQRT
C
      INTEGER           NOFAPE(MXFAPE),
     %                  N1TETS(MXSOMM),
     %                  NOTETR(8,MXTETR),
     %                  LEFACO(11,0:MXFACO),
     %                  N1FASC(NBSOMM)
C
      INTEGER           N1ARCF(0:MXARCF),
     %                  NOARCF(1:3,1:MXARCF),
     %                  NAETOI(4,MXETOI),
     %                  NOTRCF(MXTRCF),
     %                  NOSTCF(MXTRCF),
     %                  NOSTIS(MXSOMM),
     %                  NOTEFA(MXTEFA),
     %                  NTAUX(1)
C
C     DETERMINATION DU SEUIL DES ANGLES POUR L'ADMISSIBILITE DES ARETES
      EPSANG = REAL( ACOS( MAX( ABS(COSMIN), ABS(COSMAX) ) ) * 0.25D0 )
      EPSANG = MIN( EPSANG, 0.035 )
C              PRECISION = 2 DEGRES
C
C     REGROUPEMENT DES FACES PERDUES ADJACENTES COPLANAIRES
C     =====================================================
      NFP  = 0
 1    IERR = 0
      NFP = NFP + 1
      IF( NFP .LE. NBFAPE ) THEN
C
C        LE NUMERO DE LA FACE PERDUE DANS LEFACO
         NF0 = NOFAPE( NFP )
         IF( NF0 .LE. 0 ) GOTO 1
C
C        LES 2 VOLUMES DE LA FACE
         NV1 = LEFACO(4,NF0)
         NV2 = LEFACO(5,NF0)
C
C        FORMATION DU PREMIER CONTOUR FERME FORME DES
C        TRIANGLES PERDUS COPLANAIRES ADJACENTS A PARTIR DE NF0
         CALL VDR1CF( COANPL, NF0,    NBSOMM, PTXYZD, N1TETS, NOTETR,
     %                NBFAPE, NOFAPE, MXFACO, LEFACO, N1FASC, NTAUX,
     %                MXTRCF, NAETOI,
     %                MXTRCF, NBTRCF, NOTRCF,
     %                MXTRCF, NBCF,   N1ARCF, NOARCF,
     %                MXTRCF, NBSTCF, NOSTCF,
     %                MXSOMM, NBSTIS, NOSTIS,
     %                IERR )
         IF( IERR .NE. 0 ) GOTO 1
C
C        CES PLUS DE 2 FACES SONT RETIREES DES FACES PERDUES POUR NE
C        PAS ETRE RETRAITEES
         DO 4 I=1,NBTRCF
            NF = NOTRCF(I)
            DO 2 J=1,NBFAPE
               IF( NF .EQ. NOFAPE(J) ) THEN
C                 FACE PERDUE A RETIRER DE LEFACO
                  NOFAPE( J ) = 0
                  GOTO 4
               ENDIF
 2          CONTINUE
 4       CONTINUE
C
C        COMPRESSION DES FACES PERDUES
         CALL VDCFAP( NBFAPE, NOFAPE )
C
         IF( NBTRCF .LE. 1 .OR. NBTRCF .GT. 16 ) GOTO 1
C
C        SI TOUTES LES ARETES DES TRIANGLES PERDUS SONT DANS LA TETRAEDRISATION
C        NE RIEN FAIRE
         DO 7 I=1,NBTRCF
            NF = NOTRCF(I)
            DO 6 J=1,3
               IF( J .EQ. 3 ) THEN
                  J1 = 1
               ELSE
                  J1 = J + 1
               ENDIF
               CALL VDARTE( LEFACO(J,NF), LEFACO(J1,NF), N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NT )
               IF( NT .LE. 0 ) GOTO 8
 6          CONTINUE
 7       CONTINUE
C        TOUTES LES ARETES SONT DANS LA TETRAEDRISATION
         GOTO 1
C
C        DESTRUCTION DANS LEFACO DES FACES ADJACENTES COPLANAIRES PERDUES
 8       DO 9 I=1,NBTRCF
            CALL VDDSFA( NOTRCF(I), MXFACO, LEFACO, NBFACO, N1FASC )
 9       CONTINUE
C
C        LE NOMBRE DE TRIANGLES AJOUTES
         NBTRC1 = 0
         NBCF   = 1
         IF( NBTRCF .EQ. 2 ) THEN
C           CF QUADRANGULAIRE
            NBAR = 4
            GOTO 420
         ENDIF
         IF( NBSTIS .LE. 0 ) GOTO 400
C
C        ======================================================================
C        RECHERCHE DE 2 POINTS ISOLES RELIES A 2 SOMMETS DU CF POUR FORMER 2 CF
C        ======================================================================
         DO 45 I=1,NBSTIS-1
            IF( NOSTIS( I ) .LE. 0 ) GOTO 45
            DO 10 J=I+1,NBSTIS
               IF( NOSTIS( J ) .LE. 0 ) GOTO 10
C              EXISTE-T-IL DANS LA TETRAEDRISATION UNE ARETE RELIANT
C              NOSTIS(I) ET NOSTIS(J)?
               CALL VDARTE( NOSTIS(I), NOSTIS(J), N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NT )
               IF( NT .LE. 0 ) GOTO 10
C
C              BOUCLE SUR LES CF ACTUELS
               NBC = 0
 11            NBC = NBC + 1
               IF( NBC .GT. NBCF ) GOTO 10
C
C              L'ARETE NOSTIS(I)-NOSTIS(J) EST ELLE ADMISSIBLE POUR LE CF?
               CALL ARTADM( EPSANG, NOSTIS(I), NOSTIS(J),
     %                      N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 11
C
C              ICI L'ARETE NOSTIS(I)-NOSTIS(J) EST ADMISSIBLE POUR LE CF
C              EST ELLE ADMISSIBLE POUR LES POINTS ISOLES ?
               CALL SMTADM( NOSTIS(I), NOSTIS(J),
     %                      NBSTIS, NOSTIS, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 11
C
C              ICI L'ARETE NOSTIS(I)-NOSTIS(J) EST ADMISSIBLE POUR LE CF
C              ET POUR LES POINTS ISOLES
C              EXISTE T IL 2 SOMMETS DES CF RELIES A CES POINTS ET ADMISSIBLES?
C
               NA01 = N1ARCF( NBC )
               NA1  = NOARCF( 2, NA01 )
               NA00 = NA1
C
C              EXISTE-T-IL DANS LA TETRAEDRISATION UNE ARETE RELIANT
C              NOSTIS(I) ET NOARCF(1,NA1) ?
 14            CALL VDARTE( NOSTIS(I), NOARCF(1,NA1), N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NT )
               IF( NT .GT. 0 ) GOTO 20
C
C              PASSAGE A L'ARETE SUIVANTE DU CF
 15            NA01 = NA1
               NA1  = NOARCF( 2, NA01 )
               IF( NA1 .NE. NA00 ) GOTO 14
C              PASSAGE AU CF SUIVANT
               GOTO 11
C
C              L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ELLE ADMISSIBLE ?
 20            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA1),
     %                      N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 15
C
C              ICI L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ADMISSIBLE POUR LE CF
C              EST ELLE ADMISSIBLE POUR LES POINTS ISOLES ?
               CALL SMTADM( NOSTIS(I), NOARCF(1,NA1),
     %                      NBSTIS, NOSTIS, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 15
C
C              EXISTE-T-IL DANS LA TETRAEDRISATION UNE ARETE RELIANT
C              NOSTIS(J) ET NOARCF(1,NA2) ?
               NA02 = NA1
               NA2  = NOARCF( 2, NA02 )
               NAQ  = NA2
C
 24            CALL VDARTE( NOSTIS(J), NOARCF(1,NA2), N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NT )
               IF( NT .GT. 0 ) GOTO 30
C
C              PASSAGE A L'ARETE SUIVANTE DU CF
 25            NA02 = NA2
               NA2  = NOARCF( 2, NA02 )
               IF( NA2 .NE. NAQ ) GOTO 24
C              PASSAGE AU POINT ISOLE SUIVANT
               GOTO 45
C
C              L'ARETE NOSTIS(J)-NOARCF(1,NA2) EST ELLE ADMISSIBLE ?
 30            CALL ARTADM( EPSANG, NOSTIS(J), NOARCF(1,NA2),
     %                      N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 25
C
C              ICI L'ARETE NOSTIS(J)-NOARCF(1,NA2) EST ADMISSIBLE POUR LE CF
C              EST ELLE ADMISSIBLE POUR LES POINTS ISOLES ?
               CALL SMTADM( NOSTIS(J), NOARCF(1,NA2),
     %                      NBSTIS, NOSTIS, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 25
C
C              LES ARETES NOARCF(1,NA1)-NOSTIS(I)-NOSTIS(J)-NOARCF(1,NA2)
C              S'INTERSECTENT ELLES ?
C              SI OUI : CELA SIGNIFIE QU'IL EXISTE 2 ARETES DANS LA TETRAEDRISAT
C              QUI S'INTERSECTENT !   => PAS DE VERIFICATION EFFECTUEE
CCC               CALL TRARTR( NCMAGE, NOARCF(1,NA1), NOSTIS(I), PTXYZD )
CCC               CALL TRARTR( NCMAGE, NOSTIS(I),     NOSTIS(J), PTXYZD )
CCC               CALL TRARTR( NCMAGE, NOSTIS(J), NOARCF(1,NA2), PTXYZD )
CCCC
C              LE PREMIER CF
C              RECHERCHE DE 3 ARETES VIDES DE CF
               NAV1 = N1ARCF(0)
               IF( NAV1 .LE. 0 ) GOTO 997
C              LA 1-ERE ARETE VIDE EST MISE A JOUR
               N1ARCF(0) = NOARCF( 2, NAV1 )
               NAV2 = N1ARCF(0)
               IF( NAV2 .LE. 0 ) GOTO 997
               N1ARCF(0) = NOARCF( 2, NAV2 )
               NAV3 = N1ARCF(0)
               IF( NAV3 .LE. 0 ) GOTO 997
               N1ARCF(0) = NOARCF( 2, NAV3 )
C
C              AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NOSTIS(I)
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV1 ) = NAV2
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV1 ) = -1
C
C              AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NOSTIS(J)
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV2 ) = NOSTIS(I)
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV2 ) = NAV3
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV2 ) = -1
C
C              AJOUT DE L'ARETE NAV3 PARTANT DE NOSTIS(J) VERS NA2
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV3 ) = NOSTIS(J)
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV3 ) = NA2
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV3 ) = -1
C
C              SI NA1=NA2 UN CF EST INTERNE A L'AUTRE
               IF( NA1 .EQ. NA2 ) THEN
C                 LE PREMIER CF FAIT UN CYCLE INTERNE AU SECOND
                  NOARCF( 2, NAV3 ) = NAV1
               ELSE
C                 L'ARETE PRECEDENTE NA01 POINTE SUR LA NOUVELLE NAV1
                  NOARCF( 2, NA01 ) = NAV1
               ENDIF
C
C              L'ORIGINE DU CF
               N1ARCF( NBC ) = NAV2
C
C              LE SECOND CF
C              RECHERCHE DE 3 ARETES VIDES DE CF
               NAV1 = N1ARCF(0)
               IF( NAV1 .LE. 0 ) GOTO 997
C              LA 1-ERE ARETE VIDE EST MISE A JOUR
               N1ARCF(0) = NOARCF( 2, NAV1 )
               NAV2 = N1ARCF(0)
               IF( NAV2 .LE. 0 ) GOTO 997
               N1ARCF(0) = NOARCF( 2, NAV2 )
               NAV3 = N1ARCF(0)
               IF( NAV3 .LE. 0 ) GOTO 997
               N1ARCF(0) = NOARCF( 2, NAV3 )
C
C              AJOUT DE L'ARETE NAV1 PARTANT DE NA2 VERS NOSTIS(J)
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV1 ) = NOARCF( 1, NA2 )
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV1 ) = NAV2
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV1 ) = -1
C
C              AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(J) VERS NOSTIS(I)
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV2 ) = NOSTIS(J)
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV2 ) = NAV3
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV2 ) = -1
C
C              AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV3 ) = NOSTIS(I)
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV3 ) = NA1
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV3 ) = -1
C
C              L'ARETE PRECEDENTE NA02 POINTE SUR LA NOUVELLE NAV1
               NOARCF( 2, NA02 ) = NAV1
C
C              L'ORIGINE DU CF
               IF( NBCF .GE. MXTRCF ) GOTO 995
               NBCF = NBCF + 1
               N1ARCF( NBCF ) = NAV2
C
C              LES POINTS NE SONT PLUS ISOLES
               NOSTIS( I ) = 0
               NOSTIS( J ) = 0
               GOTO 45
 10         CONTINUE
 45      CONTINUE
         CALL VDCFAP( NBSTIS, NOSTIS )
C
C        ===================================================================
C        RECHERCHE D'UN POINT ISOLE RELIE A 2 SOMMETS DU CF POUR FORMER 2 CF
C        ===================================================================
C        2 POINTS ISOLES NE PEUVENT ETRE RELIES ENTRE EUX DANS LA TETRAEDRISATIO
         DO 90 I=1,NBSTIS
            IF( NOSTIS( I ) .LE. 0 ) GOTO 90
C
C           BOUCLE SUR LES CF ACTUELS SUSCEPTIBLES DE CONTENIR NOSTIS(I)
            NBC = 1
C
 55            NA01 = N1ARCF( NBC )
               NA1  = NOARCF( 2, NA01 )
               NA00 = NA1
C
C              EXISTE-T-IL DANS LA TETRAEDRISATION UNE ARETE RELIANT
C              NOSTIS(I) ET NOARCF(1,NA1) ?
 60            CALL VDARTE( NOSTIS(I), NOARCF(1,NA1), N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NT )
               IF( NT .GT. 0 ) GOTO 70
C
C              PASSAGE A L'ARETE SUIVANTE DU CF
 65            NA01 = NA1
               NA1  = NOARCF( 2, NA01 )
               IF( NA1 .NE. NA00 ) GOTO 60
C              PASSAGE AU CF SUIVANT
               IF( NBC .LT. NBCF ) THEN
                  NBC = NBC + 1
                  GOTO 55
               ENDIF
C              PASSAGE AU POINT ISOLE SUIVANT
               GOTO 90
C
C              L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ELLE ADMISSIBLE ?
 70            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA1),
     %                      N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 65
C
C              ICI L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ADMISSIBLE POUR LE CF
C              EST ELLE ADMISSIBLE POUR LES POINTS ISOLES ?
               CALL SMTADM( NOSTIS(I), NOARCF(1,NA1),
     %                      NBSTIS, NOSTIS, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 65
C
C              RECHERCHE DE LA 2-EME ARETE NOSTIS(I)-SUIVANT DU SOMMET NOARCF(1,
C              NA2 NE PEUT ETRE L'ARETE SUIVANTE DE NA1 AFIN D'ESSAYER D'EVITER
C              D'OBTENIR UN TRIANGLE AVEC UN POINT ISOLE INTERNE
               NA02 = NA1
               NA2  = NOARCF( 2, NA02 )
               IF( NA2 .EQ. NA1 ) GOTO 65
C
 80            NA02 = NA2
               NA2  = NOARCF( 2, NA02 )
               IF( NA2 .EQ. NA1 ) GOTO 65
C              L'ARETE NOSTIS(I)-NOARCF(1,NA2) EST ELLE DANS LA TETRAEDRISATION?
               CALL VDARTE( NOSTIS(I), NOARCF(1,NA2), N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NT )
               IF( NT .LE. 0 ) GOTO 80
C
C              L'ARETE EST DANS LA TETRAEDRISATION. EST ELLE ADMISSIBLE?
               CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA2),
     %                      N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 80
C
C              L'ARETE EST ADMISSIBLE POUR LE CF.
C              EST ELLE ADMISSIBLE POUR LES POINTS ISOLES ?
               CALL SMTADM( NOSTIS(I), NOARCF(1,NA2),
     %                      NBSTIS, NOSTIS, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 80
CCCC
CCCC              TRACE DES ARETES UTILISEES
CCC               CALL TRARTR( NCMAGE, NOARCF(1,NA1), NOSTIS(I), PTXYZD )
CCC               CALL TRARTR( NCMAGE, NOSTIS(I), NOARCF(1,NA2), PTXYZD )
C
C              LE PREMIER CF
C              -------------
C              RECHERCHE DE 2 ARETES VIDES DE CF
               NAV1 = N1ARCF(0)
               IF( NAV1 .LE. 0 ) GOTO 997
C              LA 1-ERE ARETE VIDE EST MISE A JOUR
               N1ARCF(0) = NOARCF( 2, NAV1 )
               NAV2 = N1ARCF(0)
               IF( NAV2 .LE. 0 ) GOTO 997
               N1ARCF(0) = NOARCF( 2, NAV2 )
C
C              AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NOSTIS(I)
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV1 ) = NAV2
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV1 ) = -1
C
C              AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA2
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV2 ) = NOSTIS(I)
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV2 ) = NA2
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV2 ) = -1
C
C              L'ARETE PRECEDENTE NA01 POINTE SUR LA NOUVELLE NAV1
               NOARCF( 2, NA01 ) = NAV1
C
C              L'ORIGINE DU CF
               N1ARCF( NBC ) = NAV2
C
C              LE SECOND CF
C              ------------
C              RECHERCHE DE 2 ARETES VIDES DE CF
               NAV1 = N1ARCF(0)
               IF( NAV1 .LE. 0 ) GOTO 997
C              LA 1-ERE ARETE VIDE EST MISE A JOUR
               N1ARCF(0) = NOARCF( 2, NAV1 )
               NAV2 = N1ARCF(0)
               IF( NAV2 .LE. 0 ) GOTO 997
               N1ARCF(0) = NOARCF( 2, NAV2 )
C
C              AJOUT DE L'ARETE NAV1 PARTANT DE NA2 VERS NOSTIS(I)
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV1 ) = NOARCF( 1, NA2 )
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV1 ) = NAV2
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV1 ) = -1
C
C              AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C              LE NUMERO DU SOMMET
               NOARCF( 1, NAV2 ) = NOSTIS(I)
C              L'ARETE SUIVANTE
               NOARCF( 2, NAV2 ) = NA1
C              LE TRIANGLE LEFACO OPPOSE EST INCONNU
               NOARCF( 3, NAV2 ) = -1
C
C              L'ARETE PRECEDENTE NA02 POINTE SUR LA NOUVELLE NAV1
               NOARCF( 2, NA02 ) = NAV1
C
C              L'ORIGINE DU CF
               IF( NBCF .GE. MXTRCF ) GOTO 995
               NBCF = NBCF + 1
               N1ARCF( NBCF ) = NAV2
C
C              LE POINT N'EST PLUS ISOLE
               NOSTIS( I ) = 0
 90      CONTINUE
         CALL VDCFAP( NBSTIS, NOSTIS )
C
C        ========================================================
C        S'IL RESTE DES POINTS ISOLES, LES RELIER AU SOMMET DU CF
C        FORMANT UNE ARETE DE LA TETRAEDRISATION ET ADMISSIBLE
C        ========================================================
         DO 190 I=1,NBSTIS
C
C           LA BOUCLE SUR LES CF ACTUELS
            NBC = 1
C
 105        NA01 = N1ARCF( NBC )
            NA1  = NOARCF( 2, NA01 )
            NA00 = NA1
            QQ   = 0
            NAQ  = 0
C
C           EXISTE-T-IL DANS LA TETRAEDRISATION UNE ARETE RELIANT
C           NOSTIS(I) ET NOARCF(1,NA1) ?
 110        CALL VDARTE( NOSTIS(I), NOARCF(1,NA1), N1TETS, NOTETR,
     %                   MXTEFA, NOTEFA, NT )
            IF( NT .LE. 0 ) GOTO 115
C
C           L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ELLE ADMISSIBLE ?
            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA1),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .NE. 0 ) GOTO 150
C
C           PASSAGE A L'ARETE SUIVANTE DU CF
 115        NA01 = NA1
            NA1  = NOARCF( 2, NA1 )
            IF( NA1 .NE. NA00 ) GOTO 110
C           PASSAGE AU CF SUIVANT
            IF( NBC .LT. NBCF ) THEN
               NBC = NBC + 1
               GOTO 105
            ENDIF
C           PASSAGE AU POINT ISOLE SUIVANT
            GOTO 190
C
C           ICI L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ADMISSIBLE
C           RECHERCHE DE LA 2-EME ARETE NOSTIS(I)-PRECEDENT DU SOMMET NOARCF(1,N
C           L'ARETE EST ELLE ADMISSIBLE?
 150        CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA01),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .NE. 0 ) THEN
C              LES 3 ARETES SONT ADMISSIBLES . LA QUALITE DU TRIANGLE A FORMER
               CALL QUTRTE( PTXYZD(1,NOSTIS(I)),
     %                      PTXYZD(1,NOARCF(1,NA01)),
     %                      PTXYZD(1,NOARCF(1,NA1)), Q )
               IF( Q .GT. QQ ) THEN
C                 MEILLEURE QUALITE
                  NAQ = NA01
                  QQ  = Q
               ENDIF
            ENDIF
C
C           RECHERCHE DE LA 2-EME ARETE NOSTIS(I)-SUIVANT DU SOMMET NOARCF(1,NA1
            NA2  = NOARCF( 2, NA1 )
C           L'ARETE EST ELLE ADMISSIBLE?
            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA2),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .NE. 0 ) THEN
C              LES 3 ARETES SONT ADMISSIBLES . LA QUALITE DU TRIANGLE A FORMER
               CALL QUTRTE( PTXYZD(1,NOSTIS(I)),
     %                      PTXYZD(1,NOARCF(1,NA1)),
     %                      PTXYZD(1,NOARCF(1,NA2)), Q )
               IF( Q .GT. QQ ) THEN
C                 MEILLEURE QUALITE
                  NAQ = NA1
                  QQ  = Q
               ENDIF
            ENDIF
C
            IF( QQ .LE. 0 ) GOTO 115
C
C           LE TRIANGLE NOSTIS(I) NOARCF(1,NAQ) NOARCF(1,NOARCF(2,NAQ))
C           EST IL COMPATIBLE SANS POINT ISOLE INTERNE ?
            NAV1 = NOARCF( 2, NAQ )
            CALL TRIADM( NOSTIS(I), NOARCF(1,NAQ), NOARCF(1,NAV1),
     %                   NBSTIS, NOSTIS, PTXYZD,  NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 115
CCCC
CCCC           TRACE DES ARETES UTILISEES
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NAQ),  NOSTIS(I), PTXYZD )
CCC            CALL TRARTR( NCMAGE, NOSTIS(I), NOARCF(1,NAV1), PTXYZD )
C
C           3 ARETES SONT ADMISSIBLES TRIANGLE SANS POINT INTERNE
C           FORMATION D'UN CF INDEPENDANT
            NA01 = NAQ
            NA1  = NOARCF( 2, NA01 )
C
C           LE PREMIER CF RECTIFIE LE CF COURANT
C           RECHERCHE DE 2 ARETES VIDES DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
            N1ARCF(0) = NOARCF( 2, NAV2 )
C
C           RECHERCHE DE L'ARETE AVANT NA01
            NAQ = N1ARCF( NBC )
 160        IF( NOARCF(2,NAQ) .NE. NA01 ) THEN
               NAQ = NOARCF(2,NAQ)
               GOTO 160
            ENDIF
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA01 VERS NOSTIS(I)
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA01 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NAV2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOSTIS(I)
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ARETE PRECEDENTE NAQ POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NAQ ) = NAV1
C
C           L'ORIGINE DU CF
            N1ARCF( NBC ) = NAV2
C
C           LE SECOND CF EST EN FAIT UN TRIANGLE AVEC L'ARETE NA01
C           RECHERCHE DE 2 ARETES VIDES DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
            N1ARCF(0) = NOARCF( 2, NAV2 )
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NOSTIS(I)
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NAV2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA01
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOSTIS(I)
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA01
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ARETE NA01 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NAV2
C
C           LE POINT N'EST PLUS ISOLE
            NOSTIS( I ) = 0
 190     CONTINUE
         CALL VDCFAP( NBSTIS, NOSTIS )
C
C        ==============================================================
C        RECHERCHE DE L'ARETE DU CF QUI JOINT AU POINT ISOLE DONNE
C        LE TRIANGLE ADMISSIBLE DE MEILLEURE QUALITE
C        ==============================================================
         DO 290 I=1,NBSTIS
C
C           LA BOUCLE SUR LES CF ACTUELS
            NBC = 1
C
 205        NA01 = N1ARCF( NBC )
            NA1  = NOARCF( 2, NA01 )
            NA00 = NA1
            QQ   = 0
            NAQ  = 0
C
C           L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ELLE ADMISSIBLE ?
 210        CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA1),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .NE. 0 ) GOTO 250
C
C           PASSAGE A L'ARETE SUIVANTE DU CF
 215        NA01 = NA1
            NA1  = NOARCF( 2, NA1 )
            IF( NA1 .NE. NA00 ) GOTO 210
C
C           EXISTE T IL DES TRIANGLES ADMISSIBLES POUR CE CF?
            IF( QQ .GT. 0 ) THEN
C              OUI: CELUI DE MEILLEURE QUALITE EST RETENU
               GOTO 280
            ENDIF
C
C           PASSAGE AU CF SUIVANT
            IF( NBC .LT. NBCF ) THEN
               NBC = NBC + 1
               GOTO 205
            ENDIF
C
C           PASSAGE AU POINT ISOLE SUIVANT
            GOTO 290
C
C           ICI L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ADMISSIBLE
C           RECHERCHE DE LA 2-EME ARETE NOSTIS(I)-SUIVANT DU SOMMET NOARCF(1,NA1
 250        NA2  = NOARCF( 2, NA1 )
C           L'ARETE EST ELLE ADMISSIBLE?
            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA2),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 215
C
C           LE TRIANGLE NOSTIS(I) NOARCF(1,NA1) NOARCF(1,NA2)
C           EST IL COMPATIBLE SANS POINT ISOLE INTERNE ?
            CALL TRIADM( NOSTIS(I), NOARCF(1,NA1), NOARCF(1,NA2),
     %                   NBSTIS, NOSTIS, PTXYZD,  NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 215
C
C           LES 3 ARETES SONT ADMISSIBLES . LA QUALITE DU TRIANGLE A FORMER
            CALL QUTRTE( PTXYZD(1,NOSTIS(I)),
     %                   PTXYZD(1,NOARCF(1,NA1)),
     %                   PTXYZD(1,NOARCF(1,NA2)), Q )
            IF( Q .GT. QQ ) THEN
C              MEILLEURE QUALITE
               NAQ = NA01
               QQ  = Q
            ENDIF
            GOTO 215
C
C           3 ARETES SONT ADMISSIBLES TRIANGLE SANS POINT INTERNE
C           FORMATION D'UN CF INDEPENDANT
 280        NA01 = NAQ
            NA1  = NOARCF( 2, NA01 )
            NA2  = NOARCF( 2, NA1  )
C
CCCC           TRACE DES ARETES UTILISEES
CCC            CALL TRARTR( NCMAGE, NOSTIS(I), NOARCF(1,NA1), PTXYZD )
CCC            CALL TRARTR( NCMAGE, NOSTIS(I), NOARCF(1,NA2), PTXYZD )
C
C           LE PREMIER CF RECTIFIE LE CF COURANT
C           RECHERCHE DE 2 ARETES VIDES DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
            N1ARCF(0) = NOARCF( 2, NAV2 )
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NOSTIS(I)
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NAV2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOSTIS(I)
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ARETE PRECEDENTE NAQ POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           L'ORIGINE DU CF
            N1ARCF( NBC ) = NAV2
C
C           LE SECOND CF EST EN FAIT UN TRIANGLE AVEC L'ARETE NA01
C           RECHERCHE DE 2 ARETES VIDES DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
            N1ARCF(0) = NOARCF( 2, NAV2 )
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA2 VERS NOSTIS(I)
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA2 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NAV2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOSTIS(I)
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ARETE NA01 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA1 ) = NAV1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NAV2
C
C           LE POINT N'EST PLUS ISOLE
            NOSTIS( I ) = 0
 290     CONTINUE
         CALL VDCFAP( NBSTIS, NOSTIS )
C
C        ==============================================================
C        RECHERCHE DU PLUS PROCHE MILIEU D'ARETE DU CF D'UN POINT ISOLE
C        ==============================================================
         DO 390 I=1,NBSTIS
            DMIN = DINFO( 'GRAND' )
            NAV  = 0
C
C           LA BOUCLE SUR LES CF ACTUELS
            NBC  = 1
C
            NA00 = N1ARCF( NBC )
            NA1  = NA00
            NA2  = NOARCF(2,NA1)
C
C           RECHERCHE DU PLUS PROCHE MILIEU D'UNE ARETE DU CF
 310        PT(1)=(PTXYZD(1,NOARCF(1,NA1))+PTXYZD(1,NOARCF(1,NA2)))*0.5
            PT(2)=(PTXYZD(2,NOARCF(1,NA1))+PTXYZD(2,NOARCF(1,NA2)))*0.5
            PT(3)=(PTXYZD(3,NOARCF(1,NA1))+PTXYZD(3,NOARCF(1,NA2)))*0.5
            D = SQRT( ( PTXYZD(1,NOSTIS(I)) - PT(1) ) ** 2
     %              + ( PTXYZD(2,NOSTIS(I)) - PT(2) ) ** 2
     %              + ( PTXYZD(3,NOSTIS(I)) - PT(3) ) ** 2 )
            IF( D .LT. DMIN ) THEN
               DMIN = D
               NAV  = NA1
            ENDIF
            NA1 = NA2
            NA2 = NOARCF( 2, NA2 )
            IF( NA2 .NE. NA00 ) GOTO 310
C
C           NAV EST LA PLUS PROCHE ARETE
            NA1 = NAV
            NA2 = NOARCF( 2, NA1 )
C           L'ARETE NOSTIS(I)-NOARCF(1,NA1) EST ELLE ADMISSIBLE?
            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA1),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 390
C           OUI   L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
            CALL ARTADM( EPSANG, NOSTIS(I), NOARCF(1,NA2),
     %                   N1ARCF(NBC), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 390
C
C           LE TRIANGLE NOSTIS(I) NOARCF(1,NA1) NOARCF(1,NA2)
C           EST IL COMPATIBLE SANS POINT ISOLE INTERNE ?
            CALL TRIADM( NOSTIS(I), NOARCF(1,NA1), NOARCF(1,NA2),
     %                   NBSTIS, NOSTIS, PTXYZD,  NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 390
CCCC
CCCC           TRACE DES ARETES UTILISEES
CCC            CALL TRARTR( NCMAGE, NOSTIS(I), NOARCF(1,NA1), PTXYZD )
CCC            CALL TRARTR( NCMAGE, NOSTIS(I), NOARCF(1,NA2), PTXYZD )
C
C           LE PREMIER CF RECTIFIE LE CF COURANT
C           RECHERCHE DE L'ARETE NA01 AVANT NA1
            NA01 = N1ARCF( NBC )
 360        IF( NOARCF(2,NA01) .NE. NA1 ) THEN
               NA01 = NOARCF(2,NA01)
               GOTO 360
            ENDIF
C
C           RECHERCHE DE 2 ARETES VIDES DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
            N1ARCF(0) = NOARCF( 2, NAV2 )
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NOSTIS(I)
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NAV2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOSTIS(I)
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ARETE PRECEDENTE NAQ POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           L'ORIGINE DU CF
            N1ARCF( NBC ) = NAV2
C
C           LE SECOND CF EST EN FAIT UN TRIANGLE AVEC L'ARETE NA01
C           RECHERCHE DE 2 ARETES VIDES DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
            N1ARCF(0) = NOARCF( 2, NAV2 )
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA2 VERS NOSTIS(I)
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA2 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NAV2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NOSTIS(I) VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOSTIS(I)
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ARETE NA1 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA1 ) = NAV1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NAV2
C
C           LE POINT N'EST PLUS ISOLE
            NOSTIS( I ) = 0
 390     CONTINUE
         CALL VDCFAP( NBSTIS, NOSTIS )
C
         IF( NBSTIS .GT. 0 ) THEN
            WRITE(IMPRIM,*) 'LIERCF: IMPOSSIBLE DE JOINDRE ',NOSTIS(I),
     &         ' AUX SOMMETS DU CF. A AMELIORER'
            CALL CLICSO
            IERR = 1
            GOTO 999
         ENDIF
C
C        ===========================================================
C               ICI PLUS AUCUN POINT ISOLE DANS LES CF ACTUELS
C        ===========================================================
 400     IF( NBCF .LE. 0 ) GOTO 800
C
C        TRAITEMENT DU CF EN SOMMET DE PILE
C
C        RECHERCHE DU NOMBRE D'ARETES DU CF
         NBAR = 0
         NA00 = N1ARCF( NBCF )
         NA1  = NA00
C
 410     NBAR = NBAR + 1
         NA1  = NOARCF( 2, NA1 )
         IF( NA1 .NE. NA00 ) GOTO 410
C
         IF( NBAR .LE. 2 ) THEN
C           PROBLEME
            CALL CLICSO
         ENDIF
C
         IF( NBAR .EQ. 3 ) THEN
C
C           LE CF EST UN TRIANGLE . IL EST AJOUTE A LEFACO
C           ==============================================
            NA00 = N1ARCF( NBCF )
            NA1  = NOARCF( 2, NA00 )
            NA2  = NOARCF( 2, NA1  )
            CALL VDAJFA( NOARCF(1,NA00), NOARCF(1,NA1), NOARCF(1,NA2),
     %                   NV1, NV2,
     %                   NOARCF(3,NA00), NOARCF(3,NA1), NOARCF(3,NA2),
     %                   MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NF )
            IF( NF .LE. 0 ) GOTO 996
CCCC
CCCC           TRACE EVENTUEL DU TRIANGLE CREE
CCC            CALL TRFATR( NCMAGE, NCBLEU, LEFACO(1,NF), PTXYZD )
C
C           AJOUT DU TRIANGLE CREE DANS SA LISTE
            IF( NBTRC1 .GE. MXTRCF ) GOTO 998
            NBTRC1 = NBTRC1 + 1
            NOTRCF( NBTRC1 ) = NF
C
C           UN CF DE MOINS
            NBCF = NBCF - 1
            GOTO 400
         ENDIF
C
 420     IF( NBAR .EQ. 4 ) THEN
C
C           LE CF EST UN QUADRANGLE
C           =======================
            NA1 = N1ARCF( NBCF )
            NA2 = NOARCF( 2, NA1 )
            NA3 = NOARCF( 2, NA2 )
            NA4 = NOARCF( 2, NA3 )
C
C           RECHERCHE DU MEILLEUR CHOIX POSSIBLE PARMI LES 2 CAS
            NS1 = NOARCF( 1, NA1 )
            NS2 = NOARCF( 1, NA2 )
            NS3 = NOARCF( 1, NA3 )
            NS4 = NOARCF( 1, NA4 )
C
C           L'ARETE 1-3 EST ELLE ADMISSIBLE?
            CALL ARTADM( EPSANG, NS1, NS3,
     %                   N1ARCF(NBCF), NOARCF, PTXYZD, NONOU3 )
C           LA QUALITE DES 2 TRIANGLES A FORMER
            CALL QUTRTE(PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3), Q1)
            CALL QUTRTE(PTXYZD(1,NS1), PTXYZD(1,NS3), PTXYZD(1,NS4), Q2)
C
C           L'AUTRE DIAGONALE
C           L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
            CALL ARTADM( EPSANG, NS2, NS4,
     %                   N1ARCF(NBCF), NOARCF, PTXYZD, NONOU4 )
C           LA QUALITE DES 2 TRIANGLES A FORMER
            CALL QUTRTE(PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS4), Q3)
            CALL QUTRTE(PTXYZD(1,NS2), PTXYZD(1,NS3), PTXYZD(1,NS4), Q4)
C
            IF( NONOU3 .EQ. 0 .AND. NONOU4 .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'LIERCF: QUADRANGLE NON ADMISSIBLE'
               WRITE(IMPRIM,*) 'SOMMETS:',NS1,NS2,NS3,NS4
               WRITE(IMPRIM,*) 'CHOIX DES SURFACES MINIMALES'
CCC               CALL CLICSO
               Q1 = SURTRD( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3))
     %            + SURTRD( PTXYZD(1,NS1), PTXYZD(1,NS3), PTXYZD(1,NS4))
               Q2 = SURTRD( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS4))
     %            + SURTRD( PTXYZD(1,NS2), PTXYZD(1,NS3), PTXYZD(1,NS4))
               IF( Q1 .LT. Q2 ) THEN
                  GOTO 430
               ELSE
                  GOTO 440
               ENDIF
            ELSE IF( NONOU3 .NE. 0 .AND. NONOU4 .NE. 0 ) THEN
C              LES 2 DIAGONALES SONT ADMISSIBLES
C              CELLE FORMANT 2 TRIANGLES DE MEILLEURE QUALITE EST CHOISIE
               IF( MIN(Q1,Q2) .LT. MIN(Q3,Q4) ) GOTO 440
            ELSE IF( NONOU4 .NE. 0 ) THEN
               GOTO 440
            ENDIF
C
C           DIAGONALE 1-3   DECALAGE D'UNE ARETE
 430        NA2 = NOARCF( 2, NA2 )
            NA3 = NOARCF( 2, NA2 )
            NA4 = NOARCF( 2, NA3 )
            NA1 = NOARCF( 2, NA4 )
C
C           DIAGONALE 2-4   CONSTRUCTION DU TRIANGLE 124
 440        CALL VDAJFA( NOARCF(1,NA1), NOARCF(1,NA2), NOARCF(1,NA4),
     %                   NV1, NV2,
     %                   NOARCF(3,NA1), -1, NOARCF(3,NA4),
     %                   MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NF124 )
            IF( NF124 .LE. 0 ) GOTO 996
CCCC           TRACE EVENTUEL DU TRIANGLE CREE
CCC            CALL TRFATR( NCMAGE, NCBLEU, LEFACO(1,NF124), PTXYZD )
C           AJOUT DU TRIANGLE CREE DANS SA LISTE
            IF( NBTRC1 .GE. MXTRCF ) GOTO 998
            NBTRC1 = NBTRC1 + 1
            NOTRCF( NBTRC1 ) = NF124
C
C           CONSTRUCTION DU TRIANGLE 234
            CALL VDAJFA( NOARCF(1,NA2), NOARCF(1,NA3), NOARCF(1,NA4),
     %                   NV1, NV2,
     %                   NOARCF(3,NA2), NOARCF(3,NA3), NF124,
     %                   MXFACO, LEFACO, NBFACO,
     %                   N1FASC, NF234 )
ccc
ccc 17/5/2008
ccc ci dessous produit une erreur 2 triangles egaux opposes sur 2 aretes
ccc            IF( NF234 .LE. 0 ) GOTO 996
cccC           MISE A JOUR DU TRIANGLE VOISIN
ccc            DO 444 KK=6,8
ccc               IF( LEFACO(KK,NF124) .EQ. -1 ) THEN
ccc                  LEFACO(KK,NF124) = NF234
ccc                  GOTO 445
ccc               ENDIF
ccc 444        CONTINUE
ccc
CCCC           TRACE EVENTUEL DU TRIANGLE CREE
CCC            CALL TRFATR( NCMAGE, NCBLEU, LEFACO(1,NF234), PTXYZD )
C
C           AJOUT DU TRIANGLE CREE DANS SA LISTE
            IF( NBTRC1 .GE. MXTRCF ) GOTO 998
            NBTRC1 = NBTRC1 + 1
            NOTRCF( NBTRC1 ) = NF234
C
C           UN CF DE MOINS
            NBCF = NBCF - 1
            GOTO 400
C
         ELSE
C
C           LE CF A AU MOINS 5 ARETES
C           =========================
C
C           RECHERCHE D'UNE ARETE TRANSVERSE MAIS NE JOIGNANT PAS 2 ARETES
C           CONSECUTIVES DU CF
            IF( NBAR .EQ. 5 ) GOTO 600
            NA01 = N1ARCF( NBCF )
            NA00 = NOARCF( 2, NA01 )
            NA1  = NA00
C
 500        NA02 = NOARCF( 2, NA1  )
            NA02 = NOARCF( 2, NA02 )
            NA2  = NOARCF( 2, NA02 )
C
C           EXISTE-T-IL DANS LA TETRAEDRISATION UNE ARETE RELIANT
C           NOARCF(1,NA1) ET NOARCF(1,NA2) ?
 510        CALL VDARTE( NOARCF(1,NA1), NOARCF(1,NA2), N1TETS, NOTETR,
     %                   MXTEFA, NOTEFA, NT )
            IF( NT .GT. 0 ) GOTO 530
C
C           PASSAGE A L'ARETE SUIVANTE DU CF
 520        NA02 = NA2
            NA2  = NOARCF( 2, NA02 )
            IF( NOARCF( 2, NOARCF(2,NA2) ) .NE. NA1 ) GOTO 510
C
C           PASSAGE A L'ARETE SUIVANTE
            NA01 = NA1
            NA1  = NOARCF( 2, NA01 )
            IF( NOARCF( 2, NOARCF(2,NA1) ) .NE. NA00 ) GOTO 500
C
C           PAS D'ARETE DANS LA TETRAEDRISATION ET ADMISSIBLE
            GOTO 600
C
C           L'ARETE NOARCF(1,NA1)-NOARCF(1,NA2) EST ELLE ADMISSIBLE ?
 530        CALL ARTADM( EPSANG, NOARCF(1,NA1), NOARCF(1,NA2),
     %                   N1ARCF(NBCF), NOARCF, PTXYZD, NONOUI )
            IF( NONOUI .EQ. 0 ) GOTO 520
C
C           ICI ARETE NOARCF(1,NA1)-NOARCF(1,NA2) ADMISSIBLE ET
C           DANS LA TETRAEDRISATION
C
CCCC           TRACE DE L'ARETE UTILISEE
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NA1), NOARCF(1,NA2), PTXYZD )
C
C           LE PREMIER CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
C
C           L'ARETE PRECEDENTE NA01 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NA2
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NA2
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           L'ORIGINE DU CF
            N1ARCF( NBCF ) = NA01
C
C           LE SECOND CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV2 )
C           L'ARETE PRECEDENTE NA02 POINTE SUR LA NOUVELLE NAV2
            NOARCF( 2, NA02 ) = NAV2
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NA2 VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOARCF( 1, NA2 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NA02
            GOTO 400
C
C           ICI PAS D'ARETE ADMISSIBLE ET DANS LA TETRAEDRISATION
C           RECHERCHE DE LA PLUS PETITE ARETE DU CF
C           RECHERCHE DU SOMMET QUI RELIE A L'ARETE DONNE LA MEILLEURE QUALITE
C           ==================================================================
C           L'ETOILE EN HAUT DE PILE A POUR PREMIERE ARETE
C           RECHERCHE DE L'ARETE LA PLUS COURTE DU CF
 600        DMIMA = DINFO( 'GRAND' )
            NAQ   = 0
            NA01  = N1ARCF( NBCF )
            NA1   = NOARCF( 2, NA01 )
            NA00  = NA1
C
C           LES 2 SOMMETS DE L'ARETE NA1 DU CF
 610        NA2 = NOARCF( 2, NA1 )
            NS1 = NOARCF( 1, NA1 )
            NS2 = NOARCF( 1, NA2 )
            D   = (PTXYZD(1,NS2)-PTXYZD(1,NS1))**2
     %          + (PTXYZD(2,NS2)-PTXYZD(2,NS1))**2
     %          + (PTXYZD(3,NS2)-PTXYZD(3,NS1))**2
            IF( D .LT. DMIMA ) THEN
               DMIMA = D
               NAQ   = NA01
            ENDIF
            NA01 = NA1
            NA1  = NA2
            IF( NA1 .NE. NA00 ) GOTO 610
C
C           NAQ EST LA PLUS COURTE ARETE DU CF
C           RECHERCHE DU SOMMET NA3 DONNANT AVEC NAQ LA MEILLEURE QUALITE
            NA01 = NAQ
            NA1  = NOARCF( 2, NA01 )
            NA2  = NOARCF( 2, NA1 )
            NA3  = NA2
            QQ   = 0
C
 615        NA03 = NA3
            NA3  = NOARCF( 2, NA03 )
            IF( NA3 .NE. NA1 ) THEN
C              L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
               CALL ARTADM( EPSANG, NOARCF(1,NA1), NOARCF(1,NA3),
     %                      N1ARCF(NBCF), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 615
C              L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
               CALL ARTADM( EPSANG, NOARCF(1,NA2), NOARCF(1,NA3),
     %                      N1ARCF(NBCF), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 615
C
C              LES 3 ARETES SONT ADMISSIBLES
C              LA QUALITE DU TRIANGLE A FORMER
               CALL QUTRTE( PTXYZD(1,NOARCF(1,NA1)),
     %                      PTXYZD(1,NOARCF(1,NA2)),
     %                      PTXYZD(1,NOARCF(1,NA3)), Q )
               IF( Q .GT. QQ ) THEN
C                 MEILLEURE QUALITE
                  NAQ = NA03
                  QQ  = Q
               ENDIF
               GOTO 615
            ENDIF
C
            IF( QQ .LE. 0 ) GOTO 650
            NA03 = NAQ
            NA3  = NOARCF( 2, NA03 )
C
C           ICI LES 3 ARETES  NOARCF(1,NA1)-NOARCF(1,NA2)-NOARCF(1,NA3) ADMISSIB
CCCC           TRACE DES ARETES UTILISEES
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NA1), NOARCF(1,NA3), PTXYZD )
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NA2), NOARCF(1,NA3), PTXYZD )
C
            IF( NA01 .EQ. NA3 ) THEN
C              L'ARETE NA3 EST L'ARETE DU CF QUI PRECEDE NA1
C              TRAITEMENT AVEC L'ARETE NOARCF(1,NA2)-NOARCF(1,NA3)
               NA01 = NA1
               NA1  = NA2
            ENDIF
C
C           LE PREMIER CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
C
C           L'ARETE PRECEDENTE NA01 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NA3
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NA3
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           L'ORIGINE DU CF
            N1ARCF( NBCF ) = NA01
C
C           LE SECOND CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV2 )
C           L'ARETE PRECEDENTE NA03 POINTE SUR LA NOUVELLE NAV2
            NOARCF( 2, NA03 ) = NAV2
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NA3 VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOARCF( 1, NA3 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NA03
            GOTO 400
C
C           ICI PAS D'ARETE ADMISSIBLE ET DANS LA TETRAEDRISATION
C           RECHERCHE DE LA PLUS GRANDE ARETE DU CF
C           RECHERCHE DU SOMMET QUI RELIE A L'ARETE DONNE LA MEILLEURE QUALITE
C           ==================================================================
C           L'ETOILE EN HAUT DE PILE A POUR PREMIERE ARETE
C           RECHERCHE DE L'ARETE LA PLUS COURTE DU CF
 650        DMIMA = 0
            NAQ   = 0
            NA01  = N1ARCF( NBCF )
            NA1   = NOARCF( 2, NA01 )
            NA00  = NA1
C
C           LES 2 SOMMETS DE L'ARETE NA1 DU CF
 660        NA2 = NOARCF( 2, NA1 )
            NS1 = NOARCF( 1, NA1 )
            NS2 = NOARCF( 1, NA2 )
            D   = (PTXYZD(1,NS2)-PTXYZD(1,NS1))**2
     %          + (PTXYZD(2,NS2)-PTXYZD(2,NS1))**2
     %          + (PTXYZD(3,NS2)-PTXYZD(3,NS1))**2
            IF( D .GT. DMIMA ) THEN
               DMIMA = D
               NAQ   = NA01
            ENDIF
            NA01 = NA1
            NA1  = NA2
            IF( NA1 .NE. NA00 ) GOTO 660
C
C           NAQ EST LA PLUS GRANDE ARETE DU CF
C           RECHERCHE DU SOMMET NA3 DONNANT AVEC NAQ LA MEILLEURE QUALITE
            NA01 = NAQ
            NA1  = NOARCF( 2, NA01 )
            NA2  = NOARCF( 2, NA1 )
            NA3  = NA2
            QQ   = 0
C
 670        NA03 = NA3
            NA3  = NOARCF( 2, NA03 )
            IF( NA3 .NE. NA1 ) THEN
C              L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
               CALL ARTADM( EPSANG, NOARCF(1,NA1), NOARCF(1,NA3),
     %                      N1ARCF(NBCF), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 670
C              L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
               CALL ARTADM( EPSANG, NOARCF(1,NA2), NOARCF(1,NA3),
     %                      N1ARCF(NBCF), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .EQ. 0 ) GOTO 670
C
C              LES 3 ARETES SONT ADMISSIBLES
C              LA QUALITE DU TRIANGLE A FORMER
               CALL QUTRTE( PTXYZD(1,NOARCF(1,NA1)),
     %                      PTXYZD(1,NOARCF(1,NA2)),
     %                      PTXYZD(1,NOARCF(1,NA3)), Q )
               IF( Q .GT. QQ ) THEN
C                 MEILLEURE QUALITE
                  NAQ = NA03
                  QQ  = Q
               ENDIF
               GOTO 670
            ENDIF
C
            IF( QQ .LE. 0 ) GOTO 700
            NA03 = NAQ
            NA3  = NOARCF( 2, NA03 )
C
C           ICI LES 3 ARETES  NOARCF(1,NA1)-NOARCF(1,NA2)-NOARCF(1,NA3) ADMISSIB
CCCC           TRACE DES ARETES UTILISEES
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NA1), NOARCF(1,NA3), PTXYZD )
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NA2), NOARCF(1,NA3), PTXYZD )
C
            IF( NA01 .EQ. NA3 ) THEN
C              L'ARETE NA3 EST L'ARETE DU CF QUI PRECEDE NA1
C              TRAITEMENT AVEC L'ARETE NOARCF(1,NA2)-NOARCF(1,NA3)
               NA01 = NA1
               NA1  = NA2
            ENDIF
C
C           LE PREMIER CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
C
C           L'ARETE PRECEDENTE NA01 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NA3
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NA3
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           L'ORIGINE DU CF
            N1ARCF( NBCF ) = NA01
C
C           LE SECOND CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV2 )
C           L'ARETE PRECEDENTE NA03 POINTE SUR LA NOUVELLE NAV2
            NOARCF( 2, NA03 ) = NAV2
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NA3 VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOARCF( 1, NA3 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NA03
            GOTO 400
C
C           TENTATIVE DE LA DERNIERE CHANCE : JONCTION DE 2 ARETES
C           CONSECUTIVES DU CF AYANT UNE ARETE DE JONCTION ADMISSIBLE
C           ET CHOIX DU TRIANGLE DE MEILLEURE QUALITE
C           =========================================================
 700        NA01 = N1ARCF( NBCF )
            NA1  = NOARCF( 2, NA01 )
            NA00 = NA1
            NA2  = NOARCF( 2, NA1 )
            QQ   = 0
C
 710        NA3  = NOARCF( 2, NA2 )
            IF( NA3 .NE. NA00 ) THEN
C              L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
               CALL ARTADM( EPSANG, NOARCF(1,NA1), NOARCF(1,NA3),
     %                      N1ARCF(NBCF), NOARCF, PTXYZD, NONOUI )
               IF( NONOUI .NE. 0 ) GOTO 730
C
C              PASSAGE AU COUPLE D'ARETES SUIVANTES
 720           NA01 = NA1
               NA1  = NA2
               NA2  = NA3
               GOTO 710
C
C              LES 3 ARETES SONT ADMISSIBLES
C              LA QUALITE DU TRIANGLE A FORMER
 730           CALL QUTRTE( PTXYZD(1,NOARCF(1,NA1)),
     %                      PTXYZD(1,NOARCF(1,NA2)),
     %                      PTXYZD(1,NOARCF(1,NA3)), Q )
               IF( Q .GT. QQ ) THEN
C                 MEILLEURE QUALITE
                  NAQ = NA01
                  QQ  = Q
               ENDIF
               GOTO 720
            ENDIF
C
            IF( QQ .LE. 0 ) THEN
C              ***********************************************************
C              2 ARETES CONSECUTIVES NE FORMENT PAS UN TRIANGLE ADMISSIBLE
C              ALGORITHME INSUFFISANT  A AMELIORER !
C              ***********************************************************
               WRITE(IMPRIM,*) 'LIERCF: 2 ARETES CONSECUTIVES NE FORMENT
     % PAS UN TRIANGLE ADMISSIBLE'
               WRITE(IMPRIM,*) 'LIERCF: ALGORITHME INSUFFISANT'
               CALL CLICSO
               IERR = 2
               GOTO 999
            ENDIF
C
C           ICI 3 ARETES ADMISSIBLES => 2 CF DONT UN TRIANGLE
            NA01 = NAQ
            NA1  = NOARCF( 2, NA01)
            NA2  = NOARCF( 2, NA1 )
            NA03 = NA2
            NA3  = NOARCF( 2, NA03 )
C
CCCC           TRACE DE L'ARETE UTILISEE
CCC            CALL TRARTR( NCMAGE, NOARCF(1,NA1), NOARCF(1,NA3), PTXYZD )
C
C           LE PREMIER CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV1 = N1ARCF(0)
            IF( NAV1 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV1 )
C
C           L'ARETE PRECEDENTE NA01 POINTE SUR LA NOUVELLE NAV1
            NOARCF( 2, NA01 ) = NAV1
C
C           AJOUT DE L'ARETE NAV1 PARTANT DE NA1 VERS NA3
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV1 ) = NOARCF( 1, NA1 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV1 ) = NA3
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV1 ) = -1
C
C           L'ORIGINE DU CF
            N1ARCF( NBCF ) = NA01
C
C           LE SECOND CF
C           RECHERCHE DE 1 ARETE VIDE DE CF
            NAV2 = N1ARCF(0)
            IF( NAV2 .LE. 0 ) GOTO 997
C           LA 1-ERE ARETE VIDE EST MISE A JOUR
            N1ARCF(0) = NOARCF( 2, NAV2 )
C           L'ARETE PRECEDENTE NA03 POINTE SUR LA NOUVELLE NAV2
            NOARCF( 2, NA03 ) = NAV2
C
C           AJOUT DE L'ARETE NAV2 PARTANT DE NA3 VERS NA1
C           LE NUMERO DU SOMMET
            NOARCF( 1, NAV2 ) = NOARCF( 1, NA3 )
C           L'ARETE SUIVANTE
            NOARCF( 2, NAV2 ) = NA1
C           LE TRIANGLE LEFACO OPPOSE EST INCONNU
            NOARCF( 3, NAV2 ) = -1
C
C           L'ORIGINE DU CF
            IF( NBCF .GE. MXTRCF ) GOTO 995
            NBCF = NBCF + 1
            N1ARCF( NBCF ) = NA03
            GOTO 400
         ENDIF
C        ICI LE CF NBCF A ETE TOTALEMENT TRAITE
         GOTO 400
C
C        MISE A JOUR DU CHAINAGE DES FACES ADJACENTES PAR LES ARETES
C        -----------------------------------------------------------
 800     DO 990 NFP0 = 1, NBTRC1
C           LE NUMERO DE LA FACE DANS LEFACO
            NF0 = NOTRCF( NFP0 )
C           BOUCLE SUR SES 3 ARETES
            DO 900 I=1,3
C              SEULE UNE ARETE SANS TRIANGLE OPPOSE EST TRAITEE
               IF( LEFACO( 5+I, NF0 ) .GT. 0 ) GOTO 900
C              LES 2 SOMMETS
               NS1 = LEFACO(I,NF0)
               IF( I .EQ. 3 ) THEN
                  NS2 = 1
               ELSE
                  NS2 = I + 1
               ENDIF
               NS2 = LEFACO(NS2,NF0)
C
C              RECHERCHE DE L'ARETE NS1-NS2 DANS LES FACES LEFACO
C              ARRET A LA PREMIERE RENCONTREE CAR FACES COPLANAIRES
C              => 2 FACES ADJACENTES AU PLUS PAR ARETE
               DO 830 NFP1=1,NBTRC1
                  NF = NOTRCF( NFP1 )
                  IF( NF .EQ. NF0 ) GOTO 830
                  DO 810 J=1,3
C                    LES SOMMETS DE L'ARETE
                     NS3 = LEFACO(J,NF)
                     IF( J .EQ. 3 ) THEN
                        NS4 = 1
                     ELSE
                        NS4 = J + 1
                     ENDIF
                     NS4 = LEFACO(NS4,NF)
                     IF( NS3 .EQ. NS1 ) THEN
                        IF( NS4 .EQ. NS2 ) GOTO 850
                     ELSE IF( NS3 .EQ. NS2 ) THEN
                        IF( NS4 .EQ. NS1 ) GOTO 850
                     ENDIF
C                    ARETE NON RETROUVEE
 810              CONTINUE
C                 PASSAGE A LA FACE SUIVANTE
 830           CONTINUE
C
C              PROBLEME : ARETE NON RETROUVEE : SURFACE OUVERTE
               WRITE(IMPRIM,*)
               NBLGRC(NRERR) = 3
               WRITE(KERR(MXLGER)( 1:10),'(I10)') NS1
               WRITE(KERR(MXLGER)(11:20),'(I10)') NS2
               KERR(1) ='LIERCF:L''ARETE' // KERR(MXLGER)(1:20)
               KERR(2) ='APPARTIENT A UN SEUL TRIANGLE'
               KERR(3) ='SURFACE OUVERTE A REFERMER'
               CALL LEREUR
               WRITE(IMPRIM,*)'FACE',NF0,':',(LEFACO(J,NF0),J=1,10)
               WRITE(IMPRIM,*)
               IERR = IERR + 1000
               GOTO 900
C
C              ARETE RETROUVEE ( NF0 , I ) <-> ( NF , J )
C              LES 2 FACES SONT COPLANAIRES ET DE MEMES SURFACES
C              => AU PLUS 2 PAR ARETE
 850           LEFACO(5+I,NF0) = NF
               LEFACO(5+J,NF ) = NF0
C
C              VERIFICATION QUE LES SURFACES NE SONT PAS DISPROPORTIONNEES
C              NF0 (NS1, NS2, NS3) ET NF(NS1 NS2 NS4) FORMENT UN QUADRANGLE
C              NS3 LE SOMMET RESTANT DE NF0
               NS3 = I - 1
               IF( NS3 .EQ. 0 ) NS3 = 3
               NS3 = LEFACO(NS3,NF0)
C
C              NS4 LE SOMMET RESTANT DE NF
               DO 860 K=1,3
                  NS4 = LEFACO(K,NF)
                  IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 870
 860           CONTINUE
C
C              LES 2 SURFACES SONT ELLES DISPROPORTIONNEES ?
 870           S0  = SURTRD(PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3))
               S1  = SURTRD(PTXYZD(1,NS2), PTXYZD(1,NS1), PTXYZD(1,NS4))
               IF( S0 .LT. 1E-2*S1 .OR. S1 .LT. 1E-2*S0 ) THEN
C                 OUI: ECHANGE DES 2 TRIANGLES
                  CALL CH2F2F(I,NF0,NF, PTXYZD, MXFACO, LEFACO, N1FASC,
     %                        NF2, NF3)
                  GOTO 900
               ENDIF
C
               CALL VDARTE( NS1, NS2, N1TETS, NOTETR,
     %                      MXTEFA, NOTEFA, NONOUI )
               IF( NONOUI .EQ. 0 ) THEN
C              L'ARETE NS1 NS2 N'EST PAS DANS LA TETRAEDRISATION
C              LA QUALITE DES 2 TRIANGLES INITIAUX < EVENTUELS ?
               CALL QUTRTE(PTXYZD(1,NS1),PTXYZD(1,NS2),PTXYZD(1,NS3),Q1)
               CALL QUTRTE(PTXYZD(1,NS1),PTXYZD(1,NS2),PTXYZD(1,NS4),Q2)
C              LA QUALITE DES 2 TRIANGLES EVENTUELS
               CALL QUTRTE(PTXYZD(1,NS3),PTXYZD(1,NS4),PTXYZD(1,NS1),Q3)
               CALL QUTRTE(PTXYZD(1,NS3),PTXYZD(1,NS4),PTXYZD(1,NS2),Q4)
               IF( MIN(Q1,Q2) .LT. MIN(Q3,Q4) ) THEN
C                 OUI: ECHANGE DES 2 TRIANGLES
                  CALL CH2F2F(I,NF0,NF, PTXYZD, MXFACO, LEFACO, N1FASC,
     %                        NF2, NF3)
               ENDIF
               ENDIF
 900        CONTINUE
 990     CONTINUE
C
C        FIN DU TRAITEMENT DE CE CF
         GOTO 1
      ENDIF
      GOTO 999
C
 995  NBLGRC(NRERR) = 1
      KERR(1) = 'LIERCF: SATURATION DU NOMBRE DE CF'
      CALL CLICSO
      IERR = 1
      GOTO 999
C
 996  NBLGRC(NRERR) = 1
      KERR(1) = 'LIERCF: SATURATION DES FACES LEFACO'
      IERR = 1
      CALL CLICSO
      GOTO 999
C
 997  NBLGRC(NRERR) = 1
      KERR(1) = 'LIERCF: SATURATION DU NOMBRE D''ARETES DE CF'
      CALL CLICSO
      GOTO 1
C
 998  NBLGRC(NRERR) = 1
      KERR(1) = 'LIERCF: SATURATION DES TRIANGLES DU CF'
      CALL LEREUR
      IERR = 1
      CALL CLICSO
C
CCCC     LES LIGNES SONT A TRACER EN BLANC
CCC 999  CALL XVCOULEUR( NCBLAN )
C
 999  RETURN
      END
