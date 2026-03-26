      SUBROUTINE TEFOAR( NARETE, NBARPI, PXYD,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, MXARTR, N1ARTR, NOARTR, NOARST,
     %                   MXARCF, N1ARCF, NOARCF, LARMIN, NOTRCF,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   FORCER L'ARETE NARETE DE NOSOAR DANS LA TRIANGULATION ACTUELLE
C -----   TRIANGULATION FRONTALE POUR LA REOBTENIR

C         ATTENTION: LE CHAINAGE LCHAIN(=6) DE NOSOAR DEVIENT ACTIF
C                    DURANT LA FORMATION DES CONTOURS FERMES (CF)

C ENTREES:
C --------
C NARETE : NUMERO NOSOAR DE L'ARETE FRONTALIERE A FORCER
C NBARPI : NUMERO DU DERNIER POINT INTERNE IMPOSE PAR L'UTILISATEUR
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE
C MXSOAR : NOMBRE MAXIMAL D'ARETES STOCKABLES DANS LE TABLEAU NOSOAR
C          ATTENTION: MXSOAR>3*MXSOMM OBLIGATOIRE!
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C MXARTR : NOMBRE MAXIMAL DE TRIANGLES STOCKABLES DANS LE TABLEAU NOARTR

C MODIFIES:
C ---------
C N1SOAR : NO DE L'EVENTUELLE PREMIERE ARETE LIBRE DANS LE TABLEAU NOSOAR
C          CHAINAGE DES VIDES SUIVANT EN 3 ET PRECEDANT EN 2 DE NOSOAR
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          HACHAGE DES ARETES = NOSOAR(1)+NOSOAR(2)*2
C          AVEC MXSOAR>=3*MXSOMM
C          UNE ARETE I DE NOSOAR EST VIDE <=> NOSOAR(1,I)=0 ET
C          NOSOAR(2,ARETE VIDE)=L'ARETE VIDE QUI PRECEDE
C          NOSOAR(3,ARETE VIDE)=L'ARETE VIDE QUI SUIT
C N1ARTR : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOARTR
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOARTR(2,.)
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1 = 0 SI TRIANGLE VIDE => ARETE2 = TRIANGLE VIDE SUIVANT
C NOARST : NOARST(I) NUMERO D'UNE ARETE DE SOMMET I

C MXARCF : NOMBRE DE VARIABLES DES TABLEAUX N1ARCF, NOARCF, LARMIN, NOTRCF

C TABLEAUX AUXILIAIRES :
C ----------------------
C N1ARCF : TABLEAU (0:MXARCF) AUXILIAIRE
C NOARCF : TABLEAU (3,MXARCF) AUXILIAIRE
C LARMIN : TABLEAU (MXARCF)   AUXILIAIRE
C NOTRCF : TABLEAU (1:MXARCF) AUXILIAIRE

C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          2 NS1 DANS AUCUN TRIANGLE
C          9 TABLEAU NOSOAR DE TAILLE INSUFFISANTE CAR TROP D'ARETES
C            A PROBLEME
C          10 UN DES TABLEAUX N1ARCF, NOARCF NOTRCF EST SATURE
C             AUGMENTER A L'APPEL MXARCF
C          11 ALGORITHME DEFAILLANT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC     MARS    1997
C MODIFS : ALAIN PERRONNET LABORATOIRE JL LIONS UPMC PARIS  OCTOBRE 2006
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2019
C....................................................................012
      INCLUDE "./incl/trvari.inc"
      PARAMETER        (MXPITR=32, mxstpe=512)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      LOGICAL           TRATRI, TRATRI0
      COMMON / DV2DCO / TRATRI
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NOARTR(MOARTR,MXARTR),
     %                  NOSOAR(MOSOAR,MXSOAR),
     %                  NOARST(*),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(3,MXARCF),
     %                  LARMIN(MXARCF),
     %                  NOTRCF(MXARCF),
     %                  nostpe(mxstpe)

      CHARACTER*80      KTITRE
      INTEGER           LAPITR(MXPITR)
      DOUBLE PRECISION  X1,Y1,X2,Y2,D12,D3,D4,X,Y,D,DMIN
      INTEGER           NOSOTR(3), NS(2)
      INTEGER           NACF(1:2), NACF1, NACF2
      EQUIVALENCE      (NACF(1),NACF1), (NACF(2),NACF2)

      IERR = 0
      TRATRI0 = TRATRI

C     TRAITEMENT DE CETTE ARETE PERDUE
      NS1 = NOSOAR( 1, NARETE )
      NS2 = NOSOAR( 2, NARETE )

CCC      WRITE(IMPRIM,*)
CCC      WRITE(IMPRIM,*) 'tefoar RECONSTRUCTION DE L''ARETE ',NS1,' ', NS2
CCC      WRITE(IMPRIM,*) 'SOMMET',NS1,' X=',PXYD(1,NS1),' Y=',PXYD(2,NS1)
CCC      WRITE(IMPRIM,*) 'SOMMET',NS2,' X=',PXYD(1,NS2),' Y=',PXYD(2,NS2)

      IF( TRATRI ) THEN
C        LES TRACES SONT DEMANDES
         CALL EFFACE
C        LE CADRE OBJET GLOBAL EN UNITES UTILISATEUR
         XX1 = REAL( MIN( PXYD(1,NS1), PXYD(1,NS2) ) )
         XX2 = REAL( MAX( PXYD(1,NS1), PXYD(1,NS2) ) )
         YY1 = REAL( MIN( PXYD(2,NS1), PXYD(2,NS2) ) )
         YY2 = REAL( MAX( PXYD(2,NS1), PXYD(2,NS2) ) )
         DD  = MAX( XX2-XX1, YY2-YY1 ) * 1.6
         IF( XX1 .GE. XX2 ) XX2 = XX1 + DD
         IF( YY1 .GE. YY2 ) YY2 = YY1 + DD*0.5
         CALL ISOFENETRE( XX1-(XX2-XX1), XX2+(XX2-XX1),
     %                    YY1-(YY2-YY1), YY2+(YY2-YY1) )
C        TRACE DE L'ARETE PERDUE
         CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
      ENDIF

C     LE SOMMET NS2 EST IL CORRECT?
      NA = NOARST( NS2 )
      IF( NA .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'tefoar: ERREUR SOMMET ',NS2,' SANS ARETE'
         IERR = 8
         GOTO 9999
      ENDIF
      IF( NOSOAR(4,NA) .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'tefoar: ERREUR SOMMET ',NS2,
     %                   ' DANS AUCUN TRIANGLE'
         IERR = 8
         GOTO 9999
      ENDIF

C     LE PREMIER PASSAGE: RECHERCHE DANS LE SENS NS1->NS2
      IPAS = 0

C     RECHERCHE DES TRIANGLES INTERSECTES PAR LE SEGMENT NS1-NS2
C     ==========================================================
 3    X1  = PXYD(1,NS1)
      Y1  = PXYD(2,NS1)
      X2  = PXYD(1,NS2)
      Y2  = PXYD(2,NS2)
      D12 = SQRT( (X2-X1)**2 + (Y2-Y1)**2 )
C
C     RECHERCHE DU TRIANGLE VOISIN DANS LE SENS INDIRECT DE ROTATION
      NSENS = -1
C
C     RECHERCHE DU NO LOCAL DU SOMMET NS1 DANS L'UN DE SES TRIANGLES
 10   NA01 = NOARST( NS1 )
      IF( NA01 .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'tefoar: SOMMET ',NS1,' SANS ARETE'
         IERR = 8
         GOTO 9999
      ENDIF
      NT0 = NOSOAR(4,NA01)
      IF( NT0 .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'tefoar: SOMMET ',NS1,' DANS AUCUN TRIANGLE'
         IERR = 8
         GOTO 9999
      ENDIF
C
C     LE NUMERO DES 3 SOMMETS DU TRIANGLE NT0 DANS LE SENS DIRECT
 20   CALL NUSOTR( NT0, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR )
      DO 22 NA00=1,3
         IF( NOSOTR(NA00) .EQ. NS1 ) GOTO 26
 22   ENDDO
C
 25   IF( IPAS .EQ. 0 ) THEN
C        LE SECOND PASSAGE: RECHERCHE DANS LE SENS NS2->NS1
C        TENTATIVE D'INVERSION DES 2 SOMMETS EXTREMITES DE L'ARETE A FORCER
         NA00 = NS1
         NS1  = NS2
         NS2  = NA00
         IPAS = 1
         GOTO 3
      ELSE
C        LES SENS NS1->NS2 ET NS2->NS1 NE DONNE PAS DE SOLUTION!
         WRITE(IMPRIM,*)'tefoar:ARETE ',NS1,' - ',NS2,' A IMPOSER'
         WRITE(IMPRIM,*)'tefoar:ANOMALIE SOMMET ',NS1,
     %   'NON DANS LE TRIANGLE DE SOMMETS ',(NOSOTR(I),I=1,3)
         IERR = 11
         GOTO 9999
      ENDIF
C
C     LE NUMERO DES ARETES SUIVANTE ET PRECEDENTE
 26   NA0 = NOSUI3( NA00 )
      NA1 = NOPRE3( NA00 )
      NS3 = NOSOTR( NA0 )
      NS4 = NOSOTR( NA1 )
C
C     TRACE DU TRIANGLE NT0 ET DE L'ARETE PERDUE
      CALL MTTRTR( PXYD, NT0, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %             NCBLAN, NCJAUN )
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
      CALL DVTRAR( PXYD, NS3, NS4, NCBLEU, NCCYAN )
C
C     POINT D'INTERSECTION DU SEGMENT NS1-NS2 AVEC L'ARETE NS3-NS4
C     ------------------------------------------------------------
      CALL INT1SD( NS1, NS2, NS3, NS4, PXYD, LINTER, X1, Y1 )
      IF( LINTER .LE. 0 ) THEN
C
C        PAS D'INTERSECTION: ROTATION AUTOUR DU POINT NS1
C        POUR TROUVER LE TRIANGLE DE L'AUTRE COTE DE L'ARETE NA01
         IF( NSENS .LT. 0 ) THEN
C           SENS INDIRECT DE ROTATION: L'ARETE DE SOMMET NS1
            NA01 = ABS( NOARTR(NA00,NT0) )
         ELSE
C           SENS DIRECT DE ROTATION: L'ARETE DE SOMMET NS1 QUI PRECEDE
            NA01 = ABS( NOARTR(NA1,NT0) )
         ENDIF
C        LE TRIANGLE DE L'AUTRE COTE DE L'ARETE NA01
         IF( NOSOAR(4,NA01) .EQ. NT0 ) THEN
            NT0 = NOSOAR(5,NA01)
         ELSE
            NT0 = NOSOAR(4,NA01)
         ENDIF
         IF( NT0 .GT. 0 ) GOTO 20
C
C        LE PARCOURS SORT DU DOMAINE
C        IL FAUT TOURNER DANS L'AUTRE SENS AUTOUR DE NS1
         IF( NSENS .LT. 0 ) THEN
            NSENS = 1
            GOTO 10
         ENDIF
C
C        DANS LES 2 SENS, PAS D'INTERSECTION => IMPOSSIBLE
C        ESSAI AVEC L'ARETE INVERSEE NS1 <-> NS2
         IF( IPAS .EQ. 0 ) GOTO 25
         WRITE(IMPRIM,*) 'tefoar: ARETE ',NS1,' ',NS2,
     %  ' SANS INTERSECTION AVEC LES TRIANGLES ACTUELS'
         WRITE(IMPRIM,*) 'REVOYEZ LES LIGNES DU CONTOUR'
         IERR = 11
         GOTO 9999
      ENDIF
C
C     IL EXISTE UNE INTERSECTION AVEC L'ARETE OPPOSEE AU SOMMET NS1
C     =============================================================
C     NBTRCF : NOMBRE DE TRIANGLES DU CF
      NBTRCF = 1
      NOTRCF( 1 ) = NT0
C
C     LE TRIANGLE OPPOSE A L'ARETE NA0 DE NT0
 30   NOAR = ABS( NOARTR(NA0,NT0) )
      IF( NOSOAR(4,NOAR) .EQ. NT0 ) THEN
         NT1 = NOSOAR(5,NOAR)
      ELSE
         NT1 = NOSOAR(4,NOAR)
      ENDIF
      IF( NT1 .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'ERREUR DANS tefoar NT1=',NT1
         READ(LECTEU,*) J
      ENDIF
C
C     TRACE DU TRIANGLE NT1 ET DE L'ARETE PERDUE
      CALL MTTRTR( PXYD, NT1, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %             NCJAUN, NCMAGE )
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
C
C     LE NUMERO DES 3 SOMMETS DU TRIANGLE NT1 DANS LE SENS DIRECT
      CALL NUSOTR( NT1, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR )
C
C     LE TRIANGLE NT1 CONTIENT IL NS2 ?
      DO 32 J=1,3
         IF( NOSOTR(J) .EQ. NS2 ) GOTO 70
 32   ENDDO
C
C     RECHERCHE DE L'ARETE NOAR, NA1 DANS NT1 QUI EST L'ARETE NA0 DE NT0
      DO 34 NA1=1,3
         IF( ABS( NOARTR(NA1,NT1) ) .EQ. NOAR ) GOTO 35
 34   ENDDO
C
C     TRACE DU TRIANGLE NT1 ET DE L'ARETE PERDUE
 35   CALL MTTRTR( PXYD, NT1, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %             NCJAUN, NCMAGE )
      CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )

C     RECHERCHE DE L'INTERSECTION DE NS1-NS2 AVEC LES 2 AUTRES ARETES DE NT1
C     ======================================================================
      NA2 = NA1
      DO 50 I1 = 1,2
C        L'ARETE SUIVANTE
         NA2 = NOSUI3(NA2)
C
C        LES 2 SOMMETS DE L'ARETE NA2 DE NT1
         NOAR = ABS( NOARTR(NA2,NT1) )
         NS3  = NOSOAR( 1, NOAR )
         NS4  = NOSOAR( 2, NOAR )
         CALL DVTRAR( PXYD, NS3, NS4, NCBLEU, NCCYAN )

C        LE POINT D'INTERSECTION DU SEGMENT NS1-NS2 AVEC L'ARETE NS3-NS4
C        ---------------------------------------------------------------
         CALL INT1SD( NS1, NS2, NS3, NS4, PXYD, LINTER, X , Y )
         IF( LINTER .GT. 0 ) THEN

C           LES 2 ARETES S'INTERSECTENT EN (X,Y)
C           DISTANCE DE (X,Y) A NS3 ET NS4
            D3 = (PXYD(1,NS3)-X)**2 + (PXYD(2,NS3)-Y)**2
            D4 = (PXYD(1,NS4)-X)**2 + (PXYD(2,NS4)-Y)**2

C           NSP EST LE POINT LE PLUS PROCHE DE (X,Y)
            IF( D3 .LT. D4 ) THEN
               NSP = NS3
               D   = D3
            ELSE
               NSP = NS4
               D   = D4
            ENDIF
            D = SQRT( D )
            IF( D .GE. 1D-5 * D12 ) GOTO 60

C           ICI LE SOMMET NSP EST PRATIQUEMENT SUR L'ARETE PERDUE NS1-NS2
C           -------------------------------------------------------------
            IF( NSP .LE. NBARPI ) THEN

C              C'EST UN POINT UTILISATEUR OU FRONTALIER DONC NON SUPPRIMABLE
C              IL EST FAIBLEMENT DEPLACER SUR L'ARETE NS3-NS4
               WRITE(IMPRIM,*)
               WRITE(IMPRIM,*) 'tefoar: Probleme au SOMMET NSP=',NSP,
     %' FRONTALIER TROP PROCHE de l''ARETE PERDUE NS1=',NS1,'-NS2=',NS2
        WRITE(IMPRIM,*) 'NSP:',NSP,' X=', PXYD(1,NSP),' Y=', PXYD(2,NSP)
        WRITE(IMPRIM,*) 'NS1:',NS1,' X=', PXYD(1,NS1),' Y=', PXYD(2,NS1)
        WRITE(IMPRIM,*) 'NS2:',NS2,' X=', PXYD(1,NS2),' Y=', PXYD(2,NS2)
            WRITE(IMPRIM,*) 'ARETE NS1-NS2 COUPE ARETE NS3-NS4 EN (X,Y)'
       WRITE(IMPRIM,*) 'NS3=',NS3,': X=', PXYD(1,NS3),' Y=', PXYD(2,NS3)
       WRITE(IMPRIM,*) 'NS4=',NS4,': X=', PXYD(1,NS4),' Y=', PXYD(2,NS4)
               WRITE(IMPRIM,*) 'INTERSECTION EN  : X=', X, ' Y=', Y
               WRITE(IMPRIM,*) 'DISTANCE NS1-NS2=', D12,
     %                        ' DISTANCE (X,Y) au PLUS PROCHE NS3 NS4=',
     %                          D,' D/D12=',D/D12,'<1D-5'

C              LES TRACES SONT DEMANDES
               TRATRI0 = TRATRI
               TRATRI = .TRUE.
C              LE CADRE OBJET GLOBAL EN UNITES UTILISATEUR
               XX1 = REAL( MIN( PXYD(1,NS1), PXYD(1,NS2) ) )
               XX1 = MIN( XX1, REAL( MIN( PXYD(1,NS3), PXYD(1,NS4) ) ) )
               XX2 = REAL( MAX( PXYD(1,NS1), PXYD(1,NS2) ) )
               XX2 = MAX( XX2, REAL( MAX( PXYD(1,NS3), PXYD(1,NS4) ) ) )

               YY1 = REAL( MIN( PXYD(2,NS1), PXYD(2,NS2) ) )
               YY1 = MIN( YY1, REAL( MIN( PXYD(2,NS3), PXYD(2,NS4) ) ) )
               YY2 = REAL( MAX( PXYD(2,NS1), PXYD(2,NS2) ) )
               YY2 = MAX( YY2, REAL( MAX( PXYD(2,NS3), PXYD(2,NS4) ) ) )

               DD  = MAX( XX2-XX1, YY2-YY1 ) * 1.6
               IF( XX1 .GE. XX2 ) XX2 = XX1 + DD
               IF( YY1 .GE. YY2 ) YY2 = YY1 + DD*0.5
               CALL ISOFENETRE( XX1-(XX2-XX1), XX2+(XX2-XX1),
     %                          YY1-(YY2-YY1), YY2+(YY2-YY1) )

C              TRACE DES TRIANGLES NT0 NT1 ET DES ARETES
               CALL MTTRTR( PXYD, NT0, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                      NCBLAN, NCMAGE )
               CALL MTTRTR( PXYD, NT1, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                      NCJAUN, NCCYAN )
C              TRACE DE L'ARETE PERDUE
               CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCNOIR )
C              TRACE DE L'ARETE NS3-NS4
               CALL DVTRAR( PXYD, NS3, NS4, NCBLEU, NCGRIS )

       KTITRE='tefoar: PB ARETE        -        PROCHE       ou        '
               WRITE(KTITRE(18:23),'(I6)') NS1
               WRITE(KTITRE(27:32),'(I6)') NS2               
               WRITE(KTITRE(41:46),'(I6)') NS3
               WRITE(KTITRE(50:55),'(I6)') NS4
               CALL SANSDBL( KTITRE, NC )
               CALL TRFINS( KTITRE(1:NC) )
               TRATRI = TRATRI0

ccc               IERR = 11
ccc               GOTO 9999

C              FAIBLE DEPLACEMENT DU SOMMET NSP SUR L'ARETE NS3-NS4
C              ----------------------------------------------------
               WRITE(IMPRIM,*) 'tefoar: le SOMMET',NSP,' XYD=',
     %                          (PXYD(K,NSP),K=1,3)
               IF( NSP .EQ. NS3 ) THEN
                  NSL = NS4
               ELSE
                  NSL = NS3
               ENDIF
               D3 = ( PXYD(1,NSL) - PXYD(1,NSP) ) ** 2
     %            + ( PXYD(2,NSL) - PXYD(2,NSP) ) ** 2
               D3 = SQRT( D3 )
               DO K=1,2
                  PXYD(K,NSP) = PXYD(K,NSP)
     %              + 1D-4 * D12 * ( PXYD(K,NSL) - PXYD(K,NSP) ) / D3
               ENDDO
               WRITE(IMPRIM,*) 'DEVIENT le SOMMET',NSP,' XYD=',
     %                          (PXYD(K,NSP),K=1,3)
               WRITE(IMPRIM,*)
               GOTO 60

            ENDIF

C           LE SOMMET INTERNE NSP EST SUPPRIME EN METTANT TOUS LES TRIANGLES
C           L'AYANT COMME SOMMET DANS LA PILE NOTRCF DES TRIANGLES A SUPPRIMER
C           ------------------------------------------------------------------
CCC            WRITE(IMPRIM,*) 'tefoar: LE SOMMET ',NSP,' EST SUPPRIME'
C           CONSTRUCTION DE LA LISTE DES TRIANGLES DE SOMMET NSP
            CALL TRP1ST( NSP,    NOARST, MOSOAR, NOSOAR,
     %                   MOARTR, MXARTR, NOARTR,
     %                   MXPITR, NBT, LAPITR )
            IF( NBT .LE. 0 ) THEN
C              LES TRIANGLES DE SOMMET NSP NE FORME PAS UNE "BOULE"
C              AVEC CE SOMMET NSP POUR "CENTRE"
               WRITE(IMPRIM,*)
     %        'tefoar: LES TRIANGLES AUTOUR DU SOMMET ',NSP,
     %        ' NE FORME PAS UNE ETOILE'
C              TRACE DES TRIANGLES DE L'ETOILE DU SOMMET NSP
               NBT = -NBT
CCC               TRATRI = .TRUE.
               CALL TRPLTR( NBT, LAPITR, PXYD,
     %                      MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                      NCROUG, NCBLAN )
CCC               TRATRI = .FALSE.
CCC               IERR = 11
CCC               GOTO 9999
               DO 36 J=1,NBT
                  NT = LAPITR(J)
C                 LE NUMERO DES 3 SOMMETS DU TRIANGLE NT
         CALL MT3STR( NT, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                NSS1, NSS2, NSS3 )
C
C        LES COORDONNEES DES 3 SOMMETS DU TRIANGLE NT
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,*)'TRIANGLE',NT,' SOMMET',NSS1,' U=',PXYD(1,NSS1)
     %                 ,' V=',PXYD(2,NSS1)
         WRITE(IMPRIM,*)'TRIANGLE',NT,' SOMMET',NSS2,' U=',PXYD(1,NSS2)
     %                 ,' V=',PXYD(2,NSS2)
         WRITE(IMPRIM,*)'TRIANGLE',NT,' SOMMET',NSS3,' U=',PXYD(1,NSS3)
     %                 ,' V=',PXYD(2,NSS3)
 36             ENDDO
            ENDIF
C
C           AJOUT DES TRIANGLES DE SOMMET NS1 A NOTRCF
            NBTRC0 = NBTRCF
            DO 38 J=1,NBT
               NT = LAPITR(J)
               DO 37 K=NBTRCF,1,-1
                  IF( NT .EQ. NOTRCF(K) ) GOTO 38
 37            ENDDO
C              TRIANGLE AJOUTE
               NBTRCF = NBTRCF + 1
               NOTRCF( NBTRCF ) = NT
               CALL MTTRTR( PXYD, NT, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                      NCJAUN, NCMAGE )
               CALL DVTRAR( PXYD, NS1, NS2, NCROUG, NCBLAN )
 38         ENDDO
C
C           CE SOMMET SUPPRIME N'APPARTIENT PLUS A AUCUN TRIANGLE
            NOARST( NSP ) = 0
C
C           NS2 EST-IL UN SOMMET DES TRIANGLES EMPILES?
C           -------------------------------------------
            DO 40 NT=NBTRC0+1,NBTRCF
C              LE TRIANGLE A SUPPRIMER NT
               NT1 = NOTRCF( NT )
C              LE NUMERO DES 3 SOMMETS DU TRIANGLE NT1 DANS LE SENS DIRECT
               CALL NUSOTR( NT1, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR)
               DO 39 K=1,3
C                 LE SOMMET K DE NT1
                  IF( NOSOTR( K ) .EQ. NS2 ) THEN
C                    BUT ATTEINT
                     GOTO 80
                  ENDIF
 39            ENDDO
 40         ENDDO
C
C           RECHERCHE DU PLUS PROCHE POINT D'INTERSECTION DE NS1-NS2
C           PAR RAPPORT A NS2 AVEC LES ARETES DES TRIANGLES AJOUTES
            NT0  = 0
            DMIN = (D12**2) * 10000
            DO 48 NT=NBTRC0+1,NBTRCF
               NT1 = NOTRCF( NT )
C              LE NUMERO DES 3 SOMMETS DU TRIANGLE NT1 DANS LE SENS DIRECT
               CALL NUSOTR( NT1, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR)
               DO 45 K=1,3
C                 LES 2 SOMMETS DE L'ARETE K DE NT
                  NS3 = NOSOTR( K )
                  NS4 = NOSOTR( NOSUI3(K) )
C
C                 POINT D'INTERSECTION DU SEGMENT NS1-NS2 AVEC L'ARETE NS3-NS4
C                 ------------------------------------------------------------
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
 45            ENDDO
 48         ENDDO
C
C           REDEMARRAGE AVEC LE TRIANGLE NT0 ET L'ARETE NA0
            IF( NT0 .GT. 0 ) GOTO 30
C
            WRITE(IMPRIM,*) 'tefoar: ALGORITHME DEFAILLANT'
            IERR = 11
            GOTO 9999
         ENDIF
 50   ENDDO
C
C     PAS D'INTERSECTION DIFFERENTE DE L'INITIALE => SOMMET SUR NS1-NS2
C     TENTATIVE D'INVERSION DES SOMMETS DE L'ARETE NS1-NS2
      IF( IPAS .EQ. 0 ) GOTO 25
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'tefoar: REVOYEZ VOS DONNEES'
      WRITE(IMPRIM,*) 'LES LIGNES FERMEES DOIVENT ETRE DISJOINTES'
      WRITE(IMPRIM,*) 'LA FONCTION TAILLE_IDEALE PAS TROP DISCONTINUE'
      IF( TRATRI ) THEN
         CALL MEMPXFENETRE
         CALL XVVOIR
      ENDIF
      IERR = 13
      GOTO 9999
C
C     CAS SANS PROBLEME : INTERSECTION DIFFERENTE DE CELLE INITIALE
C     =================   =========================================
 60   NBTRCF = NBTRCF + 1
      NOTRCF( NBTRCF ) = NT1
C     PASSAGE AU TRIANGLE SUIVANT
      NA0 = NA2
      NT0 = NT1
      GOTO 30
C
C     ----------------------------------------------------------
C     ICI TOUTES LES INTERSECTIONS DE NS1-NS2 ONT ETE PARCOURUES
C     TOUS LES TRIANGLES INTERSECTES OU ETENDUS FORMENT LES
C     NBTRCF TRIANGLES DU TABLEAU NOTRCF
C     ----------------------------------------------------------
 70   NBTRCF = NBTRCF + 1
      NOTRCF( NBTRCF ) = NT1
C
C     FORMATION DU CF DES ARETES SIMPLES DES TRIANGLES DE NOTRCF
C     ET DESTRUCTION DES NBTRCF TRIANGLES DU TABLEAU NOARTR
C     ATTENTION: LE CHAINAGE LCHAIN DU TABLEAU NOSOAR DEVIENT ACTIF
C     =============================================================
 80   IF( NBTRCF*3 .GT. MXARCF ) THEN
         WRITE(IMPRIM,*) 'SATURATION DU TABLEAU NOARCF'
         IERR = 10
         GOTO 9999
      ENDIF
C
      CALL FOCFTR( NBTRCF, NOTRCF, NBARPI, PXYD,   NOARST,
     %             MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %             MOARTR, N1ARTR, NOARTR,
     %             NBARCF, N1ARCF, NOARCF, nbstpe, nostpe,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     CHAINAGE DES ARETES VIDES DANS LE TABLEAU NOARCF
C     ------------------------------------------------
C     DECALAGE DE 2 ARETES CAR 2 ARETES SONT NECESSAIRES ENSUITE POUR
C     INTEGRER 2 FOIS L'ARETE PERDUE ET FORMER AINSI 2 CF
C     COMME NBTRCF*3 MINORE MXARCF IL EXISTE AU MOINS 2 PLACES VIDES
C     DERRIERE => PAS DE TEST DE DEBORDEMENT
      N1ARCF(0) = NBARCF+3
      MMARCF = MIN(8*NBARCF,MXARCF)
      DO 90 I=NBARCF+3,MMARCF
         NOARCF(2,I) = I+1
 90   ENDDO
      NOARCF(2,MMARCF) = 0
C
C     REPERAGE DES SOMMETS NS1 NS2 DE L'ARETE PERDUE DANS LE CF
C     ---------------------------------------------------------
      NS1   = NOSOAR( 1, NARETE )
      NS2   = NOSOAR( 2, NARETE )
      NS(1) = NS1
      NS(2) = NS2
      DO 120 I=1,2
C        LA PREMIERE ARETE DANS NOARCF DU CF
         NA0 = N1ARCF(1)
 110     IF( NOARCF(1,NA0) .NE. NS(I) ) THEN
C           PASSAGE A L'ARETE SUIVANTE
            NA0 = NOARCF( 2, NA0 )
            GOTO 110
         ENDIF
C        POSITION DANS NOARCF DU SOMMET I DE L'ARETE PERDUE
         NACF(I) = NA0
 120  ENDDO
C
C     FORMATION DES 2 CF CHACUN CONTENANT L'ARETE NS1-NS2
C     ---------------------------------------------------
C     SAUVEGARDE DE L'ARETE SUIVANTE DE CELLE DE SOMMET NS1
      NA0 = NOARCF( 2, NACF1 )
      NT1 = NOARCF( 3, NACF1 )
C
C     LE PREMIER CF
      N1ARCF( 1 ) = NACF1
C     L'ARETE SUIVANTE DANS LE PREMIER CF
      NOARCF( 2, NACF1 ) = NACF2
C     CETTE ARETE EST CELLE PERDUE
      NOARCF( 3, NACF1 ) = NARETE
C
C     LE SECOND CF
C     L'ARETE DOUBLEE
      N1 = NBARCF + 1
      N2 = NBARCF + 2
C     LE PREMIER SOMMET DE LA PREMIERE ARETE DU SECOND CF
      NOARCF( 1, N1 ) = NS2
C     L'ARETE SUIVANTE DANS LE SECOND CF
      NOARCF( 2, N1 ) = N2
C     CETTE ARETE EST CELLE PERDUE
      NOARCF( 3, N1 ) = NARETE
C     LA SECONDE ARETE DU SECOND CF
      NOARCF( 1, N2 ) = NS1
      NOARCF( 2, N2 ) = NA0
      NOARCF( 3, N2 ) = NT1
      N1ARCF( 2 ) = N1
C
C     RECHERCHE DU PRECEDENT DE NACF2
 130  NA1 = NOARCF( 2, NA0 )
      IF( NA1 .NE. NACF2 ) THEN
C        PASSAGE A L'ARETE SUIVANTE
         NA0 = NA1
         GOTO 130
      ENDIF
C     NA0 PRECEDE NACF2 => IL PRECEDE N1
      NOARCF( 2, NA0 ) = N1
C
C     DEPART AVEC 2 CF
      NBCF = 2
C
C     TRIANGULATION DIRECTE DES 2 CONTOURS FERMES
C     L'ARETE NS1-NS2 DEVIENT UNE ARETE DE LA TRIANGULATION DES 2 CF
C     ==============================================================
      CALL TRIDCF( NBCF,   nbstpe, nostpe, PXYD,   NOARST,
     %             MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %             MOARTR, N1ARTR, NOARTR,
     %             MXARCF, N1ARCF, NOARCF, LARMIN,
     %             NBTRCF, NOTRCF, IERR )

 9999 TRATRI = TRATRI0
      RETURN
      END
