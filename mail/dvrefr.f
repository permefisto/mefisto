      SUBROUTINE DVREFR( NBLFTR, NDARLF, MXSOAR, NOSOAR,
     %                   PXYD  , NLSOFR,
     %                   N1TRVI, NOTRIA, NOTRSO,
     %                   MXETRI, NAETOI, NARMIN,
     %                   MXARCF, N1ARCF, NOARCF, NOTRCF,
     %                   NBARPE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RECHERCHE DES ARETES DE LA FRONTIERE NON DANS LA TRIANGULATION
C -----   TRIANGULATION FRONTALE POUR LES REOBTENIR
C
C ENTREES:
C --------
C NBLFTR : NOMBRE DE LIGNES FERMEES LIMITANT LA SURFACE A TRIANGULER
C NDARLF : NUMERO DE LA PREMIERE ARETE DE CHAQUE LIGNE FERMEE DANS
C          LE TABLEAU NOSOAR
C MXSOAR : NOMBRE MAXIMAL D'ARETES FRONTIERES DECLARABLES
C NOSOAR : NUMERO DES 2 SOMMETS DE CHAQUE ARETE ET
C          POINTEUR SUR L'ARETE SUIVANTE
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C NLSOFR : NUMERO(1 A NBLFTR) DE LA LIGNE FERMEE DU POINT
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE
C          0 SI LE POINT EST INTERNE OU EXTERNE NON IMPOSE
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOTRIA
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOTRIA(4,.)
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
C NAETOI : TABLEAU (4,MXETRI) AUXILIAIRE D'ENTIERS
C NARMIN : TABLEAU (MXETRI)   AUXILIAIRE D'ENTIERS
C N1ARCF : TABLEAU (0:MXARCF) AUXILIAIRE D'ENTIERS
C NOARCF : TABLEAU (3,MXARCF) AUXILIAIRE D'ENTIERS
C NOTRCF : TABLEAU (  MXARCF) AUXILIAIRE D'ENTIERS
C
C SORTIE :
C --------
C NBARPE : NOMBRE D'ARETES PERDUES PUIS RETROUVEES
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          2 UN SOMMET D'ARETE PERDUE APPARTIENT A AUCUN TRIANGLE
C          9 TABLEAU NOSOAR DE TAILLE INSUFFISANTE CAR TROP D'ARETES
C            A PROBLEME
C          10 UN DES TABLEAUX N1ARCF, NOARCF NOTRCF EST SATURE
C             AUGMENTER A L'APPEL MXARCF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C....................................................................012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOSOAR(3,MXSOAR),
     %                  NDARLF(1:NBLFTR),
     %                  NOTRIA(6,*),
     %                  NLSOFR(*),
     %                  NOTRSO(*),
     %                  N1ARCF(0:MXARCF),
     %                  NOARCF(3,MXARCF)
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NAETOI(4,MXETRI),
     %                  NARMIN(MXETRI),
     %                  NOTRCF(MXARCF)
C
C     LA BOUCLE SUR LES ARETES DES LIGNES FERMEES
C     ===========================================
      NBARPE = 0
      DO 1000 LLL=1,NBLFTR
C
C        LA PREMIERE ARETE DE LA LIGNE FERMEE
         NARET1 = NDARLF( LLL )
C
C        TANT QU'IL EXISTE DES ARETES DANS CETTE LIGNE FERMEE LLL
 1       IF( NARET1 .GT. 0 ) THEN
C
C           LE NUMERO DE L'ARETE DANS NOARET EST RECHERCHE
            NS1 = NOSOAR(1,NARET1)
            NS2 = NOSOAR(2,NARET1)
C
C           PARCOURS DES TRIANGLES PAR LES ARETES DE SOMMET NS1
            NT0  = NOTRSO( NS1 )
            NT1  = NT0
C
C           REPERAGE DES SOMMETS NS2 ET NS1 DANS NT1
 5          NA1 = 0
            DO 10 I=1,3
               IF( NOTRIA(I,NT1) .EQ. NS2 ) GOTO 50
               IF( NOTRIA(I,NT1) .EQ. NS1 ) NA1 = I
 10         CONTINUE
 12         IF( NA1 .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'DVREFR:SOMMET ',NS1,
     %         ' NON DANS TRIANGLE ',NT1,' ST:',(NOTRIA(I,NT1),I=1,3)
               IERR = 2
               RETURN
            ENDIF
C
C           NS2 NON RETROUVE. L'ARETE DE GAUCHE DE NA1
            NT1 = NOTRIA(NA1+3,NT1)
            IF( NT1 .EQ. 0 ) THEN
C
C              LE PARCOURS PASSE PAR 1 DES TRIANGLES FRONTALIERS
C              LE PARCOURS EST INVERSE PAR L'ARETE DE GAUCHE
               NT1 = NT0
C              REPERAGE DES SOMMETS NS2 ET NS1 DANS NT1
 15            NA1 = 0
               DO 20 I=1,3
                  IF( NOTRIA(I,NT1) .EQ. NS2 ) GOTO 50
                  IF( NOTRIA(I,NT1) .EQ. NS1 ) NA1 = I
 20            CONTINUE
               IF( NA1 .EQ. 0 ) GOTO 12
C
C              NS2 NON RETROUVE. L'ARETE DE DROITE DE NA1
               NA1 = NA1 + 2
               IF( NA1 .EQ. 3 ) NA1 = 6
               NT1 = NOTRIA(NA1,NT1)
               IF( NT1 .EQ. 0   ) GOTO 30
               IF( NT1 .NE. NT0 ) GOTO 15
            ENDIF
            IF( NT1 .NE. NT0 ) GOTO 5
C
C           L'ARETE N : NS1-NS2 N'EST PAS RETROUVEE
C           ---------------------------------------
 30         NBARPE = NBARPE + 1
            WRITE(IMPRIM,10030) NS1,(PXYD(J,NS1),J=1,2),
     %                          NS2,(PXYD(J,NS2),J=1,2)
10030       FORMAT(' DVREFR: ARETE PERDUE',
     %             (T23,'SOMMET=',I6,' X=',G13.5,' Y=',G13.5))
C
C           TRAITEMENT DE CETTE ARETE PERDUE NS1-NS2
            CALL DVFOAR( NS1,    NS2,
     %                   PXYD  , NLSOFR,
     %                   N1TRVI, NOTRIA, NOTRSO,
     %                   MXETRI, NAETOI, NARMIN,
     %                   MXARCF, N1ARCF, NOARCF, NOTRCF,
     %                   IERR )
            IF( IERR .NE. 0 ) RETURN
C
C           FIN DU TRAITEMENT DE CETTE ARETE PERDUE ET RETROUVEE
C           ----------------------------------------------------
C
C           PASSAGE A L'ARETE SUIVANTE DE LA FRONTIERE
 50         NARET1 = NOSOAR(3,NARET1)
            GOTO 1
         ENDIF
 1000 CONTINUE
      END
