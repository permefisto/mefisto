      SUBROUTINE TRIACF( NBETOI, PXYD  , N1TRVI, NOTRIA, NOTRSO,
     %                   MXETRI, NARMIN, N1ARCF, NOARCF,
     %                   NBTRCF, NOTRCF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRIANGULATION FRONTALE POUR BOUCHER LES NBET ETOILES
C -----
C
C ENTREES:
C --------
C NBETOI : NOMBRE D'ETOILES OU CF A TRIANGULER
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : POINTE DANS NOTRIA VERS LE PREMIER TRIANGLE VIDE
C NOTRIA : LISTE DES TRIANGLES
C          SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C N1ARCF : NUMERO DE LA PREMIERE ARETE DE CHACUN DES NBETOI CF
C NOARCF : NUMERO DU SOMMET , NUMERO DE L'ARETE SUIVANTE
C          NUMERO DU TRIANGLE EXTERIEUR A L'ETOILE POUR CHAQUE ARETE
C
C AUXILIAIRES :
C -------------
C MXETRI : DIMENSION DU SECOND INDICE DES TABLEAUX AUXILIAIRES
C NARMIN : TABLEAU (MXETRI)   AUXILIAIRE
C
C SORTIE :
C --------
C NBTRCF : NOMBRE DE  TRIANGLES DES NBETOI CF
C NOTRCF : NUMERO DES TRIANGLES DES NBETOI ETOILES
C IERR   : 0 SI PAS D'ERREUR
C          3 SI CONTOUR FERME REDUIT A MOINS DE 3 ARETES
C          4 SATURATION DES TABLEAUX AUXILIAIRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1995
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,*),
     %                  NOTRSO(*),
     %                  N1ARCF(0:MXETRI),
     %                  NOARCF(3,MXETRI)
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NARMIN(MXETRI),
     %                  NOTRCF(MXETRI)
C
C     DEPART AVEC NBET CF
      NBET   = NBETOI
      NBTRCF = 0
C
C     TANT QUE LE NOMBRE DE CF CREES EST NON NUL FAIRE
C     TRIANGULATION FRONTALE DU CF
C     ================================================
 10   IF( NBET .GT. 0 ) THEN
C
C        L'ETOILE EN HAUT DE PILE A POUR PREMIERE ARETE
         NA01 = N1ARCF( NBET )
         NA1  = NOARCF( 2, NA01 )
C
C        CHOIX DU SOMMET DE L'ETOILE A RELIER A L'ARETE NA1
         CALL TRCHTD( PXYD, NA01, NA1, NOARCF,
     %                NA03, NA3, NARMIN )
         IF( NA3 .EQ. 0 ) THEN
            IERR = 3
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
     %                MXETRI, N1ARCF, NOARCF, NF )
         IF( NF .LE. 0 ) THEN
C           SATURATION DU TABLEAU NOTRIA OU NOARCF OU N1ARCF
            IERR = 4
            RETURN
         ENDIF
C
C        AJOUT DU TRIANGLE CREE A SA LISTE
         IF( NBTRCF .GE. MXETRI ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) ='SATURATION DES TRIANGLES DU CF'
            CALL LERESU
            IERR = 4
            RETURN
         ENDIF
         NBTRCF = NBTRCF + 1
         NOTRCF( NBTRCF ) = NF
         GOTO 10
      ENDIF
C
C     MISE A JOUR DU CHAINAGE DES FACES ADJACENTES PAR LES ARETES
C     -----------------------------------------------------------
      DO 700 NFP0 = 1, NBTRCF
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
            DO 530 NFP1=NFP0+1,NBTRCF
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
C                 ARETE NON RETROUVEE
 510           CONTINUE
C              PASSAGE A LA FACE SUIVANTE
 530        CONTINUE
            GOTO 600
C
C           ARETE RETROUVEE ( NF0 , I ) <-> ( NF , J )
 550        NOTRIA(3+I,NF0) = NF
            NOTRIA(3+J,NF ) = NF0
 600     CONTINUE
 700  CONTINUE
      END
