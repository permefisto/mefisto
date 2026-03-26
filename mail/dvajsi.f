      SUBROUTINE DVAJSI( NBITPI , COMIMX , AREMIN , ARETMX , ACCROI ,
     %                   MXSOMM , MXTRIA , PXYD   , NLSOFR , NOTRIA ,
     %                   DISTSO , NBSOMM , IERR )
C-------------------------------------------------------------------------
C fichier MEFISTO/doc/accroi
C Une TRIANGULATION EQUILATERALE forme le FOND du maillage
C Les ARETES FRONTALIERES peuvent etre de tailles beaucoup plus petites.
C Il est possible de passer RAPIDEMENT ou NON de la FRONTIERE
C a l'INTERIEUR en PEU ou BEAUCOUP de COUCHES de TRIANGLES
C
C Le reel ACCROI definit cette vitesse d'accroissement
C Si L est la longueur de l'arete frontaliere, alors
C L * ACCROI sera la taille de l'arete du triangle servant
C de premiere couche puis,
C L * ACCROI * ACCROI pour la suivante, et
C L * ( ACCROI ** n ) pour la n-eme couche.
C ACCROI=1.2 peut etre choisie comme VALEUR PAR DEFAUT.
C ACCROI<=1.0 impose en fait cette valeur 1.2 PAR DEFAUT.
C ACCROI=2  conduit a peu de couches pour atteindre les
C triangles equilateraux. La TAILLE DES COUCHES VARIE VITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES POINTS INTERNES
C -----    LA PREMIERE FOIS SUR UNE GRILLE REGULIERE
C          LES FOIS SUIVANTES COMME BARYCENTRES DES GRANDS TRIANGLES
C
C ENTREES:
C --------
C NBITPI : NOMBRE D'ITERATIONS SUR LA GENERATION DES POINTS INTERNES
C COMIMX : COORDONNEES MINIMALES ET MAXIMALES DES POINTS FRONTALIERS
C AREMIN : LONGUEUR MINIMALE D'UNE ARETE DE LA FRONTIERE
C ARETMX : LONGUEUR MAXIMALE DES ARETES DES TRIANGLES EQUILATERAUX
C ACCROI : FACTEUR D'ACCROISSEMENT DE L'ARETE FRONTALIERE COMME TAILLE
C          DES ARETES POUR ATTEINDRE LES TRIANGLES FRONTALIERS
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C
C ENTREES ET SORTIES :
C --------------------
C NLSOFR : NUMERO DE LA LIGNE FERMEE(1 A NBLFTR) DES POINTS
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE
C          0 SINON (POINT INTERNE NON IMPOSE PAR L'UTILISATEUR)
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE i
C
C DISTSO : DISTANCE SOUHAITEE TABLEAU AUXILIAIRE DE MXSOMM REELS
C NBSOMM : NOMBRE DE SOMMETS AVANT ET APRES AJOUT
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1991
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NLSOFR(MXSOMM)
      DOUBLE PRECISION  PXYD(3,MXSOMM), S, S1, S2, S3, SURTRI, SURTRA
      REAL              COMIMX(3,2),
     %                  DISTSO(MXSOMM)
C
C     LE NOMBRE INITIAL DE SOMMETS
      NBSOM0 = NBSOMM
C
      ACRLOG = 1.0 / LOG( ACCROI )
C     LE NOMBRE MAXIMAL DE COUCHES  POUR ATTEINDRE ARETMX A PARTIR
C     DE AREMIN ET DU FACTEUR D'ACCROISSEMENT ACCROI
      CMAX   = LOG( 1.0 + ( ARETMX*(ACCROI-1) / (AREMIN*ACCROI) ) )
      MXCOUC = NINT( CMAX * ACRLOG )
      MXCOUC = MAX( 1 , MXCOUC )
C
      IF( NBITPI .EQ. 0 ) THEN
C
C        ============================================================
C        LA PREMIERE FOIS GENERATION D'UNE GRILLE REGULIERE DE POINTS
C        ============================================================
         ARETEY = ARETMX * SQRT(3.0) * 0.5
C
C        LE NOMBRE DE COUCHES EN X ET Y
         NBCX  = NINT( ( COMIMX(1,2) - COMIMX(1,1) ) / ARETMX )
         NBCY  = NINT( ( COMIMX(2,2) - COMIMX(2,1) ) / ARETEY )
C
C        L'ORDONNEE INITIALE
C        LE RESTE/2 DU A NBCX OU NBCY ENTIER
         X = ( COMIMX(1,2) - COMIMX(1,1) - NBCX * ARETMX  ) * 0.5
         A =   COMIMX(1,1) + X
         Y = ( COMIMX(2,2) - COMIMX(2,1) - NBCY * ARETEY ) * 0.5
C        LA DERNIERE COUCHE TROP PRES DES BORDS EST SUPPRIMEE
         NBCX = NBCX - 1
         NBCY = NBCY - 1
         Y    =   COMIMX(2,1) + ARETEY + Y
         DO 20 J=1,NBCY
C           CADRAGE DU POINT POUR NE PAS ETRE TROP PRES DES BORDS
            I = MOD(J,2)
            IF( I .EQ. 0 ) THEN
               X  = A + ARETMX
               NB = NBCX
            ELSE
               X = A + ARETMX * 1.5
               NB = NBCX - 1
            ENDIF
            IF( Y .LT. COMIMX(2,2) ) THEN
               DO 10 I=1,NB
                  IF( X .LT. COMIMX(1,2) ) THEN
                     NBSOMM = NBSOMM + 1
                     IF( NBSOMM  .GT. MXSOMM ) THEN
                        NBLGRC(NRERR) = 2
                        KERR(1) = 'SATURATION DES SOMMETS'
                        WRITE(KERR(2)(1:10),'(I10)') NBSOMM
                        KERR(2)(11:NBCAER) = ' CREES '
                        CALL LEREUR
                        IERR = 1
                        GOTO 9999
                     ENDIF
                     PXYD(1,NBSOMM) = X
                     PXYD(2,NBSOMM) = Y
                     PXYD(3,NBSOMM) = ARETMX
                     NLSOFR(NBSOMM) = 0
CCCC                    LE TRACE DU POINT
CCC                     CALL MOVE2( X , Y )
CCC                     CALL ECRSYM
                  ENDIF
                  X = X + ARETMX
 10            CONTINUE
            ENDIF
            Y = Y + ARETEY
 20      CONTINUE
         NBITPI = NBITPI + 1
         IF( NBSOMM .GT. NBSOM0 ) RETURN
         NBITPI = 2
      ENDIF
C
      IF( NBITPI .LE. MXCOUC ) THEN
C
C        =========================================================
C        TRAITEMENT DES COUCHES EN MODIFIANT LES POIDS DES SOMMETS
C        DES TRIANGLES EQUILATERAUX EN CONTACT AVEC LES BORDS
C        =========================================================
C
C        PROTECTION POUR TRAITER UNE COUCHE SEULEMENT PAR ITERATION
         DO 90 I=1,NBSOMM
            DISTSO(I) = REAL( PXYD(3,I) )
 90      CONTINUE
C
         IMIN = 0
         NBSC = 0
         DO 200 NT1=1,MXTRIA
            IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 200
            CMIN = 1E38
            DO 110 I=1,3
               DETRT = DISTSO(NOTRIA(I,NT1))
               IF( CMIN .GT. DETRT ) THEN
                   CMIN = DETRT
                   IMIN = I
               ENDIF
 110        CONTINUE
C
            IF( CMIN .LT. ARETMX ) THEN
C              EXISTE T IL UN SOMMET DE DISTANCE SOUHAITEE ARETE ?
               NS2 = NOTRIA(IMIN,NT1)
               DO 120 I=1,3
                  NS1 = NOTRIA(I,NT1)
                  IF( DISTSO( NS1 ) .EQ. ARETMX ) THEN
C                     LE NOMBRE DE COUCHES NECESSAIRES POUR OBTENIR LA
C                     DISTANCE ENTRE CE SOMMET ET LE SOMMET MIN
                      D    = REAL( SQRT( (PXYD(1,NS1)-PXYD(1,NS2))**2 +
     %                                   (PXYD(2,NS1)-PXYD(2,NS2))**2 ))
                      H    = CMIN * SQRT( ACCROI*ACCROI - 0.25 )
                      CMAX = ( (ACCROI-1.0) * D ) / H
                      CMAX = LOG( 1 + CMAX ) * ACRLOG
C                     LE NOMBRE DE MAILLES
                      N    = MAX( NINT( CMAX ) , 1 ) - 1
                      CMAX = H * ( ACCROI ** N )
                      IF( CMAX .LT. ARETMX ) THEN
C                        CHANGEMENT DU POIDS DE CE SOMMET
                         PXYD(3,NS1) = CMAX
                         NBSC = NBSC + 1
                      ENDIF
                  ENDIF
 120           CONTINUE
            ENDIF
 200     CONTINUE
         IF( NBSC .LE. 0 ) NBITPI = MXCOUC + 1
C        EVITE DES TESTS INUTILES  ( BOUCLE 200 )
      ENDIF
C
C     ==================================================
C     DE NOUVEAUX SOMMETS INTERNES AUX TRIANGLES ACTUELS
C     DOIVENT ILS ETRE CREES COMME BARYCENTRES ?
C     ==================================================
C     BOUCLE SUR LES TRIANGLES
      DO 8100 NT1=1,MXTRIA
         IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 8100
C
C        UN TRIANGLE TROP ECRASE NE GENERE PAS DE BARYCENTRE
         CALL QUTR2D( PXYD(1,NOTRIA(1,NT1)),
     %                PXYD(1,NOTRIA(2,NT1)),
     %                PXYD(1,NOTRIA(3,NT1)), QUALIT )
         IF( QUALIT .LT. 0.1 ) GOTO 8100
C
C        SURFACE ACTUELLE DU TRIANGLE
         NS1 = NOTRIA(1,NT1)
         NS2 = NOTRIA(2,NT1)
         NS3 = NOTRIA(3,NT1)
         SURTRI =   ( PXYD(1,NS2)-PXYD(1,NS1) )
     %            * ( PXYD(2,NS3)-PXYD(2,NS1) )
     %            - ( PXYD(2,NS2)-PXYD(2,NS1) )
     %            * ( PXYD(1,NS3)-PXYD(1,NS1) )
C
C        MINIMUM DES SURFACES DES TRIANGLES ADJACENTS PAR UNE ARETE
         NBTA   = 0
         SURTRA = SURTRI
         DO 8010 I=1,3
C           LA SURFACE DU TRIANGLE ADJACENT PAR L'ARETE I DE NT1
            NT2   = NOTRIA(3+I,NT1)
            IF( NT2 .GT. 0 ) THEN
               NST1 = NOTRIA(1,NT2)
               NST2 = NOTRIA(2,NT2)
               NST3 = NOTRIA(3,NT2)
               S =   ( PXYD(1,NST2)-PXYD(1,NST1) )
     %             * ( PXYD(2,NST3)-PXYD(2,NST1) )
     %             - ( PXYD(2,NST2)-PXYD(2,NST1) )
     %             * ( PXYD(1,NST3)-PXYD(1,NST1) )
               SURTRA = MIN( SURTRA, S )
               NBTA = NBTA + 1
            ENDIF
 8010    CONTINUE
C
C        "SURFACE" SOUHAITEE DU TRIANGLE
         S1 = PXYD(3,NS1)
         S2 = PXYD(3,NS2)
         S3 = PXYD(3,NS3)
         S  = ( S1 + S2 + S3 ) * 0.5
         S  = S * (S-S1) * (S-S2) * (S-S3)
         IF( S .GE. 0 ) THEN
C           LE TRIANGLE SOUHAITE EST CONSTRUCTIBLE
            SURTRS = REAL( SQRT( S ) )
         ELSE
C           LE TRIANGLE SOUHAITE N'EST PAS CONSTRUCTIBLE
C           LE BARYCENTRE SANS POIDS EST AJOUTE
            SURTRS = 0
            S1     = 1D0 / 3D0
            S2     = S1
            S3     = 1D0 - S1
         ENDIF
C
         IF( SURTRI .GT. ACCROI*SURTRS ) THEN
CCC     %   .OR. (NBTA .GE. 2 .AND. SURTRI .GT. 3*SURTRA) ) THEN
C
C           UN NOUVEAU POINT NBSOMM DOIT ETRE AJOUTE DANS LE TRIANGLE NT1
C           -------------------------------------------------------------
            IF( NBSOMM .GE. MXSOMM ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES SOMMETS'
               CALL LEREUR
               IERR = 1
               GOTO 9999
            ENDIF
            NBSOMM = NBSOMM + 1
C
C           LE POINT EST INTERNE
            NLSOFR(NBSOMM) = 0
C
C           LES POIDS DU BARYCENTRE
            S1 = 1D0 / S1
            S2 = 1D0 / S2
            S3 = 1D0 / S3
C
C           L'INVERSE DE LA SOMME DES POIDS DU BARYCENTRE
            S  = 1D0 / ( S1 + S2 + S3 )
C
C           LES 2 COORDONNEES ET LA DISTANCE SOUHAITEE
            DO 8020 J=1,3
               PXYD(J,NBSOMM) = ( PXYD(J,NS1)*S1 + PXYD(J,NS2)*S2
     %                          + PXYD(J,NS3)*S3 ) * S
 8020       CONTINUE
CCCC
CCCC           TRACE DU TRIANGLE CONCERNE
CCC            CALL DVTRTR( PXYD , NOTRIA , NT1 , NCROUG , NCBLAN )
CCC            CALL DVTRTR( PXYD , NOTRIA , NOTRIA(4,NT1), NCVERT, NCBLAN )
CCC            CALL DVTRTR( PXYD , NOTRIA , NOTRIA(5,NT1), NCVERT, NCBLAN )
CCC            CALL DVTRTR( PXYD , NOTRIA , NOTRIA(6,NT1), NCVERT, NCBLAN )
CCCC           LE TRACE DU POINT
CCC            X = PXYD(1,NBSOMM)
CCC            Y = PXYD(2,NBSOMM)
CCC            CALL MOVE2(X, Y )
CCC            CALL ECRSYM
CCC            X = PXYD(1,NBSOMM)
CCC            Y = PXYD(2,NBSOMM)
CCC            CALL MOVE2( X, Y )
CCC            CALL ECRENT( NBSOMM )
         ENDIF
 8100 CONTINUE
C
 9999 NBITPI = NBITPI + 1
      END
