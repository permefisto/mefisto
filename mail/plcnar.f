      SUBROUTINE PLCNAR( NOFOTI, MXGENE, XYZ1PL, VN,
     %                   PAX1CN, RAY1CN, RAY2CN, HAUTCN, RPCONE,
     %                   NPSUIV, RPINT,  NBIP,   N1COUR,
     %                   NBARLI, NBSOM,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES LIGNES D'INTERSECTION ENTRE UN TRONC DE CONE ET
C -----    UN PLAN DEFINI PAR 3 POINTS
C
C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION TAILLE_IDEALE SI ELLE EXISTE, 0 SINON
C MXGENE : NOMBRE DE GENERATRICES DU CONE A INTERSECTER AVEC LE PLAN
C XYZ1PL : 3 COORDONNEES DANS LE REPERE DU CONE DES 1-ER POINT DE DEFINITION DU
C VN     : VECTEUR NORMAL UNITAIRE AU PLAN DES 3 POINTS
C
C PAX1CN : POINT 1 DE BASE DE L'AXE DU TRONC DE CONE
C RAY1CN : RAYON DU CERCLE AU POINT PAX1CN (PLAN ORTHOGONAL A L'AXE DU CONE )
C PAX2CN : POINT 2 DE L'AXE DU TRONC DE CONE
C RAY2CN : RAYON DU CERCLE AU POINT PAX2CN (PLAN ORTHOGONAL A L'AXE DU CONE )
C
C HAUTCN : HAUTEUR DU TRONC DE CONE
C RPCONE : REPERE PROPRE ORTHONORME DU CONE D'ORIGINE AU POINT PAX1CN
C          ORIGINE AU POINT PAX1CN
C
C NPSUIV : NUMERO DANS RPINT DU POINT D'INTERSECTION SUIVANT
C RPINT  : RPINT(1,N) ANGLE DE LA GENERATRICE DU CONE DE CE POINT N
C          EN SORTIE=>RPINT(2,N) DISTANCE AVEC LE POINT QUI PRECEDE
C NBIP   : NUMERO DANS N1COUR DU DERNIER INTERVALLE DE POINTS D'INTERSECTION
C
C MODIFIES :
C ----------
C N1COUR : DEBUT DES CHAINAGES DES PARTIES DE LA LIGNE
C          N1COUR(0) POINTE SUR LE PREMIER VIDE
C NBARLI : NOMBRE SOUHAITE D'ARETES DE LA LIGNE
C
C SORTIES:
C --------
C NBARLI : NOMBRE FINAL D'ARETES   DE LA LIGNE INTERSECTION
C NBSOM  : NOMBRE FINAL DE SOMMETS DE LA LIGNE INTERSECTION
C IERR   : =0 PAS D'ERREUR DETECTEE
C          >0 PAS DE CREATION DE LA LIGNE INTERSECTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1998
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  PAX1CN(3)
      DOUBLE PRECISION  RAY1CN, RAY2CN, HAUTCN, PARAME
      DOUBLE PRECISION  RPCONE(3,3), XYZ1PL(3), VN(3)
      DOUBLE PRECISION  RPINT(1:2,MXGENE+1)
      INTEGER           NPSUIV(MXGENE+1)
      INTEGER           N1COUR(0:6)
C
      DOUBLE PRECISION  P(3), Q(3), PT0(3), PT1(3),
     %                  PT(3), RL, RL0, RL2,
     %                  ANGLE, ANGLE0, ANGLE1, ANGLE2, DEUXPI,
     %                  TAILAR, TAILA1, TAILA2
C
C     ICI, IL EXISTE (NBIP-1)/2 INTERVALLES DE POINTS D'INTERSECTION
C     GENERATION DES EF SEGMENTS P3-HERMITE
C     NPSUIV(POINT)<0 INDIQUE UN POINT INTERMEDIAIRE NON SOMMET D'UN SEGMENT P3
C     =========================================================================
      DEUXPI = ATAN(1D0) * 8D0
      NBSOM  = 0
C
      IF( NOFOTI .LE. 0 ) THEN
C
C        PAS DE FONCTION 'TAILLE_IDEALE'
C        -------------------------------
C
C        CALCUL DE LA LONGUEUR DE LA LIGNE ACTUELLE
         RLONG = 0D0
         DO 5 N = 1, NBIP, 2
            N1 = 0
            N2 = N1COUR( N )
C           LA FIN DU CHAINAGE DE L'INTERVALLE DE POINTS D'INTERSECTION
            NF = N1COUR( N+1 )
C
C           CALCUL DU PREMIER POINT
            CALL PLCNAN( RPINT(1,N2), RAY1CN, RAY2CN, HAUTCN,
     %                   XYZ1PL, VN,
     %                   PARAME, P, Q,  IERR )
            PT0(1) = P(1) + PARAME * ( Q(1) - P(1) )
            PT0(2) = P(2) + PARAME * ( Q(2) - P(2) )
            PT0(3) = P(3) + PARAME * ( Q(3) - P(3) )
C           LA DISTANCE AVEC LE POINT QUI PRECEDE
            RPINT(2,N2) = 0D0
            GOTO 4
C
C           LES POINTS SUIVANTS
 3          N1 = N2
 4          N2 = NPSUIV(N2)
            IF( N2 .GT. 0 .AND. N1 .NE. NF ) THEN
C              CALCUL DU POINT
               CALL PLCNAN( RPINT(1,N2), RAY1CN, RAY2CN, HAUTCN,
     %                      XYZ1PL, VN,
     %                      PARAME, P, Q,  IERR )
               PT1(1) = P(1) + PARAME * ( Q(1) - P(1) )
               PT1(2) = P(2) + PARAME * ( Q(2) - P(2) )
               PT1(3) = P(3) + PARAME * ( Q(3) - P(3) )
C
C              LA DISTANCE AVEC LE POINT QUI PRECEDE
               RPINT(2,N2) = SQRT( ( PT1(1) - PT0(1) ) ** 2 +
     %                             ( PT1(2) - PT0(2) ) ** 2 +
     %                             ( PT1(3) - PT0(3) ) ** 2 )
C
C              LA LONGUEUR DE LA LIGNE
               RLONG = REAL( RLONG + RPINT(2,N2) )
               PT0(1) = PT1(1)
               PT0(2) = PT1(2)
               PT0(3) = PT1(3)
               GOTO 3
            ENDIF
 5       CONTINUE
C
C        LONGUEUR MOYENNE SOUHAITEE D'UNE ARETE
         RL = RLONG / NBARLI
C
C        LE DECOUPAGE EN ARETES DE TAILLE VOISINE DE RL
C        ==============================================
         NBARLI = 0
         DO 50 N = 1, NBIP, 2
C
C           LE DEBUT DU CHAINAGE DE L'INTERVALLE DE POINTS D'INTERSECTION
            NBA = 0
            N1  = -1
            N2  = N1COUR( N )
C           LA FIN   DU CHAINAGE DE L'INTERVALLE DE POINTS D'INTERSECTION
            NF  = N1COUR( N+1 )
C
C           LE PREMIER SOMMET DE CETTE PARTIE DE LA LIGNE
            NPSUIV( N2 ) = ABS( NPSUIV(N2) )
C
C           CALCUL DE L'ARETE PAR PARCOURS DES POINTS D'INTERSECTION
            N2DER = 0
C           NOMBRE DE SOMMETS ET ARETES DE CET INTERVALLE DE POINTS
            NBAR = 0
C
C           DEBUT D'ARC OU SEGMENT P3 DE L'INTERVALLE DE POINTS
 20         NBPARC = 0
            TAILAR = 0D0
C
 30         IF( (N2 .GT. 0 .AND. N2 .NE. NF) .OR. NBA .LE. 0 ) THEN
               NBA = NBA + 1
               N1  = N2
               N2  = ABS( NPSUIV(N2) )
C              POINT N2 A CONSIDERER
               TAILAR = TAILAR + RPINT(2,N2)
               IF( TAILAR .LT. RL ) THEN
C                 POINT A AJOUTER
                  NBPARC = NBPARC + 1
C                 SON SUIVANT EST MARQUE <0 COMME POINT INTERMEDIAIRE
C                 NPSUIV(N2)>0 INDIQUE UN SOMMET DES ARETES
                  NPSUIV( N2 ) = - ABS( NPSUIV( N2 ) )
                  GOTO 30
               ELSE
C                 DERNIER POINT DES ARETES DE CETTE PARTIE DE LA LIGNE
                  IF( NBPARC .GE. 4 ) THEN
C                    LE POINT N2 N'EST PAS PRIS. RETOUR AU PRECEDENT
                     N2 = N1
                  ENDIF
C                 LE POINT N2 EST PRIS <=> NPSUIV(N2) RESTE POSITIF
                  NPSUIV( N2 ) = ABS( NPSUIV( N2 ) )
                  NBAR = NBAR + 1
C                 SAUVEGARDE DE L'AVANT DERNIER SOMMET DE L'INTERVALLE
                  N2DER = N2
                  GOTO 20
               ENDIF
            ENDIF
C
C           FIN DE L'INTERVALLE DES POINTS D'INTERSECTION
C           LE DERNIER ARC EST IL TROP PETIT?
            IF( TAILAR .LE. 0D0 ) THEN
C              PAS D'ARC
               GOTO 45
            ENDIF
            IF( TAILAR .LT. RL*0.5D0 ) THEN
C              OUI
               IF( N2DER .GT. 0 ) THEN
C                 LE DERNIER ARC EST AJOUTE AU PRECEDENT
C                 CE QUI CORRESPOND A ELIMINER LE POINT N2
                  NPSUIV( N2DER ) = - ABS( NPSUIV( N2DER ) )
                  GOTO 45
               ENDIF
            ENDIF
C
C           NOMBRE TOTAL D'ARETES ET DE SOMMETS DE LA LIGNE FUTURE
            NBAR = NBAR + 1
 45         NBARLI = NBARLI + NBAR
            NBSOM  = NBSOM  + NBAR + 1
C           LE DERNIER SOMMET DE CETTE PARTIE DE LA LIGNE
            NPSUIV( NF ) = ABS( NPSUIV( NF ) )
C
 50      CONTINUE
C
      ELSE
C
C        EXISTENCE DE LA FONCTION 'TAILLE_IDEALE'
C        ----------------------------------------
C        LE DECOUPAGE EN SEGMENTS P3 DE TAILLE VOISINE A TAILLE_IDEALE
         NBARLI = 0
         DO 170 N = 1, NBIP, 2
C
C           LE DEBUT DU CHAINAGE DE L'INTERVALLE DE POINTS D'INTERSECTION
            N1 = -1
            N2 = N1COUR( N )
C           LA FIN   DU CHAINAGE DE L'INTERVALLE DE POINTS D'INTERSECTION
            NF = N1COUR( N+1 )
C           NOMBRE DE SOMMETS ET ARETES DE CET INTERVALLE DE POINTS
            NBAR = 0
C           NUMERO DU POINT DE L'ARETE PRECEDENTE ET AVANT PRECEDENTE
            N2DER  =  0
            N2DER0 = -1
C
C           LE PREMIER SOMMET PT0 DE CETTE PARTIE DE LA LIGNE
            NPSUIV( N2 ) = ABS( NPSUIV(N2) )
C           LA TAILLE_IDEALE AUTOUR DU SOMMET N2 INITIAL
            ANGLE0 = RPINT(1,N2)
            CALL PLCNAN( ANGLE0, RAY1CN, RAY2CN, HAUTCN, XYZ1PL, VN,
     %                   PARAME, P, Q,  IERR )
            PT(1) = P(1) + PARAME * ( Q(1) - P(1) )
            PT(2) = P(2) + PARAME * ( Q(2) - P(2) )
            PT(3) = P(3) + PARAME * ( Q(3) - P(3) )
C           RETOUR DANS LE REPERE GLOBAL INITIAL
            CALL AC3D3D( PAX1CN, RPCONE, PT, PT0 )
            CALL FONVAL( NOFOTI, 3, PT0,  NCODEV, RL0 )
            IF( NCODEV .EQ. 0 ) THEN
C               NCODEV  : 0 RL N'EST PAS INITIALISEE EN SORTIE
C                         1 RL   EST     INITIALISEE EN SORTIE
               PT1(1) = PT0(1)
               PT1(2) = PT0(2)
               PT1(3) = PT0(3)
               GOTO 9000
            ENDIF
            RL = ABS( RL0 )
C
C           ESSAI DE CREER UN ARC DE LONGUEUR CETTE TAILLE RL
C           -------------------------------------------------
C           LE POINT INITIAL DE L'ARC EST N0
 120        N0 = N2
C           LA TAILLE DE L'ARC EST NULLE
            TAILAR = 0D0
C
 130        IF( N2 .NE. NF .OR. N2DER .LE. 0 ) THEN
C
C              CE N'EST NI LE PREMIER NI LE DERNIER ARC
               N2 = ABS( NPSUIV(N2) )
C              POINT A CONSIDERER
C              CALCUL DE LA TAILLE_IDEALE AUTOUR DU SOMMET N2
               ANGLE2 = RPINT(1,N2)
               CALL PLCNAN( ANGLE2, RAY1CN, RAY2CN, HAUTCN, XYZ1PL, VN,
     %                      PARAME, P, Q,  IERR )
               PT(1) = P(1) + PARAME * ( Q(1) - P(1) )
               PT(2) = P(2) + PARAME * ( Q(2) - P(2) )
               PT(3) = P(3) + PARAME * ( Q(3) - P(3) )
C              RETOUR DANS LE REPERE GLOBAL INITIAL
               CALL AC3D3D( PAX1CN, RPCONE, PT, PT1 )
C
C              DISTANCE ENTRE PT0 ET PT1
               TAILAR = SQRT( (PT1(1)-PT0(1)) ** 2
     %                      + (PT1(2)-PT0(2)) ** 2
     %                      + (PT1(3)-PT0(3)) ** 2 )
C
C              LA TAILLE_IDEALE DE L'ARETE EN CE POINT N2
               CALL FONVAL( NOFOTI, 3, PT1,  NCODEV, RL2 )
               IF( NCODEV .EQ. 0 ) GOTO 9000
C
               IF( TAILAR .LT. RL ) THEN
C
C                 ARETE TROP COURTE: POINT N2 A ELIMINER DANS L'ARC P3
C                                    SAUF SI C'EST LE DERNIER DE LA LIGNE
                  IF( N2 .EQ. NF ) THEN
                     GOTO 150
                  ENDIF
C
C                 ELIMINATION DU POINT DANS LE CHAINAGE
                  NPSUIV( N0 ) = NPSUIV( N2 )
C                 N2 EST RENDU LIBRE
                  NPSUIV( N2 ) = N1COUR(0)
                  N1COUR(0) = N2
C                 LE SUIVANT
                  N2 = N0
                  GOTO 130
C
               ENDIF
C
C              TAILAR>=RL => LE SOMMET IDEAL EST ENTRE N0 ET N2
C              LA TAILLE_IDEALE MINIMALE ENTRE LES 2 POINTS
               RL = MIN( RL0, RL2 )
C              L'ANGLE INITIAL ET FINAL
               ANGLE1 = ANGLE0
               IF( N2 .EQ. NF ) THEN
                  IF( ABS(ANGLE2) .LE. 1D-12 ) THEN
C                    PERIODICTE DE 2 PI A PRENDRE EN COMPTE
                     ANGLE2 = ANGLE2 + DEUXPI
                  ENDIF
               ENDIF
               ANGLE = ANGLE2
C
C              DISTANCE INITIALE ENTRE LE POINT N0 ET N2
               TAILA2 = TAILAR
C              DISTANCE INITIALE ENTRE LE POINT N0 ET N1=N0
               TAILA1 = 0D0
C
C              DICHOTOMIE POUR OBTENIR UNE DISTANCE EGALE A RL
C              -----------------------------------------------
 140           IF( ABS(ANGLE2-ANGLE1) .GT. 1D-12 ) THEN
C
                  IF( ABS(TAILAR-RL) .GT. 0.01D0*TAILAR ) THEN
C
                     ANGLE = ( ANGLE1 + ANGLE2 ) / 2
                     CALL PLCNAN( ANGLE,  RAY1CN, RAY2CN, HAUTCN,
     %                            XYZ1PL, VN,
     %                            PARAME, P, Q,  IERR )
                     PT(1) = P(1) + PARAME * ( Q(1) - P(1) )
                     PT(2) = P(2) + PARAME * ( Q(2) - P(2) )
                     PT(3) = P(3) + PARAME * ( Q(3) - P(3) )
C                    RETOUR DANS LE REPERE GLOBAL INITIAL
                     CALL AC3D3D( PAX1CN, RPCONE, PT, PT1 )
C
C                    DISTANCE ENTRE PT0 ET PT1
                     TAILAR = SQRT( (PT1(1)-PT0(1)) ** 2
     %                            + (PT1(2)-PT0(2)) ** 2
     %                            + (PT1(3)-PT0(3)) ** 2 )
C
                     IF( TAILAR .GT. RL ) THEN
C
C                       LE POINT MILIEU DEVIENT LE POINT LE PLUS ELOIGNE DE PT0
                        TAILA2 = TAILAR
                        ANGLE2 = ANGLE
C
                     ELSE
C
C                       LE POINT MILIEU DEVIENT LE POINT LE MOINS ELOIGNE DE PT0
                        TAILA1 = TAILAR
                        ANGLE1 = ANGLE
C
                     ENDIF
                     GOTO 140
C
                  ENDIF
C
               ELSE
C
C                 CONVERGENCE DE L'ANGLE MAIS PT0 PT1 SONT 2 POINTS
C                 DE MEME ANGLE ET DE RACINES DIFFERENTES ET TAILAR<RL
                  CALL XVPAUSE
                  GOTO 147
C
               ENDIF
C
C              ANGLE ACCEPTE => UNE ARETE OU ARC OU SEGMENT P3 EN PLUS
C              LE NOUVEAU POINT EST INTERCALE ENTRE N0 ET N2
C              POUR EVITER DE PERDRE UNE PARTIE DE LA COURBE
 147           N1 = N1COUR(0)
               IF( N1 .LE. 0 ) GOTO 9100
               N1COUR(0) = NPSUIV(N1)
C
C              LE POINT N1 EST INTERCALE ENTRE N0 ET N2
               NPSUIV( N1 ) = N2
C              LE SUIVANT DE N0 EST DIRECTEMENT N1
               NPSUIV( N0 ) = N1
C              L'ANGLE DU POINT
               RPINT(1,N1) = ANGLE
C              LA DISTANCE DE N1 A SON PREDECESSEUR N0
               RPINT(2,N1) = TAILAR
C              LA DISTANCE DE N2 A SON PREDECESSEUR N1
               RPINT(2,N2) = RPINT(2,N2) - TAILAR
C              LA TAILLE_IDEALE DE L'ARETE EN CE POINT
               CALL FONVAL( NOFOTI, 3, PT1, NCODEV, RL0 )
               IF( NCODEV .EQ. 0 ) GOTO 9000
               RL = ABS( RL0 )
C              L'ANGLE DU SOMMET INITIAL DE LA FUTURE ARETE
               ANGLE0 = ANGLE
C              LE SOMMET INITIAL DE LA FUTURE ARETE
               PT0(1) = PT1(1)
               PT0(2) = PT1(2)
               PT0(3) = PT1(3)
C              UNE ARETE DE PLUS
               NBAR  = NBAR + 1
C              SAUVEGARDE DE L'AVANT DERNIER SOMMET DE L'INTERVALLE
               N2DER0 = N2DER
               N2DER  = N1
               N2     = N1
               GOTO 120
            ENDIF
C
C           FIN DE L'INTERVALLE DES POINTS D'INTERSECTION
C           ---------------------------------------------
C           LE DERNIER ARC EST IL TROP PETIT?
 150        IF( TAILAR .LE. 0D0 ) THEN
C              PAS D'ARC
               GOTO 160
            ENDIF
            IF( TAILAR .LT. RL*0.5D0 ) THEN
C              OUI
               IF( N2DER0 .GT. 0 ) THEN
C                 LE DERNIER ARC EST AJOUTE AU PRECEDENT
C                 ARETE TROP COURTE: POINT N2 A ELIMINER DANS L'ARC P3
                  NPSUIV(N2DER0) = NF
                  GOTO 160
               ENDIF
            ENDIF
C           UNE ARETE DE PLUS
            NBAR = NBAR + 1
C
C           NOMBRE TOTAL D'ARETES ET DE SOMMETS DE LA FUTURE LIGNE INTERSECTION
 160        NBARLI = NBARLI + NBAR
            NBSOM  = NBSOM  + NBAR + 1
C
 170     CONTINUE
      ENDIF
      GOTO 9999
C
C     ERREUR
 9000 NBLGRC(NRERR) = 2
      KERR(1) = 'FONCTION TAILLE_IDEALE INCALCULABLE AU POINT'
      KERR(2) = ' '
      WRITE(KERR(2)(1:15), '(E15.7)') PT1(1)
      WRITE(KERR(2)(18:32),'(E15.7)') PT1(2)
      WRITE(KERR(2)(35:49),'(E15.7)') PT1(3)
      CALL LEREUR
      IERR = 21
      GOTO 9999
C
C     ERREUR
 9100 NBLGRC(NRERR) = 2
      KERR(1) = 'TABLEAU N1COUR SATURE'
      KERR(2) = 'AUGMENTER NOMBRE ARETES DE LA LIGNE'
      CALL LEREUR
      IERR = 22
C
 9999 RETURN
      END
