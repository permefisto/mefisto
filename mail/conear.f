      SUBROUTINE CONEAR( NOFOTI, MXPINT,
     %                   PAX1C1, RAY1C1, RAY2C1,
     %                   PAX1C2, RAY1C2, RAY2C2,
     %                   HAUTC1, RPCON1, HAUTC2, RPCON2,
     %                   NPSUIV, NTPINT, RPINT,
     %                   NBIP1,  NBIP2,  N1COUR, RLONG,
     %                   NBARLI, NBSOM,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES LIGNES D'INTERSECTION ENTRE 2 TRONCS DE CONES
C -----
C
C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION TAILLE_IDEALE SI ELLE EXISTE, 0 SINON
C MXGENE : NOMBRE DE GENERATRICES DU CONE 1 A INTERSECTER AVEC LE CONE 2
C MXPINT : NOMBRE MAXIMAL DE POINTS D'INTERSECTION DECLARABLES (>2*MXGENE)
C
C PAX1C1 : POINT DE BASE DE L'AXE DU TRONC DE CONE 1
C RAY1C1 : RAYON DU CERCLE AU POINT PAX1C1 (PLAN ORTHOGONAL A L'AXE DU CONE )
C PAX2C1 : POINT DE L'AXE DU TRONC DE CONE 1
C RAY2C1 : RAYON DU CERCLE AU POINT PAX2C1 (PLAN ORTHOGONAL A L'AXE DU CONE )
C
C PAX1C2 : POINT DE BASE DE L'AXE DU TRONC DE CONE 2
C RAY1C2 : RAYON DU CERCLE AU POINT PAX1C2 (PLAN ORTHOGONAL A L'AXE DU CONE )
C PAX2C2 : POINT DE L'AXE DU TRONC DE CONE 2
C RAY2C2 : RAYON DU CERCLE AU POINT PAX2C2 (PLAN ORTHOGONAL A L'AXE DU CONE )
C
C HAUTC1 : HAUTEUR DU TRONC DE CONE 1
C RPCON1 : REPERE PROPRE ORTHONORME DU CONE 1 D'ORIGINE AU POINT PAX1C1
C          ORIGINE AU POINT PAX1C1
C HAUTC2 : HAUTEUR DU TRONC DE CONE 2
C RPCON2 : REPERE PROPRE ORTHONORME DU CONE 2 D'ORIGINE AU POINT PAX1C2
C          ORIGINE AU POINT PAX1C2
C NPSUIV : NUMERO DANS NTPINT, RPINT DU POINT D'INTERSECTION SUIVANT
C NTPINT : NUMERO DE LA RACINE PARMI LES 2 DISTINCTES DU POINT D'INTERSECTION
C          3 SI RACINE DOUBLE
C RPINT  : RPINT(1,N) ANGLE DE LA GENERATRICE DU CONE 1 DE CE POINT N
C          RPINT(2,N) RACINE TEL QUE PM = RACINE PQ SUR LA GENERATRICE DU CONE 1
C NBIP1  : NUMERO DANS N1COUR DU PREMIER INTERVALLE DE POINTS D'INTERSECTION
C NBIP2  : NUMERO DANS N1COUR DU DERNIER INTERVALLE DE POINTS D'INTERSECTION
C N1COUR : DEBUT DES CHAINAGES DES COURBES, PUIS DES INTERVALLES DE GENERATRICES
C          PUIS DES INTERVALLES DE POINTS D'INTERSECTION
C          N1COUR(0) POINTE SUR LE PREMIER VIDE
C RLONG  : LONGUEUR DE LA COURBE POUR CHAQUE INTERVALLE N1COUR(NBIP1,NBIP2, PAS
C NBARLI : NOMBRE SOUHAITE D'ARETES DE LA LIGNE
C
C SORTIES:
C --------
C NBARLI : NOMBRE FINAL D'ARETES   DE LA LIGNE INTERSECTION DES 2 TRONCS DE CONE
C NBSOM  : NOMBRE FINAL DE SOMMETS DE LA LIGNE INTERSECTION DES 2 TRONCS DE CONE
C IERR   : =0 PAS D'ERREUR DETECTEE
C          >0 PAS DE CREATION DE LA LIGNE INTERSECTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1998
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  PAX1C1(3)
      DOUBLE PRECISION  PAX1C2(3)
      DOUBLE PRECISION  RAY1C1, RAY2C1
      DOUBLE PRECISION  RAY1C2, RAY2C2
      DOUBLE PRECISION  HAUTC1, HAUTC2
      DOUBLE PRECISION  RPCON1(3,3), RPCON2(3,3)
      DOUBLE PRECISION  RPINT(3,1:MXPINT), RLONG(64)
      INTEGER           NPSUIV(MXPINT), NTPINT(MXPINT)
      INTEGER           N1COUR(0:64)
C
      DOUBLE PRECISION  P(3), Q(3), PT0(3), PT1(3), PT2(3), PT3(3),
     %                  PT(3), QT(3), RACIN2(2), RL, RL0, RL2,
     %                  ANGLE, ANGLE0, ANGLE1, ANGLE2, DEUXPI,
     %                  TAILAR, TAILA1, TAILA2, D20, D21, D30, D31
C
C     ICI, IL EXISTE (NBIP2-NBIP1)/2 INTERVALLES DE POINTS D'INTERSECTION
C     GENERATION DES EF SEGMENTS P3-HERMITE
C     NPSUIV(POINT)<0 INDIQUE UN POINT INTERMEDIAIRE NON SOMMET D'UN SEGMENT P3
C     =========================================================================
      DEUXPI = ATAN(1D0) * 8D0
      NBSOM  = 0
      NUR    = 0
      IF( NOFOTI .LE. 0 ) THEN
C
C        PAS DE FONCTION 'TAILLE_IDEALE'
C        -------------------------------
C        CALCUL DE LA LONGUEUR DE LA LIGNE A PARTIR DES POINTS RECENSES
         RL = 0D0
         DO 10 N = NBIP1, NBIP2, 2
C           LA LONGUEUR DE CETTE PARTIE DE LA LIGNE EST RLONG(N)
            RL = RL + RLONG(N)
 10      CONTINUE
C        LONGUEUR MOYENNE SOUHAITEE D'UN SEGMENT P3
         RL = RL / NBARLI
C
C        LE DECOUPAGE EN SEGMENTS P3 DE TAILLE VOISINE A RL
         NBARLI = 0
         DO 50 N = NBIP1, NBIP2, 2
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
C           LA LONGUEUR DE CETTE PARTIE DE LA LIGNE EST RLONG(N)
C           CALCUL DU SEGMENT P3 PAR PARCOURS DES POINTS D'INTERSECTION
            N2DER = 0
C           NOMBRE DE SOMMETS ET ARETES DE CET INTERVALLE DE POINTS
            NBAR = 0
C
C           DEBUT D'ARC OU SEGMENT P3 DE L'INTERVALLE DE POINTS
 20         NBPARC = 0
            TAILAR = 0D0
C
 30         IF( N2 .NE. NF .OR. NBA .LE. 0 ) THEN
               NBA = NBA + 1
               N1  = N2
               N2  = ABS( NPSUIV(N2) )
C              POINT N2 A CONSIDERER
               TAILAR = TAILAR + RPINT(3,N2)
               IF( TAILAR .LT. RL ) THEN
C                 POINT A AJOUTER A L'ARETE ET CONTINUER
                  NBPARC = NBPARC + 1
C                 SON SUIVANT EST MARQUE <0 COMME POINT INTERMEDIAIRE
C                 NPSUIV(N2)>0 INDIQUE UN SOMMET DES SEGMENTS P3
                  NPSUIV( N2 ) = - ABS( NPSUIV( N2 ) )
                  GOTO 30
               ELSE
C                 DERNIER POINT DE L'ARC
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
         DO 170 N = NBIP1, NBIP2, 2
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
            CALL CONERA( ANGLE0,
     %                   PAX1C1, RPCON1, RAY1C1, RAY2C1, HAUTC1,
     %                   PAX1C2, RPCON2, RAY1C2, RAY2C2, HAUTC2,
     %                   NBRAC2, RACIN2, P, Q  )
            IF( NTPINT(N2) .EQ. 3 ) THEN
C              RACINE DOUBLE
               NURAC = 1
            ELSE
C              RACINE SIMPLE
               NURAC = NTPINT(N2)
            ENDIF
            DO 110 K=1,3
               PT(K) = P(K) + RACIN2(NURAC) * (Q(K) - P(K))
 110        CONTINUE
C           RETOUR DANS LE REPERE GLOBAL INITIAL
            CALL AC3D3D( PAX1C2, RPCON2, PT, PT0 )
            CALL FONVAL( NOFOTI, 3, PT0,  NCODEV, RL0 )
            IF( NCODEV .EQ. 0 ) THEN
C               NCODEV  : 0 RL N'EST PAS INITIALISEE EN SORTIE
C                         1 RL   EST     INITIALISEE EN SORTIE
               PT1(1) = PT0(1)
               PT1(2) = PT0(2)
               PT1(3) = PT0(3)
               GOTO 9000
            ENDIF
 115        RL = ABS( RL0 )
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
               CALL CONERA( ANGLE2,
     %                      PAX1C1, RPCON1, RAY1C1, RAY2C1, HAUTC1,
     %                      PAX1C2, RPCON2, RAY1C2, RAY2C2, HAUTC2,
     %                      NBRAC2, RACIN2, P, Q  )
               NURACI = NTPINT(N2)
               IF( NURACI .EQ. 3 ) NURACI = 1
               DO 135 K=1,3
                  PT(K) = P(K) + RACIN2(NURACI) * (Q(K) - P(K))
 135           CONTINUE
C              RETOUR DANS LE REPERE GLOBAL INITIAL
               CALL AC3D3D( PAX1C2, RPCON2, PT , PT1 )
C              DISTANCE ENTRE PT0 ET PT1
               TAILAR = SQRT( (PT1(1)-PT0(1)) ** 2
     %                      + (PT1(2)-PT0(2)) ** 2
     %                      + (PT1(3)-PT0(3)) ** 2 )
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
                  IF( NTPINT(N0) .NE. NTPINT(N2) ) THEN
C                    ZONE DE RACINE DOUBLE => TRAITEMENT SPECIAL
                     IF( TAILAR .GE. RL*0.5D0 ) THEN
C                       L'ARETE PT0 PT1 EST DE TAILLE RAISONNABLE
C                       LE POINT N2 EST LE SOMMET DE L'ARETE SUIVANTE
                        RL0 = TAILAR
                     ELSE
C                       L'ARETE PT0 PT1 EST TROP PETITE
C                       LE POINT PT0 EST REPOUSSE VERS SON PRECEDENT
                        IF( N2DER0 .GT. 0 ) THEN
C                          IL EXISTE N2DER0 PRECEDENT DE N0
                           ANGLE = RPINT(1,N2) - RPINT(1,N2DER0)
                           RPINT(1,N0) = RPINT(1,N2DER0) + 0.7D0 * ANGLE
                        ELSE
C                          CHOIX DISCUTABLE
                           ANGLE = RPINT(1,N2) - RPINT(1,N0)
                           RPINT(1,N0) = RPINT(1,N0) + 0.3D0 * ANGLE
                        ENDIF
                     ENDIF
                     PT0(1) = PT1(1)
                     PT0(2) = PT1(2)
                     PT0(3) = PT1(3)
                     N2DER0 = N2DER
                     N2DER  = N2
C                    UNE ARETE DE PLUS
                     NBAR = NBAR + 1
                     GOTO 115
                  ENDIF
C                 ZONE SANS PROBLEME
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
                  IF( ABS(ANGLE2) .LE. 1D-14 ) THEN
C                    PERIODICTE DE 2 PI A PRENDRE EN COMPTE
                     ANGLE2 = ANGLE2 + DEUXPI
                  ENDIF
               ENDIF
               ANGLE = ANGLE2
C
               IF( NTPINT(N0) .NE. NTPINT(N2) ) THEN
C                 CAS A PROBLEME : N0 OU N2 EST PROCHE D'UN POINT A RACINE DOUBL
C                 LE NUMERO DE LA RACINE EST MARQUE NUL
                  NURACI = 0
               ELSE
C                 CAS A NUMERO DE RACINE SANS PB
C                 LE NUMERO DE LA RACINE DU POINT N0
                  NURACI = NTPINT(N0)
                  IF( NURACI .EQ. 3 ) NURACI = 1
               ENDIF
C              DISTANCE INITIALE ENTRE LE POINT N0 ET N2
               TAILA2 = TAILAR
C              DISTANCE INITIALE ENTRE LE POINT N0 ET N1=N0
               TAILA1 = 0D0
C
C              DICHOTOMIE POUR OBTENIR UNE DISTANCE EGALE A RL
C              -----------------------------------------------
 140           IF( ABS(ANGLE2-ANGLE1) .GT. 1D-14 ) THEN
C
                  IF( ABS(TAILAR-RL) .GT. 0.01D0*TAILAR ) THEN
C
                     ANGLE = ( ANGLE1 + ANGLE2 ) / 2
                     CALL CONERA( ANGLE,
     %                    PAX1C1, RPCON1, RAY1C1, RAY2C1, HAUTC1,
     %                    PAX1C2, RPCON2, RAY1C2, RAY2C2, HAUTC2,
     %                    NBRAC2, RACIN2, P, Q  )
C
 142                 IF( NURACI .GT. 0 ) THEN
C                       CAS LOIN D'UNE RACINE DOUBLE (SANS PB)
                        DO 145 K=1,3
                           PT(K) = P(K) + RACIN2(NURACI) * (Q(K)-P(K))
 145                    CONTINUE
C                       RETOUR DANS LE REPERE GLOBAL INITIAL
                        CALL AC3D3D( PAX1C2, RPCON2, PT , PT1 )
C                       DISTANCE ENTRE PT0 ET PT1
                        TAILAR = SQRT( (PT1(1)-PT0(1)) ** 2
     %                               + (PT1(2)-PT0(2)) ** 2
     %                               + (PT1(3)-PT0(3)) ** 2 )
                     ELSE
C                       CAS PROCHE D'UNE RACINE DOUBLE => RISQUE DE PROBLEME!
C                       IL FAUT CHOISIR PRECISEMENT LE NUMERO DE LA RACINE
C                       POUR EVITER DE PASSER DE L'AUTRE COTE DU MIN OU MAX
C                       ON CHOISIT LE POINT TEL QUI REALISE LE MINIMUM
C                       DE D(PT0,PTG)+D(PTG,PT1) POUR LE FORCER A ETRE
C                       ENTRE LES 2 POINTS PT0 ET PT1
                        DO 148 K=1,3
                           PT(K) = P(K) + RACIN2(1) * (Q(K)-P(K))
                           QT(K) = P(K) + RACIN2(2) * (Q(K)-P(K))
 148                    CONTINUE
C                       RETOUR DANS LE REPERE GLOBAL INITIAL
                        CALL AC3D3D( PAX1C2, RPCON2, PT , PT2 )
                        CALL AC3D3D( PAX1C2, RPCON2, QT , PT3 )
                        D20 = SQRT( (PT2(1)-PT0(1)) ** 2
     %                            + (PT2(2)-PT0(2)) ** 2
     %                            + (PT2(3)-PT0(3)) ** 2 )
                        D21 = SQRT( (PT2(1)-PT1(1)) ** 2
     %                            + (PT2(2)-PT1(2)) ** 2
     %                            + (PT2(3)-PT1(3)) ** 2 )
                        D30 = SQRT( (PT3(1)-PT0(1)) ** 2
     %                            + (PT3(2)-PT0(2)) ** 2
     %                            + (PT3(3)-PT0(3)) ** 2 )
                        D31 = SQRT( (PT3(1)-PT1(1)) ** 2
     %                            + (PT3(2)-PT1(2)) ** 2
     %                            + (PT3(3)-PT1(3)) ** 2 )
                        IF( D20+D21 .LT. D30+D31 ) THEN
C                          LE POINT PT2 EST ENTRE PT0 ET PT1 => IL EST CHOISI
                           PT1(1) = PT2(1)
                           PT1(2) = PT2(2)
                           PT1(3) = PT2(3)
                           TAILAR = D20
                           NUR = 1
                        ELSE
C                          LE POINT PT3 EST ENTRE PT0 ET PT1 => IL EST CHOISI
                           PT1(1) = PT3(1)
                           PT1(2) = PT3(2)
                           PT1(3) = PT3(3)
                           TAILAR = D30
                           NUR = 2
                        ENDIF
                     ENDIF
C
                     IF( TAILAR .GT. RL ) THEN
C
C                       LE POINT MILIEU DEVIENT LE POINT LE PLUS ELOIGNE DE PT0
C                       LA DISTANCE ENTRE PT0 ET PT2 DOIT DECROITRE SINON ERREUR
                        IF( TAILAR .GT. TAILA2 ) THEN
C                          LA RACINE N'EST PAS LA BONNE
                           NURACI = 3 - NURACI
                           GOTO 142
                        ENDIF
                        TAILA2 = TAILAR
                        ANGLE2 = ANGLE
C
                     ELSE
C
C                       LE POINT MILIEU DEVIENT LE POINT LE MOINS ELOIGNE DE PT0
C                       LA DISTANCE ENTRE PT0 ET PT1 DOIT CROITRE SINON ERREUR D
                        IF( TAILAR .LT. TAILA1 ) THEN
C                          LA RACINE N'EST PAS LA BONNE
                           NURACI = 3 - NURACI
                           GOTO 142
                        ENDIF
                        TAILA1 = TAILAR
                        ANGLE1 = ANGLE
C
                     ENDIF
                     GOTO 140
C
                  ENDIF
               ELSE
C                 CONVERGENCE DE L'ANGLE MAIS PT0 PT1 SONT 2 POINTS
C                 DE MEME ANGLE ET DE RACINES DIFFERENTES ET TAILAR<RL
                  WRITE(IMPRIM,*)
     %           'SP CONEAR: 2 POINTS MEME ANGLE RACINES DIFFERENTES'
               ENDIF
C
C              ANGLE ACCEPTE => UNE ARETE OU ARC OU SEGMENT P3 EN PLUS
C              LE NOUVEAU POINT EST INTERCALE ENTRE N0 ET N2
C              POUR EVITER DE PERDRE UNE PARTIE DE LA COURBE
               N1 = N1COUR(0)
               IF( N1 .LE. 0 ) GOTO 9100
               N1COUR(0) = NPSUIV(N1)
C              LE POINT N1 EST INTERCALE ENTRE N0 ET N2=NF
               NPSUIV( N1 ) = N2
C              LE NUMERO DE LA RACINE
               IF( NURACI .EQ. 0 ) NURACI = NUR
               NTPINT( N1 ) = NURACI
C              LE SUIVANT DE N0 EST DIRECTEMENT N1
               NPSUIV( N0 ) = N1
C              L'ANGLE DU POINT
               RPINT(1,N1) = ANGLE
C              LA RACINE DU POINT
               RPINT(2,N1) = RACIN2(NURACI)
C              LA DISTANCE DE N1 A SON PREDECESSEUR N0
               RPINT(3,N1) = TAILAR
C              LA DISTANCE DE N2 A SON PREDECESSEUR N1
               RPINT(3,N2) = RPINT(3,N2) - TAILAR
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
