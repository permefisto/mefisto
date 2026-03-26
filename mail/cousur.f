      SUBROUTINE COUSUR( NDIM,   NBSOM,  XYZSOM, NBEFOB, NOSOEF,
     %                   L1ARFA, L2ARFA, NARFA,  NBARXF, NAR1F,
     %                   NADAR1F, NSTCH, NEWNST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ESSAI DE COUTURE DES ARETES SIMPLES POUR FERMER LE MAILLAGE
C -----    DE TRIANGLES-QUADRANGLES D'UNE SURFACE
C
C ENTREES:
C --------
C NDIM   : DIMENSION 2 ou 3 DE L'ESPACE DE LA SURFACE
C NBSOM  : NOMBRE DE SOMMETS DE XYZSOM
C XYZSOM : LES 3 COORDONNEES DES SOMMETS DE L'OBJET
C NBSOEF : =4 NOMBRE DE SOMMETS PAR QT
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE DE FACES DU TABLEAU NARFA

C MODIFIES:
C ---------
C NBEFOB : NOMBRE D'EF TRIANGLE ou QUADRANGLE
C NOSOEF : NUMERO DES 4 SOMMETS DE CHAQUE QT
C ATTENTION: NARFA EST ICI MODIFIE  A RECONSTRUIRE EN DEHORS DE CE SP
C            LE HACHAGE SUR LES ARETES NE FONCTIONNANT PLUS EN SORTIE
C NARFA  : NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          NARFA(4:3+MXFAAR,I)= NO NOSOFA DE LA FACE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE MXFAAR(=L1ARFA-3) FACES, 
C          LE NUMERO DE FACE MXFAAR EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES FACES ADJACENTES EST INCOMPLETE
C NBARXF : NBARXF(n) NOMBRE D'ARETES APPARTENANT A n FACES n=1,2,3
C                    POUR n=3 APPARTENANT a >2 FACES
C NAR1F  : NUMERO DANS NARFA DES ARETES SIMPLES
C NADAR1F: NUMERO DANS NAR1F DE L'ARETE SIMPLE DE LA LIGNE 
C          ACTUELLE D'IDENTIFICATION DE SOMMETS OPPOSES
C NSTCH: NUMERO XYZSOM DES SOMMETS DES ARETES SIMPLES CHAINEES

C SORTIE :
C --------
C NEWNST : NOUVEAU NUMERO DES NBSOM SOMMETS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY NOVEMBRE 2015
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NOSOEF(4,NBEFOB), NEWNST(0:NBSOM)
      INTEGER           NARFA(L1ARFA,L2ARFA)
      INTEGER           NBARXF(3), NAR1F(*), NADAR1F(*), NSTCH(*)
      REAL              DARETE(3)

C     INITIALISATION DU NOUVEAU NUMERO DES SOMMETS IDENTIFIES
      DO N = 0, NBSOM
         NEWNST( N ) = N
      ENDDO

C     CONSTRUCTION DE LA LISTE DES ARETES SIMPLES DANS NAR1F
C     QUI SE SUIVENT PAR UN SOMMET COMMUN
C     L'ORDRE ENTRAINE LA CONNAISSANCE DE LA SUIVANTE
C     ------------------------------------------------------
      NBAR1F = NBARXF(1)

 10   PRINT *,'cousur:',NBAR1F,' ARETES SIMPLES'
      IF( NBAR1F .LE. 2 ) GOTO 9000

C     A PARTIR DE LA PREMIERE ARETE SIMPLE
      NBARCH = 0
      NBSTID = 0
      NBARID = 0
      NA0 = NAR1F( 1 )
      IF( NA0 .EQ. 0 ) GOTO 200
      NS1 = NEWNST( NARFA( 1, NA0 ) )
      NS2 = NEWNST( NARFA( 2, NA0 ) )
      NBARCH = 1
      NADAR1F( NBARCH ) = 1
C     ARETE SIMPLE MARQUEE
      NAR1F( 1 ) = -NA0

C     CONSTRUCTION DE LA LISTE NAR1F DES ARETES SUIVANTES DE NA0 PAR NS2
C     FORMANT UNE LIGNE POUR IDENTIFIER LES SOMMETS OPPOSES
C     ------------------------------------------------------------------
 15   DO 20 N = 2, NBAR1F
         NA1 = NAR1F( N )
         IF( NA1 .LE. 0 ) GOTO 20
         NS3 = NEWNST( NARFA( 1, NA1 ) )
         NS4 = NEWNST( NARFA( 2, NA1 ) )
         IF( NS3 .EQ. NS2 ) THEN
            NBARCH = NBARCH + 1
            NADAR1F( NBARCH ) = N
            NS2 = NS4
            NAR1F( N ) = -NA1
            GOTO 15
         ELSE IF( NS4 .EQ. NS2 ) THEN
            NBARCH = NBARCH + 1
            NADAR1F( NBARCH ) = N
            NS2 = NS3
            NAR1F( N ) = -NA1
            GOTO 15
         ENDIF
 20   ENDDO

C     NS2 N'A PAS ETE RETROUVE PARMI LES ARETES SIMPLES NON TRAITEES
      IF( NBARCH .LE. 1 ) GOTO 100

C     LES ARETES SIMPLE MARQUEES SONT DEMARQUEES
      DO N=1,NBARCH
         NA = NADAR1F( N )
         NAR1F( NA ) = ABS( NAR1F( NA ) )
      ENDDO

C     CONSTRUCTION DU TABLEAU DES SOMMETS DES ARETES SIMPLES CHAINEES
C     ---------------------------------------------------------------
      NBSTCH = 0
      DO N = 1, NBARCH
         NA  = NADAR1F( N )
         NA1 = NAR1F( NA )
         DO 22 M=1,2
            NS = NEWNST( NARFA( M, NA1 ) )
            DO K=1,NBSTCH
               IF( NS .EQ. NSTCH(K) ) GOTO 22
            ENDDO
C           NS EST AJOUTE AU TABLEAU NSTCH
            NBSTCH = NBSTCH + 1
            NSTCH(NBSTCH) = NS
 22      ENDDO
      ENDDO

C     LA LIGNE EST ELLE FERMEE?
      WRITE(IMPRIM,*)
      N   = NADAR1F( NBARCH )
      NA1 = NAR1F( N )
      NS3 = NEWNST( NARFA( 1, NA1 ) )
      NS4 = NEWNST( NARFA( 2, NA1 ) )
      IF( (NS1 .EQ. NS3 .OR. NS1 .EQ. NS4).AND.NBARCH .EQ. NBSTCH ) THEN

C        LA LIGNE EST FERMEE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'LIGNE FERMEE de',NBARCH,
     %           ' ARETES SIMPLES de',NBSTCH,' SOMMETS FORMANT une BOU
     %CLE  -----------------'
         ELSE
            WRITE(IMPRIM,*) 'CLOSED LINE of',NBARCH,
     %           ' SINGLE EDGES of',NBSTCH,' VERTICES FORMING a LOOP -
     %-----------------------'
         ENDIF

      ELSE

C        LA LIGNE N'EST PAS FERMEE
         IF( LANGAG .EQ. 0 ) THEN
           WRITE(IMPRIM,*)'LIGNE NON FERMEE de',NBARCH,' ARETES SIMPLES 
     %de',NBSTCH,' SOMMETS'
         ELSE
            WRITE(IMPRIM,*)'NOT CLOSED LINE of',NBARCH,' SINGLE EDGES of
     %',NBSTCH,' VERTICES'
         ENDIF

      ENDIF

C     TRACE DES ARETES SIMPLES ET FACES ADJACENTES DE CETTE LIGNE
C     CADRE EXTREME DES ARETES SIMPLES DE LA LIGNE
C     -----------------------------------------------------------
      DO K=1,3
         R = XYZSOM(K,NS1)
         COOEXT(K,1) = R
         COOEXT(K,2) = R
      ENDDO

      DO 30 N = 1, NBARCH
         NA1 = NAR1F( NADAR1F( N ) )
         IF( NA1 .LE. 0 ) GOTO 30
         DO M=1,2
            NS = NEWNST( NARFA( M, NA1 ) )
            DO K=1,3
               R = XYZSOM(K,NS)
               COOEXT(K,1) = MIN( COOEXT(K,1), R )
               COOEXT(K,2) = MAX( COOEXT(K,2), R )
            ENDDO
         ENDDO
 30   ENDDO

      DO K=1,3
         R = COOEXT(K,2) - COOEXT(K,1)
         COOEXT(K,1) = COOEXT(K,1) - R
         COOEXT(K,2) = COOEXT(K,2) + R
      ENDDO

C     LA VISEE SELON LES COOEXT DU COMMON / XYZEXT /
      CALL VISEE0

C     TRACER EN MAGENTA TOUTES LES ARETES APPARTENANT A  1 FACE
C            EN GRIS    TOUTES LES ARETES APPARTENANT A  2 FACES
C            EN ROUGE   TOUTES LES ARETES APPARTENANT A >2 FACES
ccc      INTERA0 = INTERA
ccc      INTERA  = 3
      CALL TRARFASU( NDIM,   NBSOM,  XYZSOM, NOSOEF,
     %               L1ARFA, L2ARFA, NARFA,  NBARXF )
ccc      INTERA = INTERA0

C     TRAITEMENT DU CAS PARTICULIER DE 3 ARETES SIMPLES EN BOUCLE
C     -----------------------------------------------------------
      IF( NBARCH .EQ. 3 .AND. NBSTCH .EQ. 3 ) THEN

C        RECHERCHE DE LA LONGUEUR DES 3 ARETES CHAINEES
         KMAX = 0
         DMAX = 0
         DO K = 1, 3

C           L'ARETE CHAINEE K
            NN  = NADAR1F( K )
            NA1 = NAR1F( NN )
            IF( NA1 .LE. 0 ) GOTO 100

C           LE SOMMET M DE L'ARETE CHAINEE KK
            NS3 = NEWNST( NARFA( 1, NA1 ) )
            NS4 = NEWNST( NARFA( 2, NA1 ) )

C           CALCUL DU CARRE DE LA DISTANCE NS3-NS4
            D = ( XYZSOM(1,NS4) - XYZSOM(1,NS3) ) ** 2
     %        + ( XYZSOM(2,NS4) - XYZSOM(2,NS3) ) ** 2
     %        + ( XYZSOM(3,NS4) - XYZSOM(3,NS3) ) ** 2
            DARETE(K) = D
            PERIM = PERIM + D

            IF( D .GT. DMAX ) THEN
               DMAX = D
               KMAX = K
            ENDIF

         ENDDO
          
C        LES 3 ARETES FORMENT ELLES UN TROU?
         D = SQRT( DMAX )
         PERIM = SQRT(DARETE(1)) + SQRT(DARETE(2)) + SQRT(DARETE(3))
         IF( PERIM-D .GT. 1.6*D ) GOTO 100

C        NON: MISE EN PREMIERE POSITION DE L'ARETE SIMPLE MAX
         IF( KMAX .EQ. 2 ) THEN
             K1 = NADAR1F( 1 )
             K2 = NADAR1F( 2 )
             K3 = NADAR1F( 3 )
             NADAR1F( 1 ) = K2
             NADAR1F( 2 ) = K3
             NADAR1F( 3 ) = K1
         ELSE IF( KMAX .EQ. 3 ) THEN
             K1 = NADAR1F( 1 )
             K2 = NADAR1F( 2 )
             K3 = NADAR1F( 3 )
             NADAR1F( 1 ) = K3
             NADAR1F( 2 ) = K1
             NADAR1F( 3 ) = K2
         ENDIF

C         ARETE SIMPLE NS1-NS2
          NA0 = NAR1F( NADAR1F( 1 ) )
          NS1 = NARFA(1,NA0)
          NS2 = NARFA(2,NA0)
C         LA FACE UNIQUE DE L'ARETE NS1-NS2
          NF  = NARFA( 4, NA0 )

C         NSM "MILIEU" DE L'ARETE NA0 NS1-NS2
          NA1 = NAR1F( NADAR1F( 2 ) )
          NSM = NARFA(1,NA1)
          IF( NSM .EQ. NS2 .OR. NSM .EQ. NS1 ) THEN
             NSM = NARFA(2,NA1)
          ENDIF
C         LA FACE UNIQUE DE L'ARETE NS1-NS2
          NF1 = NARFA( 4, NA1 )

C         DECOUPAGE DE LA FACE NF
          IF( NOSOEF(4,NF) .EQ. 0 ) THEN
             NBS = 3
          ELSE
             NBS = 4
          ENDIF

C         K1 POSITION DE NS1 DANS NF
          DO K1=1,NBS
             IF( NOSOEF(K1,NF) .EQ. NS1 ) GOTO 34
          ENDDO

C         RECHERCHE DE LA POSITION DE NS2 DANS NF
 34       DO K2=1,NBS
             IF( NOSOEF(K2,NF) .EQ. NS2 ) GOTO 36
          ENDDO

C         NSM EST ENTRE NS1-NS2 POSITION DE K2 PAR RAPPORT A K1
 36       K = K1 + 1
          IF( K .GT. NBS ) K=1

          IF( NOSOEF(4,NF) .EQ. 0 ) THEN

C            NF EST UN TRIANGLE => DECOUPAGE EN 2 TRIANGLES
C            A PARTIR DE NSM MILIEU DE L'ARETE NS1-NS2
C            RECHERCHE DES SOMMETS POUR GARDER LE SENS DE NF
             NBEFOB = NBEFOB + 1
             IF( K .EQ. K2 ) THEN

C               NS3 LE 3-EME SOMMET DE NF
                IF( K .LT. 3 ) THEN
                   K = K+1
                ELSE
                   K = 1
                ENDIF
                NS3 = NOSOEF(K,NF)

                NOSOEF(1,NBEFOB) = NS3
                NOSOEF(2,NBEFOB) = NS1
                NOSOEF(3,NBEFOB) = NSM
                NOSOEF(4,NBEFOB) = 0

C               NF EST MODIFIE
                NOSOEF(1,NF) = NS3
                NOSOEF(2,NF) = NSM
                NOSOEF(3,NF) = NS2

             ELSE

C               NS3 LE 3-EME SOMMET DE NF
                NS3 = NOSOEF(K,NF)
                NOSOEF(1,NBEFOB) = NS3
                NOSOEF(2,NBEFOB) = NS2
                NOSOEF(3,NBEFOB) = NSM
                NOSOEF(4,NBEFOB) = 0

C               NF EST MODIFIE
                NOSOEF(1,NF) = NS3
                NOSOEF(2,NF) = NSM
                NOSOEF(3,NF) = NS1

             ENDIF

         ELSE

C            NF EST UN QUADRANGLE => 3 TRIANGLES
C            A PARTIR DE NSM MILIEU DE L'ARETE NS1-NS2
             IF( K .EQ. K2 ) THEN

C               NS3 LE 3-EME SOMMET DE NF
                IF( K .LT. 4 ) THEN
                   K = K+1
                ELSE
                   K = 1
                ENDIF
                NS3 = NOSOEF(K,NF)

C               NS4 LE 4-EME SOMMET DE NF
                IF( K .LT. 4 ) THEN
                   K = K+1
                ELSE
                   K = 1
                ENDIF
                NS4 = NOSOEF(K,NF)

                NBEFOB = NBEFOB + 1
                NOSOEF(1,NBEFOB) = NSM
                NOSOEF(2,NBEFOB) = NS4
                NOSOEF(3,NBEFOB) = NS1
                NOSOEF(4,NBEFOB) = 0

                NBEFOB = NBEFOB + 1
                NOSOEF(1,NBEFOB) = NSM
                NOSOEF(2,NBEFOB) = NS3
                NOSOEF(3,NBEFOB) = NS4
                NOSOEF(4,NBEFOB) = 0

C               NF EST MODIFIE
                NOSOEF(1,NF) = NSM
                NOSOEF(2,NF) = NS2
                NOSOEF(3,NF) = NS3
                NOSOEF(4,NF) = 0

             ELSE

C               NS2 PRECEDE NS1 DANS NF. IL DEVIENT NS4
                NS20 = NS2
                NS4 = NS2

C               NS2 LE SOMMET SUIVANT DE NS1 DANS NF
                IF( K1 .LT. 4 ) THEN
                   K = K1+1
                ELSE
                   K = 1
                ENDIF
                NS2 = NOSOEF(K,NF)

C               NS3 LE 3-EME SOMMET DE NF
                IF( K .LT. 4 ) THEN
                   K = K+1
                ELSE
                   K = 1
                ENDIF
                NS3 = NOSOEF(K,NF)

                NBEFOB = NBEFOB + 1
                NOSOEF(1,NBEFOB) = NSM
                NOSOEF(2,NBEFOB) = NS1
                NOSOEF(3,NBEFOB) = NS2
                NOSOEF(4,NBEFOB) = 0

                NBEFOB = NBEFOB + 1
                NOSOEF(1,NBEFOB) = NSM
                NOSOEF(2,NBEFOB) = NS2
                NOSOEF(3,NBEFOB) = NS3
                NOSOEF(4,NBEFOB) = 0

C               NF EST MODIFIE
                NOSOEF(1,NF) = NSM
                NOSOEF(2,NF) = NS3
                NOSOEF(3,NF) = NS4
                NOSOEF(4,NF) = 0

C               RESTAURATION DE NS2 INITIAL
                NS2 = NS20

            ENDIF

         ENDIF

C        LES 3 ARETES SONT MARQUEES AVEC 2 FACES (PAS FORCEMENT LES VRAIES)
         DO KK = 1, 3

C           L'ARETE CHAINEE KK DANS NARFA
            NN  = NADAR1F( KK )
            NA1 = NAR1F( NN )
            IF( NARFA(4,NA1) .EQ. NF ) THEN
               NARFA(5,NA1) = NBEFOB
            ELSE
               NARFA(5,NA1) = NF
            ENDIF

C           ARETE TRAITEE => N'EST PLUS SIMPLE
            NAR1F( NN ) = 0

         ENDDO

C        LES 3 ARETES SIMPLES SONT TRAITEES
         NBARID = 3
         NBSTID = 1
         GOTO 100

      ENDIF

C     POUR CHAQUE SOMMET DE LA LIGNE RECHERCHE DU PLUS PROCHE SOMMET
C     POUR TENTER SON IDENTIFICATION ET COUDRE LES 2 SOMMETS
C     --------------------------------------------------------------
      DO 50 K = 1, NBSTCH

C        NUMERO XYZSOM DU SOMMET K DU CHAINAGE
         NS1 = NSTCH( K )
         IF( NS1 .LE. 0 ) GOTO 50

C        CALCUL DU CARRE DE LA DISTANCE NS1-NS2
         KMIN  = 0
         NSMIN = 0
         DMIN  = 1E28
         DS1S2 = 0

C        RECHERCHE DU PLUS PROCHE SOMMET DE NS1
         DO 38 KK = 1, NBARCH

C           L'ARETE CHAINEE KK
            NN  = NADAR1F( KK )
            NA1 = NAR1F( NN )
            IF( NA1 .LE. 0 ) GOTO 38

            DO M=1,2

C              LE SOMMET M DE L'ARETE CHAINEE KK
               NS = NEWNST( NARFA( M, NA1 ) )

C              CALCUL DU CARRE DE LA DISTANCE NS-NS1
               D = ( XYZSOM(1,NS) - XYZSOM(1,NS1) ) ** 2
     %           + ( XYZSOM(2,NS) - XYZSOM(2,NS1) ) ** 2
     %           + ( XYZSOM(3,NS) - XYZSOM(3,NS1) ) ** 2

               IF( NS .NE. NS1 ) THEN
                  IF( D .LT. DMIN ) THEN
                     KMIN  = KK
                     DMIN  = D
                     NSMIN = NS
                  ENDIF
               ELSE
                  IF( M .EQ. 1 ) THEN
                     MM = 2
                  ELSE
                     MM = 1
                  ENDIF
                  NA0 = NA1
                  NS2 = NEWNST( NARFA( MM, NA1 ) )

C                 CALCUL DU CARRE DE LA DISTANCE NS1-NS2
                  D = ( XYZSOM(1,NS2) - XYZSOM(1,NS1) ) ** 2
     %              + ( XYZSOM(2,NS2) - XYZSOM(2,NS1) ) ** 2
     %              + ( XYZSOM(3,NS2) - XYZSOM(3,NS1) ) ** 2
                  DS1S2 = MAX( DS1S2, SQRT( D ) )
               ENDIF

            ENDDO

 38      ENDDO

         IF( NSMIN .EQ. 0 ) GOTO 50
         IF( DS1S2 .EQ. 0 ) GOTO 50

         DMIN = SQRT( DMIN )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10100) DMIN, NS1, NSMIN, DS1S2, DMIN/DS1S2
            WRITE(IMPRIM,10101) NS1,  (XYZSOM(KK,NS1  ),KK=1,3)
            WRITE(IMPRIM,10101) NSMIN,(XYZSOM(KK,NSMIN),KK=1,3)
         ELSE
            WRITE(IMPRIM,20100) DMIN, NS1, NSMIN, DS1S2, DMIN/DS1S2
            WRITE(IMPRIM,20101) NS1,  (XYZSOM(KK,NS1  ),KK=1,3)
            WRITE(IMPRIM,20101) NSMIN,(XYZSOM(KK,NSMIN),KK=1,3)
         ENDIF

10100 FORMAT(/' DISTANCE D=',G13.6,' ENTRE le SOMMET',I10,
     %        ' et SON PLUS PROCHE SOMMET OPPOSE',I10,
     %        ' des ARETES SIMPLES, LONGUEUR ARETE L=',G13.6,
     %        ' D/L=',G13.6)

20100 FORMAT(/' DISTANCE D=',G13.6,' BETWEEN the VERTEX',I10,
     %        ' and its NEAREST OPPOSED VERTEX',I10,
     %        ' of SINGLE EDGES, EDGE LENGTH L=',G13.6,' D/L=',G13.6)

10101 FORMAT(' SOMMET',I10,' : X=',G14.7,' Y=',G14.7,' Z=',G14.7)
20101 FORMAT(' VERTEX',I10,' : X=',G14.7,' Y=',G14.7,' Z=',G14.7)

C        IDENTIFICATION DE NS1 NSMIN POSSIBLE?
         IF( DMIN .LE. 0.666 * DS1S2 ) THEN

C           OUI: IDENTIFICATION DE NS1 <- NSMIN?
C           -----------------------------------
C           NS1 EST IL UN SOMMET DU MAILLAGE de la SURFACE?
            CALL STDUMAIL( NS1, 4, NBEFOB, NOSOEF, NONOUI )
C           NONOUI=1 NS1 EST UN SOMMET DES EF ACTIFS DU MAILLAGE
C                 =0 SINON
            IF( NONOUI .EQ. 0 ) GOTO 50

C           NSMIN EST IL UN SOMMET DU MAILLAGE de la SURFACE?
            CALL STDUMAIL( NSMIN, 4, NBEFOB, NOSOEF, NONOUI )
C           NONOUI=1 NSMIN EST UN SOMMET DES EF ACTIFS DU MAILLAGE
C                 =0 SINON
            IF( NONOUI .EQ. 0 ) GOTO 50

            IF( NS1 .GT. NSMIN ) THEN
               NN    = NS1
               NS1   = NSMIN
               NSMIN = NN
            ENDIF

C           NS1 <- NSMIN ?
            NS = NEWNST( NSMIN )
            IF( NS .NE. NSMIN ) THEN

C              PROBLEME: SOMMET NSMIN DEJA IDENTIFIE
               GOTO 50

ccc               WRITE(IMPRIM,*) 'COUSUR: PB SOMMET',NSMIN,
ccc     %                        ' DEJA IDENTIFIE AU SOMMET',NS
ccc               WRITE(IMPRIM,*) 'COUSUR: ESSAI',NSMIN,'<->',NS
ccc               IF( NS1 .NE. NEWNST(NS1) ) THEN
ccc                  NN    = NS1
ccc                  NS1   = NSMIN
ccc                  NSMIN = NN
ccc               ELSE
ccc                  WRITE(IMPRIM,*) 'COUSUR: PB SOMMET',NS1,
ccc     %                           ' DEJA IDENTIFIE AU SOMMET',NEWNST(NS1)
ccc                  WRITE(IMPRIM,*) 'COUSUR: PROBLEME RESTANT A RESOUDRE!'
ccc                  WRITE(IMPRIM,*)
ccc               ENDIF

            ENDIF

C           NS1 <- NSMIN
            IF( LANGAG .EQ. 0 ) THEN
              WRITE(IMPRIM,*)'Le SOMMET',NSMIN,' DEVIENT le SOMMET',NS1,
     %                       ': XYZ=',(XYZSOM(NN,NS1),NN=1,3)
            ELSE
              WRITE(IMPRIM,*)'VERTEX',NSMIN,' IDENTIFIED to VERTEX',NS1,
     %                       ': XYZ=',(XYZSOM(NN,NS1),NN=1,3)
            ENDIF
            NEWNST( NSMIN ) = NS1
            NBSTID = NBSTID + 1

C           LE SOMMET K DU CHAINAGE EST MARQUE
            NSTCH( K ) = -NS1
            DO NN=1,NBSTCH
               IF( NSTCH(NN) .EQ. NSMIN ) THEN
                  NSTCH(NN) = -NSMIN
                  GOTO 40
               ENDIF
            ENDDO

C           LE SOMMET NS1 ET NSMIN DEVIENT LE MILIEU DE L'ARETE NS1-NSMIN
 40         DO NN=1,3
               XYZSOM(NN,NS1  ) = (XYZSOM(NN,NS1) + XYZSOM(NN,NSMIN)) /2
               XYZSOM(NN,NSMIN) =  XYZSOM(NN,NS1)
            ENDDO

C           IDENTIFICATION DES ARETES SIMPLES CHAINEES
C           RECHERCHE DES ARETES DE SOMMET NS1 IDENTIQUES
            DO 45 KK = 1, NBARCH
C              L'ARETE CHAINEE KK
               NA1 = NAR1F( NADAR1F( KK ) )
               IF( NA1 .LE. 0 ) GOTO 45
               DO M=1,2
C                 LE SOMMET M DE L'ARETE CHAINEE KK
                  NS = NEWNST( NARFA( M, NA1 ) )
                  IF( NS .EQ. NS1 ) THEN

                     DO 42 KKK = KK+1, NBARCH
C                       L'ARETE CHAINEE KKK
                        NA2 = NAR1F( NADAR1F( KKK ) )
                        IF( NA2 .LE. 0 ) GOTO 42
                        DO MM=1,2
                           NS = NEWNST( NARFA( MM, NA2 ) )
                           IF( NS .EQ. NS1 ) THEN
C                             LES 2 ARETES NA1 ET NA2 ONT LE SOMMET NS1
                              IF( M .EQ. 1 ) THEN
                                 L=2
                              ELSE
                                 L=1
                              ENDIF
                              NS3 = NEWNST( NARFA( L, NA1 ) )
                              IF( MM .EQ. 1 ) THEN
                                 L=2
                              ELSE
                                 L=1
                              ENDIF
                              NS4 = NEWNST( NARFA( L, NA2 ) )

                              IF( NS3 .EQ. NS4 ) THEN
C                                LES 2 ARETES NA1 ET NA2 ONT LES SOMMETS NS1-NS3
C                                LES 2 ARETES SIMPLES SONT ALORS DANS 2 FACES
                                 NAA = NAR1F( NADAR1F( KK  ) )
                                 NF1 = NARFA(4,NAA)
                                 NAAA= NAR1F( NADAR1F( KKK ) )
                                 NF2 = NARFA(4,NAAA)
                                 NARFA(5,NAA ) = NF2
                                 NARFA(5,NAAA) = NF1
C                                LES 2 ARETES SIMPLES SONT IDENTIFIEES
                                 NAR1F( NADAR1F( KK  ) ) = 0
                                 NAR1F( NADAR1F( KKK ) ) = 0
                                 NBARID = NBARID + 2
                                 GOTO 45
                              ENDIF

                           ENDIF
                        ENDDO
 42                  ENDDO
                  ENDIF
               ENDDO
 45         ENDDO

            GOTO 50

         ENDIF

 50   ENDDO

C     AFFICHAGE DU NOMBRE DE SOMMETS IDENTIFIES SUR CETTE LIGNE
C     ---------------------------------------------------------
 100  WRITE(IMPRIM,*) NBSTID,' SOMMETS IDENTIFIES',
     %                NBARID,' ARETES IDENTIFIEES pour les',
     %                NBARCH,' ARETES SIMPLES CHAINEES'
      IF( NBSTID .GT. 0 .AND. NBARID .GT. 0 ) GOTO 200

C     MARQUAGE DANS NAR1F DE TOUTES LES ARETES SIMPLES CHAINEES TRAITEES
C     ------------------------------------------------------------------
      DO K = 1, NBARCH
         NAR1F( NADAR1F( K ) ) = 0
      ENDDO

C     SUPPRESSION DE NAR1F DES ARETES SIMPLES CHAINEES TRAITEES
C     ---------------------------------------------------------  
 200  NB = 0
      DO N = 1, NBAR1F
         NA1 = NAR1F( N )
         IF( NA1 .GT. 0 ) THEN
            NB = NB + 1
            NAR1F( NB ) = NA1
         ENDIF
      ENDDO

C     NOUVEAU NOMBRE DES ARETES SIMPLES
      NBAR1F = NB
      GOTO 10

C     SUPPRESSION DES SOMMETS IDENTIFIES ET DES EF DEGENERES
C     DU MAILLAGE DE LA SURFACE
C     ------------------------------------------------------
 9000 CALL MAJSTEFSURF( NBSOM, XYZSOM, NEWNST, NBEFOB, NOSOEF )

      RETURN
      END
