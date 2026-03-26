      SUBROUTINE ADSUPO( NS,     PXYD,
     %                   N1TRVI, NOTRIA, NOTRSO,
     %                   LESUTR, NOSUTR,
     %                   MXETRI, NAETOI, NARMIN, N1ARCF, NOARCF,
     %                   NBTRCF, NOTRCF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   SUPPRIMER LE SOMMET NS D'UNE TRIANGULATION ADAPTEE
C -----   TRIANGULATION FRONTALE POUR BOUCHER LES ETOILES
C
C ENTREES:
C --------
C NS     : NUMERO DU SOMMET A DETRUIRE DE LA TRIANGULATION
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C LESUTR : 0 SI PAS DE TABLEAU NOSUTR
C          1 SI LE TABLEAU NOSUTR EXISTE ET DOIT ETRE UTILISE
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : POINTE DANS NOTRIA VERS LE PREMIER TRIANGLE VIDE
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                          ADJACENT PAR L'ARETE i
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C NOSUTR : NUMERO DE SURFACE DES TRIANGLES
C
C AUXILIAIRES :
C -------------
C MXETRI : DIMENSION DU SECOND INDICE DES TABLEAUX AUXILIAIRES
C NAETOI : TABLEAU (4,MXETRI) AUXILIAIRE
C NARMIN : TABLEAU (MXETRI)   AUXILIAIRE
C N1ARCF : TABLEAU (0:MXETRI) AUXILIAIRE
C NOARCF : TABLEAU (3,MXETRI) AUXILIAIRE
C NOTRCF : TABLEAU (MXETRI)   AUXILIAIRE
C
C SORTIE :
C --------
C NBTRCF : NOMBRE DE TRIANGLES APRES SUPPRESSION DU SOMMET NS
C NOTRCF : NUMERO DES TRIANGLES DE L'ETOILE APRES SUPPRESSION
C IERR   : 0 SI PAS D'ERREUR   NS A ETE SUPPRIME DE LA TRIANGULATION
C          1 SI LE SOMMET NE PEUT ETRE SUPPRIME
C          2 SI SOMMET APPARTENANT A AUCUN TRIANGLE
C          3 SI CONTOUR FERME REDUIT A MOINS DE 3 ARETES
C          4 SATURATION DES TABLEAUX AUXILIAIRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1994
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
C
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,*),
     %                  NOTRSO(*),
     %                  NOSUTR(*),
     %                  N1ARCF(0:MXETRI),
     %                  NOARCF(3,MXETRI)
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NAETOI(4,MXETRI),
     %                  NARMIN(MXETRI),
     %                  NOTRCF(MXETRI)
C
      INTEGER           NOSUF(2)
C
      IERR = 0
      NA2  = 0
C
C     NT0 EST LE NUMERO DANS NOTRIA DU PREMIER TRIANGLE DE SOMMET NS
      NT0 = NOTRSO( NS )
 5    IF( NT0 .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NS
         KERR(1) ='ADSUPO: SOMMET ' // KERR(MXLGER)(1:10) //
     %            ' DANS AUCUN TRIANGLE'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LES NBTRCF TRIANGLES DE SOMMET NS
C     ---------------------------------
      CALL TRC1ST( NS, NT0, NOTRIA, MXETRI, NBTRCF, NOTRCF )
      IF( NBTRCF .EQ. 0 ) THEN
         NT0 = 0
         GOTO 5
      ENDIF
      IF( NBTRCF .LT. 0 )  NBTRCF = -NBTRCF
C
C     RECHERCHE DU NOMBRE DE SURFACES DES TRIANGLES
C     ---------------------------------------------
      IF( LESUTR .EQ. 0 ) THEN
C        PAS DE TABLEAU NOSUTR
         NBSURF = 1
         NBTRS1 = NBTRCF
         NOSUF(1) = 1
         GOTO 80
      ENDIF
C
      NBSURF = 0
      DO 20 I=1,NBTRCF
         NUS = NOSUTR( NOTRCF(I) )
         DO 10 J=1,NBSURF
            IF( NUS .EQ. NARMIN(J) ) GOTO 20
 10      CONTINUE
C        AJOUT DE LA NOUVELLE SURFACE
         NBSURF = NBSURF + 1
         NARMIN(NBSURF) = NUS
 20   CONTINUE
C
      IF( NBSURF .GT. 2 ) THEN
         NBLGRC(NRERR) = 3
         WRITE(KERR(MXLGER)(1:10),'(I10)') NS
         KERR(1) ='LE SOMMET' // KERR(MXLGER)(1:10)
         KERR(2) = 'APPARTIENT A PLUS DE 2 SURFACES'
         KERR(3) = 'IL NE PEUT ETRE SUPPRIME DANS ADSUPO'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     SAUVEGARDE DU NUMERO DES 2 SURFACES
      NOSUF(1) = NARMIN(1)
      IF( NBSURF .EQ. 2 ) THEN
         NOSUF(2) = NARMIN(2)
      ELSE
         NOSUF(2) = NOSUF(1)
      ENDIF
C
C     LE NOMBRE DE TRIANGLES DE LA PREMIERE SURFACE
      NBTRS1 = NBTRCF
      IF( NBSURF .EQ. 2 ) THEN
C
C        REMISE EN ORDRE DES TRIANGLES SELON LEUR NUMERO DE SURFACE
         NS1 = 2
         NS2 = NBTRCF
 30      IF( NS1 .LT. NS2 ) THEN
            IF( NOSUTR(NOTRCF(NS1)) .EQ. NOSUF(1) ) THEN
C              PASSAGE AU SUIVANT
               NS1 = NS1 + 1
               GOTO 30
            ENDIF
C           NS1 EST UN TRIANGLE MAL PLACE
  40        IF( NOSUTR(NOTRCF(NS2)) .NE. NOSUF(1) ) THEN
C              TRIANGLE BIEN PLACE
               NS2 = NS2 - 1
               GOTO 40
            ENDIF
C           ECHANGE DES 2 TRIANGLES MAL PLACES
            J           = NOTRCF(NS1)
            NOTRCF(NS1) = NOTRCF(NS2)
            NOTRCF(NS2) = J
            NS1 = NS1 + 1
            NS2 = NS2 - 1
         ENDIF
         IF( NOSUTR(NOTRCF(NS1)) .NE. NOSUF(1) ) NS1 = NS1 - 1
C        LE NOMBRE DE TRIANGLES DE LA PREMIERE DES 2 SURFACES
         NBTRS1 = NBTRCF
      ENDIF
C
C     FORMATION DES ARETES UNIQUES DU CF DES ARETES PERIPHERIQUES
C     -----------------------------------------------------------
C     REINITIALISATION A VIDE DES ARETES DE L'ETOILE
C     FORMEE DES ARETES VUES UNE FOIS DANS LES TRIANGLES DE L'ETOILE
 80   N1AEVI = 1
      N1AEOC = 0
      DO 90 I=1,MXETRI
C        NUMERO DANS NAETOI DE L'ARETE SUIVANTE
         NAETOI(4,I) = I+1
 90   CONTINUE
      NAETOI(4,MXETRI) = 0
C
C     TRAITEMENT DES TRIANGLES DE CHACUNE DES 2 SURFACES DE L'ETOILE
C     ==============================================================
      NBARCF = 0
      N1TR   = 0
      DO 1000 NSU = 1,NBSURF
C
C        FORMATION DES ARETES DE L'ETOILE
         DO 100 I=1,NBTRCF
C           AJOUT OU RETRAIT DES 3 ARETES DU TRIANGLE NOTRCF(I) A L'ETOILE
            NT = NOTRCF(I)
C           COMPARAISON DU NUMERO DES SURFACES
            IF( LESUTR .NE. 0 ) THEN
               IF( NOSUTR(NT) .NE. NOSUF( NSU ) ) GOTO 100
            ENDIF
            CALL AJTRET( NT, NOTRIA, N1AEVI, N1AEOC, NAETOI )
 100     CONTINUE
C
C        MODIFICATION DU CONTENU DU TABLEAU NAETOI
C        LE NO DE TRIANGLE            => NO TRIANGLE AU DELA DE L'ARETE
C        LE NO LOCAL DANS LE TRIANGLE => NO 1-ER SOMMET DE L'ARETE
C        LE NO INUTILISE              => NO 2-ME SOMMET DE L'ARETE
C        --------------------------------------------------------------
         NA1 = N1AEOC
C        BOUCLE SUR LES ARETES DE LA PERIPHERIE DE L'ETOILE
 110     IF( NA1 .GT. 0 ) THEN
C           LE NO DU TRIANGLE ET LOCAL DE L'ARETE
            NT   = NAETOI(1,NA1)
            I    = ABS(NAETOI(2,NA1))
C           LE NUMERO DU TRIANGLE AU DELA DE L'ARETE
            NTOP = NOTRIA(3+I,NT)
            NAETOI(1,NA1) = NTOP
C           LE NUMERO DES 2 SOMMETS DE L'ARETE I DU TRIANGLE NT
            IF( I .EQ. 3 ) THEN
               NS4 = 1
            ELSE
               NS4 = I + 1
            ENDIF
C           NUMERO DU SOMMET 1 DE L'ARETE
            NS3 = NOTRIA(I,NT)
            NAETOI(2,NA1) = NS3
C           NUMERO DU TRIANGLE CONTENANT ENCORE NS3
            NOTRSO( NS3 ) = NTOP
C           NUMERO DU SOMMET 2 DE L'ARETE
            NS4 = NOTRIA(NS4,NT)
            NAETOI(3,NA1) = NS4
C           NUMERO DU TRIANGLE CONTENANT ENCORE NS4
            NOTRSO( NS4 ) = NTOP
C           PASSAGE A L'ARETE SUIVANTE
            NA1 = NAETOI(4,NA1)
            GOTO 110
         ENDIF
C
C        LES ARETES SONT REORDONNEES POUR FORMER UNE LIGNE FERMEE
C        --------------------------------------------------------
         NA1 = N1AEOC
C        LA PREMIERE ARETE
         NS0 = NAETOI(2,NA1)
         NS1 = NAETOI(3,NA1)
C
C        LE 1-ER SOMMET OU ARETE DU CONTOUR FERME DE LA SURFACE NSU
C        LE NOMBRE DE SOMMETS DU CONTOUR FERME DE L'ETOILE
         NBARCF = NBARCF + 1
         NBARC0 = NBARCF
         N1ARCF( NSU ) = NBARCF
C        LE PREMIER SOMMET DE L'ETOILE
         NOARCF( 1, NBARCF ) = NS0
C        LE SOMMET SUIVANT
         NOARCF( 2, NBARCF ) = NBARCF + 1
C        LE NUMERO DU TRIANGLE DE L'AUTRE COTE DE CETTE ARETE
         NOARCF( 3, NBARCF ) = NAETOI(1,NA1)
C
C        TRACE DE L'ARETE
         IF( TRATRI ) THEN
            X0 = REAL( PXYD(1,NS0) )
            Y0 = REAL( PXYD(2,NS0) )
            X  = REAL( PXYD(1,NS1) )
            Y  = REAL( PXYD(2,NS1) )
            CALL TRAIT2D(  NCVERT, X0, Y0, X, Y )
CCC            CALL ENTIER2D( NCVERT, X0, Y0, NS0 )
CCC            CALL ENTIER2D( NCVERT,  X,  Y, NS1 )
         ENDIF
C
C        L'ARETE SUIVANTE
         NA1    = NAETOI(4,NA1)
         N1AEOC = NA1
C
 120     IF( N1AEOC .GT. 0 ) THEN
C
C           RECHERCHE DE L'ARETE DE 1-ER SOMMET NS1
            NA0 = N1AEOC
            NA1 = NA0
 160        IF( NA1 .GT. 0 ) THEN
C
C              LE NUMERO DU PREMIER SOMMET DE L'ARETE
               IF ( NS1 .NE. NAETOI(2,NA1) .AND.
     %              NS1 .NE. NAETOI(3,NA1) ) THEN
C                 PASSAGE A L'ARETE SUIVANTE
                  NA0 = NA1
                  NA1 = NAETOI(4,NA1)
                  GOTO 160
               ENDIF
C
C              ARETE PERIPHERIQUE RETROUVEE
               NBARCF = NBARCF + 1
C              LE NUMERO DES 2 SOMMETS DE L'ARETE
               IF( NS1 .EQ. NAETOI(2,NA1) ) THEN
                  NS1 = NAETOI(3,NA1)
                  NS2 = 2
               ELSE
                  NS1 = NAETOI(2,NA1)
                  NS2 = 3
               ENDIF
C              LE NUMERO DU SOMMET DE L'ARETE
               NOARCF( 1, NBARCF ) = NAETOI(NS2,NA1)
C              L'ARETE SUIVANTE
               NOARCF( 2, NBARCF ) = NBARCF + 1
C              LE NUMERO DU TRIANGLE DE L'AUTRE COTE
               NOARCF( 3, NBARCF ) = NAETOI(1,NA1)
C
C              TRACE DE L'ARETE
               IF( TRATRI ) THEN
                  NS2 = NAETOI(2,NA1)
                  X0  = REAL( PXYD(1,NS2) )
                  Y0  = REAL( PXYD(2,NS2) )
                  NS2 = NAETOI(3,NA1)
                  X   = REAL( PXYD(1,NS2) )
                  Y   = REAL( PXYD(2,NS2) )
                  CALL TRAIT2D(  NCVERT, X0, Y0, X, Y )
                  CALL ENTIER2D( NCVERT, X0, Y0, NAETOI(2,NA1) )
                  CALL ENTIER2D( NCVERT,  X,  Y, NS2 )
               ENDIF
C
C              SUPPRESSION DE L'ARETE DANS LE CHAINAGE SUR L'ARETE SUIVANTE
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
C           ARETE NON RETROUVEE : L'ETOILE NE SE REFERME PAS
            WRITE(IMPRIM,*) 'ADSUPO: ANOMALIE ADSUPO A CORRIGER'
            WRITE(IMPRIM,*) 'ADSUPO: CF NON REFERME'
            CALL XVPAUSE
         ENDIF
C
C        LE SOMMET SUIVANT DU DERNIER SOMMET EST LE PREMIER
C        CHAINAGE CIRCULAIRE
         NOARCF( 2, NBARCF ) = NBARC0
C
C        MISE A JOUR DU NUMERO DU TRIANGLE QUI PRECEDE
         N1TR = NBTRS1
 1000 CONTINUE
C
C     CHAINAGE DES ARETES DU TABLEAU ARCF VIDES ( 2 SONT NECESSAIRES ENSUITE )
C     -----------------------------------------
      N1ARCF(0) = NBARCF+3
      DO 1050 I=NBARCF+3,MXETRI
         NOARCF(2,I) = I+1
 1050 CONTINUE
      NOARCF(2,MXETRI) = 0
C
C     SI LE SOMMET NS EST UN SOMMET DU CF IL DOIT DISPARAITRE
C     => IL FAUT TRANSFORMER LES 2 ARETES DE SOMMET NS
C        EN UNE SEULE ARETE ET CELA DANS CHACUN DES CF
C     -------------------------------------------------------
      DO 1200 NSU=1,NBSURF
C
C        RECHERCHE DU SOMMET NS DANS LES ARETES DU CF POUR LE SUPPRIMER
C        LA PREMIERE ARETE DU CF
         NA00 = N1ARCF(NSU)
         NA0  = NA00
C        L'ARETE SUIVANTE
         NA1 = NOARCF(2,NA0)
 1150    IF ( NS .NE. NOARCF(1,NA1) ) THEN
C           PASSAGE A L'ARETE SUIVANTE
            NA0 = NA1
            NA1 = NOARCF(2,NA1)
            IF( NA1 .NE. NA00 ) GOTO 1150
C           LE SOMMET NS NE FAIT PAS PARTIE DES ARETES DU CF
            GOTO 1210
         ENDIF
C
C        LE SOMMET NS EST COMMUN AUX ARETES NA0 ET NA1
C        L'ARETE NA1 EST SUPPRIMEE
         IF( N1ARCF(NSU) .EQ. NA1 ) N1ARCF(NSU) = NA0
         NOARCF(2,NA0) = NOARCF(2,NA1)
C        LE TRIANGLE OPPOSE VA DISPARAITRE DONC INCONNU (CODE 0)
         NOARCF(3,NA0) = 0
C        L'ARETE EST RENDUE AUX ARETES VIDES
         NOARCF(2,NA1) = N1ARCF(0)
         N1ARCF(0) = NA1
C
 1200 CONTINUE
C
C     DESTRUCTION DES TRIANGLES DES CF
C     --------------------------------
 1210 DO 1250 I=1,NBTRCF
         NT0 = NOTRCF( I )
C        MISE DE NT0 DANS LES TRIANGLES VIDES
         NOTRIA( 1, NT0 ) = 0
         NOTRIA( 4, NT0 ) = N1TRVI
         N1TRVI = NT0
 1250 CONTINUE
C
C     TRIANGULATION FRONTALE DES CF DE L'ETOILE DU SOMMET NS
C     ======================================================
      IF( NBSURF .EQ. 2 ) THEN
C        SAUVEGARDE DU CF DU SECOND CF
         NA2 = N1ARCF(2)
      ENDIF
C
C     TRIANGULATION FRONTALE DU CF NUS
      IERR = 0
      NUS  = 1
 1300 CALL TRIACF( 1,      PXYD  , N1TRVI, NOTRIA, NOTRSO,
     %             MXETRI, NARMIN, N1ARCF, NOARCF,
     %             NBTRCF, NOTRCF, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     MISE A JOUR DU NUMERO DE SURFACE DES TRIANGLES
      DO 1310 I=1,NBTRCF
C        LE NUMERO DE SURFACE EST MIS A JOUR
         IF( LESUTR .NE. 0 ) NOSUTR( NOTRCF( I ) ) = NOSUF( NUS )
C        LE TRACE DU TRIANGLE
         CALL DVTRTR( PXYD, NOTRIA, NOTRCF(I), NCVERT, NCBLAN )
 1310 CONTINUE
C
      IF( NBSURF .EQ. 2 ) THEN
C        TRIANGULATION DU SECOND CF
         N1ARCF(1) = NA2
         NUS = 2
         GOTO 1300
      ENDIF
C
C     LE SOMMET NS NE FAIT PLUS PARTIE DE LA TRIANGULATION
      NOTRSO( NS ) = 0

      RETURN
      END
