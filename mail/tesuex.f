      SUBROUTINE TESUEX( NBLFTR, NULFTR,
     %                   NDTRI0, NBSOMM, PXYD, NSLIGN,
     %                   MOSOAR, MXSOAR, NOSOAR,
     %                   MOARTR, MXARTR, N1ARTR, NOARTR, NOARST,
     %                   NBTRIA, LETRSU, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER DU TABLEAU NOARTR LES TRIANGLES EXTERNES AU DOMAINE
C -----    EN ANNULANT LE NUMERO DE LEUR 1-ERE ARETE DANS NOARTR
C          ET EN LES CHAINANT COMME TRIANGLES VIDES
C
C ENTREES:
C --------
C NBLFTR : NOMBRE DE  LIGNES FERMEES DEFINISSANT LA SURFACE
C NULFTR : NUMERO DES LIGNES FERMEES DEFINISSANT LA SURFACE
C NDTRI0 : PLUS GRAND NUMERO DANS NOARTR D'UN TRIANGLE
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C NSLIGN : TABLEAU DU NUMERO DE SOMMET DANS SA LIGNE POUR CHAQUE
C          SOMMET FRONTALIER
C          NUMERO DU POINT DANS LE LEXIQUE POINT SI INTERNE IMPOSE
C          0 SI LE POINT EST INTERNE NON IMPOSE PAR L'UTILISATEUR
C         -1 SI LE SOMMET EST EXTERNE AU DOMAINE
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE
C MXSOAR : NOMBRE MAXIMAL D'ARETES STOCKABLES DANS LE TABLEAU NOSOAR
C          ATTENTION: MXSOAR>3*MXSOMM OBLIGATOIRE!
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          HACHAGE DES ARETES = NOSOAR(1)+NOSOAR(2)*2
C          AVEC MXSOAR>=3*MXSOMM
C          UNE ARETE I DE NOSOAR EST VIDE <=> NOSOAR(1,I)=0 ET
C          NOSOAR(2,ARETE VIDE)=L'ARETE VIDE QUI PRECEDE
C          NOSOAR(3,ARETE VIDE)=L'ARETE VIDE QUI SUIT
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C MXARTR : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES
C N1ARTR : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOARTR
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOARTR(2,.)
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1 = 0 SI TRIANGLE VIDE => ARETE2 = TRIANGLE VIDE SUIVANT
C NOARST : NOARST(I) NUMERO NOSOAR D'UNE ARETE DE SOMMET I
C
C SORTIES:
C --------
C NBTRIA : NOMBRE DE TRIANGLES INTERNES AU DOMAINE
C LETRSU : LETRSU(NT)=NUMERO DU TRIANGLE INTERNE, 0 SINON
C NOARST : NOARST(I) NUMERO NOSOAR D'UNE ARETE DU SOMMET I (MODIFI'E)
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
CC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1999
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  JANVIER 2008
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NULFTR(NBLFTR),NSLIGN(NBSOMM),
     %                  NOSOAR(MOSOAR,MXSOAR),
     %                  NOARTR(MOARTR,MXARTR),
     %                  NOARST(*)
      INTEGER           LETRSU(1:NDTRI0)
      DOUBLE PRECISION  DMIN
C
C     LES TRIANGLES SONT A PRIORI NON MARQUES
      DO 5 NT=1,NDTRI0
         LETRSU(NT) = 0
 5    CONTINUE
C
C     LES ARETES SONT MARQUEES NON CHAINEES
      DO 10 NOAR1=1,MXSOAR
         NOSOAR(6,NOAR1) = -2
 10   CONTINUE
C
C     RECHERCHE DU SOMMET DE LA TRIANGULATION DE PLUS PETITE ABSCISSE
C     ===============================================================
      NTMIN = 0
      DMIN  = 1D38
      DO 20 I=1,NBSOMM
         IF( PXYD(1,I) .LT. DMIN ) THEN
C           LE NOUVEAU MINIMUM
            NOAR1 = NOARST(I)
            IF( NOAR1 .GT. 0 ) THEN
C              LE SOMMET APPARTIENT A UNE ARETE DE TRIANGLE
               IF( NOSOAR(4,NOAR1) .GT. 0 ) THEN
C                 LE NOUVEAU MINIMUM
                  DMIN  = PXYD(1,I)
                  NTMIN = I
               ENDIF
            ENDIF
         ENDIF
 20   CONTINUE
C
C     UNE ARETE DE SOMMET NTMIN
      NOAR1 = NOARST( NTMIN )
C     UN TRIANGLE D'ARETE NOAR1
      NTMIN = NOSOAR( 4, NOAR1 )
      IF( NTMIN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PAS DE TRIANGLE D''ABSCISSE MINIMALE'
         ELSE
            KERR(1) = 'NO TRIANGLE WTH MINIMUM ABSCISSAE'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9990
      ENDIF
C
C     29 janvier 2008  Correction d'une erreur sur chainage des 3 aretes de ntmi
C     ..........................................................................
C     DEPART AVEC UNE ARETE DU TRIANGLE NTMIN NON SUR UNE LIGNE FERMEE
C     ================================================================
      DO 25 I=1,3
         NOAR1 = ABS( NOARTR(I,NTMIN) )
C        LE NUMERO DE 1 A NBLFTR DANS NULFTR DE LA LIGNE DE L'ARETE
         IF( NOSOAR(3,NOAR1) .LE. 0 ) THEN
C           CETTE ARETE N'A PAS DE SUIVANTE POUR L'INSTANT
            NOSOAR(6,NOAR1) = 0
            GOTO 30
         ENDIF
 25   CONTINUE
      NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRIANGLE D''ABSCISSE MINIMALE AVEC SES 3 ARETES'
         KERR(2) = 'SUR LES LIGNES FERMEES DU CONTOUR'
      ELSE
         KERR(1) = 'The TRIANGLE WTH MINIMUM ABSCISSAE HAS ITS 3 EDGES'
         KERR(2) = 'ON THE CLOSED LINES of the CONTOUR'
      ENDIF
      CALL LEREUR
      IERR = 3
      GOTO 9990
C     ..........................................................................
C
C     LE TRIANGLE NTMIN EST A L'EXTERIEUR DU DOMAINE
C     TOUS LES TRIANGLES EXTERNES SONT MARQUES -123 456 789
C     LES TRIANGLES DE L'AUTRE COTE D'UNE ARETE SUR UNE LIGNE
C     SONT MARQUES: NO DE LA LIGNE DE L'ARETE * SIGNE OPPOSE
C     =======================================================
 30   LIGNE0 = 0
      LIGNE  = -123 456 789
C
 40   IF( NOAR1 .NE. 0 ) THEN
C
C        L'ARETE NOAR1 DU TABLEAU NOSOAR EST A TRAITER
C        ---------------------------------------------
         NOAR = NOAR1
C        L'ARETE SUIVANTE DEVIENT LA PREMIERE A TRAITER ENSUITE
         NOAR1 = NOSOAR(6,NOAR1)
C        L'ARETE NOAR EST TRAITEE
         NOSOAR(6,NOAR) = -3
C
         DO 60 I=4,5
C
C           L'UN DES 2 TRIANGLES DE L'ARETE
            NT = NOSOAR(I,NOAR)
            IF( NT .GT. 0 ) THEN
C
C              TRIANGLE DEJA TRAITE POUR UNE LIGNE ANTERIEURE?
               IF(     LETRSU(NT)  .NE. 0      .AND.
     %             ABS(LETRSU(NT)) .NE. LIGNE ) GOTO 60
C
C              TRACE DU TRIANGLE NT EN COULEUR LIGNE0
               CALL MTTRTR( PXYD,   NT, MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                      LIGNE0, NCNOIR )
C
C              LE TRIANGLE EST MARQUE AVEC LA VALEUR DE LIGNE
               LETRSU(NT) = LIGNE
C
C              CHAINAGE EVENTUEL DES AUTRES ARETES DE CE TRIANGLE
C              SI CE N'EST PAS ENCORE FAIT
               DO 50 J=1,3
C
C                 LE NUMERO NA DE L'ARETE J DU TRIANGLE NT DANS NOSOAR
                  NA = ABS( NOARTR(J,NT) )
                  IF( NOSOAR(6,NA) .NE. -2 ) GOTO 50
C
C                 LE NUMERO DE 1 A NBLFTR DANS NULFTR DE LA LIGNE DE L'ARETE
                  NL = NOSOAR(3,NA)
C
C                 SI L'ARETE EST SUR UNE LIGNE FERMEE DIFFERENTE DE CELLE ENVELO
C                 ET NON MARQUEE ALORS EXAMEN DU TRIANGLE OPPOSE
                  IF( NL .GT. 0 ) THEN
C
                     IF( NL .EQ. LIGNE0 ) GOTO 50
C
C                    ARETE FRONTALIERE DE LIGNE NON TRAITEE
C                    => PASSAGE DE L'AUTRE COTE DE LA LIGNE
C                    LE TRIANGLE DE L'AUTRE COTE DE LA LIGNE EST RECHERCHE
                     IF( NT .EQ. ABS( NOSOAR(4,NA) ) ) THEN
                        NT2 = 5
                     ELSE
                        NT2 = 4
                     ENDIF
                     NT2 = ABS( NOSOAR(NT2,NA) )
                     IF( NT2 .GT. 0 ) THEN
C
C                       LE TRIANGLE NT2 DE L'AUTRE COTE EST MARQUE
C                       AVEC LE SIGNE OPPOSE DE CELUI DE LIGNE
                        IF( LIGNE .GE. 0 ) THEN
                           LSIGNE = -1
                        ELSE
                           LSIGNE =  1
                        ENDIF
                        LETRSU(NT2) = LSIGNE * NL
C
C                       TEMOIN DE LIGNE A TRAITER ENSUITE DANS NULFTR
                        NULFTR(NL) = -ABS( NULFTR(NL) )
C
C                       TRACE DU TRIANGLE NT2 EN JAUNE BORDE DE MAGENTA
                        CALL MTTRTR( PXYD,NT2,
     %                               MOARTR,NOARTR,MOSOAR,NOSOAR,
     %                               NCJAUN, NCMAGE )
C
C                       L'ARETE EST TRAITEE
                        NOSOAR(6,NA) = -3
C
                     ENDIF
C
C                    L'ARETE EST TRAITEE
                     GOTO 50
C
                  ENDIF
C
C                 ARETE NON TRAITEE => ELLE EST CHAINEE
                  NOSOAR(6,NA) = NOAR1
                  NOAR1 = NA
C
 50            CONTINUE
C
            ENDIF
 60      CONTINUE
         GOTO 40
      ENDIF
C     LES TRIANGLES DE LA LIGNE FERMEE ONT TOUS ETE MARQUES
C     PLUS D'ARETE CHAINEE
C
C     RECHERCHE D'UNE NOUVELLE LIGNE FERMEE A TRAITER
C     ===============================================
 65   DO 70 NL=1,NBLFTR
         IF( NULFTR(NL) .LT. 0 ) GOTO 80
 70   CONTINUE
C     PLUS DE LIGNE FERMEE A TRAITER
      GOTO 110
C
C     TOUS LES TRIANGLES DE CETTE COMPOSANTE CONNEXE
C     ENTRE LIGNE ET LIGNE0 VONT ETRE MARQUES
C     ==============================================
C     REMISE EN ETAT DU NUMERO DE LIGNE
C     NL EST LE NUMERO DE LA LIGNE DANS NULFTR A TRAITER
 80   NULFTR(NL) = -NULFTR(NL)
      DO 90 NT2=1,NDTRI0
         IF( ABS(LETRSU(NT2)) .EQ. NL ) GOTO 92
 90   CONTINUE
C
C     RECHERCHE DE L'ARETE J DU TRIANGLE NT2 AVEC CE NUMERO DE LIGNE NL
 92   DO 95 J=1,3
C
C        LE NUMERO DE L'ARETE J DU TRIANGLE DANS NOSOAR
         NOAR1 = 0
         NA0   = ABS( NOARTR(J,NT2) )
         IF( NL .EQ. NOSOAR(3,NA0) ) THEN
C
C           NA0 EST L'ARETE DE LIGNE NL
C           L'ARETE SUIVANTE DU TRIANGLE NT2
            I   = MOD(J,3) + 1
C           LE NUMERO DANS NOSOAR DE L'ARETE I DE NT2
            NA1 = ABS( NOARTR(I,NT2) )
            IF( NOSOAR(6,NA1) .EQ. -2 ) THEN
C              ARETE NON TRAITEE => ELLE EST LA PREMIERE DU CHAINAGE
               NOAR1 = NA1
C              PAS DE SUIVANTE DANS CE CHAINAGE
               NOSOAR(6,NA1) = 0
            ELSE
               NA1 = 0
            ENDIF
C
C           L'EVENTUELLE SECONDE ARETE SUIVANTE
            I  = MOD(I,3) + 1
            NA = ABS( NOARTR(I,NT2) )
            IF( NOSOAR(6,NA) .EQ. -2 ) THEN
               IF( NA1 .EQ. 0 ) THEN
C                 1 ARETE NON TRAITEE ET SEULE A CHAINER
                  NOAR1 = NA
                  NOSOAR(6,NA) = 0
               ELSE
C                 2 ARETES A CHAINER
                  NOAR1 = NA
                  NOSOAR(6,NA) = NA1
               ENDIF
            ENDIF
C
            IF( NOAR1 .GT. 0 ) THEN
C
C              IL EXISTE AU MOINS UNE ARETE A VISITER POUR LIGNE
C              MARQUAGE DES TRIANGLES INTERNES A LA LIGNE NL
               LIGNE  = LETRSU(NT2)
               LIGNE0 = NL
               GOTO 40
C
            ELSE
C
C              NT2 EST LE SEUL TRIANGLE DE LA LIGNE FERMEE
               GOTO 65
C
            ENDIF
         ENDIF
 95   CONTINUE
C
C     REPERAGE DES SOMMETS INTERNES OU EXTERNES DANS NSLIGN
C     NSLIGN(SOMMET EXTERNE AU DOMAINE)=-1
C     NSLIGN(SOMMET INTERNE AU DOMAINE)= 0
C     =====================================================
 110  DO 170 NS1=1,NBSOMM
C        TOUT SOMMET NON SUR LA FRONTIERE OU INTERNE IMPOSE
C        EST SUPPOSE EXTERNE
         IF( NSLIGN(NS1) .EQ. 0 ) NSLIGN(NS1) = -1
 170  CONTINUE
C
C     LES TRIANGLES EXTERNES SONT MARQUES VIDES DANS LE TABLEAU NOARTR
C     ================================================================
      NBTRIA = 0
      DO 200 NT=1,NDTRI0
C
         IF( LETRSU(NT) .LE. 0 ) THEN
C
C           TRIANGLE NT EXTERNE
            IF( NOARTR(1,NT) .NE. 0 ) THEN
C              LA PREMIERE ARETE EST ANNULEE
               NOARTR(1,NT) = 0
C              LE TRIANGLE NT EST CONSIDERE COMME ETANT VIDE
               NOARTR(2,NT) = N1ARTR
               N1ARTR = NT
            ENDIF
C
         ELSE
C
C           TRIANGLE NT INTERNE
            NBTRIA = NBTRIA + 1
            LETRSU(NT) = NBTRIA
C
C           MARQUAGE DES 3 SOMMETS DU TRIANGLE NT
            DO 190 I=1,3
C              LE NUMERO NOSOAR DE L'ARETE I DU TRIANGLE NT
               NOAR = ABS( NOARTR(I,NT) )
C              LE NUMERO DES 2 SOMMETS
               NS1 = NOSOAR(1,NOAR)
               NS2 = NOSOAR(2,NOAR)
C              MISE A JOUR DU NUMERO D'UNE ARETE DES 2 SOMMETS DE L'ARETE
               NOARST( NS1 ) = NOAR
               NOARST( NS2 ) = NOAR
C              NS1 ET NS2 SONT DES SOMMETS DE LA TRIANGULATION DU DOMAINE
               IF( NSLIGN(NS1) .LT. 0 ) NSLIGN(NS1)=0
               IF( NSLIGN(NS2) .LT. 0 ) NSLIGN(NS2)=0
 190        CONTINUE
C
         ENDIF
C
 200  CONTINUE
C     ICI TOUT SOMMET EXTERNE NS VERIFIE NSLIGN(NS)=-1
C
C     LES TRIANGLES EXTERNES SONT MIS A ZERO DANS NOSOAR
C     ==================================================
      DO 300 NOAR=1,MXSOAR
C
         IF( NOSOAR(1,NOAR) .GT. 0 ) THEN
C
C           LE SECOND TRIANGLE DE L'ARETE NOAR
            NT = NOSOAR(5,NOAR)
            IF( NT .GT. 0 ) THEN
C              SI LE TRIANGLE NT EST EXTERNE
C              ALORS IL EST SUPPRIME POUR L'ARETE NOAR
               IF( LETRSU(NT) .LE. 0 ) NOSOAR(5,NOAR)=0
            ENDIF
C
C           LE PREMIER TRIANGLE DE L'ARETE NOAR
            NT = NOSOAR(4,NOAR)
            IF( NT .GT. 0 ) THEN
               IF( LETRSU(NT) .LE. 0 ) THEN
C                 SI LE TRIANGLE NT EST EXTERNE
C                 ALORS IL EST SUPPRIME POUR L'ARETE NOAR
C                 ET L'EVENTUEL TRIANGLE OPPOSE PREND SA PLACE
C                 EN POSITION 4 DE NOSOAR
                  IF( NOSOAR(5,NOAR) .GT. 0 ) THEN
                     NOSOAR(4,NOAR)=NOSOAR(5,NOAR)
                     NOSOAR(5,NOAR)=0
                  ELSE
                     NOSOAR(4,NOAR)=0
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
 300  CONTINUE
C
C     REMISE EN ETAT POUR EVITER LES MODIFICATIONS DE LADEFI
 9990 DO 9991 NL=1,NBLFTR
         IF( NULFTR(NL) .LT. 0 ) NULFTR(NL)=-NULFTR(NL)
 9991 CONTINUE
      RETURN
      END
