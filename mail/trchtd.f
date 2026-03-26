      SUBROUTINE TRCHTD( PXYD,   NAR00, NAR0,  NOARCF,
     %                   NAMIN0, NAMIN, LARMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHE DANS LE CONTOUR FERME DU SOMMET QUI JOINT A LA PLUS
C -----    COURTE ARETE NAR00 DONNE LE TRIANGLE SANS INTERSECTION
C          AVEC LE CONTOUR FERME DE MEILLEURE QUALITE
C
C ENTREES:
C --------
C PXYD   : TABLEAU DES COORDONNEES DES SOMMETS ET DISTANCE_SOUHAITEE
C
C ENTREES ET SORTIES:
C -------------------
C NAR00  : NUMERO DANS NOARCF DE L'ARETE AVANT NAR0
C NAR0   : NUMERO DANS NOARCF DE LA PLUS PETITE ARETE DU CONTOUR FERME
C          A JOINDRE A NOARCF(1,NAMIN) POUR FORMER LE TRIANGLE IDEAL
C NOARCF : NUMERO DU SOMMET , NUMERO DE L'ARETE SUIVANTE,
C          NUMERO DE L'ARETE DANS LE TABLEAU NOSOAR
C
C SORTIE :
C --------
C NAMIN0 : NUMERO DANS NOARCF DE L'ARETE AVANT NAMIN
C NAMIN  : NUMERO DANS NOARCF DU SOMMET CHOISI
C          0 SI CONTOUR FERME REDUIT A MOINS DE 3 ARETES
C            OU PAS DE SOMMET DU CF POUR FORMER UN TRIANGLE AVEC LA
C            PLUS PETITE ARETE DU CF
C LARMIN : TABLEAU AUXILIAIRE POUR STOCKER LA LISTE DES NUMEROS DES
C          ARETES DE MEILLEURE QUALITE POUR FAIRE LE CHOIX FINAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE    UPMC PARIS  FEVRIER 1992
C MODIFS : ALAIN PERRONNET LABORATOIRE JL LIONS UPMC PARIS  OCTOBRE 2006
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  PXYD(1:3,1:*)
      INTEGER           NOARCF(1:3,1:*),
     %                  LARMIN(1:*)
      DOUBLE PRECISION  DINFO, DMAXIM, DD, DMIMA,
     %                  UNPEPS, RAYON, SURTD2
      LOGICAL           OUI
      DOUBLE PRECISION  CENTRE(3)
      INTEGER           NOSOTR(3)
      EQUIVALENCE      (NOSOTR(1),NS1),(NOSOTR(2),NS2),(NOSOTR(3),NS3)
C
C     INITIALISATIONS
C     DMAXIM : LE PLUS GRAND REEL MACHINE
      DMAXIM = DINFO( 'GRAND' )
      UNPEPS = 1D0 + 100D0 * DINFO( 'PRECISION' )
C     SUR DEC-ALPHA LA PRECISION EST DE 10**-14 SEULEMENT
C
C     RECHERCHE DE LA PLUS COURTE ARETE DU CONTOUR FERME
      NBMIN = 0
      NA00  = NAR00
      DMIMA = DMAXIM
      NBAR  = 0
C
 2    NA0  = NOARCF( 2, NA00 )
      NA1  = NOARCF( 2, NA0  )
      NBAR = NBAR + 1
C     LES 2 SOMMETS DE L'ARETE NA0 DU CF
      NS1  = NOARCF( 1, NA0 )
      NS2  = NOARCF( 1, NA1 )
      DD   = (PXYD(1,NS2)-PXYD(1,NS1))**2 + (PXYD(2,NS2)-PXYD(2,NS1))**2
      IF( DD .LT. DMIMA ) THEN
         DMIMA = DD
         LARMIN(1) = NA00
      ENDIF
      NA00 = NA0
      IF( NA00 .NE. NAR00 ) THEN
C        DERNIERE ARETE NON ATTEINTE
         GOTO 2
      ENDIF
C
      IF( NBAR .EQ. 3 ) THEN
C
C        CONTOUR FERME REDUIT A UN TRIANGLE
C        ----------------------------------
         NAMIN  = NAR00
         NAR0   = NOARCF( 2, NAR00 )
         NAMIN0 = NOARCF( 2, NAR0  )
         RETURN
C
      ELSE IF( NBAR .LE. 2 ) THEN
         WRITE(IMPRIM,*) 'ERREUR TRCHTD: CF<3 ARETES'
         CALL XVPAUSE
         NAMIN  = 0
         NAMIN0 = 0
         RETURN
      ENDIF
C
C     CF NON REDUIT A UN TRIANGLE
C     LA PLUS PETITE ARETE EST NAR0 DANS NOARCF
      NAR00 = LARMIN( 1 )
      NAR0  = NOARCF( 2, NAR00 )
      NAR   = NOARCF( 2, NAR0  )
C
      NS1   = NOARCF( 1, NAR0 )
      NS2   = NOARCF( 1, NAR  )
C
C     RECHERCHE DANS CETTE ETOILE DU SOMMET OFFRANT LA MEILLEURE QUALITE
C     DU TRIANGLE NS1-NS2 NS3 SANS INTERSECTION AVEC LE CONTOUR FERME
C     ==================================================================
      NAR3  = NAR
      QMIMA = -1
C
C     PARCOURS DES SOMMETS POSSIBLES NS3
 10   NAR3  = NOARCF( 2, NAR3 )
      IF( NAR3 .NE. NAR0 ) THEN
C
C        IL EXISTE UN SOMMET NS3 DIFFERENT DE NS1 ET NS2
         NS3 = NOARCF( 1, NAR3 )
C
C        LES ARETES NS1-NS3 et NS2-NS3 INTERSECTENT-ELLES UNE ARETE DU CONTOUR F
C        -----------------------------------------------------------------------
C        INTERSECTION DE L'ARETE NS2-NS3 ET DES ARETES DU CF JUSQU'AU SOMMET NS3
         NAR1 = NOARCF( 2, NAR )
 12      IF( NAR1 .NE. NAR3 .AND. NOARCF( 2, NAR1 ) .NE. NAR3 ) THEN
C           L'ARETE SUIVANTE
            NAR2 = NOARCF( 2, NAR1 )
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NP1  = NOARCF( 1, NAR1 )
            NP2  = NOARCF( 1, NAR2 )
            CALL INT2AR( PXYD(1,NS2), PXYD(1,NS3),
     %                   PXYD(1,NP1), PXYD(1,NP2), OUI )
            IF( OUI ) GOTO 10
C           LES 2 ARETES NE S'INTERSECTENT PAS ENTRE LEURS SOMMETS
            NAR1 = NAR2
            GOTO 12
         ENDIF
C
C        INTERSECTION DE L'ARETE NS1-NS3 ET DES ARETES DU CF JUSQU'AU SOMMET NS3
         NAR1 = NAR
 14      IF( NAR1 .NE. NAR3.AND. NOARCF( 2, NAR1 ) .NE. NAR3 ) THEN
C           L'ARETE SUIVANTE
            NAR2 = NOARCF( 2, NAR1 )
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NP1  = NOARCF( 1, NAR1 )
            NP2  = NOARCF( 1, NAR2 )
            CALL INT2AR( PXYD(1,NS1), PXYD(1,NS3),
     %                   PXYD(1,NP1), PXYD(1,NP2), OUI )
            IF( OUI ) GOTO 10
C           LES 2 ARETES NE S'INTERSECTENT PAS ENTRE LEURS SOMMETS
            NAR1 = NAR2
            GOTO 14
         ENDIF
C
C        INTERSECTION DE L'ARETE NS1-NS3 ET DES ARETES DU CF JUSQU'AU SOMMET D'A
         NAR1 = NOARCF( 2, NAR3 )
 16      IF( NAR1 .NE. NAR0 .AND. NOARCF( 2, NAR1 ) .NE. NAR0 ) THEN
C           L'ARETE SUIVANTE
            NAR2 = NOARCF( 2, NAR1 )
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NP1  = NOARCF( 1, NAR1 )
            NP2  = NOARCF( 1, NAR2 )
            CALL INT2AR( PXYD(1,NS1), PXYD(1,NS3),
     %                   PXYD(1,NP1), PXYD(1,NP2), OUI )
            IF( OUI ) GOTO 10
C           LES 2 ARETES NE S'INTERSECTENT PAS ENTRE LEURS SOMMETS
            NAR1 = NAR2
            GOTO 16
         ENDIF
C
C        INTERSECTION DE L'ARETE NS2-NS3 ET DES ARETES DU CF JUSQU'AU SOMMET D'A
         NAR1 = NOARCF( 2, NAR3 )
 18      IF( NAR1 .NE. NAR0 .AND. NOARCF( 2, NAR1 ) .NE. NAR0) THEN
C           L'ARETE SUIVANTE
            NAR2 = NOARCF( 2, NAR1 )
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NP1  = NOARCF( 1, NAR1 )
            NP2  = NOARCF( 1, NAR2 )
            CALL INT2AR( PXYD(1,NS2), PXYD(1,NS3),
     %                   PXYD(1,NP1), PXYD(1,NP2), OUI )
            IF( OUI ) GOTO 10
C           LES 2 ARETES NE S'INTERSECTENT PAS ENTRE LEURS SOMMETS
            NAR1 = NAR2
            GOTO 18
         ENDIF
C
C        ICI  LE TRIANGLE NS1-NS2-NS3 N'INTERSECTE PAS UNE ARETE DU CONTOUR FERM
C        MAIS LE TRIANGLE NS1-NS2-NS3 CONTIENT IL UN SOMMET DU CONTOUR FERME?
C        -----------------------------------------------------------------------
C        SOMMETS DES ARETES DU CF AU DELA DE NS2 JUSQU'AU SOMMET NS3
         NAR1 = NOARCF( 2, NAR )
 22      IF( NAR1 .NE. NAR3 ) THEN
C           LE NUMERO DU SOMMET DE L'ARETE NAR1
            NP1 = NOARCF( 1, NAR1 )
            CALL  PTDATR( PXYD(1,NP1), PXYD, NOSOTR, NSIGNE )
C           NSIGNE : >0 SI LE POINT NP1 EST DANS LE TRIANGLE OU SUR UNE DES 3 AR
C                    =0 SI LE TRIANGLE EST DEGENERE OU INDIRECT OU NE CONTIENT
C                       PAS LE POINT
            IF( NSIGNE .GT. 0 ) GOTO 10
C           L'ARETE SUIVANTE
            NAR1 = NOARCF( 2, NAR1 )
            GOTO 22
         ENDIF
C
C        SOMMETS DES ARETES DU CF AU DELA DE NS3 JUSQU'AU SOMMET NS1
         NAR1 = NOARCF( 2, NAR3 )
 24      IF( NAR1 .NE. NAR0 ) THEN
C           LE NUMERO DU SOMMET DE L'ARETE NAR1
            NP1 = NOARCF( 1, NAR1 )
            CALL  PTDATR( PXYD(1,NP1), PXYD, NOSOTR, NSIGNE )
C           NSIGNE : >0 SI LE POINT NP1 EST DANS LE TRIANGLE OU SUR UNE DES 3 AR
C                    =0 SI LE TRIANGLE EST DEGENERE OU INDIRECT OU NE CONTIENT
C                       PAS LE POINT
            IF( NSIGNE .GT. 0 ) GOTO 10
C           L'ARETE SUIVANTE
            NAR1 = NOARCF( 2, NAR1 )
            GOTO 24
         ENDIF
C
C        ICI LE TRIANGLE NS1-NS2-NS3 N'INTERSECTE PAS UNE ARETE DU CONTOUR FERME
C        ET  LE TRIANGLE NS1-NS2-NS3 NE CONTIENT PAS DE SOMMETS DU CONTOUR FERME
C        -----------------------------------------------------------------------
C
C        LE CALCUL DE LA SURFACE DU TRIANGLE NS1-NS2-NS3
         DD = SURTD2( PXYD(1,NS1), PXYD(1,NS2), PXYD(1,NS3) )
         IF( DD .LE. 0D0 ) THEN
C           SURFACE NEGATIVE => TRIANGLE A REJETER
            Q = 0
         ELSE
C           CALCUL DE LA QUALITE DU  TRIANGLE  NS1-NS2-NS3
            CALL QUTR2D( PXYD(1,NS1), PXYD(1,NS2), PXYD(1,NS3), Q )
         ENDIF
C
         IF( Q .GE. QMIMA*1.00001 ) THEN
C           Q EST UN VRAI MAXIMUM DE LA QUALITE
            QMIMA = Q
            NBMIN = 1
            LARMIN(1) = NAR3
         ELSE IF( Q .GE. QMIMA*0.999998 ) THEN
C           Q EST VOISIN DE QMIMA
C           IL EST EMPILE
            NBMIN = NBMIN + 1
            LARMIN( NBMIN ) = NAR3
         ENDIF
         GOTO 10
      ENDIF
C
C     BILAN
C     =====
      IF( NBMIN .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'TRCHTD: PAS DE SOMMET DU CF POUR FORMER UN TRI
     %ANGLE AVEC L''ARETE',NS1,NS2
         NAMIN = 0
         RETURN
      ENDIF
C
C     EXISTE T IL PLUSIEURS SOMMETS DE MEME QUALITE?
C     ===============================================
      IF( NBMIN .GT. 1 ) THEN
C
C        OUI:RECHERCHE DE CEUX DE CERCLE NE CONTENANT PAS D'AUTRES SOMMETS
         DO 80 I=1,NBMIN
C           LE SOMMET
            NAR = LARMIN( I )
            IF( NAR .LE. 0 ) GOTO 80
            NS3 = NOARCF(1,NAR)
C           LES COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT
C           ET SON RAYON
            IER = -1
            CALL CENCED( PXYD(1,NS1), PXYD(1,NS2), PXYD(1,NS3),
     %                   CENTRE, IER )
            IF( IER .NE. 0 ) THEN
C              LE SOMMET NS3 NE CONVIENT PAS
               LARMIN( I ) = 0
               GOTO 80
            ENDIF
            RAYON = CENTRE(3) * UNPEPS
            DO 70 J=1,NBMIN
               IF( J .NE. I ) THEN
C                 L'AUTRE SOMMET
                  NAR1 = LARMIN(J)
                  IF( NAR1 .LE. 0 ) GOTO 70
                  NS4 = NOARCF(1,NAR1)
C                 APPARTIENT T IL AU CERCLE NS1 NS2 NS3 ?
                  DD = (CENTRE(1)-PXYD(1,NS4))**2 +
     %                 (CENTRE(2)-PXYD(2,NS4))**2
                  IF( DD .LE. RAYON ) THEN
C                    NS4 EST DANS LE CERCLE CIRCONSCRIT  NS1 NS2 NS3
C                    LE SOMMET NS3 NE CONVIENT PAS
                     LARMIN( I ) = 0
                     GOTO 80
                  ENDIF
               ENDIF
 70         CONTINUE
 80      CONTINUE
C
C        EXISTE T IL PLUSIEURS SOMMETS ?
         J = 0
         DO 90 I=1,NBMIN
            IF( LARMIN( I ) .GT. 0 ) THEN
C              COMPACTAGE DES MIN
               J = J + 1
               LARMIN(J) = LARMIN(I)
            ENDIF
 90      CONTINUE
C
         IF( J .GT. 1 ) THEN
C           OUI : CHOIX DU PLUS PETIT RAYON DE CERCLE CIRCONSCRIT
            DMIMA = DMAXIM
            DO 120 I=1,NBMIN
               NS3 = NOARCF(1,LARMIN(I))
C
C              LES COORDONNEES DU CENTRE DE CERCLE CIRCONSCRIT
C              AU TRIANGLE NT ET SON RAYON
               IER = -1
               CALL CENCED( PXYD(1,NS1), PXYD(1,NS2), PXYD(1,NS3),
     %                      CENTRE, IER )
               IF( IER .NE. 0 ) THEN
C                 LE SOMMET NS3 NE CONVIENT PAS
                  GOTO 120
               ENDIF
               RAYON = SQRT( CENTRE(3) )
               IF( RAYON .LT. DMIMA ) THEN
                  DMIMA = RAYON
                  LARMIN(1) = LARMIN(I)
               ENDIF
 120        CONTINUE
         ENDIF
      ENDIF
C
C     LE CHOIX FINAL
C     ==============
      NAMIN = LARMIN(1)
C
C     RECHERCHE DE L'ARETE AVANT NAMIN ( NAR0 <> NAMIN )
C     ==================================================
      NAR1 = NAR0
 200  IF( NAR1 .NE. NAMIN ) THEN
         NAMIN0 = NAR1
         NAR1   = NOARCF( 2, NAR1 )
         GOTO 200
      ENDIF
      END
