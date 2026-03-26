      SUBROUTINE OTCHTR( PTXYZD, N1ARCF, NOARCF,
     %                   NAR00,  NAR0,   NAMIN0, NAMIN, NARMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHE DANS LE CONTOUR FERME DU SOMMET QUI JOINT A LA PLUS
C -----    COURTE ARETE NAR00 DONNE LE TRIANGLE SANS INTERSECTION
C          AVEC LE CONTOUR FERME DE MEILLEURE QUALITE
C
C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES SOMMETS ET DISTANCE_SOUHAITEE
C N1ARCF : POINTEUR DANS NOARCF DE LA PREMIERE ARETE DU CF
C NOARCF : NUMERO DES ARETES DE LA LIGNE DU CONTOUR FERME SELON UN SENS
C
C ENTREES ET SORTIES:
C -------------------
C NAR00  : NUMERO DANS NOARCF DE L'ARETE AVANT NAR0
C NAR0   : NUMERO DANS NOARCF DE LA PLUS PETITE ARETE DU CONTOUR FERME
C          A JOINDRE A NOARCF(1,NAMIN) POUR FORMER LE TRIANGLE IDEAL
C NOARCF : NUMERO DU SOMMET , NUMERO DE L'ARETE SUIVANTE
C          NUMERO DU TRIANGLE EXTERIEUR A L'ETOILE
C
C SORTIE :
C --------
C NAMIN0 : NUMERO DANS NOARCF DE L'ARETE AVANT NAMIN
C NAMIN  : NUMERO DANS NOARCF DU SOMMET CHOISI
C          0 SI CONTOUR FERME REDUIT A MOINS DE 3 ARETES
C NARMIN : TABLEAU AUXILIAIRE POUR STOCKER LES ARETES A EGALES DISTANCES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1993
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  PTXYZD(1:4,1:*)
      INTEGER           NOARCF(1:3,1:*),
     %                  NARMIN(1:*)
      DOUBLE PRECISION  RAYON,
     %                  CENTRE(3),
     %                  XYZNOR(3),
     %                  XY2D(2,3),
     %                  XY2DAR(2,2)
      DOUBLE PRECISION  DD,D,DMIMA,A,B,C,XP, SURTRD
      DOUBLE PRECISION  D2D3(3,3)
      LOGICAL           OUI
C
C     INITIALISATIONS
C     DMAXIM : LE PLUS GRAND REEL MACHINE
      DMAXIM = RINFO( 'GRAND' )
      UNPEPS = 1. + 5. * RINFO( 'PRECISION' )
C
      NAR00 = N1ARCF
      NAR0  = NOARCF( 2, NAR00 )
C
C     RECHERCHE DE LA PLUS COURTE ARETE DU CONTOUR FERME
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
      DD   = 0
      DO 5 I=1,3
         D  = PTXYZD(I,NS2) - PTXYZD(I,NS1)
         DD = DD + D * D
 5    CONTINUE
      IF( DD .LT. DMIMA ) THEN
         DMIMA = DD
         NARMIN(1) = NA00
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
         WRITE(IMPRIM,*) 'ERREUR OTCHTR:  CF<3 ARETES'
         CALL XVPAUSE
         NAMIN  = 0
         NAMIN0 = 0
         RETURN
C
      ELSE IF( NBAR .EQ. 4 ) THEN
C
C        CONTOUR FERME REDUIT A UN QUADRANGLE  => 2 TRIANGLES AU MIEUX
C        ------------------------------------
         NA1 = N1ARCF
         NA2 = NOARCF( 2, NA1 )
         NA3 = NOARCF( 2, NA2 )
         NA4 = NOARCF( 2, NA3 )
C
         NS1 = NOARCF( 1, NA1 )
         NS2 = NOARCF( 1, NA2 )
         NS3 = NOARCF( 1, NA3 )
         NS4 = NOARCF( 1, NA4 )
C        L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
         CALL ARTADM( NS1, NS3,
     %                N1ARCF, NOARCF, PTXYZD, NONOU3 )
C        LA QUALITE DES 2 TRIANGLES A FORMER
         CALL QUTRTE( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3), Q1 )
         CALL QUTRTE( PTXYZD(1,NS1), PTXYZD(1,NS3), PTXYZD(1,NS4), Q2 )
C
C        L'AUTRE DIAGONALE
C        L'ARETE AVEC LE SOMMET SUIVANT EST ELLE ADMISSIBLE?
         CALL ARTADM( NS2, NS4,
     %                N1ARCF, NOARCF, PTXYZD, NONOU4 )
C        LA QUALITE DES 2 TRIANGLES A FORMER
         CALL QUTRTE( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS4), Q3 )
         CALL QUTRTE( PTXYZD(1,NS2), PTXYZD(1,NS3), PTXYZD(1,NS4), Q4 )
C
         IF( NONOU3 .EQ. 0 .AND. NONOU4 .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'OTCHTR: QUADRANGLE NON ADMISSIBLE'
            WRITE(IMPRIM,*) 'CHOIX DES SURFACES MINIMALES'
            CALL XVPAUSE
            Q1 = REAL(
     %           SURTRD(  PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3) )
     %         + SURTRD(  PTXYZD(1,NS1), PTXYZD(1,NS3), PTXYZD(1,NS4) ))
            Q2 = REAL(
     %           SURTRD(  PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS4) )
     %         + SURTRD(  PTXYZD(1,NS2), PTXYZD(1,NS3), PTXYZD(1,NS4) ))
            IF( Q1 .LT. Q2 ) THEN
               GOTO 6
            ELSE
               GOTO 8
            ENDIF
         ELSE IF( NONOU3 .NE. 0 .AND. NONOU4 .NE. 0 ) THEN
C           LES 2 DIAGONALES SONT ADMISSIBLES
C           CELLE FORMANT 2 TRIANGLES DE MEILLEURE QUALITE EST CHOISIE
            IF( MIN(Q1,Q2) .LT. MIN(Q3,Q4) ) GOTO 8
         ELSE IF( NONOU4 .NE. 0 ) THEN
            GOTO 8
         ENDIF
C
C        DIAGONALE 1-3
  6      NAR00  = NA4
         NAR0   = NA1
         NAMIN0 = NA2
         NAMIN  = NA3
         RETURN
C
C        DIAGONALE 2-4
 8       NAR00  = NA1
         NAR0   = NA2
         NAMIN0 = NA3
         NAMIN  = NA4
         RETURN
      ENDIF
C
C     CF NON REDUIT A UN TRIANGLE
C     LA PLUS PETITE ARETE EST NAR0 DANS NOARCF
      NAR00 = NARMIN( 1 )
      NAR0  = NOARCF( 2, NAR00 )
      NAR   = NOARCF( 2, NAR0  )
C
C     PASSAGE DANS LE PLAN NS1-NS2 SUIVANT DE NS2
      NS1   = NOARCF( 1, NAR0 )
      NS2   = NOARCF( 1, NAR  )
      NAR3  = NAR
C     LA MATRICE DE PASSAGE  DE 3D EN 2D
 9    NS3 = NOARCF( 1, NOARCF( 2, NAR3 ) )
      CALL DI3D2D( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %             D2D3, IERR )
      IF( IERR .NE. 0 ) THEN
C        NS3  ALIGNE AVEC NS1-NS2
         NAR3 = NOARCF( 2, NAR3 )
         GOTO 9
      ENDIF
C     LA NOUVELLE ORIGINE
      NS0 = NS1
C
C     LES COORDONNEES DANS LE PLAN
      CALL CI3D2D( PTXYZD(1,NS0) , D2D3 , PTXYZD(1,NS1), XY2D(1,1) )
      CALL CI3D2D( PTXYZD(1,NS0) , D2D3 , PTXYZD(1,NS2), XY2D(1,2) )
C
C     RECHERCHE DANS CETTE ETOILE DU SOMMET OFFRANT LA MEILLEURE QUALITE
C     DU TRIANGLE NS1-NS2 NS3 SANS INTERSECTION AVEC LE CONTOUR FERME
C     ==================================================================
      NBMIN = 0
      NAR3  = NAR
      DMIMA = -1
C
C     PARCOURS DES SOMMETS POSSIBLES NS3
 10   NAR3  = NOARCF( 2, NAR3 )
      IF( NAR3 .NE. NAR0 ) THEN
C
C        IL EXISTE UN SOMMET NS3 DIFFERENT DE NS1 ET NS2
         NS3 = NOARCF( 1, NAR3 )
C
C        LES ARETES NS1-NS3 ET NS2-NS3 INTERSECTENT-ELLES UNE ARETE
C        DU CONTOUR FERME ?
C        ----------------------------------------------------------
C        LES COORDONNEES DU SOMMET NS3 DANS LE PLAN
         CALL CI3D2D( PTXYZD(1,NS0), D2D3, PTXYZD(1,NS3), XY2D(1,3) )
C
C        INTERSECTION DE L'ARETE NS2-NS3 ET DES ARETES DU CF
C        JUSQU'AU SOMMET NS3
         NAR1 = NOARCF( 2, NAR )
C
 15      IF( NAR1 .NE. NAR3 .AND. NOARCF( 2, NAR1 ) .NE. NAR3 ) THEN
C           L'ARETE SUIVANTE
            NAR2 = NOARCF( 2, NAR1 )
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NP1  = NOARCF( 1, NAR1 )
            NP2  = NOARCF( 1, NAR2 )
            CALL CI3D2D(PTXYZD(1,NS0), D2D3, PTXYZD(1,NP1), XY2DAR(1,1))
            CALL CI3D2D(PTXYZD(1,NS0), D2D3, PTXYZD(1,NP2), XY2DAR(1,2))
            CALL INT2AR( XY2D(1,2)  , XY2D(1,3),
     %                   XY2DAR(1,1), XY2DAR(1,2), OUI )
            IF( OUI ) GOTO 10
C           LES 2 ARETES NE S'INTERSECTENT PAS ENTRE LEURS SOMMETS
            NAR1 = NAR2
            GOTO 15
         ENDIF
C
C        INTERSECTION DE L'ARETE NS3-NS1 ET DES ARETES DU CF
C        JUSQU'AU SOMMET DE L'ARETE NAR0
         NAR1 = NOARCF( 2, NAR3 )
C
 18      IF( NAR1 .NE. NAR0 .AND. NOARCF( 2, NAR1 ) .NE. NAR0 ) THEN
C           L'ARETE SUIVANTE
            NAR2 = NOARCF( 2, NAR1 )
C           LE NUMERO DES 2 SOMMETS DE L'ARETE
            NP1  = NOARCF( 1, NAR1 )
            NP2  = NOARCF( 1, NAR2 )
            CALL CI3D2D(PTXYZD(1,NS0), D2D3, PTXYZD(1,NP1), XY2DAR(1,1))
            CALL CI3D2D(PTXYZD(1,NS0), D2D3, PTXYZD(1,NP2), XY2DAR(1,2))
            CALL INT2AR( XY2D(1,1)  , XY2D(1,3),
     %                   XY2DAR(1,1), XY2DAR(1,2), OUI )
            IF( OUI ) GOTO 10
C           LES 2 ARETES NE S'INTERSECTENT PAS ENTRE LEURS SOMMETS
            NAR1 = NAR2
            GOTO 18
         ENDIF
C
C        LE TRIANGLE NS1-NS2-NS3 N'INTERSECTE PAS UNE ARETE DU CONTOUR FERME
C        CALCUL DE LA QUALITE DU  TRIANGLE  NS1-NS2 NS3:
C         = 2 * RAY_INSCRIT / RAY_CIRCONSCRIT
         A = 0
         B = 0
         C = 0
         DO 20 I=1,3
            D = PTXYZD(I,NS2)-PTXYZD(I,NS1)
            A = A + D * D
            D = PTXYZD(I,NS3)-PTXYZD(I,NS2)
            B = B + D * D
            D = PTXYZD(I,NS1)-PTXYZD(I,NS3)
            C = C + D * D
20       CONTINUE
         A  = SQRT ( A )
         B  = SQRT ( B )
         C  = SQRT ( C )
         XP = (A+B+C) * 0.5D0
C        QUALITE: 2 * RAYON_INSCRIT / RAYON_CIRCONSCRIT
         IF ( (A*B*C) .NE. 0 ) THEN
           D = 8D0 * (XP-A)*(XP-B)*(XP-C) / (A*B*C)
         ELSE
           D = 0
         ENDIF
C
         IF( D .GE. DMIMA*1.00001 ) THEN
C           D EST UN VRAI MAXIMUM DE LA QUALITE
            DMIMA = D
            NBMIN = 1
            NARMIN(1) = NAR3
         ELSE IF( D .GE. DMIMA*0.999998 ) THEN
C           D EST VOISIN DE DMIMA
C           IL EST EMPILE
            NBMIN = NBMIN + 1
            NARMIN( NBMIN ) = NAR3
         ENDIF
         GOTO 10
      ENDIF
C
C     BILAN : EXISTE T IL PLUSIEURS SOMMETS DE MEME QUALITE?
C     ======================================================
      IF( NBMIN .GT. 1 ) THEN
C
C        OUI:RECHERCHE DE CEUX DE CERCLE NE CONTENANT PAS D'AUTRES SOMMETS
         DO 80 I=1,NBMIN
C           LE SOMMET
            NAR = NARMIN( I )
            IF( NAR .LE. 0 ) GOTO 80
            NS3 = NOARCF(1,NAR)
C           LES COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT
C           ET SON RAYON
            CALL VDCE2T( NS1, NS2, NS3, PTXYZD,
     %                   CENTRE, RAYON, XYZNOR, IERR )
            IF( IERR .NE. 0 ) THEN
C              LE SOMMET NS3 NE CONVIENT PAS
               NARMIN( I ) = 0
               GOTO 80
            ENDIF
            RAYON = RAYON * RAYON * UNPEPS
            DO 70 J=1,NBMIN
               IF( J .NE. I ) THEN
C                 L'AUTRE SOMMET
                  NAR1 = NARMIN(J)
                  IF( NAR1 .LE. 0 ) GOTO 70
                  NS4 = NOARCF(1,NAR1)
C                 APPARTIENT T IL AU CERCLE NS1 NS2 NS3 ?
                  XG = REAL( CENTRE(1) - PTXYZD(1,NS4) )
                  YG = REAL( CENTRE(2) - PTXYZD(2,NS4) )
                  ZG = REAL( CENTRE(3) - PTXYZD(3,NS4) )
                  IF( XG*XG + YG*YG + ZG*ZG .LE. RAYON ) THEN
C                    NS4 EST DANS LE CERCLE CIRCONSCRIT  NS1 NS2 NS3
C                    LE SOMMET NS3 NE CONVIENT PAS
                     NARMIN( I ) = 0
                     GOTO 80
                  ENDIF
               ENDIF
 70         CONTINUE
 80      CONTINUE
C
C        EXISTE T IL PLUSIEURS SOMMETS ?
         J = 0
         DO 90 I=1,NBMIN
            IF( NARMIN( I ) .GT. 0 ) THEN
C              COMPACTAGE DES MIN
               J = J + 1
               NARMIN(J) = NARMIN(I)
            ENDIF
 90      CONTINUE
C
         IF( J .GT. 1 ) THEN
C           OUI : CHOIX DU PLUS PETIT RAYON DE CERCLE CIRCONSCRIT
            DMIMA = DMAXIM
            DO 120 I=1,NBMIN
               NS3 = NOARCF(1,NARMIN(I))
C
C              LES COORDONNEES DU CENTRE DE CERCLE CIRCONSCRIT
C              AU TRIANGLE NT ET SON RAYON
               CALL VDCE2T( NS1, NS2, NS3, PTXYZD,
     %                      CENTRE, RAYON, XYZNOR, IERR )
               IF( IERR .NE. 0 ) THEN
C                 LE SOMMET NS3 NE CONVIENT PAS
                  GOTO 120
               ENDIF
               IF( RAYON .LT. DMIMA ) THEN
                  DMIMA = RAYON
                  NARMIN(1) = NARMIN(I)
               ENDIF
 120        CONTINUE
         ENDIF
      ENDIF
C
C     LE CHOIX FINAL
C     ==============
      NAMIN = NARMIN(1)
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
