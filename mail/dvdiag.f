      SUBROUTINE DVDIAG( MXTRIA, PXYD, NOTRIA, NOTRSO, NOSUTR, NUISOP,
     %                   LISOGD, NBCHGT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECHANGER LES DIAGONALES POUR RETROUVER LA TRIANGULATION
C -----    DE TYPE DELAUNAY.
C          LE TRIANGLE NT0 EST SUPPOSE DANS LA TRIANGULATION
C          LORSQU'UNE ARETE EST TRAITEE, LE NUMERO DU TRIANGLE
C          OPPOSE EST NEGATIF DANS LES 2 TRIANGLES
C          EN FIN DE TRAITEMENT DE CE SP, LE SIGNE REDEVIENT POSITIF
C
C ENTREES:
C --------
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C
C ENTREES ET SORTIES :
C --------------------
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE i
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C NOSUTR : NUMERO DE SURFACE DE CHAQUE TRIANGLE
C NUISOP : NUMERO DE L'ISO DU POINT ET 0 SINON
C LISOGD : 1 SI LES ARETES ISO    DOIVENT     ETRE GARDEES
C          0 SI LES ARETES ISO NE DOIVENT PAS ETRE GARDEES
C
C SORTIES:
C --------
C NBCHGT : NOMBRE DE CHANGEMENTS DE DIAGONALES
C IERR   : 0 SI PAS D'ERREUR
C          4 SATURATION DES ARETES DE LA PILE
C          5 ANOMALIE DANS LA TOPOLOGIE DES TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1995
C....................................................................012
      DOUBLE PRECISION  EPSCER
      PARAMETER        (EPSCER=1D-8)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NOTRSO(*),
     %                  NOSUTR(MXTRIA),
     %                  NUISOP(*)
      DOUBLE PRECISION  PXYD(3,*),
     %                  CETRI1(3),CETRI2(3)
      DOUBLE PRECISION  D1, D2, SURTD2, SUR123, SUR234
C
C     ECHANGES DES DIAGONALES DES QUADRANGLES POUR OBTENIR
C     UNE TRIANGULATION DELAUNAY
C     ====================================================
      NBCHGT = 0
      DO 1000 NT=1,MXTRIA
         IF( NOTRIA(1,NT) .LE. 0 ) GOTO 1000
         DO 900 NAR=1,3
C           LE TRIANGLE OPPOSE
            NTOP = NOTRIA(3+NAR,NT)
            IF( NTOP .LE. 0 ) GOTO 900
C
C           SI    LES 2 TRIANGLES SONT SUR DES SURFACES DIFFERENTES
C           ALORS ILS NE PEUVENT ETRE MODIFIES
            IF( NOSUTR(NT) .NE. NOSUTR(NTOP) ) GOTO 900
C
C           LES SOMMETS DE L'ARETE COMMUNE
            NS1 = NOTRIA(NAR,NT)
            IF( NAR .LT. 3 ) THEN
               NA1 = NAR + 1
            ELSE
               NA1 = 1
            ENDIF
            IF( NA1 .LT. 3 ) THEN
               NA2 = NA1 + 1
            ELSE
               NA2 = 1
            ENDIF
            NS2 = NOTRIA(NA1,NT)
C
            IF( LISOGD .NE. 0 ) THEN
               IF( NUISOP(NS1) .GT. 0 .AND. NUISOP(NS2) .GT. 0 .AND.
     %             NUISOP(NS1) .EQ. NUISOP(NS2) ) GOTO 900
            ENDIF
C
C           LE SOMMET OPPOSE A L'ARETE NAR DANS NT
            NS3 = NOTRIA(NA2,NT)
C
C           RECHERCHE DE L'ARETE NS2-NS1 DANS NTOP
            DO 110 N=1,3
               IF( NOTRIA(3+N,NTOP) .EQ. NT ) GOTO 115
 110        CONTINUE
C
C           ERREUR DANS LE VOISINAGE DES TRIANGLES NT ET NTOP
            WRITE(IMPRIM,*) 'DVDIAG: ARETE ',NS1, NS2,
     %        ' DU TRIANGLE ',NTOP, ' SANS LE TRIANGLE OPPOSE ',NT
            WRITE(IMPRIM,*)'TRIANGLE',NTOP,': ',(NOTRIA(N1,NTOP),N1=1,6)
            WRITE(IMPRIM,*)'TRIANGLE ',NT  ,': ',(NOTRIA(N1,NT),N1=1,6)
            NBLGRC(NRERR) = 1
            KERR(1) ='ANOMALIE DANS TOPOLOGIE DES TRIANGLES'
            CALL LEREUR
            IERR = 5
            GOTO 9999
C
 115        IF( N .LT. 3 ) THEN
               N1 = N + 1
            ELSE
               N1 = 1
            ENDIF
            IF( N1 .LT. 3 ) THEN
               N2 = N1 + 1
            ELSE
               N2 = 1
            ENDIF
C           LE SOMMET OPPOSE A L'ARETE NAR DANS NTOP
            NS4 = NOTRIA(N2,NTOP)
C
C           L'ARETE NAR DE NT EST MARQUEE DANS LES 2 TRIANGLES NT ET NTOP
            NOTRIA(3+NAR,NT) = NTOP
            NOTRIA(3+N,NTOP) = NT
C
C           NUMERO MINIMUM ET MAXIMUM DES TRIANGLES TRAITES
            IF( NTOP .LT. NTMIN ) NTMIN = NTOP
            IF( NTOP .GT. NTMAX ) NTMAX = NTOP
C
C           CALCUL DU CENTRE DU CERCLE CIRCONSCRIT A NTOP ET DU CARRE DU RAYON
            IER = -1
            CALL CENCED( PXYD(1,NOTRIA(1,NTOP)), PXYD(1,NOTRIA(2,NTOP)),
     %                   PXYD(1,NOTRIA(3,NTOP)), CETRI1, IER )
            IF( IER .NE. 0 ) THEN
C              TRIANGLE DEGENERE
               GOTO 120
            ENDIF
C
C           DISTANCE AU CARRE DE NS3 AU CENTRE DU CERCLE DE NTOP
            D1 = (PXYD(1,NS3)-CETRI1(1)) **2 +
     %           (PXYD(2,NS3)-CETRI1(2)) **2
C
C           DISTANCE AU CARRE DE NS4 AU CENTRE DU CERCLE DE NT
C           CALCUL DU CENTRE DU CERCLE CIRCONSCRIT ET DU CARRE DU RAYON
            IER = -1
            CALL CENCED( PXYD(1,NOTRIA(1,NT)), PXYD(1,NOTRIA(2,NT)),
     %                   PXYD(1,NOTRIA(3,NT)), CETRI2, IER )
            IF( IER .NE. 0 ) THEN
C              TRIANGLE DEGENERE
               GOTO 120
            ENDIF
            D2 = (PXYD(1,NS4)-CETRI2(1)) **2 +
     %           (PXYD(2,NS4)-CETRI2(2)) **2
C
C           NS3 ET NS4 SONT ILS DANS LE CERCLE DU TRIANGLE OPPOSE?
            IF( (CETRI1(3)-D1) .LT. EPSCER*D1  .OR.
     %          (CETRI2(3)-D2) .LT. EPSCER*D2 ) GOTO 900
C
C           LES POINTS NS3 ET NS4 SONT INTERNES AUX CERCLES
C           ===============================================
C           LES NOUVEAUX TRIANGLES SERAIENT ILS DEGENERES?
C           CALCUL DU CENTRE DU CERCLE CIRCONSCRIT ET DU CARRE DU RAYON
 120        IER = -1
            CALL CENCED( PXYD(1,NS2), PXYD(1,NS3), PXYD(1,NS4),
     %                CETRI1, IER )
            IF( IER .NE. 0 ) THEN
C              TRIANGLE DEGENERE REFUS DE L'ECHANGE
               GOTO 900
            ENDIF
            IER = -1
            CALL CENCED( PXYD(1,NS3), PXYD(1,NS1), PXYD(1,NS4),
     %                   CETRI1, IER )
            IF( IER .NE. 0 ) THEN
C              TRIANGLE DEGENERE REFUS DE L'ECHANGE
               GOTO 900
            ENDIF
C
C           LE QUADRANGLE EST IL CONVEXE ?
            SUR123 = ABS( SURTD2(PXYD(1,NS1),PXYD(1,NS2),PXYD(1,NS3)) )
     %             + ABS( SURTD2(PXYD(1,NS1),PXYD(1,NS4),PXYD(1,NS2)) )
            SUR234 = ABS( SURTD2(PXYD(1,NS2),PXYD(1,NS3),PXYD(1,NS4)) )
     %             + ABS( SURTD2(PXYD(1,NS3),PXYD(1,NS1),PXYD(1,NS4)) )
            IF( ABS(SUR123-SUR234) .GT. 1D-4 * SUR123 ) GOTO 900
C
C           LE QUADRANGLE FORME DES TRIANGLES 123+142 EST CONVEXE => ECHANGE POS
C
C           ECHANGE EFFECTIF DES TRIANGLES 123-214 EN 143-234
C           =================================================
C           LE NUMERO D'UN TRIANGLE CONTENANT CHACUN DES 4 SOMMETS
            NOTRSO( NS1 )   = NTOP
            NOTRSO( NS2 )   = NT
            NOTRSO( NS3 )   = NT
            NOTRSO( NS4 )   = NT
C
C           ECHANGE DES SOMMETS
C           NOTRIA(NA1,NT ) = NS2
C           NOTRIA(NA2,NT ) = NS3
            NOTRIA(NAR,NT ) = NS4
C
            NOTRIA(N ,NTOP) = NS3
C           NOTRIA(N1,NTOP) = NS1
C           NOTRIA(N2,NTOP) = NS4
C
C           LES TRIANGLES OPPOSES AUX 4 COTES DU QUADRANGLE
            NTT23 = NOTRIA(NA1+3,NT)
            NTT31 = NOTRIA(NA2+3,NT)
            NTT14 = NOTRIA(N1+3, NTOP)
            NTT42 = NOTRIA(N2+3, NTOP)
C
            NOTRIA(NAR+3,NT ) = NTT42
            NOTRIA(NA1+3,NT ) = NTT23
C           LA DIAGONALE COMMUNE EST TRAITEE
            NOTRIA(NA2+3,NT ) = NTOP
C
            NOTRIA(N +3,NTOP) = NTT31
            NOTRIA(N1+3,NTOP) = NTT14
C           LA DIAGONALE COMMUNE EST TRAITEE
            NOTRIA(N2+3,NTOP) = NT
C
C           CORRECTION NT NTOP ?
            CALL CORTRI( NT,   NOTRIA )
            CALL CORTRI( NTOP, NOTRIA )
C
C           LE NUMERO D'ARETE DE NS4-NS2 DANS NTT42
            IF( NTT42 .GT. 0 ) THEN
               IF( NOTRIA(1,NTT42) .EQ. NS2 ) THEN
                  I = 4
               ELSE IF( NOTRIA(2,NTT42) .EQ. NS2 ) THEN
                  I = 5
               ELSE
                  I = 6
               ENDIF
               NOTRIA(I,NTT42) = NT
            ENDIF
C
C           LE NUMERO D'ARETE DE NS3-NS2 DANS NTT23
            IF( NTT23 .GT. 0 ) THEN
               IF( NOTRIA(1,NTT23) .EQ. NS3 ) THEN
                  I = 4
               ELSE IF( NOTRIA(2,NTT23) .EQ. NS3 ) THEN
                  I = 5
               ELSE
                  I = 6
               ENDIF
C              L'ARETE EST DEMARQUEE
               NOTRIA(I,NTT23) = NT
            ENDIF
C
C           LE NUMERO D'ARETE DE NS3-NS1 DANS NTT31
            IF( NTT31 .GT. 0 ) THEN
               IF( NOTRIA(1,NTT31) .EQ. NS1 ) THEN
                  I = 4
               ELSE IF( NOTRIA(2,NTT31) .EQ. NS1 ) THEN
                  I = 5
               ELSE
                  I = 6
               ENDIF
               NOTRIA(I,NTT31) = NTOP
            ENDIF
C
C           LE NUMERO D'ARETE DE NS1-NS4 DANS NTT14
            IF( NTT14 .GT. 0 ) THEN
               IF( NOTRIA(1,NTT14) .EQ. NS4 ) THEN
                  I = 4
               ELSE IF( NOTRIA(2,NTT14) .EQ. NS4 ) THEN
                  I = 5
               ELSE
                  I = 6
               ENDIF
C              L'ARETE EST DEMARQUEE
               NOTRIA(I,NTT14) = NTOP
            ENDIF
C
C           CORRECTION ?
            CALL CORTRI( NTT42, NOTRIA )
            CALL CORTRI( NTT23, NOTRIA )
            CALL CORTRI( NTT31, NOTRIA )
            CALL CORTRI( NTT14, NOTRIA )
C
C           TRACE DES 2 TRIANGLES NT ET NTOP
            CALL DVTRTR( PXYD, NOTRIA, NT,  NCCYAN, NCMAGE )
            CALL DVTRTR( PXYD, NOTRIA, NTOP,NCCYAN, NCMAGE )
            NBCHGT = NBCHGT + 1
C
 900     CONTINUE
 1000 CONTINUE
C
 9999 RETURN
      END
