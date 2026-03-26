      SUBROUTINE TRVMTR22( CONTMN, CONTMX,
     %                     NBPIEX, CRITER,
     %                     NBELEM, NBPOEL, NUPTEL, NBPOIT, XYZPOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN 2D LE CRITERE DES CONTRAINTES DE VON MISES TRESCA
C -----    PAR ZONES DE COULEURS SUR LES EF 2D TRIANGLE OU QUADRANGLE
C          (CF $MEFISTO/ther/trzon2.f)
C
C ENTREES :
C ---------
C CONTMN : CONTRAINTE MIN
C CONTMX : CONTRAINTE MAX
C NDIM   : ESPACE DE TRAVAIL 2 OU 3
C NBPIEX : NOMBRE DE NOEUDS SUR LESQUELS LE CRITERE A ETE INTERPOLE
C CRITER : LES NBPIEX*NBELEM VALEURS DU CRITERE
C NBELEM : NOMBRE D'ELEMENTS DE CE TYPE
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C23456---------------------------------------------------------------012
      PARAMETER     ( LIGCON=0 )
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NUPTEL(NBELEM,NBPOEL)
      REAL              XYZPOI(3,NBPOIT)
      DOUBLE PRECISION  CRITER(NBPIEX,NBELEM)
      REAL              CS(2,4), CV(4)
      INTRINSIC         REAL
C
C     PALETTE ARC-EN-CIEL
      CALL PALCDE( 11 )
C
C     LE TRACE
C     --------
C     LIGNES EPAISSIES ET CONTINUES
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
C     LE NOMBRE DE COULEURS DISPONIBLES DANS LA PALETTE
      NBCOUL = NDCOUL - N1COUL
C
      DO 20 K=1,NBELEM
C
C        LES NBPIEX SOMMETS DE L'EF POUR LE CRITERE
         DO 10 I=1,NBPIEX
            NOSI = NUPTEL(K,I)
            CS(1,I) = XYZPOI(1,NOSI)
            CS(2,I) = XYZPOI(2,NOSI)
            COULEUR = REAL( (CRITER(I,K)-CONTMN) / (CONTMX-CONTMN) )
C           PROJECTION SI DEPASSEMENT AU DELA DU MIN OU MAX
            IF( COULEUR .GT. 1.0 ) COULEUR = 1.0
            IF( COULEUR .LT. 0.0 ) COULEUR = 0.0
            CV(I) = N1COUL + NBCOUL * COULEUR
 10      CONTINUE
C
         IF( NBPIEX .EQ. 3 ) THEN
C           TRACE DU TRIANGLE SELON LA COULEUR AUX 3 SOMMETS
            CALL TRIACOUL2D( CS, CV )
         ELSE
C           TRACE DU QUADRANGLE SELON LA COULEUR AUX 4 SOMMETS
            CALL QUADCOUL2D( CS, CV )
         ENDIF
C
 20   CONTINUE
C
      IF( NCOUAF .GE. 0 ) THEN
C
C        LE TRACE DES ARETES
C        -------------------
C        LIGNES NON EPAISSIES ET POINTILLEES
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( NTLAFR )
C
         DO 50 K=1,NBELEM
            DO 40 I=1,NARET
C              LE NUMERO DU 1-ER SOMMET DE L'ARETE
               NS0 = NUPTEL( K, NOSOAR(1,I) )
C              LE NUMERO DU 2-EME SOMMET DE L'ARETE
               NS  = NUPTEL( K, NOSOAR(2,I) )
               IF( NBNOAR(I) .LE. 0 ) THEN
C                 TRACE DE L'ARETE
                  CALL TRAIT2D( NCOUAF, XYZPOI(1,NS0), XYZPOI(2,NS0),
     %                                  XYZPOI(1,NS ), XYZPOI(2,NS ) )
               ELSE
C                 TRACE DU SOMMET1 AU MILIEU
                  NM = NUPTEL( K, NONOAR(1,I) )
                  CALL TRAIT2D( NCOUAF, XYZPOI(1,NS0), XYZPOI(2,NS0),
     %                                  XYZPOI(1,NM ), XYZPOI(2,NM ) )
C                 TRACE DU MILIEU AU SOMMET 2
                  CALL TRAIT2D( NCOUAF, XYZPOI(1,NM), XYZPOI(2,NM),
     %                                  XYZPOI(1,NS), XYZPOI(2,NS) )
               ENDIF
 40         CONTINUE
 50      CONTINUE
C
      ENDIF
C
C     LES LIGNES SONT REMISES A LEUR EPAISSEUR NORMALE
C     ------------------------------------------------
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
