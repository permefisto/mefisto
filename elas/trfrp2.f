      SUBROUTINE TRFRP2( AMPLDE, NCAS,   NDIMES,
     %                   NBELEM, NBNOEL, NUNDEL, XYZNOE,
     %                   NTDL  , NBVECT, DEPLAC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    POUR UN TYPE D'EF 2D TRACER LES EF ET LES EF DEFORMES AVEC
C -----    UNE AMPLIFICATION DES DEPLACEMENTS ET UN
C          MOUVEMENT SELON LA FREQUENCE PROPRE DE L'OBJET
C
C ENTREES:
C --------
C AMPLDE : FACTEUR D'AMPLIFICATION DES DEPLACEMENTS
C NCAS   : NUMERO DU CAS A TRAITER
C NDIMES : ESPACE DE TRAVAIL 2 OU 3
C NBELEM : NOMBRE D'ELEMENTS DE CE TYPE
C NUNDEL : NUMERO DES NOEUDS PAR ELEMENTS
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE
C NTDL   : NOMBRE TOTAL DE COMPOSANTES DES DEPLACEMENTS
C NBVECT : NOMBRE TOTAL DE VECTEURS DEPLACEMENTS
C DEPLAC : LES COMPOSANTES DES DEPLACEMENTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & VEULETTES sur MER AVRIL 2009
C23456---------------------------------------------------------------012
      PARAMETER        (LIGCON=0, LIGTIR=1 )
      include"./incl/langue.inc"
      include"./incl/ponoel.inc"
      include"./incl/mecoit.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NUNDEL(NBELEM,NBNOEL)
      REAL              XYZNOE(3,*)
      DOUBLE PRECISION  DEPLAC(NTDL,NBVECT)
      INTRINSIC         REAL
C
      IF( NCOUAF .GE. 0 ) THEN
C
C        LE TRACE DES ARETES SANS DEPLACEMENT
C        ------------------------------------
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( NTLAFR )
C
         DO 100 K=1,NBELEM
            DO 20 I=1,NARET
C              LE NUMERO DES 2 SOMMETS DE L'ARETE I DE L'EF K
               NS1 = NUNDEL( K, NOSOAR(1,I) )
               NS2 = NUNDEL( K, NOSOAR(2,I) )
               IF( NBNOAR(I) .LE. 0 ) THEN
C                 TRACE DE L'ARETE
                  CALL TRAIT2D( NCOUAF, XYZNOE(1,NS1), XYZNOE(2,NS1),
     %                                  XYZNOE(1,NS2), XYZNOE(2,NS2) )
               ELSE
C                 TRACE DU SOMMET1 AU MILIEU
                  NM = NUNDEL( K, NONOAR(1,I) )
                  CALL TRAIT2D( NCOUAF, XYZNOE(1,NS1), XYZNOE(2,NS1),
     %                                  XYZNOE(1,NM ), XYZNOE(2,NM ) )
C                 TRACE DU MILIEU AU SOMMET 2
                  CALL TRAIT2D( NCOUAF, XYZNOE(1,NM ), XYZNOE(2,NM ),
     %                                  XYZNOE(1,NS2), XYZNOE(2,NS2) )
               ENDIF
 20         CONTINUE
 100     CONTINUE
C
      ENDIF
C
C     LE TRACE DES ARETES DEFORMEES PAR LE DEPLACEMENT
C     ------------------------------------------------
      CALL XVEPAISSEUR( 2 )
      CALL XVTYPETRAIT( LIGCON )
C
      DO 200 K=1,NBELEM
         DO 120 I=1,NARET
C           LE NUMERO DES 2 SOMMETS DE L'ARETE I DE L'EF K
            NS1 = NUNDEL( K, NOSOAR(1,I) )
            XS1 = REAL( XYZNOE(1,NS1)
     %          + AMPLDE * DEPLAC(NDIMES*(NS1-1)+1,NCAS) )
            YS1 = REAL( XYZNOE(2,NS1)
     %          + AMPLDE * DEPLAC(NDIMES*(NS1-1)+2,NCAS) )
            NS2 = NUNDEL( K, NOSOAR(2,I) )
            XS2 = REAL( XYZNOE(1,NS2)
     %          + AMPLDE * DEPLAC(NDIMES*(NS2-1)+1,NCAS) )
            YS2 = REAL( XYZNOE(2,NS2)
     %          + AMPLDE * DEPLAC(NDIMES*(NS2-1)+2,NCAS) )
C
            IF( NBNOAR(I) .LE. 0 ) THEN
C              TRACE DE L'ARETE I DE L'EF K
               CALL TRAIT2D( NCOUAD, XS1, YS1, XS2, YS2 )
            ELSE
C              TRACE DU SOMMET1 AU MILIEU DE L'ARETE I DE L'EF K
               NM = NUNDEL( K, NONOAR(1,I) )
               XM = REAL( XYZNOE(1,NM)
     %            + AMPLDE * DEPLAC(NDIMES*(NM-1)+1,NCAS) )
               YM = REAL( XYZNOE(2,NM)
     %            + AMPLDE * DEPLAC(NDIMES*(NM-1)+2,NCAS) )
               CALL TRAIT2D( NCOUAD, XS1, YS1, XM, YM )
C              TRACE DU MILIEU AU SOMMET 2 DE L'ARETE I DE L'EF K
               CALL TRAIT2D( NCOUAD, XM, YM, XS2, YS2 )
            ENDIF
C
 120     CONTINUE
 200  CONTINUE
C
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
