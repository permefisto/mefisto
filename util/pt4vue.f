      SUBROUTINE PT4VUE( NC, XYZPT, KNOMPT, NUPOIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER UN POINT DANS LE MODE 4 VUES
C -----    VERSION xvue

C ENTREES :
C ---------
C NC     : NUMERO DE LA COULEUR DE TRACE DU POINT
C XYZPT  : LES 3 COORDONNEES DU POINT
C KNOMPT : LE NOM DU POINT
C NUPOIN : LE NUMERO DU POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       JUIN 1994
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/cm4vue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              XYZPT(3)
      CHARACTER*(*)     KNOMPT
      CHARACTER*25      CARSYM

C     LE TRACE DANS LE PLAN XZ
      CALL PT1VUE( NC, 1, XYZPT, KNOMPT, NUPOIN )

C     LE TRACE DANS LE PLAN YZ
      CALL PT1VUE( NC, 2, XYZPT, KNOMPT, NUPOIN )

C     LE TRACE DANS LE PLAN XY
      CALL PT1VUE( NC, 4, XYZPT, KNOMPT, NUPOIN )

C     LE TRACE EN AXONOMETRIE
      CALL AXOYZ( XYZPT(1), XYZPT(2), XYZPT(3),  YAX, ZAX )

C     TRANSFORMATION DES COORDONNEES EN CM ENTRE AXX1,AXX2 ET AXY1,AXY2
      YAXCM=((YAXMAX-YAX) * AXX1 + (YAX-YAXMIN) * AXX2 )
     %         / (YAXMAX-YAXMIN)
      ZAXCM=((ZAXMAX-ZAX) * AXY1 + (ZAX-ZAXMIN) * AXY2 )
     %         / (ZAXMAX-ZAXMIN)

C     TRACE DU POINT EN PROJECTION DANS LE PLAN N1 N2
C     ===============================================
      IF( YAXCM .GE. AXX1 .AND. YAXCM .LE. AXX2 .AND.
     %    ZAXCM .GE. AXY1 .AND. ZAXCM .LE. AXY2 ) THEN

C         LE POINT VISIBLE EST MIS DANS LE TABLEAU MNITEP
          IF( MCN(MNITEP+2) .GE. MCN(MNITEP+1) ) THEN
C            TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
             CALL ITEMAU( MNITEP )
          ENDIF

C         LE NOMBRE D'ITEMS EST AUGMENTE DE 1
          MCN(MNITEP+2) = MCN(MNITEP+2) + 1
C         L'ADRESSE MCN DE L'ITEM
          MNI = MNITEP + MCN(MNITEP) * MCN(MNITEP+2)
C         LES COORDONNEES PIXELS DE CET ITEM, SON CODE
          MCN( MNI     ) = NUPXEX( YAXCM )
          MCN( MNI + 1 ) = NUPXEY( ZAXCM )
C         NUMERO DU POINT DANS SON LEXIQUE
          MCN( MNI + 2 ) = NUPOIN

C         TRACE EFFECTIF DU POINT
          CARSYM = '+' // KNOMPT
          CALL SYMBOLE2D( NC, YAXCM, ZAXCM, CARSYM )

      ENDIF
      END
