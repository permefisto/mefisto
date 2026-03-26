      SUBROUTINE PT1VUE( NC, NOVUE, XYZPT, KNOMPT, NUPOIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER UN POINT DANS LE MODE VUE pour NOVUE= 1 A 3
C -----    VERSION xvue

C ENTREES :
C ---------
C NC     : NUMERO DE LA COULEUR DE TRACE DU POINT
C NOVUE  : 1 => XZ, 2 => YZ, 3 => XY
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
      REAL              XYZPT(3),XYZCM(3)
      CHARACTER*(*)     KNOMPT
      CHARACTER*25      CARSYM

C     LE NUMERO DES COORDONNEES SELON LA VUE
      IF( NOVUE .EQ. 1 ) THEN
         N1 = 1
         N2 = 3
         XCM1 = XZX2
         XCM2 = XZX1
         YCM1 = XZY1
         YCM2 = XZY2
      ELSE IF( NOVUE .EQ. 2 ) THEN
         N1 = 2
         N2 = 3
         XCM1 = YZX1
         XCM2 = YZX2
         YCM1 = YZY1
         YCM2 = YZY2
      ELSE IF( NOVUE .EQ. 4 ) THEN
         N1 = 2
         N2 = 1
         XCM1 = XYX1
         XCM2 = XYX2
         YCM1 = XYY2
         YCM2 = XYY1
      ELSE
         RETURN
      ENDIF

C     TRANSFORMATION DES COORDONNEES EN CM ENTRE XCM1 XCM2 ET YCM1 YCM2
      XYZCM(N1)=((CADRSA(N1,2)-XYZPT(N1)) * XCM1 +
     %           (XYZPT(N1)-CADRSA(N1,1)) * XCM2 )
     %         / (CADRSA(N1,2)-CADRSA(N1,1))
      XYZCM(N2)=((CADRSA(N2,2)-XYZPT(N2)) * YCM1 +
     %           (XYZPT(N2)-CADRSA(N2,1)) * YCM2 )
     %         / (CADRSA(N2,2)-CADRSA(N2,1))

C     TRACE DU POINT EN PROJECTION DANS LE PLAN N1 N2
C     ===============================================
      IF( XCM1 .GT. XCM2 ) THEN
          A    = XCM1
          XCM1 = XCM2
          XCM2 = A
      ENDIF
      IF( YCM1 .GT. YCM2 ) THEN
          A    = YCM1
          YCM1 = YCM2
          YCM2 = A
      ENDIF
      IF( XYZCM(N1) .GE. XCM1 .AND. XYZCM(N1) .LE. XCM2 .AND.
     %    XYZCM(N2) .GE. YCM1 .AND. XYZCM(N2) .LE. YCM2 ) THEN

C         LE POINT VISIBLE EST TRACE
C         LE POINT VISIBLE EST MIS DANS LE TABLEAU MNITEP
          IF( MCN(MNITEP+2) .GE. MCN(MNITEP+1) ) THEN
C            TABLEAU TROP PETIT: LA TAILLE DU TABLEAU EST AUGMENTEE
             CALL ITEMAU( MNITEP )
          ENDIF

C         LE NOMBRE D'ITEMS EST AUGMENTE DE 1
          MCN(MNITEP+2) = MCN(MNITEP+2) + 1
C         L'ADRESSE MCN DE L'ITEM
          MNI  = MNITEP + MCN(MNITEP) * MCN(MNITEP+2)
C         LES COORDONNEES IMAGE DE CET ITEM, SON CODE
          MCN( MNI     ) = NUPXEX( XYZCM(N1) )
          MCN( MNI + 1 ) = NUPXEY( XYZCM(N2) )
C         NUMERO DU POINT DANS SON LEXIQUE
          MCN( MNI + 2 ) = NUPOIN

C         TRACE EFFECTIF DU POINT
          CARSYM = '+' // KNOMPT
          CALL SYMBOLE2D( NC, XYZCM(N1), XYZCM(N2), CARSYM )

      ENDIF
      END
