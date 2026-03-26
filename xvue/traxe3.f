      SUBROUTINE TRAXE3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER DES 3 AXES
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS      MAI 1994
C AUTEUR:ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      include"./incl/mecoit.inc"
      include"./incl/traaxe.inc"
      COMMON /UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL             XYZ1(3), XYZ2(3), XYZ3(3)
      CHARACTER*16     KFORM
      CHARACTER*24     KVAL0, KVAL
C
      IF( INTERA .LE. 0 ) RETURN
C
C     PROTECTION CONTRE UNE DIVISION PAR ZERO
      DO 5 I=1,4
         IF( XYZAMPLI(I) .EQ. 0 ) XYZAMPLI(I)=1.0
 5    CONTINUE
C
      IF( NTRAXE .LE. 0 ) GOTO 9900
ccc      IF( NETAXE .GT. 0 ) GOTO 9900
C
C     LE TRACE DES 3 AXES
      CALL XVEPAISSEUR( 3 )
C
      IF( NTRAXE .EQ. 1 ) THEN
C
C        TRACE DU PLAN XY ET AXE Z EN REDUCTION SANS GRADUATION
C        ======================================================
C        LES COORDONNEES DU CENTRE DES AXES
         ECART   = MAX( ABS( XOBMAX-XOBMIN ),
     %                  ABS( YOBMAX-YOBMIN ) ) * 0.09
         XYZ1(1) = XOBMIN + ABS( XOBMAX-XOBMIN ) * 0.1
         XYZ1(2) = YOBMIN + ABS( YOBMAX-YOBMIN ) * 0.1
         XYZ1(3) = 0.0
C
         XYZ2(1) = XYZ1(1)
         XYZ2(2) = XYZ1(2)
         XYZ2(3) = ECART
C
C        LES COORDONNEES XYZ DU POINT
         CALL AXOXYZ( XYZ1, XYZ1 )
         CALL AXOXYZ( XYZ2, XYZ2 )
         ECART = SQRT( (XYZ2(1)-XYZ1(1))**2 +
     %                 (XYZ2(2)-XYZ1(2))**2 +
     %                 (XYZ2(3)-XYZ1(3))**2 )
C
C        LA TRONCATURE PAR LES PLANS ARRIERE ET AVANT EST ANNIHILEE
         AXAR = AXOARR
         AXAV = AXOAVA
         AXOARR = -1E20
         AXOAVA =  1E20
C
C        LE POINT MAX EN X Y A Z=0
         XYZ3(1) = XYZ1(1) + ECART
         XYZ3(2) = XYZ1(2) + ECART
         XYZ3(3) = XYZ1(3)
C
C        L'AXE DES X
         XYZ2(1) = XYZ1(1) + ECART
         XYZ2(2) = XYZ1(2)
         XYZ2(3) = XYZ1(3)
         CALL TRAIT3D( NCVERT, XYZ3, XYZ2 )
         CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )
C        TRACE DU NOM DE L'AXE
         CALL TEXTE3D( NCNOIR, XYZ2, '.X' )
C
C        L'AXE DES Y
         XYZ2(1) = XYZ1(1)
         XYZ2(2) = XYZ1(2) + ECART
         XYZ2(3) = XYZ1(3)
         CALL TRAIT3D( NCROUG, XYZ2, XYZ3 )
         CALL TRAIT3D( NCVERT, XYZ1, XYZ2 )
C        TRACE DU NOM DE L'AXE
         CALL TEXTE3D( NCNOIR, XYZ2, '.Y' )
C
C        L'AXE DES Z
         XYZ2(1) = XYZ1(1)
         XYZ2(2) = XYZ1(2)
         XYZ2(3) = XYZ1(3) + ECART
         CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )
C        TRACE DU NOM DE L'AXE
         CALL TEXTE3D( NCNOIR, XYZ2, '.Z' )
C
C        LES PLANS ARRIERE ET AVANT SONT REACTIVES
         AXOARR = AXAR
         AXOAVA = AXAV
C
      ELSE IF( NTRAXE .GE. 2 ) THEN
C
C        TRACE DES 3 AXES MIN-MAX de COOEXT AVEC GRADUATIONS
C        ===================================================
C        ECART EN % DEFINIT LA MARGE ENTRE LES AXES ET LE BORD ECRAN
         ECARTX = (COOEXT(1,2)-COOEXT(1,1)) * 0.016
         ECARTY = (COOEXT(2,2)-COOEXT(2,1)) * 0.016
         ECARTZ = (COOEXT(3,2)-COOEXT(3,1)) * 0.023
C
C        L'AXE DES X
C        ===========
         XYZ2(1) = COOEXT(1,2)
         XYZ2(2) = COOEXT(2,1)
         XYZ2(3) = COOEXT(3,1)
         CALL TRAIT3D( NCROUG, COOEXT(1,1), XYZ2 )
C
C        LE TRACE DE LA GRADUATION SUR L'AXE X
C        -------------------------------------
         CALL GRADUA( XYZAMPLI(1), COOEXT(1,1), COOEXT(1,2),
     %                PREMIER, PAS, NBPAS, NBDFOR )
C        LE FORMAT D'ECRITURE DE CHAQUE VALEUR
         WRITE( KFORM, 10000 ) NBDFOR+8,NBDFOR
10000    FORMAT('(G',I2,'.',I1,') ')
ccc         print *,'traxe3 X: KFORM=',KFORM,' NBDFOR=',NBDFOR
C
         DO 10 I=0,NBPAS
C           L'ABSCISSE DU NOMBRE A TRACER
            X = PREMIER + I * PAS
            IF( X .LT. COOEXT(1,1) ) GOTO 10
            IF( X .GT. COOEXT(1,2) ) GOTO 15
C           TRACE D'UN TRAIT SELON Y
            XYZ1(1) = X
            XYZ1(2) = COOEXT(2,1)
            XYZ1(3) = COOEXT(3,1)
            XYZ2(1) = X
            XYZ2(2) = COOEXT(2,1) + ECARTY
            XYZ2(3) = COOEXT(3,1)
            CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )
            XYZ2(1) = X
            XYZ2(2) = COOEXT(2,1)
            XYZ2(3) = COOEXT(3,1) + ECARTZ
            CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )
C           COORDONNEES DU POINT
            W = X / XYZAMPLI(1)
C           LA VALEUR W A TRACER EN EVITANT LES 1E-7 ...
            IF( ABS(W) .LT. 1E-3 * ABS(PAS/XYZAMPLI(1)) ) W=0
            WRITE( KVAL0, KFORM ) W
            CALL REELCA( KVAL0, KVAL )
            XYZ1(2) = XYZ1(2) - ECARTY * 1.1
            XYZ1(3) = XYZ1(3) - ECARTZ * 0.4
            CALL SANSBL( KVAL, NBC )
            CALL TEXTE3D( NCNOIR, XYZ1, KVAL(1:NBC) )
 10      CONTINUE
C
C        // A L'AXE DES Y
C        ================
 15      XYZ1(1) = COOEXT(1,2)
         XYZ1(2) = COOEXT(2,1)
         XYZ1(3) = COOEXT(3,1)
         XYZ2(1) = COOEXT(1,2)
         XYZ2(2) = COOEXT(2,2)
         XYZ2(3) = COOEXT(3,1)
         CALL TRAIT3D( NCVERT, XYZ1, XYZ2 )
C
C        LE TRACE DE LA GRADUATION SUR L'AXE Y
C        -------------------------------------
         CALL GRADUA( XYZAMPLI(2),COOEXT(2,1),COOEXT(2,2),
     %                PREMIER,PAS,NBPAS,NBDFOR)
C        LE FORMAT D'ECRITURE DE CHAQUE VALEUR
         WRITE( KFORM,10000 ) NBDFOR+8,NBDFOR
ccc         print *,'traxe3 Y: KFORM=',KFORM,' NBDFOR=',NBDFOR
C
         DO 20 I=0,NBPAS
C           L'ABSCISSE DU NOMBRE A TRACER
            Y = PREMIER + I * PAS
            IF( Y .LT. COOEXT(2,1) ) GOTO 20
            IF( Y .GT. COOEXT(2,2) ) GOTO 25
C           TRACE D'UN TRAIT SELON X
            XYZ1(1) = COOEXT(1,2)
            XYZ1(2) = Y
            XYZ1(3) = COOEXT(3,1)
            XYZ2(1) = COOEXT(1,2) - ECARTX
            XYZ2(2) = Y
            XYZ2(3) = COOEXT(3,1)
            CALL TRAIT3D( NCVERT, XYZ1, XYZ2 )
            XYZ2(1) = COOEXT(1,2)
            XYZ2(2) = Y
            XYZ2(3) = COOEXT(3,1) + ECARTZ
            CALL TRAIT3D( NCVERT, XYZ1, XYZ2 )
C           COORDONNEES DU POINT
            W = Y / XYZAMPLI(2)
C           LA VALEUR W A TRACER EN EVITANT LES 1E-7 ...
            IF( ABS(W) .LT. 1E-3 * ABS(PAS/XYZAMPLI(2)) ) W=0
            if( abs(w) .gt. 1E20 ) then
               KVAL = '???'
               GOTO 17
            ENDIF
ccc            print *,'traxe3: kform=',kform,' w=',w
            WRITE( KVAL0, KFORM ) W
            CALL REELCA( KVAL0, KVAL )
 17         XYZ1(1) = XYZ1(1) + ECARTX * 1.1
            XYZ1(3) = XYZ1(3) - ECARTZ * 0.3
            CALL SANSBL( KVAL, NBC )
            CALL TEXTE3D( NCNOIR, XYZ1, KVAL(1:NBC) )
 20      CONTINUE
C
C        L'AXE DES Z
C        ===========
 25      XYZ2(1) = COOEXT(1,1)
         XYZ2(2) = COOEXT(2,1)
         XYZ2(3) = COOEXT(3,2)
         CALL TRAIT3D( NCBLEU, COOEXT(1,1), XYZ2 )
C
C        LE TRACE DE LA GRADUATION SUR L'AXE Z
C        -------------------------------------
         IF( NTRAXZ .EQ. 0 ) THEN
            ZMIN = COOEXT(3,1)
            ZMAX = COOEXT(3,2)
            ECH  = 1.0
         ELSE
            ZMIN = ZMIAXZ
            ZMAX = ZMXAXZ
            ECH  = (COOEXT(3,2)-COOEXT(3,1)) / (ZMAX-ZMIN)
         ENDIF
         CALL GRADUA( XYZAMPLI(3),ZMIN,ZMAX,
     %                PREMIER,PAS,NBPAS,NBDFOR)
C        LE FORMAT D'ECRITURE DE CHAQUE VALEUR
         WRITE( KFORM,10000 ) NBDFOR+8,NBDFOR
ccc         print *,'traxe3 Z: KFORM=',KFORM,' NBDFOR=',NBDFOR
C
         DO 30 I=0,NBPAS
C           LA COTE DU NOMBRE A TRACER
            IF( NTRAXZ .EQ. 0 ) THEN
               Z = PREMIER + I * PAS
               V = Z
            ELSE
               V = PREMIER + I * PAS
               Z = COOEXT(3,1) + (PREMIER-ZMIN)*ECH + ECH * I * PAS
            ENDIF
            IF( Z .LT. COOEXT(3,1) ) GOTO 30
            IF( Z .GT. COOEXT(3,2) ) GOTO 35
C           TRACE D'UN TRAIT SELON X
            XYZ1(1) = COOEXT(1,1)
            XYZ1(2) = COOEXT(2,1)
            XYZ1(3) = Z
            XYZ2(1) = COOEXT(1,1) + ECARTX
            XYZ2(2) = COOEXT(2,1)
            XYZ2(3) = Z
            CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )
C           TRACE D'UN TRAIT SELON Y
            XYZ2(1) = COOEXT(1,1)
            XYZ2(2) = COOEXT(2,1) + ECARTY
            XYZ2(3) = Z
            CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )
C           TRACE D'UN TRAIT SELON Y
            XYZ2(1) = COOEXT(1,1) + ECARTX
            XYZ2(2) = COOEXT(2,1)
            XYZ2(3) = Z
            CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )
C           LA VALEUR W A TRACER EN EVITANT LES 1E-7 ...
            W = V / XYZAMPLI(3)
            IF( ABS(W) .LT. 1E-3 * ABS(PAS/XYZAMPLI(3)) ) W=0
            WRITE( KVAL0, KFORM ) W
            CALL REELCA( KVAL0, KVAL )
            XYZ1(2) = XYZ1(2) - ECARTY * 1.1
            XYZ1(3) = XYZ1(3) - ECARTZ * 0.3
            CALL SANSBL( KVAL, NBC )
            CALL TEXTE3D( NCNOIR, XYZ1, KVAL(1:NBC) )
 30      CONTINUE
      ENDIF
C
 35   IF( NTRAXE .EQ. 3 ) THEN
C
C        AJOUT DE LA PARALLELE A L'AXE DES X
         XYZ1(1) = COOEXT(1,1)
         XYZ1(2) = COOEXT(2,2)
         XYZ1(3) = COOEXT(3,1)
C
         XYZ2(1) = COOEXT(1,2)
         XYZ2(2) = COOEXT(2,2)
         XYZ2(3) = COOEXT(3,1)
         CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )
C
C        L'AXE DES Y
         CALL TRAIT3D( NCVERT, COOEXT(1,1), XYZ1 )
      ENDIF
C
C     TRACE SUR LA FENETRE
CCC 9900 CALL MEMPXFENETRE   ici CCC EVITE LE SCINTILLEMENT
C
 9900 CALL XVEPAISSEUR( 0 )
ccc      NETAXE = 3  ne sert pas ailleurs...
      RETURN
      END
