      SUBROUTINE IDSTMA( NBSOEF, NBEF, NOSOEF, NBSTIN, XYZSIN,
     %                   NBSTFI, XYZSFI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER LES SOMMETS AYANT LES MEMES 3 COORDONNEES
C -----    MISE A JOUR DU NO DES SOMMETS DES EF
C
C ENTREES:
C --------
C NBSOEF : NOMBRE DE SOMMETS PAR EF
C NBEF   : NOMBRE D'ELEMENTS FINIS
C NBSTIN : NOMBRE DE SOMMETS INITIAUX AVANT IDENTIFICATION
C XYZSIN : 3 COORDONNEES DES NBSTIN SOMMETS INITIAUX
C
C MODIFIE:
C --------
C NOSOEF : NO DES NBSOEF SOMMETS DES NBEF EF AVANT ET APRES IDENTIFICATION
C
C SORTIES:
C --------
C NBSTFI : NOMBRE DE SOMMETS APRES IDENTIFICATION
C XYZSFI : 3 COORDONNEES DES NBSTFI SOMMETS APRES IDENTIFICATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C MODIFS : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray Novembre 2015
C2345X7...............................................................12
      PARAMETER        (NBIPAX=63, NBIPX1=64,
     %                  NBIPAY=63, NBIPY1=64,
     %                  NBIPAZ=63, NBIPZ1=64)
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
C
      REAL         XYZ( 3 )
      EQUIVALENCE (XYZ(1),X), (XYZ(2),Y), (XYZ(3),Z)

      REAL         XYZEXT( 3, 2 )
      EQUIVALENCE (XYZEXT(1,1),XMINS),(XYZEXT(2,1),YMINS),
     %            (XYZEXT(3,1),ZMINS)
      EQUIVALENCE (XYZEXT(1,2),XMAXS),(XYZEXT(2,2),YMAXS),
     %            (XYZEXT(3,2),ZMAXS)
C
      REAL         XYZSIN( 3, NBSTIN ), XYZSFI( 3, NBSTIN )
      INTEGER      NOSOEF( NBSOEF, NBEF )
      INTRINSIC    NINT

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10001) EPZERO, EPSXYZ
      ELSE
         WRITE(IMPRIM,20001) EPZERO, EPSXYZ
      ENDIF
10001 FORMAT(' idstma: SEUIL IDENTIFICATION XYZ AUTOUR DE ZERO=',G14.6/
     %       '         SEUIL IDENTIFICATION XYZ AILLEURS      =',G14.6,
     %               ' SUR CHACUNE DES 3 COMPOSANTES X Y Z' )
20001 FORMAT(' idstma: XYZ IDENTIFICATION TRESHOLD AROUND ZERO=',G14.6/
     %       '         XYZ IDENTIFICATION TRESHOLD ELSEWHERE  =',G14.6,
     %               ' ON EACH OF 3 COMPONENTS X Y Z' )

C     INITIALISATION DU PAVAGE
C     MAXIMUM DES EPSILON
      EPS    = MAX( EPZERO, EPSXYZ )
      UNMEPS = 1.0 - EPS
      EPS2   = 2 * EPS

C     DECLARATION ET OUVERTURE EN MC DU TABLEAU DES
C     NOUVEAU NUMERO DES SOMMETS APRES IDENTIFICATION
      NBSTFI = 0
      MNNOUV = 0
      CALL TNMCDC( 'ENTIER', 1+NBSTIN, MNNOUV )
      CALL AZEROI( 1+NBSTIN, MCN(MNNOUV) )

C     DECLARATION ET OUVERTURE EN MC DU TABLEAU DES PAVES
      MNPAVE = 0
      MOPAVE = 1 + (NBIPX1+1) * (NBIPY1+1) * (NBIPZ1+1)
      CALL TNMCDC( 'ENTIER', MOPAVE, MNPAVE )
C     MISE A ZERO DU CHAINAGE SUIVANT DES PAVES
      CALL AZEROI( MOPAVE, MCN( MNPAVE) )
C
C     LE TABLEAU DES CHAINAGES SUIVANT DES SOMMETS FINAUX DANS LES PAVES
      MNSUIV = 0
      CALL TNMCDC( 'ENTIER', 1+NBSTIN, MNSUIV )
      CALL AZEROI( 1+NBSTIN, MCN(MNSUIV) )

C     COORDONNEES EXTREMES DES SOMMETS
      DO K=1,3
         X = XYZSIN(K,1)
         XYZEXT(K,1) = X
         XYZEXT(K,2) = X
      ENDDO
      DO I=2,NBSTIN
         DO K=1,3
            X = XYZSIN(K,I)
            IF( X .LT. XYZEXT(K,1) ) XYZEXT(K,1)=X
            IF( X .GT. XYZEXT(K,2) ) XYZEXT(K,2)=X
         ENDDO
      ENDDO

C     CALCUL DES 'ECHELLES' DANS LES 3 DIRECTIONS
      ECHPAX = XMAXS - XMINS
      IF( ECHPAX .LE. EPS ) THEN
         XMAXS  = XMAXS + 1.
         ECHPAX = 1.
      ENDIF
      XMINS  = XMINS - EPS2 * ECHPAX
      XMAXS  = XMAXS + EPS2 * ECHPAX
      ECHPAX = NBIPX1 / ( XMAXS - XMINS )
C
      ECHPAY = YMAXS - YMINS
      IF( ECHPAY .LE. EPS ) THEN
         YMAXS  = YMAXS + 1.
         ECHPAY = 1.
      ENDIF
      YMINS  = YMINS - EPS2 * ECHPAY
      YMAXS  = YMAXS + EPS2 * ECHPAY
      ECHPAY = NBIPY1 / ( YMAXS - YMINS )
C
      ECHPAZ = ZMAXS - ZMINS
      IF( ECHPAZ .LE. EPS ) THEN
         ZMAXS  = ZMAXS + 1.
         ECHPAZ = 1.
      ENDIF
      ZMINS  = ZMINS - EPS2 * ECHPAZ
      ZMAXS  = ZMAXS + EPS2 * ECHPAZ
      ECHPAZ = NBIPZ1 / ( ZMAXS - ZMINS )
C
C     IDENTIFICATION DES SOMMETS INITIAUX
C     -----------------------------------
      DO 50 I=1,NBSTIN

C        LE NUMERO DES PAVES SUSCEPTIBLES DE CONTENIR CE SOMMET I
C        RECHERCHE DES INDICES DU PAVAGE EN X POUVANT CONTENIR LE SOMMET
         X  = XYZSIN(1,I)
         Y  = XYZSIN(2,I)
         Z  = XYZSIN(3,I)

C        LE NUMERO DES PAVES SUSCEPTIBLES DE CONTENIR CE SOMMET I
C        RECHERCHE DES INDICES DU PAVAGE EN X POUVANT CONTENIR LE SOMMET
         D  = ( X - XMINS ) * ECHPAX
         LX = NINT( D )
         D  = D - LX
         IF( D .LE. EPS ) THEN
            LX2 = LX
            LX1 = MAX( 0, LX-1 )
         ELSE IF( D .GE. UNMEPS ) THEN
            LX1 = LX
            LX2 = MIN( NBIPAX, LX+1 )
         ELSE
            LX1 = LX
            LX2 = LX
         ENDIF
C
C        RECHERCHE DES INDICES DU PAVAGE EN Y POUVANT CONTENIR LE SOMMET
         D  = ( Y - YMINS ) * ECHPAY
         LY = NINT( D )
         D  = D - LY
         IF( D .LE. EPS ) THEN
            LY2 = LY
            LY1 = MAX( 0, LY-1 )
         ELSE IF( D .GE. UNMEPS ) THEN
            LY1 = LY
            LY2 = MIN( NBIPAY, LY+1 )
         ELSE
            LY1 = LY
            LY2 = LY
         ENDIF
C
C        RECHERCHE DES INDICES DU PAVAGE EN Z POUVANT CONTENIR LE SOMMET
         D  = ( Z - ZMINS ) * ECHPAZ
         LZ = NINT( D )
         D  = D - LZ
         IF( D .LE. EPS ) THEN
            LZ2 = LZ
            LZ1 = MAX( 0, LZ-1 )
         ELSE IF( D .GE. UNMEPS ) THEN
            LZ1 = LZ
            LZ2 = MIN( NBIPAZ, LZ+1 )
         ELSE
            LZ1 = LZ
            LZ2 = LZ
         ENDIF

C        PARCOURS DES PAVES SUSCEPTIBLES DE CONTENIR XYZ
         DO NZ = LZ1, LZ2
            DO NY = LY1, LY2
               DO NX = LX1, LX2
                  NUPAVE = NX + NBIPX1 * ( NY + NBIPY1 * NZ )
C
C                 LE PREMIER SOMMET DU PAVE NUPAVE
                  NOSOMM = MCN( MNPAVE + NUPAVE )

 10               IF( NOSOMM .GT. 0 ) THEN

C                    LE SOMMET NOSOMM EXISTE : COMPARAISON DES 3 COORDONNEES
                     CALL XYZIDE( XYZSFI(1,NOSOMM), XYZ, IDENTQ )

                     IF( IDENTQ .EQ. 0 ) THEN
C                       LES 2 SOMMETS NE SONT PAS IDENTIQUES
                        GOTO 20
                     ENDIF
C
C                    LES 2 SOMMETS I et NOSOMM ET SONT IDENTIQUES
C                    --------------------------------------------
C                    LE SOMMET I EST MARQUE IDENTIFIE A NOSOMM
                     MCN( MNNOUV + I ) = NOSOMM

C                    AFFICHAGE SI LE 5-EME CHIFFRE SIGNIFICATIF EST DIFFERENT
                     IF( ABS(XYZSFI(1,NOSOMM)-X) .GT. 1E-5*ABS(X) .OR.
     %                   ABS(XYZSFI(2,NOSOMM)-Y) .GT. 1E-5*ABS(Y) .OR.
     %                   ABS(XYZSFI(3,NOSOMM)-Z) .GT. 1E-5*ABS(Z) ) THEN
                        WRITE(IMPRIM,*)
                        IF( LANGAG .EQ. 0 ) THEN
                           WRITE(IMPRIM,10020) I, (XYZ(K),K=1,3),
     %                                  NOSOMM, (XYZSFI(K,NOSOMM),K=1,3)
                        ELSE
                           WRITE(IMPRIM,20020) I, (XYZ(K),K=1,3),
     %                                  NOSOMM, (XYZSFI(K,NOSOMM),K=1,3)
                        ENDIF
                     ENDIF

                     GOTO 50
C
C                    LE SOMMET NOSOMM N'EST PAS XYZ
C                    ------------------------------
C                    PASSAGE AU SOMMET SUIVANT
 20                  NOSOMM = MCN( MNSUIV + NOSOMM )
                     GOTO 10

                  ENDIF

               ENDDO
            ENDDO
         ENDDO

10020 FORMAT(' REMARQUE  :  SOMMET',I8,':',3G15.7/
     %       ' IDENTIFIE AU SOMMET',I8,':',3G15.7)
20020 FORMAT(' REMARK:       VERTEX',I8,':',3G15.7/
     %       ' IDENTIFIED to VERTEX',I8,':',3G15.7)
C
C        SOMMET NON IDENTIFIE . AJOUT DE CE SOMMET
C        -----------------------------------------
         NBSTFI = NBSTFI + 1
         NOSOMM = NBSTFI

C        LE NOUVEAU NUMERO DU SOMMET I INITIAL NON IDENTIFIE
         MCN( MNNOUV + I ) = NOSOMM

C        CHAINAGE SUIVANT DU NOUVEAU SOMMET NOSOMM XYZ EN DEBUT DE SON PAVE
         NUPAVE = LX + NBIPX1 * ( LY + NBIPY1 * LZ )
         MCN( MNSUIV+NOSOMM ) = MCN( MNPAVE+NUPAVE )
         MCN( MNPAVE+NUPAVE ) = NOSOMM

C        LES 3 COORDONNEES DU NOUVEAU SOMMET
         XYZSFI(1,NOSOMM) = XYZ(1)
         XYZSFI(2,NOSOMM) = XYZ(2)
         XYZSFI(3,NOSOMM) = XYZ(3)

 50   ENDDO
C
C     RENUMEROTATION DES SOMMETS DES EF APRES IDENTIFICATION
C     TRANSFORMER TOUT QUADRANGLE AYANT 2 SOMMETS IDENTIQUES EN 1 TRIANGLE
C     SUPPRIMER LES TRIANGLES AYANT 2 SOMMETS IDENTIQUES
ccc   CF CALL MAJSTEFSURF( NBSTFI, XYZSFI, MCN(MNNOUV), NBEF, NOSOEF )
      NEWNBEF = 0
      DO 70 NEF=1,NBEF

C        EF NEF
         DO L=1,NBSOEF
C           LE NO DU SOMMET APRES IDENTIFICATION
            NOSOEF(L,NEF) = MCN( MNNOUV + NOSOEF(L,NEF) )
         ENDDO

C        NOMBRE DE SOMMETS DE NUMERO NUL CREES
         IF( NOSOEF(4,NEF) .EQ. 0 ) THEN
C           TRIANGLE
            NB0 = 1
            NBS = 3
         ELSE
C           QUADRANGLE
            NB0 = 0
            NBS = 4
         ENDIF

C        RECHERCHE DE LA MULTIPLICITE DES SOMMETS DE L'EF NEF
         DO L=1,NBS-1

C           SOMMET L DU QT NEF
            NS0 = NOSOEF(L,NEF)
            IF( NS0 .GT. 0 ) THEN

C              MULTIPLICITE DU SOMMET NS0
 60            NB = 1
               DO M = L+1, NBS
                  IF( NOSOEF(M,NEF) .EQ. NS0 ) THEN
                     NB = NB + 1
                     M2 = M
                  ENDIF
               ENDDO

               IF( NB .GE. 2 ) THEN
C                 NOSOEF(L,NEF) = NOSOEF(M2,NEF)
C                 LE SOMMET DOUBLE EST SUPPRIME
                  NOSOEF(M2,NEF) = 0
                  NB0 = NB0 + 1
                  GOTO 60
               ENDIF

            ENDIF

         ENDDO

         IF( NB0 .EQ. 1 ) THEN

C           LE ZERO EST IL EN POSITION 4?
            IF( NOSOEF(4,NEF) .NE. 0 ) THEN
C              NON: MISE EN 4-EME POSITION
               DO M=1,3
                  IF( NOSOEF(M,NEF) .EQ. 0 ) THEN
                     DO N=M+1,4
                        NOSOEF(N-1,NEF)=NOSOEF(N,NEF)
                     ENDDO
                     NOSOEF(4,NEF) = 0
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF( NB0 .GE. 2 ) THEN

C           EF A SUPPRIMER
            GOTO 70

         ENDIF

C        MISE EN POSITION DU QT
         NEWNBEF = NEWNBEF + 1
         DO L=1,4
            NOSOEF(L,NEWNBEF) = NOSOEF(L,NEF)
         ENDDO

 70   ENDDO

C     NOMBRE FINAL DE QT
      NBEF = NEWNBEF
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
      CALL TNMCDS( 'ENTIER', 1+NBSTIN, MNSUIV )
      CALL TNMCDS( 'ENTIER',   MOPAVE, MNPAVE )
      CALL TNMCDS( 'ENTIER', 1+NBSTIN, MNNOUV )

      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'idstma: NOMBRE INITIAL DE SOMMETS=',NBSTIN
         PRINT *,'idstma: NOMBRE FINAL   DE SOMMETS=',NBSTFI
         PRINT *,'idstma: NOMBRE SOMMETS IDENTIFIES=',NBSTIN-NBSTFI
      ELSE
         PRINT *,'idstma: INITIAL NUMBER of VERTICES=',NBSTIN
         PRINT *,'idstma: FINAL   NUMBER of VERTICES=',NBSTFI
         PRINT *,'idstma: IDENTIFIED VERTICES NUMBER=',NBSTIN-NBSTFI
      ENDIF

      RETURN
      END
