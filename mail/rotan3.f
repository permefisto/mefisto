      SUBROUTINE ROTAN3( NBANGL, ANGLES, PT1,    PT2,    MNSOSU,
     %                   FERME,  NBSOSU, SOMSUR, NOSOAX,NOSOAX2,
     %                   NUTYMS, NBFASU,
     %                   NBSOVO, SOMVOL, NUTYMV, NSCUVO,
     %                   NBPENT1, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERATION DU MAILLAGE DU VOLUME PAR ROTATION D'UNE SURFACE
C -----
C ENTREES :
C ---------
C NBANGL  : NOMBRE DES ANGLES DE ROTATION
C ANGLES  : LES NBANGL ANGLES DE ROTATION( EN DEGRES )
C PT1 PT2 : LES 3 COORDONNEES DES 2 POINTS DE DEFINITION DE L'AXE
C MNSOSU  : ADRESSE MCN DU TABLEAU NSEF DE LA SURFACE INITIALE
C FERME   : VRAI SI LA SOMME DES ANGLES VAUT 360 DEGRES, FAUX SINON
C NBSOSU  : LE NOMBRE DE SOMMETS DE LA SURFACE
C SOMSUR  : LES 3 COORDONNEES DES NBSOSU SOMMETS DE LA SURFACE
C NBSOAX  : LE NOMBRE DE SOMMETS DE LA SURFACE  SUR L'AXE
C NOSOAX  : +LE NUMERO DU SOMMET S'IL N'EST PAS SUR L'AXE PT1-PT2
C           -LE NUMERO DU SOMMET S'IL   EST     SUR L'AXE PT1-PT2
C NOSOAX2 : +LE NUMERO DU SOMMET S'IL N'EST PAS SUR L'AXE PT1-PT2
C           DE LA COUCHE SUIVANTE
C NUTYMS  : NUMERO TYPE DU MAILLAGE DE LA SURFACE
C NBFASU  : NOMBRE DE FACES DE LA SURFACE
C NBSOVO  : LE NOMBRE DE SOMMETS DU VOLUME MAILLE
C NUTYMV  : LE NUMERO TYPE DU MAILLAGE DU VOLUME
C           0      MAILLAGE NON STRUCTURE => TABLEAU NSCUVO     GENERE
C           6 OU 7 MAILLAGE     STRUCTURE => TABLEAU NSCUVO NON GENERE
C
C SORTIES :
C ---------
C SOMVOL  : LES 3 COORDONNEES DES NBSOVO SOMMETS DU VOLUME MAILLE
C NSCUVO  : LES 8 SOMMETS DES CUBES DU VOLUME
C NBPENT1 : NOMBRE DE PENTAEDRES DUS AUX QUADRANGLES AVEC 1 SOMMET SUR L'AXE
C IERR    : 1 SI UN SOMMET EST SUR L'AXE
C           0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1990
C MODIFS : PERRONNET ALAIN Labo J.L. LIONS UPMC PARIS       FEVRIER 2007
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      REAL     ANGLES(NBANGL),PT1(3),PT2(3),
     %         SOMSUR(3,NBSOSU),SOMVOL(3,NBSOVO),XYZ(3)
      INTEGER  NOSOAX(NBSOSU), NOSOAX2(NBSOSU),
     %         NSCUVO(8,NBFASU,NBANGL)
      INTEGER  NOSOEL(64)
      LOGICAL  FERME
      DOUBLE PRECISION D2D3(3,3), A, ARAD, PI, COSA, SINA
C
      PI = ATAN(1D0) * 4D0
C
C     CHANGEMENT DE REPERE PT1: PT1-PT2, PT1-1-ER POINT DE LA SURFACE
C                                            NON COLINEAIRE A L'AXE
      DO 10 I=1,NBSOSU
         IF( NOSOAX(I) .LT. 0 ) GOTO 10
C        LE SOMMET I N'EST PAS COLINEAIRE A L'AXE
         CALL DF3D2D( PT1 , PT2 , SOMSUR(1,I), D2D3 , IERR )
         IF( IERR .EQ. 0 ) GOTO 20
 10   CONTINUE
      RETURN
C
C     LES PARAMETRES DES FACES DE LA SURFACE
 20   CALL NSEFPA( MCN(MNSOSU) ,
     %             NUTYMS , NBSOEL , NBSOEF , NBTGEF,
     %             LDAPEF , LDNGEF , LDTGEF , NBEFOB,
     %             NX     , NY     , NZ     ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LA SURFACE N'EST PAS COLINEAIRE A PT1-PT2
C     LES COORDONNEES DES SOMMETS DE LA 1-ERE SURFACE DANS LE VOLUME
      CALL TRTATA( SOMSUR , SOMVOL , 3*NBSOSU )
C
C     LES NBANGL ROTATIONS
      A    = 0D0
      NBSS = NBSOSU
C
C     NOMBRE DE PENTAEDRES POUR QUADRANGLES AVEC UN SOMMET SUR L'AXE
      NBTETR2 = 0
      NBPYRAT = 0
      NBPYRAQ = 0
      NBPENT1 = 0
      NBPENT2 = 0
      NBPENT0 = 0
      NBHEXA0 = 0
C
      DO 100 N=1,NBANGL
C
C        L'ANGLE ACTUEL DE ROTATION
         A = A + ANGLES(N)
C        LE COSINUS ET SINUS
         ARAD = A * PI / 180D0
         COSA = COS( ARAD )
         SINA = SIN( ARAD )
C
C        LES COORDONNEES DES SOMMETS DE LA SURFACE DE ROTATION
         IF( FERME .AND. N.EQ.NBANGL ) THEN
C
C           LA DERNIERE COUCHE EST PERIODIQUE
            DO 30 I=1,NBSOSU
               IF( NOSOAX(I) .GT. 0 ) THEN
C                 SOMMET NON SUR L'AXE = SOMMET INITIAL
                  NOSOAX2(I) = I
               ELSE
C                 SOMMET SUR L'AXE
                  NOSOAX2(I) = NOSOAX(I)
               ENDIF
 30         CONTINUE
C
         ELSE
C
C           UNE COUCHE DE PLUS
            DO 50 I=1,NBSOSU
               IF( NOSOAX(I) .GT. 0 ) THEN
C                 SOMMET A AJOUTER CAR NON SUR L'AXE
                  NBSS = NBSS + 1
                  NOSOAX2(I) = NBSS
C                 LE SOMMET DE LA SURFACE EST RAMENE DANS LE NOUVEAU REPERE
                  CALL CH3D3D( PT1, D2D3, SOMSUR(1,I), SOMVOL(1,NBSS) )
C                 LA ROTATION AUTOUR DE L'AXE PT1-PT2 CAD X
                  S2 = SOMVOL(2,NBSS)
                  S3 = SOMVOL(3,NBSS)
                  XYZ(1) = SOMVOL(1,NBSS)
                  XYZ(2) = REAL( COSA * S2 - SINA * S3 )
                  XYZ(3) = REAL( SINA * S2 + COSA * S3 )
C                 RETOUR AU REPERE INITIAL
                  CALL CH3D3R( PT1, D2D3, XYZ, SOMVOL(1,NBSS) )
               ELSE
C                 SOMMET SUR L'AXE
                  NOSOAX2(I) = NOSOAX(I)
               ENDIF
 50         CONTINUE
         ENDIF
C
C        LES CUBES NON STRUCTURES ENGENDRES PAR CETTE ROTATION
         IF( NUTYMV .NE. 6 .AND. NUTYMV .NE. 7 ) THEN
C
C           LA BOUCLE SUR LES TRIANGLES OU QUADRANGLES DE LA SURFACE INITIALE
            DO 90 I=1,NBEFOB
C              LE NUMERO DES 4 SOMMETS DE CETTE FACE I
               CALL NSEFNS( I      , NUTYMS , NBSOEF , NBTGEF,
     %                      LDAPEF , LDNGEF , LDTGEF ,
     %                      MNSOSU , NX , NY , NZ ,
     %                      NCOGEL , NUGEEF , NUEFTG, NOSOEL, IERR )
C              LE NOMBRE DE SOMMETS DE CET ELEMENT FINI TRIANGLE OU QUAD
               NBSO = NBSOME( NCOGEL )
C
C              RECHERCHE DU NOMBRE DE SOMMETS DE CET EF SUR L'AXE
               NBSSAX = 0
               DO 62 J=1,NBSO
                  IF( NOSOAX( NOSOEL(J) ) .LT. 0 ) THEN
                     NBSSAX = NBSSAX + 1
                  ENDIF
 62            CONTINUE
C
               IF( NBSO .EQ. 3) THEN
                  IF( NBSSAX .EQ. 0 ) THEN
C
C                    PENTAEDRE A PARTIR D'UN TRIANGLE AVEC 0 SOMMET SUR AXE
                     DO 70 J=1,NBSO
                        NOST = NOSOEL(J)
                        NSCUVO(J  ,I,N) = NOSOAX( NOSOEL(J) )
                        NSCUVO(J+3,I,N) = NOSOAX2(NOSOEL(J) )
 70                  CONTINUE
                     NBPENT0 = NBPENT0 + 1
                     GOTO 78
C
                  ELSE IF( NBSSAX .EQ. 1 ) THEN
C
C                    PYRAMIDE A PARTIR D'UN TRIANGLE AVEC 1 SOMMET SUR AXE
                     DO 72 J=1,NBSO
                        NOST = NOSOEL(J)
                        IF( NOSOAX(NOST) .LT. 0 ) THEN
C                          LE SOMMET SUR L'AXE
                           NSCUVO(5,I,N) = NOST
C                          LE SOMMET SUIVANT
                           IF( J .EQ. NBSO ) THEN
                              NS = NOSOEL(1)
                           ELSE
                              NS = NOSOEL(J+1)
                           ENDIF
                           NSCUVO(1,I,N) = NOSOAX(NS)
                           NSCUVO(4,I,N) = NOSOAX2(NS)
C                          LE SOMMET PRECEDANT
                           IF( J .EQ. 1 ) THEN
                              NS = NOSOEL(NBSO)
                           ELSE
                              NS = NOSOEL(J-1)
                           ENDIF
                           NSCUVO(2,I,N) = NOSOAX(NS)
                           NSCUVO(3,I,N) = NOSOAX2(NS)
                           NSCUVO(6,I,N) = 0
                           NBPYRAT = NBPYRAT + 1
                           GOTO 78
                        ENDIF
 72                  CONTINUE
C
                  ELSE IF( NBSSAX .EQ. 2 ) THEN
C
C                    TETRAEDRE A PARTIR D'UN TRIANGLE AVEC 2 SOMMETS SUR AXE
                     DO 76 J=1,NBSO
                        NOST = NOSOEL(J)
                        IF( NOSOAX(NOST) .GT. 0 ) THEN
C                          LE SOMMET NON SUR L'AXE EST NOST
                           NS1 = ABS( NOSOAX( NOSOEL(1) ) )
                           NS2 = ABS( NOSOAX( NOSOEL(2) ) )
                           NS3 = ABS( NOSOAX( NOSOEL(3) ) )
                           NS4 = NOSOAX2(NOST)
                           NSCUVO(1,I,N) = NS1
                           NSCUVO(2,I,N) = NS2
                           NSCUVO(3,I,N) = NS3
                           NSCUVO(4,I,N) = NS4
                           NSCUVO(5,I,N) = 0
                           NSCUVO(6,I,N) = 0
                           NBTETR2 = NBTETR2 + 1
                           GOTO 78
                        ENDIF
 76                  CONTINUE
C
                  ENDIF
C
C                 COMPLETION DES 2 DERNIERS NO DE SOMMETS
 78               NSCUVO(7,I,N) = 0
                  NSCUVO(8,I,N) = 0
C
               ELSE
C
C                 EF QUADRANGLE DE LA SURFACE
                  IF( NBSSAX .EQ. 0 ) THEN
C
C                    HEXAEDRE A PARTIR D'UN QUADRANGLE AVEC 0 SOMMET SUR AXE
                     DO 80 J=1,NBSO
                        NOST = NOSOEL(J)
                        NSCUVO(J  ,I,N) = NOSOAX(NOST)
                        NSCUVO(J+4,I,N) = NOSOAX2(NOST)
 80                  CONTINUE
                     NBHEXA0 = NBHEXA0 + 1
C
                  ELSE IF( NBSSAX .EQ. 1 ) THEN
C
C                    PYRAMIDE + PENTAEDRE A PARTIR D'UN QUADRANGLE AVEC
C                    AVEC 1 SOMMET SUR L'AXE
                     DO 82 J=1,NBSO
                        NOST = NOSOEL(J)
                        IF( NOSOAX(NOST) .LT. 0 ) THEN
C                          LE SOMMET SUR L'AXE EST NOST
C
C                          LA PYRAMIDE
                           NSCUVO(5,I,N) = NOST
C                          LE SOMMET SUIVANT
                           IF( J .EQ. NBSO ) THEN
                              K = 1
                           ELSE
                              K = J+1
                           ENDIF
                           NS1 = NOSOEL( K )
                           NSCUVO(1,I,N) = NOSOAX(NS1)
                           NSCUVO(4,I,N) = NOSOAX2(NS1)
C                          LE SOMMET PRECEDANT
                           IF( J .EQ. 1 ) THEN
                              NS = NOSOEL(NBSO)
                           ELSE
                              NS = NOSOEL(J-1)
                           ENDIF
                           NSCUVO(2,I,N) = NOSOAX(NS)
                           NSCUVO(3,I,N) = NOSOAX2(NS)
                           NSCUVO(6,I,N) = 0
                           NSCUVO(7,I,N) = 0
                           NSCUVO(8,I,N) = 0
                           NBPYRAQ = NBPYRAQ + 1
C
C                          LE PENTAEDRE AU DELA DES EF 3D
C                          LE SOMMET SUIVANT LE SUIVANT K
                           IF( K .EQ. NBSO ) THEN
                              NS2 = NOSOEL(1)
                           ELSE
                              NS2 = NOSOEL(K+1)
                           ENDIF
                           NBPENT1 = NBPENT1 + 1
                           NSCUVO(1,NBEFOB+NBPENT1,NBANGL)= NOSOAX(NS1)
                           NSCUVO(2,NBEFOB+NBPENT1,NBANGL)= NOSOAX(NS)
                           NSCUVO(3,NBEFOB+NBPENT1,NBANGL)= NOSOAX(NS2)
                           NSCUVO(4,NBEFOB+NBPENT1,NBANGL)= NOSOAX2(NS1)
                           NSCUVO(5,NBEFOB+NBPENT1,NBANGL)= NOSOAX2(NS)
                           NSCUVO(6,NBEFOB+NBPENT1,NBANGL)= NOSOAX2(NS2)
                           NSCUVO(7,NBEFOB+NBPENT1,NBANGL)= 0
                           NSCUVO(8,NBEFOB+NBPENT1,NBANGL)= 0
                           GOTO 90
                        ENDIF
 82                  CONTINUE
C
                  ELSE IF( NBSSAX .EQ. 2 ) THEN
C
C                    PENTAEDRE A PARTIR D'UN QUADRANGLE AVEC 2 SOMMETS SUR AXE
                     NS  = 0
                     NS1 = 0
                     NS2 = 0
                     NOST= 0
                     DO 84 J=1,NBSO
                        NOST = NOSOEL(J)
                        IF( NOSOAX(NOST) .LT. 0 ) THEN
C                          LE 1ER  SOMMET SUR L'AXE
C                          LE 2EME SOMMET SUR L'AXE EST LE SUIVANT DE J
C                                                    OU LE DERNIER NBSO
                           K  = J+1
                           NS = NOSOEL(K)
                           IF( NOSOAX(NS) .LT. 0 ) THEN
C                             LES 2 AUTRES SOMMETS DU QUADRANGLE NON SUR L'AXE
C                             SONT LES 2 SUIVANTS DE J+1
                              K = K + 1
                              IF( K .GT. NBSO ) K=1
                              NS1 = NOSOEL(K)
                              K = K+1
                              IF( K .GT. NBSO ) K=1
                              NS2 = NOSOEL(K)
                           ELSE
C                             LES 2 AUTRES SOMMETS DU QUADRANGLE NON SUR L'AXE
C                             SONT LES NUMEROS 2 ET 3
                              NS =  NOSOEL(4)
                              NS1 = NOSOEL(2)
                              NS2 = NOSOEL(3)
C                             PERMUTATION DES 2 SOMMETS SUR L'AXE
                              K    = NOST
                              NOST = NS
                              NS   = K
                           ENDIF
                           GOTO 85
                        ENDIF
 84                  CONTINUE
C                    NOST-NS SONT SUR L'AXE, NS1 NS2 SUIVENT
 85                  NBPENT2 = NBPENT2 + 1
                     NSCUVO(1,I,N) = NOST
                     NSCUVO(2,I,N) = NOSOAX2(NS2)
                     NSCUVO(3,I,N) = NOSOAX(NS2)
                     NSCUVO(4,I,N) = NS
                     NSCUVO(5,I,N) = NOSOAX2(NS1)
                     NSCUVO(6,I,N) = NOSOAX(NS1)
                     NSCUVO(7,I,N) = 0
                     NSCUVO(8,I,N) = 0
                  ENDIF
               ENDIF
 90         CONTINUE
C
         ENDIF
C
C        MISE A JOUR DU NO DE SOMMET DE LA DERNIERE COUCHE => 1-ERE COUCHE ENSUI
         DO 95 I=1,NBSOSU
            NOSOAX(I) = NOSOAX2(I)
 95      CONTINUE
C
 100  CONTINUE
C
      IF(LANGAG .EQ. 0 ) THEN
      WRITE(IMPRIM,*) NBTETR2,
     %' TETRAEDRES DUS AUX TRIANGLES   AVEC 2 SOMMETS SUR L''AXE'
      WRITE(IMPRIM,*) NBPYRAT,
     %' PYRAMIDES  DUS AUX TRIANGLES   AVEC 1 SOMMET  SUR L''AXE'
      WRITE(IMPRIM,*) NBPYRAQ,
     %' PYRAMIDES  DUS AUX QUADRANGLES AVEC 1 SOMMET  SUR L''AXE'
      WRITE(IMPRIM,*) NBPENT0,
     %' PENTAEDRES DUS AUX QUADRANGLES AVEC 0 SOMMET  SUR L''AXE'
      WRITE(IMPRIM,*) NBPENT1,
     %' PENTAEDRES DUS AUX QUADRANGLES AVEC 1 SOMMET  SUR L''AXE'
      WRITE(IMPRIM,*) NBPENT2,
     %' PENTAEDRES DUS AUX QUADRANGLES AVEC 2 SOMMETS SUR L''AXE'
      WRITE(IMPRIM,*) NBHEXA0,
     %' HEXAEDRES  DUS AUX QUADRANGLES AVEC 0 SOMMET  SUR L''AXE'
      ELSE
      WRITE(IMPRIM,*) NBTETR2,
     %' TETRAHEDRA DUE to TRIANGLES   WITH 2 VERTICES on the AXIS'
      WRITE(IMPRIM,*) NBPYRAT,
     %' PYRAMIDS   DUE to TRIANGLES   WITH 1 VERTEX   on the AXIS'
      WRITE(IMPRIM,*) NBPYRAQ,
     %' PYRAMIDS   DUE to QUADRANGLES WITH 1 VERTEX   on the AXIS'
      WRITE(IMPRIM,*) NBPENT0,
     %' PENTAHEDRA DUE to QUADRANGLES WITH 0 VERTEX   on the AXIS'
      WRITE(IMPRIM,*) NBPENT1,
     %' PENTAHEDRA DUE to QUADRANGLES WITH 1 VERTEX   on the AXIS'
      WRITE(IMPRIM,*) NBPENT2,
     %' PENTAHEDRA DUE to QUADRANGLES WITH 2 VERTICES on the AXIS'
      WRITE(IMPRIM,*) NBHEXA0,
     %' HEXAHEDRA  DUE to QUADRANGLES WITH 0 VERTEX   on the AXIS'
      ENDIF
C
      RETURN
      END
