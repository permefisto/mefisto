      SUBROUTINE NSEFNS( NUELEM, NUTYMA, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNNSEF, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETROUVER LES NUMEROS DES SOMMETS DE L'EF NUELEM
C -----     QUE LE MAILLAGE SOIT STRUCTURE OU NON
C           CF ~/td/d/a___nsef
C
C ENTREES:
C ---------
C NUELEM :  NUMERO DE L'EF
C NUTYMA :  NUMERO DE TYPE DU MAILLAGE
C           0 : 'NON STRUCTURE'
C           2 : 'SEGMENT   STRUCTURE'  ,
C           3 : 'TRIANGLE  STRUCTURE'  , 4 : 'QUADRANGLE STRUCTURE' ,
C           5 : 'TETRAEDRE STRUCTURE'  , 6 : 'PENTAEDRE  STRUCTURE' ,
C           7 : 'HEXAEDRE  STRUCTURE'  , 8 : '6-CUBE     STRUCTURE'
C NBSOEF :  NOMBRE DE SOMMETS   STOCKES PAR EF (0 POUR COMPLETER)
C NBTGEF :  NOMBRE DE TANGENTES STOCKES PAR EF (0 POUR COMPLETER)
C           ( TRIANGLE STRUCTURE SANS TG => NBSOEF=4 NBTGEF=0
C             TRIANGLE STRUCTURE AVEC TG => NBSOEF=4 NBTGEF=8 ... )
C LDAPEF :  NOMBRE DE MOTS DE DECALAGE POUR POINTER SUR LPEFAP(1)
C           ADRESSE MCN = MNNSEF + LDAPEF
C LDNGEF :  NOMBRE DE MOTS DE DECALAGE POUR POINTER SUR NGEFTG(1)
C           ADRESSE MCN = MNNSEF + LDNGEF
C LDTGEF :  NOMBRE DE MOTS DE DECALAGE POUR POINTER SUR NUTGEF(1,1)
C           ADRESSE MCN = MNNSEF + LDTGEF
C
C MNNSEF :  ADRESSE MCN DU TABLEAU 'NSEF'
C NX,NY,NZ: LE NOMBRE D'ARETES DANS LES DIRECTION X Y Z SI STRUCTURE
C
C SORTIES:
C --------
C NCOGEL :  NUMERO DU CODE GEOMETRIQUE DE L'EF NUELEM
C           1:POINT , 2:SEGMENT , 3:TRIANGLE , 4:QUADRANGLE ,
C           5:TETRAEDRE , 6:PENTAEDRE , 7:HEXAEDRE, 8:6-CUBE, 9:PYRAMIDE
C NUGEEF :  >=0 NUMERO DE DEFINITION GEOMETRIQUE DE L'EF
C           CF ~/td/d/a___nsef
C           ...  A COMPLETER AU FUR ET A MESURE
C NUEFTG : >0 NUMERO DE CET EF NUELEM PARMI LES EF A TG
C          =0 SI L'EF NUELEM N'A PAS DE TG
C
C NOSOEL :  NUMERO DES NBSOEF SOMMETS
C           SUIVI EVENTUELLEMENT DES NBTGEF +-TANGENTES DE L'EF NUELEM
C           LE - INDIQUE QUE LES COMPOSANTES DE LA TANGENTE SONT A INVERSER
C
C           SI PAS DE TANGENTES (NBTGEF=0):
C              POINT     : NO SOMMET1
C              ARETE     : NO SOMMET1 , NS2
C              TRIANGLE  : NO SOMMET1 , NS2 , NS3, 0
C              QUADRANGLE: NO SOMMET1 , NS2 , NS3, NS4
C              TETRAEDRE : NO SOMMET1 , NS2 , NS3, NS4, 0,   0,   0,   0
C              PENTAEDRE : NO SOMMET1 , NS2 , NS3, NS4, NS5, NS6, 0,   0
C              HEXAEDRE  : NO SOMMET1 , NS2 , NS3, NS4, NS5, NS6, NS7, NS8
C
C           S'IL EXISTE DES TANGENTES (NBTGEF>0) ELLES SONT RANGEES PAR SOMMET:
C           ATTENTION: SUR LA FACE SUPERIEURE DU PENTAEDRE ET HEXAEDRE
C                      LE TRIEDRE EST INDIRECT CAR Tk:Si Si+1, Tk+1:Si Si-1, Tk+
C              POINT     : NO SOMMET1,
C                          NO TANGENTE1
C              ARETE     : NO SOMMET1, NS2,
C                          NO TANGENTE1(S1S2), NT2(S2S1)
C              TRIANGLE  : NO  SOMMET1, NS2 , NS3, 0,
C                          NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1)
C                          NO TANGENTE5(S3S1), NT6(S3S2),   0,         0
C              QUADRANGLE: NO  SOMMET1, NS2 , NS3, NS4,
C                          NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1)
C                          NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C              TETRAEDRE : NO SOMMET1, NS2 , NS3, NS4, 0 , 0 , 0 , 0 ,
C                          NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                          NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S4),
C                          NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S4),
C                          NO TANGENT10(S4S1), NT11(S4S2), NT12(S4S3),
C                              0,               0,          0,
C                              0,               0,          0,
C                              0,               0,          0,
C                              0,               0,          0
C              PENTAEDRE : NO SOMMET1 , NS2 , NS3, NS4, NS5, NS6, 0 , 0 ,
C                          NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                          NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S5),
C                          NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S6),
C                          NO TANGENT10(S4S5), NT11(S4S6), NT12(S4S1),
C                          NO TANGENT13(S5S6), NT14(S5S4), NT15(S5S2),
C                          NO TANGENT16(S6S4), NT17(S6S5), NT18(S6S3),
C                              0,               0,          0,
C                              0,               0,          0
C              HEXAEDRE  : NO SOMMET1 , NS2 , NS3, NS4, NS5, NS6, NS7 , NS8
C                          NO TANGENTE1(S1S2), NT2(S1S4),  NT3(S1S5),
C                          NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S6),
C                          NO TANGENTE7(S3S4), NT8(S3S2),  NT9(S3S7),
C                          NO TANGENT10(S4S1), NT11(S4S3), NT12(S4S8),
C                          NO TANGENT13(S5S6), NT14(S5S8), NT15(S5S1),
C                          NO TANGENT16(S6S7), NT17(S6S5), NT18(S6S2),
C                          NO TANGENT19(S7S8), NT20(S7S6), NT21(S7S3),
C                          NO TANGENT22(S8S5), NT23(S8S7), NT24(S8S4)
C
C              6-CUBE    : NO SOMMET1 , NS2 , ... , NS NX*NY
C
C           CE CHOIX PERMET UNE BOUCLE SUR LES TANGENTES PAR LES SOMMETS
C           EN PRESENCE DE TANGENTES LE TABLEAU NUTGEF(1:NBTGEF,1:NBEFOB)
C           EST SITUE A LA FIN DU TMS ~/td/d/a___nsef
C
C IERR   :  CODE D'ERREUR 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1989
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS        MAI 1996
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/a___nsef.inc"

C     8 SOMMETS AU PLUS ET 32 AVEC 3 TANGENTES EN LES 8 SOMMETS
C     LES NUMEROS DES SOMMETS PUIS EVENTUELLEMENT LES NO DES TANGENTES
      INTEGER           NOSOEL(1:64), NSOM(6)
      INTRINSIC         INT, SQRT

      IF( NUELEM .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'nsefns:',NUELEM,' NUMERO d''EF INTERDIT'
         ELSE
            PRINT*,'nsefns:',NUELEM,' FORBIDDEN FINITE ELEMENT NUMBER'
         ENDIF
C        MISE A ZERO DU NO DES NBSOEF SOMMETS
         DO I=1,NBSOEF
            NOSOEL(I) = 0
         ENDDO
         IERR = 1
         RETURN
      ENDIF

      IERR = 0

C     RECHERCHE SELON LE TYPE DU MAILLAGE
      IF( NUTYMA .EQ. 0 ) THEN

C        MAILLAGE NON STRUCTURE
C        ======================
         MN = MNNSEF + WUSOEF + NBSOEF * (NUELEM-1) - 1
         DO 1 I=1,NBSOEF
            NOSOEL(I) = MCN( MN + I )
  1      ENDDO

C        NO DU TYPE DE L'ELEMENT FINI (SURFACE ou VOLUME)
         NCOGEL = MCN( MNNSEF + WUTYOB )
         GOTO( 100 , 100 , 3 , 4 , 4 ) , NCOGEL

C        POINT OU SEGMENT
C 2      GOTO 100

C        TRIANGLE OU QUADRANGLE
  3      IF( NOSOEL(4) .EQ. 0 ) THEN
            NCOGEL = 3
         ELSE
            NCOGEL = 4
         ENDIF
         GOTO 100

C        6-CUBE HEXAEDRE PENTAEDRE TETRAEDRE
  4      IF( NBSOEF .EQ. 64 ) THEN
            NCOGEL = 8
         ELSE IF( NOSOEL(8) .GT. 0 ) THEN
            NCOGEL = 7
         ELSE IF( NOSOEL(6) .GT. 0 ) THEN
            NCOGEL = 6
         ELSE IF( NOSOEL(5) .GT. 0 ) THEN
            NCOGEL = 9
         ELSE IF( NOSOEL(4) .GT. 0 ) THEN
            NCOGEL = 5
         ENDIF
         GOTO 100

      ELSE IF ( NUTYMA .GT. 0 ) THEN

C        MAILLAGE STRUCTURE
C        ==================
C        LE CODE GEOMETRIQUE DE L'ELEMENT
         NCOGEL = NUTYMA
         GOTO ( 10, 20, 30, 40, 50, 60, 70, 80 ) , NUTYMA
C
C        NOEUDSOMMET STRUCTURE
C        =====================
  10     NOSOEL(1) = NUELEM
         GOTO 100

C        SEGMENT STRUCTURE
C        =================
  20     NOSOEL(1) = NUELEM
         NOSOEL(2) = NUELEM + 1
         GOTO 100

C        TRIANGLE STRUCTURE
C        ==================
  30     NBAR = NX
         NBELEM = NBAR * NBAR
         IF (NUELEM.GT.NBELEM) GOTO 9900
C        RECHERCHE DU NUMERO DE L'ELEMENT
         N1B = INT( SQRT(NUELEM*1D0) )
         DO N1=N1B,N1B+1
            NEL  = (N1-1)*(N1-1)
            NEL1 = NEL - 1
            I1   = (N1*N1-N1)/2 + 1
            I2   = (N1*N1+N1)/2 + 1
            DO N2=1,N1
               NEL1 = NEL1 + 2
               IF ( NEL1 .EQ. NUELEM ) THEN
                  NOSOEL(1) = I1
                  NOSOEL(2) = I2 + 1
                  NOSOEL(3) = I2
                  NOSOEL(4) = 0
                  GOTO 100
               END IF
               I1 = I1 + 1
               I2 = I2 + 1
            ENDDO
            NEL2 = NEL
            I1   = (N1*N1+N1)/2 + 1
            I2   = (N1*N1-N1)/2 + 1
            DO N2=1,N1-1
               NEL2 = NEL2 + 2
               IF ( NEL2 .EQ. NUELEM ) THEN
                  NOSOEL(1) = I1 + 1
                  NOSOEL(2) = I2
                  NOSOEL(3) = I2 + 1
                  NOSOEL(4) = 0
                  GOTO 100
               END IF
               I1 = I1 + 1
               I2 = I2 + 1
            ENDDO
         ENDDO
         GOTO 100

C        QUADRANGLE STRUCTURE
C        ====================
 40      NBELEM = NX * NY
         IF (NUELEM.GT.NBELEM) GOTO 9900
         NPY = NUELEM / NX
         NPX = NUELEM - NPY * NX
         IF( NPX .EQ. 0 ) THEN
            NPX = NX
            NPY = NPY - 1
         ENDIF
         NX1 = NX + 1
C        NUMERO DU SOMMET GAUCHE INFERIEUR DU SOUSOBJET NUELEM
         NOSO      = NPX  + NPY * NX1
         NOSOEL(1) = NOSO
         NOSOEL(2) = NOSO + 1
         NOSOEL(3) = NOSO + 1 + NX1
         NOSOEL(4) = NOSO + NX1
         GOTO 100

C        TETRAEDRE STRUCTURE
C        ===================
  50     NBELEM = NX * NX * NX
         NB0    = 0
         IF (NUELEM.GT.NBELEM) GOTO 9900
         DO NRC=1,NX
            NCUBE = NRC*NRC*NRC
            IF(NCUBE.GE.NUELEM) THEN
               NB0 = NRC
               GO TO 501
            ENDIF
         ENDDO

 501     NBVO = (NB0-1)*(NB0-1)*(NB0-1)
         NSPH = (NB0-1)*NB0*(NB0+1) / 6
         NSPB = NB0*(NB0+1)*(NB0+2) / 6
C        RECHERCHE DU NUMERO DE L'ELEMENT
C        ** 1) LES ELEMENTS AYANT LA BASE EN BAS
         DO NB1=1,NB0
            DO NB2=1,NB1
C              LE NUMERO DU VOLUME
               NBVO = NBVO + 1
               IF(NUELEM.EQ.NBVO) THEN
C                 LES SOMMETS
                  NOSOEL(1) = NB1*(NB1-1)/2+NB2+NSPH
                  NOSOEL(2) = NB1*(NB1-1)/2+NB2+NSPB
                  NOSOEL(3) = NB1*(NB1+1)/2+NB2+NSPB
                  NOSOEL(4) = NB1*(NB1+1)/2+NB2+1+NSPB
                  NOSOEL(5) = 0
                  NOSOEL(6) = 0
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
            ENDDO
         ENDDO
C        ** 2) LES ELEMENTS INTERMEDIAIRES (BASE EN HAUT OU EN BAS)
         DO NB1=1,NB0-1
            DO NB2=1,NB1
C              LES SOMMETS DU HAUT
               NSOM(1) = NB1*(NB1-1)/2+NB2+NSPH
               NSOM(2) = NB1*(NB1+1)/2+NB2+NSPH
               NSOM(3) = NB1*(NB1+1)/2+NB2+1+NSPH
C              LES SOMMETS DU BAS
               NSOM(4) = NB1*(NB1+1)/2+NB2+NSPB
               NSOM(5) = (NB1+1)*(NB1+2)/2+NB2+1+NSPB
               NSOM(6) = NB1*(NB1+1)/2+NB2+1+NSPB
C              LES 4 TETRAEDRES
               NBVO = NBVO + 1
               IF(NUELEM.EQ.NBVO) THEN
                  NOSOEL(1) = NSOM(1)
                  NOSOEL(2) = NSOM(2)
                  NOSOEL(3) = NSOM(3)
                  NOSOEL(4) = NSOM(4)
                  NOSOEL(5) = 0
                  NOSOEL(6) = 0
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
               NBVO = NBVO + 1
               IF(NUELEM.EQ.NBVO) THEN
                  NOSOEL(1) = NSOM(2)
                  NOSOEL(2) = NSOM(3)
                  NOSOEL(3) = NSOM(4)
                  NOSOEL(4) = NSOM(5)
                  NOSOEL(5) = 0
                  NOSOEL(6) = 0
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
               NBVO = NBVO + 1
               IF(NUELEM.EQ.NBVO) THEN
                  NOSOEL(1) = NSOM(3)
                  NOSOEL(2) = NSOM(4)
                  NOSOEL(3) = NSOM(5)
                  NOSOEL(4) = NSOM(6)
                  NOSOEL(5) = 0
                  NOSOEL(6) = 0
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
               NBVO = NBVO + 1
               IF(NUELEM.EQ.NBVO) THEN
                  NOSOEL(1) = NSOM(4)
                  NOSOEL(2) = NSOM(6)
                  NOSOEL(3) = NSOM(3)
                  NOSOEL(4) = NSOM(1)
                  NOSOEL(5) = 0
                  NOSOEL(6) = 0
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
            ENDDO
         ENDDO
C        ** 3) LES ELEMENTS AYANT LA BASE EN HAUT
         DO  NB1=1,NB0-2
            DO NB2=1,NB1
C              LE NUMERO DU VOLUME
               NBVO = NBVO + 1
               IF(NUELEM.EQ.NBVO) THEN
C              LES SOMMETS
                  NOSOEL(1) = NB1*(NB1+1)/2+NB2+NSPH
                  NOSOEL(2) = (NB1+1)*(NB1+2)/2+NB2+1+NSPH
                  NOSOEL(3) = NB1*(NB1+1)/2+NB2+1+NSPH
                  NOSOEL(4) = (NB1+1)*(NB1+2)/2+NB2+1+NSPB
                  NOSOEL(5) = 0
                  NOSOEL(6) = 0
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
            ENDDO
         ENDDO
         GOTO 100

C        PENTAEDRE STRUCTURE
C        ===================
  60     NBAT = NX
         NBAH = NY
         NBTR = NBAT * NBAT
         NSCH = (NBAT+1)*(NBAT+2)/2
         NBELEM = NBTR * NBAH
         IF (NUELEM.GT.NBELEM) GOTO 9900
         NUCO = NUELEM / NBTR
         NUEL = NUELEM - NUCO * NBTR
         IF (NUEL.EQ.0) THEN
            NUCO = NUCO - 1
            NUEL = NBTR
         END IF
         N1B  = INT( SQRT(NUEL*1D0) )
         NSPH = NSCH * NUCO
         NSPB = NSPH + NSCH
C        RECHERCHE DU NUMERO DE L'ELEMENT FINI
         DO N1=N1B,N1B+1
            NEL  = (N1-1)*(N1-1)
            NEL1 = NEL - 1
            I1   = (N1*N1-N1)/2 + 1
            I2   = (N1*N1+N1)/2 + 1
            DO N2=1,N1
               NEL1 = NEL1 + 2
               IF ( NEL1 .EQ. NUEL ) THEN
                  NOSOEL(1) = I1     + NSPH
                  NOSOEL(2) = I2     + NSPH
                  NOSOEL(3) = I2 + 1 + NSPH
                  NOSOEL(4) = I1     + NSPB
                  NOSOEL(5) = I2     + NSPB
                  NOSOEL(6) = I2 + 1 + NSPB
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
               I1 = I1 + 1
               I2 = I2 + 1
            ENDDO
            NEL2 = NEL
            I1   = (N1*N1+N1)/2 + 1
            I2   = (N1*N1-N1)/2 + 1
            DO N2=1,N1-1
               NEL2 = NEL2 + 2
               IF ( NEL2 .EQ. NUEL ) THEN
                  NOSOEL(1) = I1 + 1 + NSPH
                  NOSOEL(2) = I2 + 1 + NSPH
                  NOSOEL(3) = I2     + NSPH
                  NOSOEL(4) = I1 + 1 + NSPB
                  NOSOEL(5) = I2 + 1 + NSPB
                  NOSOEL(6) = I2     + NSPB
                  NOSOEL(7) = 0
                  NOSOEL(8) = 0
                  GOTO 100
               END IF
               I1 = I1 + 1
               I2 = I2 + 1
            ENDDO
         ENDDO
         GOTO 100

C        HEXAEDRE STRUCTURE
C        ==================
  70     NBELEM = NX * NY * NZ
         IF( NUELEM .GT. NBELEM ) GOTO 9900
         NXY = NX * NY
         NUE = NUELEM - 1
         IF( NX .GT. 1 .AND. NY .GT. 1 ) THEN
            NPZ = NUE / NXY
            NPY = ( NUE - NPZ * NXY ) / NX
            NPX = NUE - NPZ * NXY - NPY * NX
         ELSE IF( NX .EQ. 1 ) THEN
C           LA BASE DOIT ETRE >1 POUR QUE LE CALCUL PRECEDENT SOIT CORRECT
            NPZ = NUE / NY
            NPY = NUE - NPZ * NY
            NPX = 0
         ELSE IF( NY .EQ. 1 ) THEN
            NPZ = NUE / NX
            NPY = 0
            NPX = NUE - NPZ * NX
         ELSE
            NPZ = 0
            NPY = 0
            NPX = 0
         ENDIF
C        NUMERO DU SOMMET GAUCHE INFERIEUR DU SOUSOBJET NUELEM
         NX1       = NX + 1
         NXY1      = NX1 * ( NY + 1 )
         NOSO      = NPX + 1 + NPY * NX1 + NPZ * NXY1
         NOSOEL(1) = NOSO
         NOSOEL(2) = NOSO + 1
         NOSOEL(3) = NOSO + 1 + NX1
         NOSOEL(4) = NOSO + NX1
C        LA COUCHE DU DESSUS
         NOSO      = NOSO + NXY1
         NOSOEL(5) = NOSO
         NOSOEL(6) = NOSO + 1
         NOSOEL(7) = NOSO + 1 + NX1
         NOSOEL(8) = NOSO + NX1
         GOTO 100

C        6-CUBE STRUCTURE
C        ================
C        NX      NOMBRE D'ARETES DANS CHAQUE COORDONNEE
C        NY=n=6  NOMBRE DE COORDONNEES DU 6-CUBE
 80      NBELEM = NX ** NY
         IF (NUELEM .GT. NBELEM) THEN
            IERR = 1
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:10),'(I10)') NUELEM
            KERR(1) ='NUMERO D''EF '//KERR(MXLGER)(1:10)
            WRITE(KERR(MXLGER)(1:10),'(I10)') NBELEM
            KERR(2) =' AU DELA DE '//KERR(MXLGER)(1:10)
            CALL LEREUR
            RETURN
         ENDIF

C        CALCUL DU NO DE L'EF NUELEM SELON LE NIVEAU DES COORDONNEES
C        NUELEM = NIVO(6) * NX**5 + NIVO(5) * NX**4 + ... + NIVO(1) * NX**0
C        NOSO = NO DU PREMIER SOMMET DU 6-CUBE NUELEM
C        NOSO = NIVO(6) * NX1**5 + NIVO(5) * NX1**4 + ... + NIVO(1)* NX1**0
         NOSO = 0
         N    = NUELEM - 1
         NX1  = NX + 1
         NX12 = NX1  * NX1
         NX13 = NX12 * NX1
         NX14 = NX13 * NX1
         NX15 = NX14 * NX1
         N21  = NX15 * NX1
         N2   = NX  ** 6
         DO I=1,5
            N2   = N2  / NX
            NIVO = N   / N2
            N21  = N21 / NX1
            NOSO = NOSO + NIVO * N21
            N    = N    - NIVO * N2
         ENDDO
         NOSO = NOSO + N + 1

C        NOSO = NO DU PREMIER SOMMET DU 6-CUBE NUELEM
C        LES AUTRES NO DE SOMMETS SE DEDUISENT PAR TRANSLATION DE NX1**j
         NU   = 0
         DO N=0,1
            DO M=0,1
               DO L=0,1
                  DO K=0,1
                     DO J=0,1
                        DO I=0,1
                           NU = NU + 1
                           NOSOEL(NU) = NOSO
     %                                + I        + J * NX1  + K * NX12
     %                                + L * NX13 + M * NX14 + N * NX15
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         GOTO 100

      ELSE

C        TYPE INCORRECT DU MAILLAGE
         IERR = 1
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NUTYMA
         KERR(1) ='TYPE MAILLAGE NUTYMA INCORRECT '//KERR(MXLGER)(1:10)
         CALL LEREUR
         RETURN

      ENDIF

C     LES NUMEROS DES NBTGEF TANGENTES AUX NBSOEF SOMMETS
C     ===================================================
C     LE NUMERO DU CODE GEOMETRIQUE DE CET EVENTUEL EF A TG
 100  NUGEEF = 0
      NUEFTG = 0
      IF( NBTGEF .GT. 0 ) THEN
C        CET EF NUELEM EST IL UN EF A TG DE NUMERO NUEFTG?
         NUEFTG = MCN( MNNSEF + LDAPEF + NUELEM - 1 )
         IF( NUEFTG .GT. 0 ) THEN
C           OUI: RECUPERATION DES NUMEROS DES TANGENTES PAR SOMMETS
            MN = MNNSEF + LDTGEF + NBTGEF * (NUEFTG-1) - 1
            DO I=1,NBTGEF
               NOSOEL( NBSOEF + I ) = MCN( MN + I )
            ENDDO
C           OUI: RECUPERATION DU NUMERO DU CODE GEOMETRIQUE DE CET EF A TG
            NUGEEF = MCN( MNNSEF + LDNGEF - 1 + NUEFTG )
         ELSE
C           PAS DE TG => NUMEROS NULS
            DO I=1,NBTGEF
               NOSOEL( NBSOEF + I ) = 0
            ENDDO
         ENDIF
      ENDIF
      RETURN

C     ERREUR DETECTEE
 9900 IERR = 1
      NBLGRC(NRERR) = 4
      WRITE(KERR(MXLGER)( 1:10),'(I10)') NUTYMA
      WRITE(KERR(MXLGER)(11:20),'(I10)') NUELEM
      WRITE(KERR(MXLGER)(21:30),'(I10)') NBELEM
      WRITE(KERR(MXLGER)(31:40),'(I10)') NX
      WRITE(KERR(MXLGER)(41:50),'(I10)') NY
      WRITE(KERR(MXLGER)(51:60),'(I10)') NZ
      IF( LANGAG .EQ. 0 ) THEN
        KERR(1)='NSEFNS: MAILLAGE STRUCTURE DE TYPE'//KERR(MXLGER)(1:10)
         KERR(2)='LE NUMERO D''EF '//KERR(MXLGER)(11:20)
       KERR(3)='EST AU DELA DU NOMBRE TOTAL D''EF '//KERR(MXLGER)(21:30)
      ELSE
         KERR(1)='NSEFNS: STRUCTURED MESH OF TYPE '//KERR(MXLGER)(1:10)
         KERR(2)='THE FINITE ELEMENT NUMBER '//KERR(MXLGER)(11:20)
         KERR(3)='IS OVER THE TOTAL NUMBER OF FE '//KERR(MXLGER)(21:30)
      ENDIF
      KERR(4)='NX='//KERR(MXLGER)(31:40)//'  NY='//KERR(MXLGER)(41:50)
     %    //'  NZ='//KERR(MXLGER)(51:60)
      CALL LEREUR

      RETURN
      END
