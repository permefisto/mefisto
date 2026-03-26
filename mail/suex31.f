        SUBROUTINE SUEX31( NTLXSU, LADEFI,
     %                     NTFATR, MNFATR, NTSOTR, MNSOTR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     GENERER LA TRIANGULATION D'UNE QUADRANGULATION
C -----    (LES QUADRANGLES SONT DECOUPES EN 2 TRIANGLES ET
C           LES EVENTUELS TRIANGLES PRESENTS DANS LA QUADRANGULATION
C           RESTENT INCHANGES)
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DE LA SURFACE
C LADEFI : TABLEAU DE DEFINITION DE LA SURFACE PARTITIONNEE
C          CF '~td/d/a_surface__definition'
C
C SORTIES:
C --------
C NTFATR : NUMERO      DU TMS 'NSEF' DES NUMEROS DES TRIANGLES
C MNFATR : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES TRIANGLES
C          CF '~TD/D/A___NSEF'
C NTSOTR : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOTR : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~TD/D/A___XYZSOMMET'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LABORATOIRE ANALYSE NUMERIQUE UPMC  NOVEMBRE 1996
C2345X7.................................................................
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*), NOSOEL(12),NOS(3)
      REAL              XYZ12(3), XYZ13(3)
      EQUIVALENCE      (XYZ12(1),NOS(1))
      CHARACTER*24      KNOMQU
C
      IERR = 0
C
C     NUMERO DE LA SURFACE QUADRANGULEE
C     ==================================
      NUSUQU = LADEFI(WUSUQU)
      IF( NUSUQU .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'suex31: SURFACE INCONNUE'
         ELSE
            KERR(1) = 'suex31: UNKNOWN SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     RECUPERATION DES TABLEAUX SOMMETS ET NSEF DE LA QUADRANGULATION
C     ===============================================================
      CALL LXNLOU( NTSURF, NUSUQU, NTLXQU, MNLXQU )
      IF( NTLXQU .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  'QUADRANGULATION INCONNUE'
         ELSE
            KERR(1) =  'UNKNOWN QUADRANGULATION'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE NOM DE LA SURFACE QUADRANGULEE
      CALL NMOBNU( 'SURFACE', NUSUQU, KNOMQU )
C
C     LE TABLEAU 'XYZSOMMET'
      CALL LXTSOU( NTLXQU, 'XYZSOMMET', NTSOQU, MNSOQU )
      IF( NTSOQU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'SURFACE ' // KNOMQU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE SANS TMS XYZSOMMET'
         ELSE
            KERR(2) = 'SURFACE WITHOUT XYZSOMMET TMS'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS
      NBSOM = MCN( MNSOQU + WNBSOM )
C     LE NOMBRE DE TANGENTES
      NBTGS = MCN( MNSOQU + WNBTGS )
C
C     LE TABLEAU 'NSEF'
      CALL LXTSOU( NTLXQU, 'NSEF', NTSSQU, MNFAQU )
      IF( NTSSQU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'SURFACE ' // KNOMQU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE SANS TMS NSEF'
         ELSE
            KERR(2) = 'SURFACE WITHOUT NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LES CARACTERISTIQUES DES NSEF DE CETTE QUADRANGULATION
      CALL NSEFPA( MCN(MNFAQU),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBSSQU,
     %             NX    , NY    , NZ    ,
     %             IERR   )
C     NUTYMA : 'NUMERO DE TYPE DU MAILLAGE'    ENTIER
C              0 : 'NON STRUCTURE'      , 2 : 'SEGMENT    STRUCTURE',
C              3 : 'TRIANGLE  STRUCTURE', 4 : 'QUADRANGLE STRUCTURE',
C              5 : 'TETRAEDRE STRUCTURE', 6 : 'PENTAEDRE  STRUCTURE',
C              7 : 'HEXAEDRE  STRUCTURE'
C     NBSOEL : NOMBRE DE SOMMETS DES NSEF
C              0 SI MAILLAGE NON STRUCTURE
C     NBSOEF : NOMBRE DE SOMMETS DE STOCKAGE DES NSEF
C              ( TRIANGLE NBSOEL=3  NBSOEF=4 )
C     NBSSQU : NOMBRE DE NSEF DU MAILLAGE ICI QUADRANGLES
C     NX, NY, NZ : LE NOMBRE D'ARETES DANS LES DIRECTION X Y Z
C                CF LE TMS ~/TD/D/A___NSEF
C
      IF( NUTYMA .NE. 0 .AND. NUTYMA .NE. 4 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'SURFACE ' // KNOMQU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'MAILLAGE DIFFERENT D''UNE QUADRANGULATION'
         ELSE
            KERR(2) = 'THE MESH IS NOT A QUADRANGULATION'
         ENDIF
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
C     CALCUL DU NOMBRE DE TRIANGLES ET QUADRANGLES AVEC TANGENTES
      NBTRIA = 0
      NBQUAD = 0
      NBTRTG = 0
      NBQUTG = 0
      DO 2 NUELEM = 1, NBSSQU
C        LE NUMERO DES SOMMETS DU QUADRANGLE NUELEM
         CALL NSEFNS( NUELEM, NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNFAQU, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( IERR .NE. 0 ) RETURN
         IF( NCOGEL .NE. 3 .AND. NCOGEL .NE. 4 ) THEN
C           ERREUR
            NBLGRC(NRERR) = 2
            KERR(1) = 'SURFACE ' // KNOMQU
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'EF DIFFERENT DE TRIANGLE OU QUADRANGLE'
            ELSE
               KERR(2) = 'FE DIFFERENT of  TRIANGLE or QUADRANGLE'
            ENDIF
            CALL LEREUR
            IERR = 5
            RETURN
         ENDIF
         IF( NCOGEL .EQ. 4 ) THEN
C           QUADRANGLE
            NBQUAD = NBQUAD + 1
            IF( NUEFTG .GT. 0 ) NBQUTG = NBQUTG + 1
         ELSE IF( NCOGEL .EQ. 3 ) THEN
C           TRIANGLE
            NBTRIA = NBTRIA + 1
            IF( NUEFTG .GT. 0 ) NBTRTG = NBTRTG + 1
         ENDIF
 2    CONTINUE
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE LA TRIANGULATION
C     -------------------------------------------------------
      NBTGST = NBTGS + NBQUTG * 2
      L      = WYZSOM + 3 * ( NBSOM + NBTGST )
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', L )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOTR , MNSOTR )
C     COPIE DU TABLEAU 'XYZSOMMET' DE LA QUADRANGULATION DANS CELUI
C     DE LA TRIANGULATION
      CALL TRTATA( MCN(MNSOQU), MCN(MNSOTR), WYZSOM+3*(NBSOM+NBTGS) )
C     LES EVENTUELLES TANGENTES SUPPLEMENTAIRES SONT AJOUTEES PLUS LOIN
      MNXYZS = MNSOTR + WYZSOM
      MNXYZT = MNXYZS + 3 * NBSOM
C
C     CONSTRUCTION DU TABLEAU 'NSEF'
C     ------------------------------
      NBEFOB = NBTRIA + NBQUAD * 2
      NBEFTG = NBTRTG + NBQUTG * 2
      IF( NBEFTG .GT. 0 ) THEN
         NBEFAP = NBEFOB
      ELSE
         NBEFAP = 0
      ENDIF
      CALL LXTNDC( NTLXSU, 'NSEF', 'MOTS',
     S             WUSOEF + NBEFOB * 4 + NBEFAP + NBEFTG*(1+NBTGEF) )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFATR , MNFATR )
C
C     CALCUL EFFECTIF DES 2 TRIANGLES DE CHAQUE QUADRANGLE
C     ----------------------------------------------------
C     LE NOMBRE D'EF A TG
      NBEFT = 0
C     LES NOUVELLES TANGENTES
      NBTG  = NBTGS
      MNTR  = MNFATR + WUSOEF - 1
      MNAP  = MNTR + 4 * NBEFOB
      MNCG  = MNAP + NBEFAP
      MNTG  = MNCG + NBEFTG
      DO 100 NUELEM = 1, NBSSQU
C
C        LE NUMERO DES SOMMETS DU QUADRANGLE NUELEM
         CALL NSEFNS( NUELEM, NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNFAQU, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C
C        QUALITE DE L'EF
         CALL QUALEF( NCOGEL, NOSOEL, NBSOM,
     %                RMCN(MNXYZS), SURFVOLU, QUALIT, IERR )
         IF( QUALIT .LE. 0 ) GOTO 9000
C
         IF( NOSOEL(4) .EQ. 0 .AND. NOSOEL(3) .GT. 0 ) THEN
C
C           TRIANGLE RECOPIE INTEGRALEMENT
C           ------------------------------
            DO 5 I=1,4
               MCN( MNTR + I ) = NOSOEL( I )
 5          CONTINUE
            MNTR = MNTR + 4
C
            IF( NBTGS .GT. 0 ) THEN
C              LE MAILLAGE A DES TANGENTES
C              LE POINTEUR SUR LES EF A TG
               MNAP = MNAP + 1
               IF( NUEFTG .GT. 0 ) THEN
C                 EF A TG
                  NBEFT = NBEFT + 1
                  MCN( MNAP ) = NBEFT
                  MNCG  = MNCG  + 1
                  MCN( MNCG ) = NUGEEF
                  DO 6 I=1,8
                     MCN( MNTG + I ) = NOSOEL( 4 + I )
 6                CONTINUE
                  MNTG = MNTG + 8
               ELSE
C                 EF SANS TG
                  MCN( MNAP ) = 0
               ENDIF
            ENDIF
            GOTO 100
         ENDIF
C
C        DECOUPAGE SELON LE CODE DE DECOUPAGE
         GOTO ( 10, 20, 30, 40 ), LADEFI( WDECOU )
C
C        TRIANGLE 123 ET 134
C        -------------------
 10      MCN( MNTR + 1 ) = NOSOEL( 1 )
         MCN( MNTR + 2 ) = NOSOEL( 2 )
         MCN( MNTR + 3 ) = NOSOEL( 3 )
         MCN( MNTR + 4 ) = 0
C
         MCN( MNTR + 5 ) = NOSOEL( 1 )
         MCN( MNTR + 6 ) = NOSOEL( 3 )
         MCN( MNTR + 7 ) = NOSOEL( 4 )
         MCN( MNTR + 8 ) = 0
         MNTR = MNTR + 8
C
         IF( NBTGS .GT. 0 ) THEN
C           LE MAILLAGE A DES TANGENTES
C           LE POINTEUR SUR LES EF A TG
            MNAP = MNAP + 2
            IF( NUEFTG .GT. 0 ) THEN
C              EF A TG
               NBEFT = NBEFT + 2
C              LE POINTEUR SUR LES EF A TG
               MCN( MNAP-1 ) = NBEFT - 1
               MCN( MNAP   ) = NBEFT
C              LE CODE GEOMETRIQUE DE L'EF A TANGENTE
               MNCG  = MNCG  + 2
               MCN( MNCG-1 ) = NUGEEF
               MCN( MNCG   ) = NUGEEF
C
C              LE CALCUL DES 3 COMPOSANTES DE LA TANGENTE EN S1
C              POUR LA DIRECTION S1S3
               CALL DEDVSR( MNXYZS, MNXYZT, NOSOEL,
     %                      0.0, 0.0, 1.0, 1.0,
     %                      RMCN(MNXYZT+3*NBTG) )
               NBTG  = NBTG + 1
               NBTG1 = NBTG
C              LE CALCUL DES 3 COMPOSANTES DE LA TANGENTE EN S3
C              POUR LA DIRECTION S3S1
               CALL DEDVSR( MNXYZS, MNXYZT, NOSOEL,
     %                      1.0, 1.0, -1.0, -1.0,
     %                      RMCN(MNXYZT+3*NBTG) )
               NBTG = NBTG + 1
C              LES NUMEROS DES TANGENTES DU TRIANGLE 123
               MCN( MNTG + 1 ) =  NOSOEL( 5 )
               MCN( MNTG + 2 ) =  NBTG1
               MCN( MNTG + 3 ) =  NOSOEL( 7 )
               MCN( MNTG + 4 ) =  NOSOEL( 8 )
               MCN( MNTG + 5 ) =  NBTG
               MCN( MNTG + 6 ) =  NOSOEL( 10 )
               MCN( MNTG + 7 ) =  0
               MCN( MNTG + 8 ) =  0
               MNTG = MNTG + 8
C              LES NUMEROS DES TANGENTES DU TRIANGLE 134
               MCN( MNTG + 1 ) =  NBTG1
               MCN( MNTG + 2 ) =  NOSOEL( 6 )
               MCN( MNTG + 3 ) =  NOSOEL( 9 )
               MCN( MNTG + 4 ) =  NBTG
               MCN( MNTG + 5 ) =  NOSOEL( 11 )
               MCN( MNTG + 6 ) =  NOSOEL( 12 )
               MCN( MNTG + 7 ) =  0
               MCN( MNTG + 8 ) =  0
               MNTG = MNTG + 8
            ELSE
C              EF SANS TG
               MCN( MNAP-1 ) = 0
               MCN( MNAP   ) = 0
            ENDIF
         ENDIF
         GOTO 100
C
C        TRIANGLE 124 ET 234
C        -------------------
 20      MCN( MNTR + 1 ) = NOSOEL( 1 )
         MCN( MNTR + 2 ) = NOSOEL( 2 )
         MCN( MNTR + 3 ) = NOSOEL( 4 )
         MCN( MNTR + 4 ) = 0
C
         MCN( MNTR + 5 ) = NOSOEL( 2 )
         MCN( MNTR + 6 ) = NOSOEL( 3 )
         MCN( MNTR + 7 ) = NOSOEL( 4 )
         MCN( MNTR + 8 ) = 0
         MNTR = MNTR + 8
C
         IF( NBTGS .GT. 0 ) THEN
C           LE MAILLAGE A DES TANGENTES
C           LE POINTEUR SUR LES EF A TG
            MNAP = MNAP + 2
            IF( NUEFTG .GT. 0 ) THEN
C              EF A TG
               NBEFT = NBEFT + 2
C              LE POINTEUR SUR LES EF A TG
               MCN( MNAP-1 ) = NBEFT - 1
               MCN( MNAP   ) = NBEFT
C              LE CODE GEOMETRIQUE DE L'EF A TANGENTE
               MNCG  = MNCG  + 2
               MCN( MNCG-1 ) = NUGEEF
               MCN( MNCG   ) = NUGEEF
C
C              LE CALCUL DES 3 COMPOSANTES DE LA TANGENTE EN S2
C              POUR LA DIRECTION S2S4
               CALL DEDVSR( MNXYZS, MNXYZT, NOSOEL,
     %                      1.0, 0.0, -1.0, 1.0,
     %                      RMCN(MNXYZT+3*NBTG) )
               NBTG  = NBTG + 1
               NBTG1 = NBTG
C              LE CALCUL DES 3 COMPOSANTES DE LA TANGENTE EN S4
C              POUR LA DIRECTION S4S2
               CALL DEDVSR( MNXYZS, MNXYZT, NOSOEL,
     %                      0.0, 1.0, 1.0, -1.0,
     %                      RMCN(MNXYZT+3*NBTG) )
               NBTG = NBTG + 1
C              LES NUMEROS DES TANGENTES DU TRIANGLE 124
               MCN( MNTG + 1 ) =  NOSOEL( 5 )
               MCN( MNTG + 2 ) =  NOSOEL( 6 )
               MCN( MNTG + 3 ) =  NBTG1
               MCN( MNTG + 4 ) =  NOSOEL( 8 )
               MCN( MNTG + 5 ) =  NOSOEL( 11 )
               MCN( MNTG + 6 ) =  NBTG
               MCN( MNTG + 7 ) =  0
               MCN( MNTG + 8 ) =  0
               MNTG = MNTG + 8
C              LES NUMEROS DES TANGENTES DU TRIANGLE 234
               MCN( MNTG + 1 ) =  NOSOEL( 7 )
               MCN( MNTG + 2 ) =  NBTG1
               MCN( MNTG + 3 ) =  NOSOEL( 9 )
               MCN( MNTG + 4 ) =  NOSOEL( 10 )
               MCN( MNTG + 5 ) =  NBTG
               MCN( MNTG + 6 ) =  NOSOEL( 12 )
               MCN( MNTG + 7 ) =  0
               MCN( MNTG + 8 ) =  0
               MNTG = MNTG + 8
            ELSE
C              EF SANS TG
               MCN( MNAP-1 ) = 0
               MCN( MNAP   ) = 0
            ENDIF
         ENDIF
         GOTO 100
C
C        DECOUPAGE SELON LE PLUS PETIT COSINUS DE L'ANGLE FORME
C        ENTRE 3 POINTS CONSECUTIFS  ( COS PEUT ETRE <0 )
 30      COSMIN = 2.0
         IMIN   = 0
         DO 38 I=1,4
C           LE POINT SUIVANT
            IF( I .LT. 4 ) THEN
               IP1 = I + 1
            ELSE
               IP1 = 1
            ENDIF
C           LE POINT PRECEDENT
            IF( I .GT. 1 ) THEN
               IM1 = I - 1
            ELSE
               IM1 = 4
            ENDIF
            MNXY  = MNSOTR + WYZSOM - 4
            MNS1  = MNXY + NOSOEL( I ) * 3
            MNS2  = MNXY + NOSOEL(IP1) * 3
            MNS3  = MNXY + NOSOEL(IM1) * 3
            DO 35 J=1,3
               XYZ12(J) = RMCN(MNS2+J) - RMCN(MNS1+J)
               XYZ13(J) = RMCN(MNS3+J) - RMCN(MNS1+J)
 35         CONTINUE
C           LE COSINUS DE L'ANGLE (I,I+1) (I,I-1)
            A =  XYZ12(1)**2 + XYZ12(2)**2 + XYZ12(3)**2
            B =  XYZ13(1)**2 + XYZ13(2)**2 + XYZ13(3)**2
            IF( A .LE. 0.0 ) GOTO 9000
            IF( B .LE. 0.0 ) GOTO 9000
            COSM = ( XYZ12(1) * XYZ13(1)
     %           +   XYZ12(2) * XYZ13(2)
     %           +   XYZ12(3) * XYZ13(3) ) / SQRT( A * B )
            IF( COSM .LT. COSMIN ) THEN
               IMIN   = I
               COSMIN = COSM
            ENDIF
 38      CONTINUE
         IF( IMIN .EQ. 1 .OR. IMIN .EQ. 3 ) THEN
C           LES ANGLES 1 ET 3 SONT DECOUPES
            GOTO 10
         ELSE
C           LES ANGLES 2 ET 4 SONT DECOUPES
            GOTO 20
         ENDIF
C
C        MAX DES QUALITES MINIMALES DES TRIANGLES
C        QUALITE DU TRIANGLE 123
 40      CALL QUALEF( 3, NOSOEL, NBSOM,
     %                RMCN(MNXYZS), S123, Q123, IERR )
C        QUALITE DU TRIANGLE 234
         CALL QUALEF( 3, NOSOEL(2), NBSOM,
     %                RMCN(MNXYZS), S234, Q234, IERR )
C        QUALITE DU TRIANGLE 134
         NOS(1) = NOSOEL(1)
         NOS(2) = NOSOEL(3)
         NOS(3) = NOSOEL(4)
         CALL QUALEF( 3, NOS, NBSOM,
     %                RMCN(MNXYZS), S134, Q134, IERR )
C        QUALITE DU TRIANGLE 124
         NOS(1) = NOSOEL(1)
         NOS(2) = NOSOEL(2)
         NOS(3) = NOSOEL(4)
         CALL QUALEF( 3, NOS, NBSOM,
     %                RMCN(MNXYZS), S124, Q124, IERR )
C
C        LE MAXIMUM DES MINIMUM DES QUALITES DES TRIANGLES EST RETENU
         IF( MIN(Q123,Q134) .GE. MIN(Q124,Q234) ) THEN
C           LES TRIANGLES 123 134 SONT FORMES
            GOTO 10
         ELSE
C           LES TRIANGLES 124 234 SONT FORMES
            GOTO 20
         ENDIF
C
 100  CONTINUE
C
C     MISE A JOUR DU TABLEAU 'NSEF' DE CETTE TRIANGULATION
C     ----------------------------------------------------
C     TYPE DE L'OBJET : SURFACE
      MCN( MNFATR + WUTYOB ) = 3
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNFATR + WUTFMA ) = MCN( MNFAQU + WUTFMA )
C     LES EF A TANGENTES STOCKEES
      MCN( MNFATR + WBTGEF ) = NBTGEF
      MCN( MNFATR + WBEFAP ) = NBEFAP
      MCN( MNFATR + WBEFTG ) = NBEFTG
C     NUMERO DU TYPE DU MAILLAGE : NON STRUCTURE
      MCN( MNFATR + WUTYMA ) = 0
C     NOMBRE DE SOMMETS PAR NSEF
      MCN( MNFATR + WBSOEF ) = 4
C     NOMBRE DE TRIANGLES
      MCN( MNFATR + WBEFOB ) = NBEFOB
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFATR) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFATR + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     MISE A JOUR DU TABLEAU 'XYZSOMMET'
C     ----------------------------------
C     LE NOMBRE DE SOMMETS DE LA TRIANGULATION
      MCN( MNSOTR + WNBSOM ) = NBSOM
C     LE NOMBRE DE TANGENTES DE LA TRIANGULATION
      MCN( MNSOTR + WNBTGS ) = NBTG
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOTR + WBCOOR ) = 3
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOTR) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOTR + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
      RETURN
C
C     ERREUR
 9000 NBLGRC(NRERR) = 3
      KERR(1) = 'SURFACE ' // KNOMQU
      WRITE(KERR(MXLGER)(1:6),'(I6)') NUELEM
      IF( LANGAG .EQ. 0 ) THEN
         IF( NOSOEL(4) .NE. 0 ) THEN
            KERR(2) = 'QUADRANGLE' // KERR(MXLGER)(1:6) // ' DEGENERE'
         ELSE
            KERR(2) = 'TRIANGLE' // KERR(MXLGER)(1:6) // ' DEGENERE'
         ENDIF
      ELSE
         IF( NOSOEL(4) .NE. 0 ) THEN
            KERR(2) = 'QUADRANGLE' // KERR(MXLGER)(1:6) //' DEGENERATED'
         ELSE
            KERR(2) = 'TRIANGLE' // KERR(MXLGER)(1:6) //' DEGENERATED'
         ENDIF
      ENDIF
      DO 9010 I=1,4
         WRITE(KERR(MXLGER)(6*I-5:6*I),'(I6)') NOSOEL(I)
 9010 CONTINUE
      IF( LANGAG .EQ. 0 ) THEN
         KERR(3) = 'SOMMETS :' // KERR(MXLGER)(1:24)
      ELSE
         KERR(3) = 'VERTICES :' // KERR(MXLGER)(1:24)
      ENDIF
      CALL LEREUR
      IERR = 1
      END
