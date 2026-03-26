        SUBROUTINE SUEX32( NTLXSU, LADEFI,
     %                     NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   GENERER LA TRIANGULATION D'UNE SURFACE PAR JONCTION D'UN POINT
c -----   AUX 2 EXTREMITES DES ARETES D'UNE LIGNE
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE A CREER
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF '~/td/d/a_surface__definition'
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF '~/td/d/a___nsef'
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~/td/d/a___xyzsommet'
C IERR   : = 0 SI PAS D'ERREUR
C          <>0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1996
C23456...............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*), NOSOEL(4)
      REAL              XYZP(3),XYZ(3)
C
      IERR = 0
C
C     LE NOMBRE D'ARETES DANS LE SENS POINT LIGNE
C     -------------------------------------------
      NBARPL = LADEFI( WBARPL )
      IF( NBARPL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE D''ARETES <=0 INCORRECT'
         CALL LEREUR
         IERR = 1
         GOTO 9900
      ENDIF
C
C     RESTAURATION DU POINT A JOINDRE
C     -------------------------------
      CALL LXNLOU( NTPOIN, LADEFI(WUPTJN), NTPTJN, MN )
      IF( NTPTJN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'POINT INCONNU'
         CALL LEREUR
         IERR = 2
         GOTO 9900
      ENDIF
C     RESTAURATION DES TABLEAUX SOMMETS ET NSEF
      CALL LXTSOU( NTPTJN, 'XYZSOMMET', NTSOPT, MNSOPT )
      IF( NTSOPT .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'POINT SANS SOMMET'
         CALL LEREUR
         IERR = 3
         GOTO 9900
      ENDIF
C
C     RESTAURATION DU MAILLAGE DE LA LIGNE INITIALE
C     ---------------------------------------------
      CALL LXNLOU( NTLIGN, LADEFI(WULGJN), NTLGJN, MN )
      IF( NTLGJN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE INCONNUE'
         CALL LEREUR
         IERR = 4
         GOTO 9900
      ENDIF
C     RESTAURATION DES TABLEAUX SOMMETS ET NSEF
      CALL LXTSOU( NTLGJN, 'XYZSOMMET', NTSOLI, MNSOLI )
      IF( NTSOLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS SOMMETS'
         CALL LEREUR
         IERR = 5
         GOTO 9900
      ENDIF
      CALL LXTSOU( NTLGJN, 'NSEF', NTARLI, MNARLI )
      IF( NTARLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS NSEF'
         CALL LEREUR
         IERR = 6
         GOTO 9900
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DE LA LIGNE A DEPLACER
      CALL NSEFPA( MCN(MNARLI),
     %             NUTYML, NBSOEL, NBSOEF, NBTGEL,
     %             LDAPEF, LDNGEF, LDTGEF, NBARLI,
     %             NX  , NY  , NZ  ,
     %             IERR   )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE A RELIER AU POINT DE MAILLAGE INCORRECT'
         CALL LEREUR
         IERR = 7
         RETURN
      ENDIF
C
C     LE NOMBRE D'EF DE LA SURFACE
      NBEFOB = NBARPL * NBARLI
C
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBSOLI = MCN( MNSOLI + WNBSOM )
C     LE NOMBRE DE SOMMETS DE LA SURFACE
      NBSOSU = 1 + NBARPL * NBSOLI
C
C     LE NOMBRE DE TANGENTES DE LA LIGNE
      NBTGLI = MCN( MNSOLI + WNBTGS )
      NBTGSU = NBTGLI * NBARPL
      IF( NBTGLI .GT. 0 ) THEN
         NBTGEF = 8
         NBEFAP = NBEFOB
         NBEFTG = NBEFOB
      ELSE
         NBTGEF = 0
         NBEFAP = 0
         NBEFTG = 0
      ENDIF
C
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
C     -----------------------------------------------
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER',
     %             WUSOEF + 4*NBEFOB + NBEFAP + NBEFTG*(1+NBTGEF) )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU, MNFASU )
C
      MNF   = MNFASU + WUSOEF
      MNAP  = MNF  + 4 * NBEFOB
      MNCG  = MNAP + NBEFAP
      MNTG  = MNCG + NBEFTG - 1
      NS1   = 1 - NBSOLI
      NBEFT = 0
      NBTG1 =-NBTGLI
      DO 30 I=1,NBARPL
C
C        LA TRANCHE SUIVANTE
         NS0 = NS1
         NS1 = NS1 + NBSOLI
C
C        LE NOMBRE DE TANGENTES AVANT ET APRES LA TRANCHE I
         NBTG0 = NBTG1
         NBTG1 = NBTG1 + NBTGLI
C
C        LA BOUCLE SUR LES ARETES DU MAILLAGE DE LA LIGNE
         DO 10 N=1,NBARLI
C           LE NUMERO DES 2 SOMMETS DE L'ARETE N
            CALL NSEFNS( N    , NUTYML, NBSOEF, NBTGEL,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNARLI, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
            IF( NCOGEL .NE. 2 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
               KERR(1) = 'LIGNE AVEC UN ELEMENT FINI DE CODE'
     %                 // KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 8
               CALL LXTSDS( NTLXSU, 'NSEF' )
               RETURN
            ENDIF
C
            IF( I .EQ. 1 ) THEN
C
C              CREATION D'UN TRIANGLE
C              ----------------------
C              LES 2 SOMMETS DE L'ARETE PUIS LE POINT => LES TRIANGLES
               MCN( MNF     ) = 1
               MCN( MNF + 1 ) = 1 + NOSOEL(1)
               MCN( MNF + 2 ) = 1 + NOSOEL(2)
               MCN( MNF + 3 ) = 0
C
C              LES EVENTUELLES TANGENTES
               IF( NBTGLI .GT. 0 ) THEN
                  IF( NUEFTG .GT. 0 ) THEN
C                    ARETE A TG
C                    LE POINTEUR SUR L'EF A TG
                     NBEFT = NBEFT + 1
                     MCN( MNAP ) = NBEFT
                     MNAP = MNAP + 1
C                    LE CODE GEOMETRIQUE P3 STANDARD
                     MCN( MNCG ) = 0
                     MNCG = MNCG + 1
C                    LE NUMERO DES 8 TANGENTES DU QUADRANGLE
                     MCN( MNTG + 1 ) = 0
                     MCN( MNTG + 2 ) = 0
                     MCN( MNTG + 3 ) = NOSOEL(3)
                     MCN( MNTG + 4 ) = 0
                     MCN( MNTG + 5 ) = 0
                     MCN( MNTG + 6 ) = NOSOEL(4)
                     MCN( MNTG + 7 ) = 0
                     MCN( MNTG + 8 ) = 0
                     MNTG = MNTG + 8
                  ENDIF
               ENDIF
            ELSE
C
C              CREATION D'UN QUADRANGLE
C              ------------------------
C              LES 2 SOMMETS DE L'ARETE PUIS L'ARETE OPPOSEE
C              => LES QUADRANGLES
               MCN( MNF     ) = NOSOEL(1) + NS0
               MCN( MNF + 1 ) = NOSOEL(1) + NS1
               MCN( MNF + 2 ) = NOSOEL(2) + NS1
               MCN( MNF + 3 ) = NOSOEL(2) + NS0
C
C              LES EVENTUELLES TANGENTES
               IF( NBTGLI .GT. 0 ) THEN
                  IF( NUEFTG .GT. 0 ) THEN
C                    ARETE A TG
C                    LE POINTEUR SUR L'EF A TG
                     NBEFT = NBEFT + 1
                     MCN( MNAP ) = NBEFT
                     MNAP = MNAP + 1
C                    LE CODE GEOMETRIQUE P3 STANDARD
                     MCN( MNCG ) = 0
                     MNCG = MNCG + 1
C                    LE NUMERO DES 8 TANGENTES DU QUADRANGLE
                     IF( NOSOEL(3) .GE. 0 ) THEN
                        LESIG3 = 1
                     ELSE
                        LESIG3 = -1
                     ENDIF
                     IF( NOSOEL(4) .GE. 0 ) THEN
                        LESIG4 = 1
                     ELSE
                        LESIG4 = -1
                     ENDIF
                     MCN( MNTG + 1 ) = 0
                     MCN( MNTG + 2 ) = NOSOEL(3) + LESIG3 * NBTG0
                     MCN( MNTG + 3 ) = NOSOEL(3) + LESIG3 * NBTG1
                     MCN( MNTG + 4 ) = 0
                     MCN( MNTG + 5 ) = 0
                     MCN( MNTG + 6 ) = NOSOEL(4) + LESIG4 * NBTG1
                     MCN( MNTG + 7 ) = NOSOEL(4) + LESIG4 * NBTG0
                     MCN( MNTG + 8 ) = 0
                     MNTG = MNTG + 8
                  ENDIF
               ENDIF
            ENDIF
            MNF = MNF + 4
 10      CONTINUE
 30   CONTINUE
C
C     TYPE DE L'OBJET : SURFACE TRIANGULEE NON STRUCTUREE
      MCN( MNFASU + WUTYOB ) = 3
      MCN( MNFASU + WUTYMA ) = 0
      MCN( MNFASU + WBSOEF ) = 4
      MCN( MNFASU + WBEFOB ) = NBEFOB
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = 0
C     LES TANGENTES STOCKEES
      MCN( MNFASU + WBTGEF ) = NBTGEF
      MCN( MNFASU + WBEFAP ) = NBEFAP
      MCN( MNFASU + WBEFTG ) = NBEFTG
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
C     ----------------------------------------------------
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS',
     %             WYZSOM + 3 * ( NBSOSU + NBTGSU ) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOFA, MNSOFA )
C     AJOUT DES 3 COORDONNEES DU POINT A JOINDRE ( NUMERO 1 )
      XYZP(1) = RMCN( MNSOPT + WYZSOM )
      XYZP(2) = RMCN( MNSOPT + WYZSOM + 1 )
      XYZP(3) = RMCN( MNSOPT + WYZSOM + 2 )
      MNF = MNSOFA + WYZSOM
      RMCN( MNF     ) = XYZP(1)
      RMCN( MNF + 1 ) = XYZP(2)
      RMCN( MNF + 2 ) = XYZP(3)
C
C     CALCUL DES COORDONNEES DES SOMMETS SUR LES DROITES : POINT ---> SOMMET
      MNF = MNSOLI + WYZSOM
      DO 50 N=1,NBSOLI
C        LES COORDONNEES DU SOMMET DE LA LIGNE
         XYZ(1) = RMCN( MNF )
         XYZ(2) = RMCN( MNF + 1 )
         XYZ(3) = RMCN( MNF + 2 )
         DO 40 I=1,NBARPL-1
            MNA = MNSOFA + WYZSOM - 1 + 3 * ( NBSOLI * ( I - 1 ) + N )
            DO 45 K=1,3
               RMCN(MNA+K) = XYZP(K) + I * (XYZ(K) - XYZP(K)) / NBARPL
 45         CONTINUE
 40      CONTINUE
         MNF = MNF + 3
 50   CONTINUE
C
C     LA LIGNE ELLE MEME: LES NBSOLI SOMMETS
      MNA = MNSOFA + WYZSOM + 3 * ( 1 + NBSOLI * ( NBARPL - 1 ) )
      CALL TRTATA( MCN(MNSOLI+WYZSOM), MCN(MNA), 3*NBSOLI )
C
C     LES EVENTUELLES TANGENTES DE LA SURFACE
      IF( NBTGLI .GT. 0 ) THEN
         MNTG = MNSOFA + WYZSOM + 3 * NBSOSU
         DO 80 I=1, NBARPL
C           LE RAPPORT D'HOMOTHETIE
            RAPPOR = REAL( I ) / NBARPL
C           ADRESSE MCN DE LA PREMIERE COMPOSANTE DES TGS DE LA LIGNE
            MNF    = MNSOLI + WYZSOM + 3 * NBSOLI
            DO 70 N=1,NBTGLI
               DO 60 K=1,3
                  RMCN(MNTG  ) = RMCN(MNF  ) * RAPPOR
                  RMCN(MNTG+1) = RMCN(MNF+1) * RAPPOR
                  RMCN(MNTG+2) = RMCN(MNF+2) * RAPPOR
 60            CONTINUE
               MNTG = MNTG + 3
               MNF  = MNF  + 3
 70         CONTINUE
 80      CONTINUE
      ENDIF
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM) = NBSOSU
C     LE NOMBRE DE TANGENTES
      MCN( MNSOFA + WNBTGS) = NBTGSU
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOFA + WBCOOR ) = 3
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
 9900 RETURN
      END
