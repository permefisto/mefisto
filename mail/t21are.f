      SUBROUTINE T21ARE( NMOBJT, NUOBJT, MNNSEF, MNSOMM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DU TABLEAU 'NSEF' D'ADRESSE MNNSEF
C -----    EN 1D ou 2D
C
C ENTREES:
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DE L'OBJET
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOMM : ADRESSE MCN DU TABLEAU 'XYZSOMMET'    A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1988
C ......................................................................
      PARAMETER        (CMDECA=0.3,EFFACE=0.03)
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*(*)     NMOBJT
      CHARACTER*10      NMSOMM
      CHARACTER*11      KSYMBO
      CHARACTER*40      KNLIGN
      REAL              XYZ(1:2,1:2),XYTG(1:2,1:2),XYZP(1:2)
C
C     LE TYPE DE L'OBJET
      N = MCN( MNNSEF + WUTYOB )
      IF( N .NE. 2 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = NMOBJT// 'N''EST PAS UNE LIGNE'
         ELSE
            KERR(1) = NMOBJT// 'IS NOT A LINE'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TYPE DU MAILLAGE
      NUTYMA = MCN( MNNSEF + WUTYMA )
C
C     LE NOMBRE DE TGS PAR EF
      NBTGEF = MCN( MNNSEF + WBTGEF )
C
C     LE NOMBRE D'EF A TG
      NBEFTG = MCN( MNNSEF + WBEFTG )
C
C     LE NOMBRE DE POINTEURS SUR LES EF A TG
      NBEFAP = MCN( MNNSEF + WBEFAP )
C
C     ADRESSE DU NUMERO DU 1-ER SOMMET DE LA 1-ERE ARETE
      MNSS = MNNSEF + WUSOEF
C
      IF( NUTYMA .EQ. 0 ) THEN
C
C        LIGNE NON STRUCTUREE( LA LIGNE PEUT ETRE EN PLUSIEURS MORCEAUX)
C        NOMBRE D'ARETES DE LA LIGNE
         NBARLI = MCN( MNNSEF + WBEFOB )
C        ADRESSE DU NUMERO DU PREMIER POINTEUR SUR LES EF A TG
         MNTT = MNSS + 2 * NBARLI
C
      ELSE IF( NUTYMA .EQ. 2 ) THEN
C
C        SEGMENT STRUCTURE
C        NOMBRE D'ARETES
         NBARLI = MCN( MNNSEF + WBARSE )
C        ADRESSE DU NUMERO DU PREMIER POINTEUR SUR LES EF A TG
         MNTT = MNNSEF + WBARSE + 1
C
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYMA
         KERR(1) ='T21ARE:TYPE INCORRECT DE NSEF'
     %           //KERR(MXLGER)(1:4)
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE NOMBRE DE TGS STOCKEES DANS XYZSOMMET
      NBTGS = MCN( MNSOMM + WNBTGS )
      IF( NBEFTG .GT. 0 .AND. NBTGS .GT. 0 .AND. NBTGEF .GT. 0 ) THEN
C        LES TANGENTES SONT A PRENDRE EN COMPTE
         INDTGS = 1
      ELSE
C        LES TANGENTES NE SONT PAS A PRENDRE EN COMPTE
         INDTGS = 0
      ENDIF
C
C     LE DECALAGE POUR ATTEINDRE LES EF A TG
      LDEFTG = MNTT + NBEFAP + NBEFTG
C
C     REDUCTION DES ARETES
      REDUCA = PREDUA * 0.005
C
C     VARIABLE AUXILIAIRE POUR POSITIONNER LA POIGNEE DE LA LIGNE
      NBARL2 = NBARLI / 2  + 1
C
C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS  = MNSOMM + WYZSOM - 3
C     ADRESSE-3 DE LA 1-ERE COMPOSANTE DE LA 1-ERE TG DU TMS 'XYZSOMMET'
      MNTG = MNS + 3 * MCN( MNSOMM + WNBSOM )
C
C     ORBITE OU NON?
      IF( LORBITE .GT. 0 ) THEN
C        INITIALISATION DU ZOOM DEPLACEMENT
         CALL ZOOM2D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
      ENDIF
C
C     TRACE DES AXES 2D
 10   CALL TRAXE2
C
C     BOUCLE SUR LES NBARLI ARETES DE LA LIGNE TRACEES EN CYAN OU BLANC
      MNSNS = MNSS
      DO 100 I=1,NBARLI
C
C        LE NUMERO DES 2 SOMMETS DE L'ARETE
         IF( NUTYMA .EQ. 0 ) THEN
C           LIGNE NON STRUCTUREE
            NS1 = MCN( MNSNS )
            NS2 = MCN( MNSNS + 1 )
            MNSNS = MNSNS + 2
         ELSE
            NS1 = I
            NS2 = I + 1
         ENDIF
C
C        LES 3 COORDONNEES DES 2 EXTREMITES DE LA I-EME ARETE
         MN = MNS + 3 * NS1
         XYZ(1,1) = RMCN( MN     )
         XYZ(2,1) = RMCN( MN + 1 )
C
         MN = MNS + 3 * NS2
         XYZ(1,2) = RMCN( MN     )
         XYZ(2,2) = RMCN( MN + 1 )
C
         IF( INDTGS .GT. 0 ) THEN
C           EF A TG?
            NUEFTG = MCN( MNTT - 1 + I )
            IF( NUEFTG .LE. 0 ) GOTO 40
C           OUI : NUMERO DES 2 TGS DE L'ARETE
            MN   = LDEFTG - 2 + 2 * NUEFTG
            NTG1 = MCN( MN   )
            NTG2 = MCN( MN+1 )
C           TRACE AVEC LES 2 TANGENTES
            IF( NTG1 .NE. 0 ) THEN
               MN = MNTG + 3 * ABS(NTG1)
               XYTG(1,1) = RMCN( MN     )
               XYTG(2,1) = RMCN( MN + 1 )
               IF( NTG1 .LT. 0 ) THEN
C                 TANGENTE A INVERSER
                  XYTG(1,1) = -XYTG(1,1)
                  XYTG(2,1) = -XYTG(2,1)
               ENDIF
            ELSE
               XYTG(1,1) = XYZ(1,2) - XYZ(1,1)
               XYTG(2,1) = XYZ(2,2) - XYZ(2,1)
            ENDIF
            IF( NTG2 .NE. 0 ) THEN
               MN = MNTG + 3 * ABS(NTG2)
               XYTG(1,2) = RMCN( MN     )
               XYTG(2,2) = RMCN( MN + 1 )
               IF( NTG2 .LT. 0 ) THEN
C                 TANGENTE A INVERSER
                  XYTG(1,2) = -XYTG(1,2)
                  XYTG(2,2) = -XYTG(2,2)
               ENDIF
            ELSE
               XYTG(1,2) = XYZ(1,1) - XYZ(1,2)
               XYTG(2,2) = XYZ(2,1) - XYZ(2,2)
            ENDIF
C
C           LE TRACE DE L'ARETE P3 ** 2
            CALL XVEPAISSEUR( NEPARL )
            CALL TRAR2D( NCOUAL, PREDUA,
     %                   XYZ(1,1),  XYZ(2,1),
     %                   XYZ(1,2),  XYZ(2,2),
     %                   XYTG(1,1), XYTG(2,1),
     %                   XYTG(1,2), XYTG(2,2) )
            GOTO 45
C
         ENDIF
C
C        LE TRACE DE L'ARETE DROITE
C        REDUCTION DE L'ARETE
 40      DIRX = ( XYZ(1,2) - XYZ(1,1) ) * REDUCA
         DIRY = ( XYZ(2,2) - XYZ(2,1) ) * REDUCA
         CALL XVEPAISSEUR( NEPARL )
         CALL TRAIT2D( NCOUAL, XYZ(1,1)+DIRX, XYZ(2,1)+DIRY,
     %                         XYZ(1,2)-DIRX, XYZ(2,2)-DIRY )
C
C        TRACE EVENTUEL DU NO DE L'EF
 45       IF( IAVNEF .NE. 0 ) THEN
            WRITE( NMSOMM, '(I7)' ) I
            KSYMBO = '.' // NMSOMM
            CALL SANSBL( KSYMBO, L )
            IF( INDTGS .GT. 0 ) THEN
C              LE POINT MILIEU EN COORDONNEES OBJET 2D
               X = VALP3H( 0.5, XYZ(1,1),XYZ(1,2),XYTG(1,1),XYTG(1,2) )
               Y = VALP3H( 0.5, XYZ(2,1),XYZ(2,2),XYTG(2,1),XYTG(2,2) )
            ELSE
               X = (XYZ(1,1)+XYZ(1,2)) * 0.5
               Y = (XYZ(2,1)+XYZ(2,2)) * 0.5
            ENDIF
            CALL SYMBOLE2D( NCONEF, X, Y, KSYMBO(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            WRITE( NMSOMM, '(I7)' ) NS1
            KSYMBO = '.' // NMSOMM
            CALL SANSBL( KSYMBO, L )
            CALL SYMBOLE2D( NCONSO, XYZ(1,1), XYZ(2,1), KSYMBO(1:L) )
            WRITE( NMSOMM, '(I7)' ) NS2
            KSYMBO = '.' // NMSOMM
            CALL SANSBL( KSYMBO, L )
            CALL SYMBOLE2D( NCONSO, XYZ(1,2), XYZ(2,2), KSYMBO(1:L) )
         ENDIF
C
C        LE TRACE DE LA 'POIGNEE DE LA LIGNE'
         IF( I .EQ. NBARL2 ) THEN
            IF( INDTGS .GT. 0 ) THEN
C              LE POINT MILIEU EN COORDONNEES OBJET 2D
               XYZP(1)=VALP3H(0.5,XYZ(1,1),XYZ(1,2),XYTG(1,1),XYTG(1,2))
               XYZP(2)=VALP3H(0.5,XYZ(2,1),XYZ(2,2),XYTG(2,1),XYTG(2,2))
            ELSE
               XYZP(1) = (XYZ(1,1)+XYZ(1,2)) * 0.5
               XYZP(2) = (XYZ(2,1)+XYZ(2,2)) * 0.5
            ENDIF
         ENDIF
 100  CONTINUE
C
C     TRACE DE LA POIGNEE SI ELLE VISIBLE SUR L'IMAGE
      CALL ITEML2( XYZP, NMOBJT, NUOBJT )
C
C     REPRISE DE TRANSLATION ZOOM SI LORBITE>0
      IF( LORBITE .GT. 0 ) THEN
C        ZOOM OU TRANSLATION ACTIFS
         IF( LANGAG .EQ. 0 ) THEN
            KNLIGN = 'LIGNE ' // NMOBJT
         ELSE
            KNLIGN = 'LINE ' // NMOBJT
         ENDIF
         CALL TRFINS( KNLIGN )
         CALL ZOOM2D1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
         GOTO 10
      ENDIF
C
 9000 CALL XVEPAISSEUR( 0 )
      RETURN
      END
