      SUBROUTINE T31ARE( NMOBJT, NUOBJT, MNNSEF, MNSOMM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES D'UNE LIGNE A PARTIR DU TMS 'NSEF'
C -----    D'ADRESSE MNNSEF ET DU TMS 'XYZSOMMET' D'ADRESSE MNSOMM
C
C ENTREES:
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DE L'OBJET
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOMM : ADRESSE MCN DU TABLEAU 'XYZSOMMET'    A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C ......................................................................
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
      CHARACTER*40      KNLIGN
      REAL              XYZ(1:3,1:2),XYZP(1:3),XYZN(1:3),
     %                  XYZTG(1:3,1:2)
C
C     LE TYPE DE L'OBJET
      N = MCN( MNNSEF + WUTYOB )
      IF( N .NE. 2 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = NMOBJT // ' N''EST PAS UNE LIGNE'
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
      MNSS   = MNNSEF + WUSOEF
C
      IF( NUTYMA .EQ. 0 ) THEN
C
C        LIGNE NON STRUCTUREE( LA LIGNE PEUT ETRE EN PLUSIEURS MORCEAUX)
C        NOMBRE D'ARETES DE LA LIGNE
         NBARLI = MCN( MNNSEF + WBEFOB )
C        ADRESSE DU NUMERO DE LA 1-ERE COMPOSANTE DE LA 1-ERE TANGENTE
         MNTT = MNSS + 2 * NBARLI
C
      ELSE IF( NUTYMA .EQ. 2 ) THEN
C
C        SEGMENT STRUCTURE
C        NOMBRE D'ARETES
         NBARLI = MCN( MNNSEF + WBARSE )
C        ADRESSE DU NUMERO DE LA 1-ERE COMPOSANTE DE LA 1-ERE TANGENTE
         MNTT = MNNSEF + WBARSE + 1
C
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYMA
         KERR(1) = 'T31ARE:TYPE STRUCTURE INCORRECT DANS NSEF '
     %           // KERR(MXLGER)(1:4)
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE NOMBRE DE TGS STOCKEES DANS XYZSOMMET
      NBTGS = MCN( MNSOMM + WNBTGS )
      IF( NBEFTG .GT. 0 .AND. NBTGS .GT. 0 .AND. NBTGEF  .GT. 0 ) THEN
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
      MNS = MNSOMM + WYZSOM - 3
C     ADRESSE-3 DE LA 1-ERE COMPOSANTE DE LA 1-ERE TG DU TMS 'XYZSOMMET'
      MNTG = MNS + 3 * MCN( MNSOMM + WNBSOM )
C
C     ORBITE OU NON?
      IF( LORBITE .GT. 0 ) THEN
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
      ENDIF
C
C     TRACE DES AXES 3D
 10   CALL TRAXE3
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
         XYZ(3,1) = RMCN( MN + 2 )
C
         MN = MNS + 3 * NS2
         XYZ(1,2) = RMCN( MN     )
         XYZ(2,2) = RMCN( MN + 1 )
         XYZ(3,2) = RMCN( MN + 2 )
C
         IF( INDTGS .GT. 0 ) THEN
C
C           EF A TG?
            NUEFTG = MCN( MNTT - 1 + I )
            IF( NUEFTG .LE. 0 ) GOTO 40
C
C           OUI : NUMERO DES 2 TGS DE L'ARETE
            MN   = LDEFTG - 2 + 2 * NUEFTG
            NTG1 = MCN( MN   )
            NTG2 = MCN( MN+1 )
C
C           LES 3 COMPOSANTES DE LA TANGENTE 1
            IF( NTG1 .NE. 0 ) THEN
               MN = MNTG + 3 * ABS(NTG1)
               IF( NTG1 .GT. 0 ) THEN
                  XYZTG(1,1) = RMCN( MN     )
                  XYZTG(2,1) = RMCN( MN + 1 )
                  XYZTG(3,1) = RMCN( MN + 2 )
               ELSE
C                 TANGENTE A INVERSER
                  XYZTG(1,1) =-RMCN( MN     )
                  XYZTG(2,1) =-RMCN( MN + 1 )
                  XYZTG(3,1) =-RMCN( MN + 2 )
               ENDIF
            ELSE
               XYZTG(1,1) = XYZ(1,2) - XYZ(1,1)
               XYZTG(2,1) = XYZ(2,2) - XYZ(2,1)
               XYZTG(3,1) = XYZ(3,2) - XYZ(3,1)
            ENDIF
C           LES 3 COMPOSANTES DE LA TANGENTE 2
            IF( NTG2 .NE. 0 ) THEN
               MN = MNTG + 3 * ABS(NTG2)
               IF( NTG2 .GT. 0 ) THEN
                  XYZTG(1,2) = RMCN( MN     )
                  XYZTG(2,2) = RMCN( MN + 1 )
                  XYZTG(3,2) = RMCN( MN + 2 )
               ELSE
C                 TANGENTE A INVERSER
                  XYZTG(1,2) =-RMCN( MN     )
                  XYZTG(2,2) =-RMCN( MN + 1 )
                  XYZTG(3,2) =-RMCN( MN + 2 )
               ENDIF
            ELSE
               XYZTG(1,2) = XYZ(1,1) - XYZ(1,2)
               XYZTG(2,2) = XYZ(2,1) - XYZ(2,2)
               XYZTG(3,2) = XYZ(3,1) - XYZ(3,2)
            ENDIF
C
C           LE TRACE DE L'ARETE P3 ** 3
            CALL XVEPAISSEUR( NEPARL )
            CALL TRAR3D( NCOUAL, PREDUA,
     %                   XYZ(1,1),   XYZ(1,2),
     %                   XYZTG(1,1), XYZTG(1,2) )
            GOTO 45
         ENDIF
C
C        LE TRACE DE L'ARETE I DROITE EN OUBLIANT LES EXTREMITES
 40      DIRX = ( XYZ(1,2) - XYZ(1,1) ) * REDUCA
         DIRY = ( XYZ(2,2) - XYZ(2,1) ) * REDUCA
         DIRZ = ( XYZ(3,2) - XYZ(3,1) ) * REDUCA
         XYZ(1,1) = XYZ(1,1) + DIRX
         XYZ(2,1) = XYZ(2,1) + DIRY
         XYZ(3,1) = XYZ(3,1) + DIRZ
         XYZ(1,2) = XYZ(1,2) - DIRX
         XYZ(2,2) = XYZ(2,2) - DIRY
         XYZ(3,2) = XYZ(3,2) - DIRZ
         CALL XVEPAISSEUR( NEPARL )
         CALL TRAIT3D( NCOUAL, XYZ(1,1), XYZ(1,2) )
C
C        TRACE EVENTUEL DU NO DES SOMMETS
 45      IF( IAVNSO .NE. 0 ) THEN
            WRITE( NMSOMM, '(I8)' ) NS1
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONSO, XYZ(1,1), NMSOMM(1:L) )
            WRITE( NMSOMM, '(I8)' ) NS2
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONSO, XYZ(1,2), NMSOMM(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DE L'EF
         IF( IAVNEF .NE. 0 ) THEN
            WRITE( NMSOMM, '(I8)' ) I
            CALL SANSBL( NMSOMM, L )
            IF( INDTGS .GT. 0 ) THEN
C              LE POINT MILIEU EN COORDONNEES OBJET 3D
               DO 50 J=1,3
                  XYZN(J)=VALP3H( 0.5, XYZ(J,1),   XYZ(J,2),
     %                                 XYZTG(J,1), XYZTG(J,2) )
 50            CONTINUE
            ELSE
               DO 55 J=1,3
                  XYZN(J) = ( XYZ(J,1) + XYZ(J,2) ) * 0.5
 55            CONTINUE
            ENDIF
            CALL TEXTE3D( NCONEF, XYZN, NMSOMM(1:L) )
         ENDIF
C
C        LE POINT DE TRACE DE LA 'POIGNEE DE LA LIGNE'
         IF( I .EQ. NBARL2 ) THEN
            IF( INDTGS .GT. 0 ) THEN
C              LE POINT MILIEU EN COORDONNEES OBJET 3D
               DO 65 J=1,3
                  XYZP(J)=VALP3H( 0.5, XYZ(J,1),   XYZ(J,2),
     %                                 XYZTG(J,1), XYZTG(J,2) )
 65            CONTINUE
            ELSE
C              SES COORDONNEES
               DO 70 J=1,3
                  XYZP(J) = ( XYZ(J,1) + XYZ(J,2) ) * 0.5
 70            CONTINUE
            ENDIF
         ENDIF
 100  CONTINUE
C
C     TRACE DE LA POIGNEE SI ELLE VISIBLE DANS LA FENETRE
      CALL ITEML3( XYZP, NMOBJT, NUOBJT )
C
C     REPRISE DE L'ORBITE SI LORBITE>0
      IF( LORBITE .GT. 0 ) THEN
C        ORBITE OU ZOOM OU TRANSLATION ACTIFS
         IF( LANGAG .EQ. 0 ) THEN
            KNLIGN = 'LIGNE ' // NMOBJT
         ELSE
            KNLIGN = 'LINE ' // NMOBJT
         ENDIF
         CALL TRFINS( KNLIGN )
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
         GOTO 10
      ENDIF
C
 9000 CALL XVEPAISSEUR( 0 )
      RETURN
      END
