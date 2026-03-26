      SUBROUTINE SUEX03( NTLXSU, LADEFI, RADEFI,
     %                   NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LES FACES DE LA SURFACE B-SPLINE UNIFORME
C ----- D'INTERPOLATION A PARTIR D'UNE GRILLE DE POINTS D'INTERPOLATION
C       OU SOMMETS DE LIGNES

C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA SURFACE
C          CES 2 TABLEAUX LADEFI ET RADEFI ONT MEME ADRESSE A L'APPEL

C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C NTSOSU : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1997
C2345X7..............................................................012
      PARAMETER        (NBCOMP=3)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_ligne__bspline.inc"
      include"./incl/a___tableau1r.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_surface__bspline.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
C
C     LES PARAMETRES NECESSAIRES POUR EXECUTER XYZBST DANS LE SP
C     TEHOTE APPELE PAR SUEX09 SI LA FONCTION TAILLE_IDEALE EXISTE
      COMMON / S09S03 / LDEXSB, LRX, MNRX,
     %                  LDEYSB, LRY, MNRY, MNSPQB
C     LA LONGUEUR DES ARETES DE LA LIGNE ENVELOPPE DANS R**3
      COMMON / PEL1R3 / PERMT3
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      INTEGER           LADEFL(0:WULFTR)
      REAL              RADEFL(0:WULFTR)
      EQUIVALENCE      (LADEFL(0),RADEFL(0))
C
      CHARACTER*24      KNOMLI
      INTEGER           NTLXLI(5), NTARLI(5), MNARLI(5),
     %                  NTSOLI(5), MNSOLI(5)
C
C     PROTECTIONS
      MNMATR = 0
      MNXYZL = 0
      MNBX   = 0
      MNBY   = 0
      MNAX   = 0
      MNAY   = 0
      MNFACM = 0
      MNSPQB = 0
      MNUX   = 0
      MNUY   = 0
      MNTX   = 0
      MNTY   = 0
      MNRX   = 0
      MNRY   = 0
      MNNORX = 0
      MNNORY = 0
      MNSX   = 0
      MNPTSB = 0
      MNLADE = 0
      MNTGSO = 0
      KX = 0
      KY = 0
C
C     VERIFICATIONS DES DONNEES
C     =========================
      IERR   = 0
      NUTYSU = LADEFI(WUTYSU)
C
C     LE NOMBRE D'ARETES DANS LA DIRECTION X ET Y
      NBAXQB = LADEFI(WBAXQB)
      NBAYQB = LADEFI(WBAYQB)
      IF( NBAXQB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE INCORRECT D''ARETES DANS LA DIRECTION X'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
      IF( NBAYQB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE INCORRECT D''ARETES DANS LA DIRECTION Y'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
C     LA RAISON DE LA PROGRESSION GEOMETRIQUE D'ESPACEMENT EN X ET Y
      RAGXQB = RADEFI(WAGXQB)
      RAGYQB = RADEFI(WAGYQB)
      IF( RAGXQB .LE. 0 .OR. RAGYQB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAISON GEOMETRIQUE=<0'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     LE DEGRE EN X ET Y DES POLYNOMES
      LDEXSB = LADEFI( WDEXSB )
      LDEYSB = LADEFI( WDEYSB )
C
      IF( NUTYSU .EQ. 3 ) THEN
C
C        B-SPLINE A L'AIDE DE POINTS
C        ===========================
C        LE NOMBRE DE POINTS D'INTERPOLATION
         NBPIEX = LADEFI( WBPIEX )
         NBPIEY = LADEFI( WBPIEY )
         IF( NBPIEX .LE. 1 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'LE NOMBRE DE POINTS EN X DOIT ETRE >1'
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
         IF( NBPIEY .LE. 1 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'LE NOMBRE DE POINTS EN Y DOIT ETRE >1'
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
         IF( LDEXSB .LE. 0 .OR. LDEXSB .GE. NBPIEX ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:5),'(I5)') NBPIEX
            KERR(1) = 'LE DEGRE EN X DOIT ETRE >0 ET <'
     %               //KERR(MXLGER)(1:5)
            CALL LEREUR
            IERR = 6
            GOTO 9999
         ENDIF
         IF( LDEYSB .LE. 0 .OR. LDEYSB .GE. NBPIEY ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:5),'(I5)') NBPIEY
            KERR(1) = 'LE DEGRE EN Y DOIT ETRE >0 ET <'
     %                //KERR(MXLGER)(1:5)
            CALL LEREUR
            IERR = 6
            GOTO 9999
         ENDIF
C
C        LECTURE DES 3 COORDONNEES DES POINTS D'INTERPOLATION
C        TABLEAU POINSB(1:NBPIEX,1:NBPIEY,1:3)
         NBINCS = NBPIEX * NBPIEY
         CALL TNMCDC( 'REEL', 3*NBINCS, MNPTSB )
         MNPX = MNPTSB
         MNPY = MNPX + NBINCS
         MNPZ = MNPY + NBINCS
         NUP  = 0
         DO 20 J=1,NBPIEY
            DO 10 I=1,NBPIEX
C              LE NUMERO DU POINT
               N = LADEFI( WOINSB - 1 + I + NBPIEX * (J-1) )
C              OUVERTURE DU LEXIQUE DU POINT
               CALL LXNLOU( NTPOIN, N, NTLXPO, MN )
               IF( NTLXPO .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:4),'(I4)') I
                  WRITE(KERR(MXLGER)(5:8),'(I4)') J
                  KERR(1) = 'POINT INCONNU ' // KERR(MXLGER)(1:4)
     %                      // KERR(MXLGER)(5:8)
                  CALL LEREUR
                  IERR = 7
                  GOTO 10
               ENDIF
C              OUVERTURE DU TABLEAU XYZSOMMET DU POINT
               CALL LXTSOU( NTLXPO, 'XYZSOMMET', NT, MNS )
               IF( NT .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:4),'(I4)') I
                  WRITE(KERR(MXLGER)(5:8),'(I4)') J
                  KERR(1) = 'POINT ' // KERR(MXLGER)(1:4)
     %                      // KERR(MXLGER)(5:8) //
     %                      ' SANS COORDONNEES'
                  CALL LEREUR
                  IERR = 7
                  GOTO 10
               ENDIF
C              LES 3 COORDONNEES DU POINT
               MNS = MNS + WYZSOM
               RMCN( MNPX + NUP ) = RMCN( MNS )
               RMCN( MNPY + NUP ) = RMCN( MNS+1 )
               RMCN( MNPZ + NUP ) = RMCN( MNS+2 )
               NUP = NUP + 1
 10         CONTINUE
 20      CONTINUE
         IF( IERR .NE. 0 ) GOTO 9900
C
C        ICI LES B-SPLINES DES BORDS SERONT OUVERTES ET D'INTERPOLATION
         NUTYLI = 11
C
      ELSE IF( NUTYSU .EQ. 4 ) THEN
C
C        B-SPLINE D'INTERPOLATION A L'AIDE DES SOMMETS DE LIGNES QUELCONQUES
C        LES SOMMETS DE CHACUNE SERVENT DE POINTS D'INTERPOLATION EN U
C        ATTENTION: LES SOMMETS DOIVENT PROGRESSER DANS LE MEME SENS
C        POUR QUE LE QUADRANGLE SPLINE NE SOIT PAS TROP MOUVEMENTE
C        ===================================================================
C        LE NOMBRE DE LIGNES D'INTERPOLATION
         NBLISB = LADEFI(WBLISB)
         IF( LDEYSB .GE. NBLISB ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:5),'(I5)') NBLISB
            KERR(1) = 'LE DEGRE EN Y DOIT ETRE >0 ET <'//
     %                 KERR(MXLGER)(1:5)
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
         NBPIEY = NBLISB
C
C        LES LIGNES INTERPOLATION DE LA SURFACE B-SPLINE
C        LE TABLEAU DES ADRESSES MCN DES TABLEAUX 'SPLINE'  DES LIGNES B-SPLI
         CALL TNMCDC( 'REEL', NBLISB, MNXYZL )
         NBPIEX = 0
         NBX    = 0
         DO 3 I=0,NBLISB-1
            N = LADEFI( WIGNSB + I )
C           RECUPERATION DE LA LIGNE
            CALL LXNLOU( NTLIGN, N, NTLXLI(1), MN )
            IF( NTLXLI(1) .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LIGNE INCONNUE: ' // KERR(MXLGER)(1:4)
               ELSE
                  KERR(1) = 'UNKNOWN LINE: ' // KERR(MXLGER)(1:4)
               ENDIF
               CALL LEREUR
               IERR = 6
               GOTO 3
            ENDIF
C
C           RECUPERATION DES COORDONNEES DES SOMMETS DES LIGNES
            CALL LXTSOU( NTLXLI(1), 'XYZSOMMET', NTXYZL, MCN(MNXYZL+I) )
            IF( NTXYZL .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               IF( LANGAG .EQ. 0 ) THEN
                 KERR(1)='LIGNE '//KERR(MXLGER)(1:4)//' SANS XYZSOMMET'
              ELSE
               KERR(1)='LINE '//KERR(MXLGER)(1:4)//' WITHOUT XYZSOMMET'
              ENDIF
               CALL LEREUR
               IERR = 7
               GOTO 3
            ENDIF
C
C           RECUPERATION DES NO SOMMETS DES EF 'NSEF' DES LIGNES
            CALL LXTSOU( NTLXLI(1), 'NSEF', NTNSEF, MNNSEF )
            IF( NTNSEF .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') I+1
               IF( LANGAG .EQ. 0 ) THEN
                 KERR(1)='LIGNE '//KERR(MXLGER)(1:4)//' SANS XYZSOMMET'
              ELSE
               KERR(1)='LINE '//KERR(MXLGER)(1:4)//' WITHOUT XYZSOMMET'
              ENDIF
               CALL LEREUR
               IERR = 7
               GOTO 3
            ENDIF
C
C           LIGNE OUVERTE OU FERMEE?
            NUTFMA = MCN( MNNSEF + WUTFMA )
C           ( -1 : 'unknown' , 0 : 'not-closed' , 1 : 'closed line or surface' )
C
C           LE NOMBRE DE SOMMETS DE CETTE LIGNE EST IL EGAL
C           A CELUI DES AUTRES LIGNES ?
            MNXYZ  = MCN(MNXYZL+I)
            NBSOML = MCN( MNXYZ + WNBSOM )
            IF( NBPIEX .EQ. 0 ) THEN
C              PREMIERE LIGNE
               NBPIEX = NBSOML
               NUTYL0 = NUTFMA
            ELSE
C              UNE LIGNE SUIVANTE
               IF( NBSOML .NE. NBPIEX ) THEN
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NBSOML
                  WRITE(KERR(MXLGER)(5:8),'(I4)') I+1
                  WRITE(KERR(MXLGER)(9:12),'(I4)') NBPIEX
                  IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LE NOMBRE ' // KERR(MXLGER)(1:4) //
     %                      ' DE SOMMETS DE LA LIGNE ' //
     %                       KERR(MXLGER)(5:8)
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NBPIEX
                  KERR(2) = 'EST DIFFERENT ' // KERR(MXLGER)(9:12) //
     %                      ' DES PRECEDENTES LIGNES'
                  ELSE
                  KERR(1) = 'The VERTICES NUMBER '// KERR(MXLGER)(1:4)//
     %                      ' of LINE ' //
     %                       KERR(MXLGER)(5:8)
                  WRITE(KERR(MXLGER)(1:4),'(I4)') NBPIEX
                  KERR(2) = 'IS NOT EQUAL TO' // KERR(MXLGER)(9:12) //
     %                      ' of PREVIOUS LINES'
                  ENDIF
                  CALL LEREUR
                  IERR = 11
                  GOTO 3
               ENDIF
            ENDIF
            IF( LDEXSB .LE. 0 .OR. LDEXSB .GE. NBPIEX ) THEN
               NBLGRC(NRERR) = 2
               WRITE(KERR(MXLGER)(1:5),'(I5)') NBPIEX
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'Le DEGRE en U DOIT ETRE >0 ET <'
     %                       //KERR(MXLGER)(1:5)
                  KERR(2) = 'NOMBRE de SOMMETS de la LIGNE'
               ELSE
                  KERR(1) = 'The U-DEGREE MUST BE >0 and <'
     %                       //KERR(MXLGER)(1:5)
                  KERR(2) = 'the LINE VERTICES NUMBER'
               ENDIF
               CALL LEREUR
               IERR = 6
               GOTO 9900
            ENDIF
C
    3    CONTINUE
C
C        ICI LES B-SPLINES DES BORDS SERONT OUVERTES ET D'INTERPOLATION
         NUTYLI = 11
         IF( IERR .NE. 0 ) GOTO 9900
      ENDIF
C
C     LES CONSTANTES EN X POUR UX TX RX DE LA SURFACE B-SPLINE
      IF( NUTYLI .GE. 21 ) THEN
         LUX = NBPIEX + 1
      ELSE
         LUX = NBPIEX
      ENDIF
      LTX = LUX
      LRX = LTX - LDEXSB
C
C     LES CONSTANTES EN Y POUR UY TY RY DE LA SURFACE B-SPLINE
      LUY = NBPIEY
      LTY = LUY
      LRY = LTY - LDEYSB
C
C     LE NOMBRE D'INCONNUES DU SYSTEME
      NBINCS = LUX * LUY
C
C     LE TABLEAU DES B(J,M) EN X ET EN Y
      KX = LDEXSB + 1
      KY = LDEYSB + 1
      CALL TNMCDC( 'REEL', KX * (LDEXSB+2), MNBX )
      CALL TNMCDC( 'REEL', KY * (LDEYSB+2), MNBY )
C
C     LE TABLEAU DES A(J,M) EN X ET EN Y
      CALL TNMCDC( 'REEL', KX * KX, MNAX )
      CALL TNMCDC( 'REEL', KY * KY, MNAY )
C
C     LES 1/FACTORIELLES DE M
      CALL TNMCDC( 'REEL', MAX(KX,KY), MNFACM )
C
C     RESERVATION DE LA MATRICE A INVERSER ET DES 3 SECONDS MEMBRES
      CALL TNMCDC( 'REEL2', NBINCS*(NBINCS+NBCOMP), MNMATR )
C
C     LE TABLEAU SSPLINE( 0:LDEXSB, 0:LDEYSB, 1:LRX, 1:LRY, 1:NBCOMP )
      MOSSPL = (LDEXSB+1)*(LDEYSB+1)*LRX*LRY*NBCOMP
      CALL TNMCDC( 'REEL', MOSSPL, MNSPQB )
C
C     LE TABLEAU SX( 0:LDEXSB, 1:NBCOMP, 0:LDEYSB )
      MOSX = (LDEXSB+1)*(LDEYSB+1)*NBCOMP
      CALL TNMCDC( 'REEL', MOSX, MNSX )
C
C     DECLARATION DU TABLEAU UX PARAMETRE UNIFORME EN X DE LA SURFACE
      MOTUX = LUX + 1
      CALL TNMCDC( 'REEL', MOTUX, MNUX )
C
C     DECLARATION DU TABLEAU UY PARAMETRE UNIFORME EN Y DE LA SURFACE
      MOTUY = LUY + 1
      CALL TNMCDC( 'REEL', MOTUY, MNUY )
C
C     LE TABLEAU DES NOEUDS TX
      MOTTX = LTX + KX + 1
      CALL TNMCDC( 'REEL', MOTTX, MNTX )
C
C     LE TABLEAU DES NOEUDS TY
      MOTTY = LTY + KY + 1
      CALL TNMCDC( 'REEL', MOTTY, MNTY )
C
C     LE TABLEAU DU PARAMETRE RX
      CALL TNMCDC( 'REEL', LRX+1, MNRX )
C     LE TABLEAU DU PARAMETRE RY
      CALL TNMCDC( 'REEL', LRY+1, MNRY )
C
C     LE NUMERO DU DERNIER NOEUD TX EN CHAQUE POINT RX
      CALL TNMCDC( 'ENTIER', LRX+1, MNNORX )
C     LE NUMERO DU DERNIER NOEUD TY EN CHAQUE POINT RY
      CALL TNMCDC( 'ENTIER', LRY+1, MNNORY )
C
C     LA CONSTRUCTION DES POLYNOMES SUR CHAQUE INTERVALLE DE LA B-SPLINE
C     ==================================================================
      IF( NUTYSU .EQ. 3 ) THEN
         MNXYZL = MNPTSB
      ELSE IF( NUTYSU .EQ. 4 ) THEN
         MNPTSB = MNXYZL
      ENDIF
      CALL BSPLSI( NUTYSU,
     %             LDEXSB, LUX, LTX, LRX,
     %             RMCN(MNUX)  , RMCN(MNTX)  , RMCN(MNRX),
     %             LDEYSB, LUY, LTY, LRY,
     %             RMCN(MNUY)  , RMCN(MNTY)  , RMCN(MNRY),
     %             NBINCS      , RMCN(MNPTSB), MCN(MNXYZL),
     %             RMCN(MNAX)  , RMCN(MNBX)  ,
     %             RMCN(MNAY)  , RMCN(MNBY)  , RMCN(MNFACM),
     %             RMCN(MNMATR), MCN(MNNORX) , MCN(MNNORY),
     %             RMCN(MNSX)  , RMCN(MNSPQB), IERR )
      IF( NUTYSU .EQ. 3 ) THEN
         MNXYZL = 0
      ELSE IF( NUTYSU .EQ. 4 ) THEN
         MNPTSB = 0
      ENDIF
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
      CALL TNMCDS( 'REEL2', NBINCS*(NBINCS+NBCOMP), MNMATR )
      CALL TNMCDS( 'REEL', MOSX, MNSX )
      IF( IERR .NE. 0 ) GOTO 9900
C
C     CONSTRUCTION DU TABLEAU BSPLINE DE LA SURFACE
C     ---------------------------------------------
      CALL LXTSOU( NTLXSU, 'BSPLINE',  NTBSPL, MNBSPL )
      IF( NTBSPL .GT. 0 ) THEN
         CALL LXTSDS( NTLXSU, 'BSPLINE' )
      ENDIF
      I = WARXQB + 1 + LRX + 1 + LRY + MOSSPL
      CALL LXTNDC( NTLXSU, 'BSPLINE', 'MOTS' , I )
      CALL LXTSOU( NTLXSU, 'BSPLINE',  NTBSPL, MNBSPL )
C     DEGRBS 'degre en X des polynomes de la BSPLINE'
      MCN( MNBSPL + WEGXQB ) = LDEXSB
C     DEGRBS 'degre en Y des polynomes de la BSPLINE'
      MCN( MNBSPL + WEGYQB ) = LDEYSB
C     NBCOBS 'nombre de composantes de la BSPLINE 3 ou 4'
      MCN( MNBSPL + WBCOQB ) = NBCOMP
C     LRX 'nombre d'intervalles en X'
      MCN( MNBSPL + WBIXQB ) = LRX
C     LRY 'nombre d'intervalles en Y'
      MCN( MNBSPL + WBIYQB ) = LRY
C     PARXQB(0..LRX) 'valeurs du parametre X de la BSPLINE'
      CALL TRTATA( MCN(MNRX), MCN(MNBSPL+WARXQB), 1+LRX )
C     PARYQB(0..LRY) 'valeurs du parametre Y de la BSPLINE'
      MN = MNBSPL+WARXQB+1+LRX
      CALL TRTATA( MCN(MNRY), MCN(MN), 1+LRY )
C     LE TABLEAU SPLIQB( 0:LDEXSB, 0:LDEYSB, 1:LRX, 1:LRY, 1:NBCOMP )
      MN = MN + 1 + LRY
      CALL TRTATA( MCN(MNSPQB), MCN(MN),  MOSSPL )
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNBSPL) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNBSPL + MOTVAR(6) ) = NONMTD( '~>SURFACE>>BSPLINE' )
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE' DES ARETES AUTOUR DU POINT
C     ICI LA CARTE EST SUPPOSEE ISOTROPE
      NOFOTI = NOFOTIEL()
C
C     ATTENTION: SI SURFACE B-SPLINE INTERPOLATION DE LIGNES ALORS
C     LA FONCTION 'TAILLE_IDEALE' N'EST PAS PRISE EN COMPTE
C     CAR LES ARETES ONT DEJA ETE CALCULEES SELON CETTE FONCTION
      IF( NUTYSU .EQ. 4 ) NOFOTI=0
C
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NOFOTI .GT. 0 ) GOTO 200
C
C     ===============================
C     PAS DE FONCTION 'TAILLE_IDEALE'
C     ===============================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     -----------------------------------
      NBSOSU = ( NBAXQB + 1 ) * ( NBAYQB + 1 )
      NBEFOB = NBAXQB * NBAYQB
      IF( LDEXSB .GT. 1 .OR. LDEYSB .GT. 1 ) THEN
C        CALCUL DES TANGENTES : 8 PAR QUADRANGLE
         NBTGEF = 8
         NBTGS  = NBTGEF * NBEFOB
         NBEFAP = NBEFOB
         NBEFTG = NBEFOB
      ELSE
         NBTGEF = 0
         NBTGS  = 0
         NBEFAP = 0
         NBEFTG = 0
      ENDIF
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'ENTIER',
     %             WYZSOM + 3 * (NBSOSU + NBTGS) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOSU , MNSOSU )
C
C     LE NOMBRE DE COORDONNEES D'UN SOMMET
      MCN( MNSOSU + WBCOOR ) = 3
C     LE NOMBRE DE SOMMETS
      MCN( MNSOSU + WNBSOM ) = NBSOSU
C     LE NOMBRE DE TANGENTES CALCULEES
      MCN( MNSOSU + WNBTGS ) = NBTGS
C
C     TABLEAU AUXILIAIRE DES TANGENTES EN CHAQUE SOMMET
      CALL TNMCDC( 'REEL', NBSOSU * 6, MNTGSO )
C
C     CALCUL DES 3 COORDONNEES DES SOMMETS ET DES TANGENTES DE CHAQUE QUADRANGLE
      IF( RAGXQB .EQ. 1 .AND. MOD(NBAXQB,LUX-1) .EQ. 0 .AND.
     %    RAGYQB .EQ. 1 .AND. MOD(NBAYQB,LUY-1) .EQ. 0  ) THEN
C
C        RESPECT DES POINTS D'INTERPOLATION EN X ET Y DANS LE MAILLAGE
         CALL BSPLAC( NBCOMP,
     %                LDEXSB, LRX, LUX, RMCN(MNRX), RMCN(MNUX), NBAXQB,
     %                LDEYSB, LRY, LUY, RMCN(MNRY), RMCN(MNUY), NBAYQB,
     %                RMCN(MNSPQB), RMCN(MNSOSU+WYZSOM),
     %                NBTGS, RMCN(MNSOSU+WYZSOM+3*NBSOSU), RMCN(MNTGSO))
C
      ELSE IF( RAGXQB .EQ. 1 .AND. MOD(NBAXQB,LUX-1) .EQ. 0 ) THEN
C
C        RESPECT DES POINTS D'INTERPOLATION EN X DANS LE MAILLAGE
C        POUR LES 2 EXTREMITES DE Y
         CALL BSPLAX( NBCOMP,
     %                LDEXSB, LRX, RMCN(MNRX), LUX, RMCN(MNUX), NBAXQB,
     %                LDEYSB, LRY, RMCN(MNRY), NBAYQB, RAGYQB,
     %                RMCN(MNSPQB), RMCN(MNSOSU+WYZSOM),
     %                NBTGS, RMCN(MNSOSU+WYZSOM+3*NBSOSU), RMCN(MNTGSO))
C
      ELSE IF( RAGYQB .EQ. 1 .AND. MOD(NBAYQB,LUY-1) .EQ. 0  ) THEN
C
C        RESPECT DES POINTS D'INTERPOLATION EN Y DANS LE MAILLAGE
C        POUR LES 2 EXTREMITES DE X
         CALL BSPLAY( NBCOMP,
     %                LDEXSB, LRX, RMCN(MNRX), NBAXQB, RAGXQB,
     %                LDEYSB, LRY, RMCN(MNRY), LUY, RMCN(MNUY), NBAYQB,
     %                RMCN(MNSPQB), RMCN(MNSOSU+WYZSOM),
     %                NBTGS, RMCN(MNSOSU+WYZSOM+3*NBSOSU), RMCN(MNTGSO))
C
      ELSE
C
C        LE CALCUL DES COORDONNEES DES SOMMETS DES FACES DE CETTE SURFACE
         CALL BSPLSC( NBCOMP,
     %                LDEXSB, LRX, RMCN(MNRX), NBAXQB, RAGXQB,
     %                LDEYSB, LRY, RMCN(MNRY), NBAYQB, RAGYQB,
     %                RMCN(MNSPQB), RMCN(MNSOSU+WYZSOM),
     %                NBTGS, RMCN(MNSOSU+WYZSOM+3*NBSOSU), RMCN(MNTGSO))
C
      ENDIF
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE
      CALL TNMCDS( 'REEL', NBSOSU * 6, MNTGSO )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' SURFACE STRUCTUREE
C     -------------------------------------------------
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER', 1+WBARYQ+NBEFAP+9*NBEFTG )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU, MNFASU )
C     LE TYPE DE L'OBJET : ICI SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = 0
C     LE NOMBRE DE SOMMETS PAR FACE
      MCN ( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C0
      MCN ( MNFASU + WBTGEF ) = NBTGEF
      MCN ( MNFASU + WBEFAP ) = NBEFAP
      MCN ( MNFASU + WBEFTG ) = NBEFTG
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNFASU + WBEFOB ) = NBEFOB
C     LE TYPE DU MAILLAGE : ICI QUADRANGLE STRUCTURE
      MCN( MNFASU + WUTYMA ) = 4
C     LE NOMBRE D'ARETES EN X DU QUADRANGLE STRUCTURE
      MCN( MNFASU + WBARXQ ) = NBAXQB
C     LE NOMBRE D'ARETES EN Y DU QUADRANGLE STRUCTURE
      MCN( MNFASU + WBARYQ ) = NBAYQB
C
      IF( NBTGS .GT. 0 ) THEN
C
C        LE POINTEUR SUR LES EF A TG
         MN = MNFASU + WBARYQ
         DO 120 J=1,NBEFOB
            MCN(MN+J) = J
 120     CONTINUE
C
C        LE CODE GEOMETRIQUE DES EF A TG
         MN = MN + NBEFTG
         DO 130 J=1,NBEFTG
C           B-SPLINE
            MCN(MN+J) = 3
 130     CONTINUE
C
C        LE NUMERO DES 8 TANGENTES DE CHAQUE QUADRANGLE
         MN = MN + NBEFTG
         DO 140 J=1,NBTGS
            MCN(MN+J) = J
 140     CONTINUE
C
      ENDIF
      GOTO 9000
C
C     ========================================
C     EXISTENCE DE LA FONCTION 'TAILLE_IDEALE'
C     ========================================
C     LES VARIABLES DE LA SURFACE BSPLINE
C     variable NBAXQB 'nombre d''aretes en X' entier ;
C     variable NBAYQB 'nombre d''aretes en Y' entier ;
C     variable RAGXQB 'raison geometrique d''espacement en X des sommets' reel ;
C     variable RAGYQB 'raison geometrique d''espacement en Y des sommets' reel ;
C     variable LDEXSB 'degre en X des polynomes' entier ;
C     variable LDEYSB 'degre en Y des polynomes' entier ;
C     variable NBPIEX 'nombre de points d''interpolation en X' entier;
C     variable NBPIEY 'nombre de points d''interpolation en Y' entier;
C     tableau  POINSB(1..NBPIEX,1..NBPIEY)
C             'nom des points d''interpolation en X puis en Y' ^~>POINT ;
C
C     CONSTRUCTION DES 4 LIGNES B-SPLINES D'INTERPOLATION DANS R**3
C     variable NBARLI 'nombre d''aretes de la ligne B-spline' entier ;
C     variable RAGEBL 'raison geometrique espacant les sommets' reel ;
C     variable NDEGBL 'degre des polynomes de la B-spline' entier ;
C     variable NBPCBL 'nombre de points d''interpolation (>degre)' entier ;
C     tableau  NUPCBL(1..NBPCBL) 'nom des points d''interpolation'^~>POINT ;
C
 200  NBPIEX = LADEFI(WBPIEX)
      NBPIEY = LADEFI(WBPIEY)
C
C     RESERVATION D'UN TABLEAU DES NUMEROS DES POINTS DES 4 COTES
      CALL TNMCDC( 'ENTIER', MAX(NBPIEX,NBPIEY), MNPOCB )
C
C     LES LIGNES B-SPLINES DE POINTS D'INTERPOLATION DES COTES 1 A 4
C     ==============================================================
      PERMT3 = 0.0
      NBST = 0
      N    = 0
      DO 248 I=1,4
C
C        LE NOM DE LA LIGNE DU COTE I EST FORME A PARTIR DES NOMS
C        DES POINTS D'INTERPOLATION DE LA COURBE B-SPLINE
         GOTO( 211, 212, 213, 214 ),I
C
C        LES POINTS DE LA B-SPLINE DU COTE 1
 211     DO 221 J=0,NBPIEX-1
            MCN(MNPOCB+J) = LADEFI(WOINSB+J)
 221     CONTINUE
         NBP = NBPIEX
         GOTO 230
C
C        LES POINTS DE LA B-SPLINE DU COTE 3
 213     N = WOINSB + NBPIEX * ( NBPIEY - 1 )
         DO 223 J=0,NBPIEX-1
            MCN(MNPOCB+J) = LADEFI(N+J)
 223     CONTINUE
         NBP = NBPIEX
         GOTO 230
C
C        LES POINTS DE LA B-SPLINE DU COTE 2
 212     DO 222 J=1,NBPIEY
            MCN(MNPOCB-1+J) = LADEFI(WOINSB-1+J*NBPIEX)
 222     CONTINUE
         NBP = NBPIEY
         GOTO 230
C
C        LES POINTS DE LA B-SPLINE DU COTE 4
 214     DO 224 J=0,NBPIEY-1
            MCN(MNPOCB+J) = LADEFI(WOINSB+J*NBPIEX)
 224     CONTINUE
         NBP = NBPIEY
C
C        LE NOM DE LA LIGNE EN SENS RETROGRADE
 230     CALL NMLINM( 'POINT', NBP, MCN(MNPOCB), -1, KNOMLI )
         CALL LXLXOU( NTLIGN, KNOMLI, NTLXLI(I), MNLXLI )
         IF( NTLXLI(I) .GT. 0 ) THEN
C           LA LIGNE EXISTE DEJA MAIS PAS DANS LE BON SENS
            CALL LXTSOU( NTLXLI(I), 'NSEF',      NTARLI(I), MNARLI(I) )
            CALL LXTSOU( NTLXLI(I), 'XYZSOMMET', NTSOLI(I), MNSOLI(I) )
C           L'ADRESSE MCN EST INVERSEE POUR MARQUER LE SENS RETROGRADE
            MNSOLI(I) = -MNSOLI(I)
            GOTO 247
         ENDIF
C
C        LE NOM DE LA LIGNE DANS LE SENS DIRECT
         CALL NMLINM( 'POINT', NBP, MCN(MNPOCB), 1, KNOMLI )
C
C        CREATION DU LEXIQUE DE LA LIGNE
         CALL LXLXOU( NTLIGN, KNOMLI, NTLXLI(I), MNLXLI )
C
         IF( NTLXLI(I) .GT. 0 ) THEN
C           LA LIGNE EXISTE DANS LE BON SENS
            CALL LXTSOU( NTLXLI(I), 'NSEF',      NTARLI(I), MNARLI(I) )
            CALL LXTSOU( NTLXLI(I), 'XYZSOMMET', NTSOLI(I), MNSOLI(I) )
            GOTO 247
         ENDIF
C
C        LA LIGNE N'EXISTE PAS => ELLE EST CREEE
         CALL LXLXDC( NTLIGN, KNOMLI, 24, 8 )
         CALL LXLXOU( NTLIGN, KNOMLI, NTLXLI(I), MNLXLI )
C
C        LE TMS DEFINITION DE LA LIGNE
         MXLADE = WUPCBL + 1
         IF( MOD(I,2) .EQ. 0 ) THEN
            MXLADE = MXLADE + NBPIEY
         ELSE
            MXLADE = MXLADE + NBPIEX
         ENDIF
C
         CALL LXTNDC( NTLXLI(I), 'DEFINITION', 'ENTIER', MXLADE )
         CALL LXTSOU( NTLXLI(I), 'DEFINITION',  NTLADE , MNLADE )
C
C        LA CONSTRUCTION DU TABLEAU DEFINITION DE LA LIGNE
         IF( MOD(I,2) .EQ. 0 ) THEN
C
C           COTES 2 ET 4
             MCN(MNLADE+WBARLI) = LADEFI(WBAYQB)
            RMCN(MNLADE+WAGEBL) = RADEFI(WAGYQB)
             MCN(MNLADE+WDEGBL) = LADEFI(WDEYSB)
             MCN(MNLADE+WBPCBL) = NBPIEY
            IF( I .EQ. 2 ) THEN
               N = 0
            ELSE
               N = 1 - NBPIEX
            ENDIF
            DO 244 J = 1, NBPIEY
               MCN(MNLADE+WUPCBL-1+J) = LADEFI(WOINSB-1+J*NBPIEX+N)
 244        CONTINUE
C
         ELSE
C
C           COTES 1 ET 3
             MCN(MNLADE+WBARLI) = LADEFI(WBAXQB)
            RMCN(MNLADE+WAGEBL) = RADEFI(WAGXQB)
             MCN(MNLADE+WDEGBL) = LADEFI(WDEXSB)
             MCN(MNLADE+WBPCBL) = NBPIEX
            IF( I .EQ. 3 ) N = NBPIEX * ( NBPIEY - 1 )
            DO 246 J = 1, NBPIEX
               MCN(MNLADE+WUPCBL-1+J) = LADEFI(WOINSB-1+J+N)
 246        CONTINUE
C
         ENDIF
C
C        LA TRANSFORMATION
         MCN( MNLADE + WTYTRL ) = 1
C        LE TYPE DE LA LIGNE
         MCN( MNLADE + WUTYLI ) = 11
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNLADE) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNLADE + MOTVAR(6) ) = NONMTD( '~>LIGNE>>DEFINITION' )
C
C        LE MAILLAGE DE LA LIGNE B-SPLINE D'INTERPOLATION DU COTE I
C        CREATION DES TMS XYZSOMMET NSEF BSPLINE
         CALL LIEX11( NTLXLI(I), MCN(MNLADE), RMCN(MNLADE),
     %                NTARLI(I), MNARLI(I), NTSOLI(I), MNSOLI(I), IERR )
         IF( IERR .NE. 0 ) GOTO 9900
C
C        PERIMETRE DU QUADRANGLE DANS R**3
 247     CALL LISTLO( ABS(MNSOLI(I)), RLG )
         PERMT3 = PERMT3 + RLG
C
C        LA SOMME DES SOMMETS DES 4 LIGNES
         NBST = NBST + MCN( ABS(MNSOLI(I)) + WNBSOM ) - 1
 248  CONTINUE
C
C     DESTRUCTION DU TABLEAU MCN
      CALL TNMCDS( 'ENTIER', MAX(NBPIEX,NBPIEY), MNPOCB )
C
C     CONSTRUCTION DE LA LIGNE CONTOUR DU RECTANGLE PLAN
C     FORMEE DES ARETES DES 4 COTES DU RECTANGLE DANS
C     L'ESPACE DES PARAMETRES
C
C     CONSTRUCTION DE LA LIGNE FERMEE CONTOUR DU RECTANGLE DANS
C     L'ESPACE DES PARAMETRES C'EST A DIRE DANS LE PLAN
C     ---------------------------------------------------------
      KNOMLI = 'CONTOUR_RECTANGLE_B-SPLW'
 250  CALL LXLXOU( NTLIGN, KNOMLI, NTLXLI(5), MNLXLI )
C
C     S'IL N'EXISTE PAS IL EST CREE
      IF( NTLXLI(5) .LE. 0 ) THEN
         CALL LXLXDC( NTLIGN, KNOMLI, 24, 8 )
         GOTO 250
      ENDIF
C
C     SI LE TABLEAU XYZSOMMET EXISTE ALORS IL EST DETRUIT
      CALL LXTSOU( NTLXLI(5), 'XYZSOMMET',  NTSOLI(5) , MNSOLI(5) )
      IF( NTSOLI(5) .GT. 0 ) THEN
         CALL LXTSDS( NTLXLI(5), 'XYZSOMMET' )
      ENDIF
      CALL LXTNDC( NTLXLI(5), 'XYZSOMMET', 'ENTIER', WYZSOM+3*NBST )
      CALL LXTSOU( NTLXLI(5), 'XYZSOMMET',  NTSOLI(5) , MNSOLI(5) )
C
C     LE COTE 1 SUR L'AXE DES X DU RECTANGLE PLAN
C     -------------------------------------------
      CALL LXTSOU( NTLXLI(1), 'TABLEAU1R"PARAMST', NTPARA, MNPARA )
      IF( NTPARA .LE. 0 ) THEN
C        LA LIGNE N'EST PAS UNE B-SPLINE
         NBLGRC(NRERR) = 1
         KERR(1) = 'COTE 1 NON LIGNE B-SPLINE D''INTERPOLATION'
         CALL LEREUR
         IERR = 7
         GOTO 9900
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      N = MCN( ABS(MNSOLI(1)) + WNBSOM )
C     LE 1-ER SOMMET DU RECTANGLE COURBE
      MNS = MNSOLI(5) + WYZSOM
      IF( MNSOLI(1) .LT. 0 ) THEN
C        SENS RETROGRADE
         MNSOLI(1) = - MNSOLI(1)
         RL        = RMCN(MNPARA+WABL1R-1+N) - RMCN(MNPARA+WABL1R)
         DO 251 J = N, 1, -1
            RMCN(MNS  ) = RL - RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+1) = 0.0
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 251     CONTINUE
      ELSE
C        SENS DIRECT
         DO 252 J = 1, N
            RMCN(MNS  ) = RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+1) = 0.0
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 252     CONTINUE
      ENDIF
C
C     LE COTE 2 : DU SOMMET 2 AU SOMMET 3 DU RECTANGLE PLAN
C     ------------------------------------------------------
      CALL LXTSOU( NTLXLI(2), 'TABLEAU1R"PARAMST', NTPARA, MNPARA )
      IF( NTPARA .LE. 0 ) THEN
C        LA LIGNE N'EST PAS UNE B-SPLINE
         NBLGRC(NRERR) = 1
         KERR(1) = 'COTE 2 NON LIGNE B-SPLINE D''INTERPOLATION'
         CALL LEREUR
         IERR = 7
         GOTO 9900
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      N   = MCN( ABS(MNSOLI(2)) + WNBSOM )
      MNS = MNS - 3
      IF( MNSOLI(2) .LT. 0 ) THEN
C        SENS RETROGRADE
         MNSOLI(2) = - MNSOLI(2)
         RL        = RMCN(MNPARA+WABL1R-1+N) - RMCN(MNPARA+WABL1R)
         DO 253 J = N, 1, -1
            RMCN(MNS  ) = RMCN(MNRX+LRX)
            RMCN(MNS+1) = RL - RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 253     CONTINUE
      ELSE
C        SENS DIRECT
         DO 254 J = 1, N
            RMCN(MNS  ) = RMCN(MNRX+LRX)
            RMCN(MNS+1) = RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 254     CONTINUE
      ENDIF
C
C     LE COTE 3 : DU SOMMET 3 AU SOMMET 4 DU RECTANGLE PLAN
C     ------------------------------------------------------
      CALL LXTSOU( NTLXLI(3), 'TABLEAU1R"PARAMST', NTPARA, MNPARA )
      IF( NTPARA .LE. 0 ) THEN
C        LA LIGNE N'EST PAS UNE B-SPLINE
         NBLGRC(NRERR) = 1
         KERR(1) = 'COTE 3 NON LIGNE B-SPLINE D''INTERPOLATION'
         CALL LEREUR
         IERR = 7
         GOTO 9900
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      N   = MCN( ABS(MNSOLI(3)) + WNBSOM )
      MNS = MNS - 3
      IF( MNSOLI(3) .LT. 0 ) THEN
C        LIGNE EN SENS RETROGRADE A RENVERSER : - - => +
         MNSOLI(3) = - MNSOLI(3)
         RL        = RMCN(MNPARA+WABL1R-1+N) - RMCN(MNPARA+WABL1R)
         DO 255 J = 1, N
            RMCN(MNS  ) = RL - RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+1) = RMCN(MNRY+LRY)
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 255     CONTINUE
      ELSE
C        SENS DIRECT
         DO 256 J = N, 1, -1
            RMCN(MNS  ) = RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+1) = RMCN(MNRY+LRY)
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 256     CONTINUE
      ENDIF
C
C     LE COTE 4 : DU SOMMET 4 AU SOMMET 1 DU RECTANGLE PLAN
C     ------------------------------------------------------
      CALL LXTSOU( NTLXLI(4), 'TABLEAU1R"PARAMST', NTPARA, MNPARA )
      IF( NTPARA .LE. 0 ) THEN
C        LA LIGNE N'EST PAS UNE B-SPLINE
         NBLGRC(NRERR) = 1
         KERR(1) = 'COTE 4 NON LIGNE B-SPLINE D''INTERPOLATION'
         CALL LEREUR
         IERR = 7
         GOTO 9900
      ENDIF
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      N   = MCN( ABS(MNSOLI(4)) + WNBSOM )
      MNS = MNS - 3
      IF( MNSOLI(4) .LT. 0 ) THEN
C        LIGNE EN SENS RETROGRADE A RENVERSER : - - => +
         MNSOLI(4) = - MNSOLI(4)
         RL        = RMCN(MNPARA+WABL1R-1+N) - RMCN(MNPARA+WABL1R)
         DO 257 J = 1, N-1
            RMCN(MNS  ) = 0.0
            RMCN(MNS+1) = RL - RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 257     CONTINUE
      ELSE
C        SENS DIRECT
         DO 258 J = N, 2, -1
            RMCN(MNS  ) = 0.0
            RMCN(MNS+1) = RMCN(MNPARA+WABL1R-1+J)
            RMCN(MNS+2) = 0.0
            MNS = MNS + 3
 258     CONTINUE
      ENDIF
C
C     FIN DE REMPLISSAGE DU TABLEAUX XYZSOMMET DE LA LIGNE CONTOUR
      MCN( MNSOLI(5) + WNBSOM ) = NBST
      MCN( MNSOLI(5) + WBCOOR ) = 3
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOLI(5) + WNBTGS ) = 0
      MCN( MNSOLI(5) + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     SI LE TABLEAU NSEF EXISTE ALORS IL EST DETRUIT
      CALL LXTSOU( NTLXLI(5), 'NSEF',  NTARLI(5) , MNARLI(5) )
      IF( NTARLI(5) .GT. 0 ) THEN
         CALL LXTSDS( NTLXLI(5), 'NSEF' )
      ENDIF
C     CREATION DU TABLEAU NSEF
      CALL LXTNDC( NTLXLI(5), 'NSEF', 'ENTIER', WUSOEF + 2 * NBST )
      CALL LXTSOU( NTLXLI(5), 'NSEF',  NTARLI(5), MNARLI(5) )
      MCN( MNARLI(5) + WUTYOB ) = 2
      MCN( MNARLI(5) + WUTYMA ) = 0
      MCN( MNARLI(5) + WBSOEF ) = 2
      MCN( MNARLI(5) + WBEFOB ) = NBST
      MNL0 = MNARLI(5) + WUSOEF
      DO 260 J=1,NBST
         MCN( MNL0     ) = J
         MCN( MNL0 + 1 ) = J+1
         MNL0 = MNL0 + 2
 260  CONTINUE
      MCN(MNL0-1) = 1
C     LE TYPE DE FERMETURE: LIGNE FERMEE
      MCN( MNARLI(5) + WUTFMA ) = 1
C     PAS DE TANGENTES STOCKEES
      MCN( MNARLI(5) + WBTGEF ) = 0
      MCN( MNARLI(5) + WBEFAP ) = 0
      MCN( MNARLI(5) + WBEFTG ) = 0
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARLI(5) + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     SI LE TABLEAU DEFINITION EXISTE ALORS IL EST DETRUIT
      CALL LXTSOU( NTLXLI(5), 'DEFINITION', NTDFLI, MNDFLI )
      IF( NTDFLI .GT. 0 ) THEN
         CALL LXTSDS( NTLXLI(5), 'DEFINITION' )
      ENDIF
C     LA LIGNE EST DEFINIE PAR SES TABLEAUX XYZSOMMET ET NSEF OPTION 10
      CALL LXTNDC( NTLXLI(5), 'DEFINITION', 'ENTIER', WUTSSL+1 )
      CALL LXTSOU( NTLXLI(5), 'DEFINITION',  NTDFLI , MNDFLI )
C     LA TRANSFORMATION
      MCN( MNDFLI + WTYTRL ) = 1
C     LE TYPE DE LA LIGNE
      MCN( MNDFLI + WUTYLI ) = 10
C     TABLEAU XYZSOMMET DES SOMMETS
      MCN( MNDFLI + WUTSOL ) = NTSOLI(5)
C     TABLEAU NSEF NO DES SOMMETS DES ARETES
      MCN( MNDFLI + WUTSSL ) = NTARLI(5)
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNDFLI) )
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOLI(5)) )
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNARLI(5)) )
C
C     ------------------------------------------
C     TRIANGULATION DU RECTANGLE PLAN PAR SUEX09
C     ------------------------------------------
C     SIMULATION DES DONNEES DU TABLEAU LADEFL POUR LA CONSTRUCTION DE LA TRIANG
C     LA TRANSFORMATION
      LADEFL( WTYTRL ) = 1
C     LE TYPE DE LA SURFACE
      LADEFL( WUTYSU ) = 9
C     variable ARETMX 'taille max des aretes des triangles equilateraux'
      RADEFL(WRETMX) = 2.6*( RMCN(MNRX+LRX) + RMCN(MNRY+LRY) ) / NBST
C     variable NBLFTR 'nombre de lignes fermees contour de la surface'
      LADEFL(WBLFTR) = 1
C     variable NBPTIT 'nombre de points internes futurs sommets'
      LADEFL(WBPTIT) = 0
C     RECHERCHE DU NUMERO DE LA LIGNE CONTOUR
C     tableau NULFTR(1..NBLFTR) 'nom des lignes fermees(enveloppe en premier)'
      CALL LXNMNO( NTLIGN , KNOMLI , I, MN )
      LADEFL(WULFTR) = I
C     TRIANGULATION EFFECTIVE
      CALL SUEX09( 3, NTLXSU, LADEFL, RADEFL,
     %             NTFASU, MNFASU, NTSOSU, MNSOSU, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     DESTRUCTION DE LA LIGNE CONTOUR DU RECTANGLE DE LA B-SPLINE
      CALL LXLXDS( NTLIGN, KNOMLI )
C
C     LE NOMBRE DE SOMMETS DE LA TRIANGULATION DU RECTANGLE PLAN
      NBSOSU = MCN( MNSOSU + WNBSOM )
C     LE NOMBRE DE TRIANGLES DE LA TRIANGULATION DU RECTANGLE PLAN
      NBEFOB = MCN( MNFASU + WBEFOB )
C
C     TRANSFORMATION DES TMS NSEF  XYZSOMMET D'APRES LA TRANSFORMATION
C     RECTANGLE_PLAN -> SURFACE B-SPLINE D'INTERPOLATION DANS R**3
C     ================================================================
      IF( LDEXSB .GT. 1 .OR. LDEYSB .GT. 1 ) THEN
C        CALCUL DES TANGENTES : 8 PAR QUADRANGLE
         NBTGEF = 8
         NBTGS  = 6 * NBEFOB
         NBEFAP = NBEFOB
         NBEFTG = NBEFOB
C        IL FAUT AUGMENTER LES TMS XYZSOMMET ET NSEF
         CALL TAMSAU( NTSOSU, WYZSOM + 3 * (NBSOSU+NBTGS) )
         CALL TAMSAU( NTFASU, WUSOEF+14*NBEFOB )
C        RE-OUVERTURE CAR L'ADRESSE MCN A PU CHANGER
         CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOSU, MNSOSU )
         CALL LXTSOU( NTLXSU, 'NSEF',      NTFASU, MNFASU )
      ELSE
         NBTGEF = 0
         NBTGS  = 0
         NBEFAP = 0
         NBEFTG = 0
      ENDIF
C
C     TABLEAU AUXILIAIRE DES GRADIENTS EN CHAQUE SOMMET DU RECTANGLE
      CALL TNMCDC( 'REEL', NBSOSU * 6, MNTGSO )
C
C     LE CALCUL DES COORDONNEES DES SOMMETS ET DES TANGENTES
C     AUX ARETES DES TRIANGLES DE CETTE SURFACE
      CALL BSPLRT( LDEXSB, LRX, RMCN(MNRX),
     %             LDEYSB, LRY, RMCN(MNRY), RMCN(MNSPQB),
     %             NBSOSU, RMCN(MNSOSU+WYZSOM),
     %             NBTGS,  RMCN(MNSOSU+WYZSOM+3*NBSOSU), RMCN(MNTGSO),
     %             NBEFOB, MCN(MNFASU+WUSOEF),
     %             MCN(MNFASU+WUSOEF+6*NBEFOB) )
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE DES GRADIENTS
      CALL TNMCDS( 'REEL', NBSOSU * 6, MNTGSO )
C
      IF( NBTGS .GT. 0 ) THEN
C        AJOUT DU POINTEUR SUR LES EF A TG ET LE CODE GEOMETRIQUE
         J = MNFASU + WUSOEF + 4*NBEFOB - 1
         DO 270 I=1,NBEFOB
C           LE POINTEUR SUR L'EF A TG
            MCN(J+I) = I
C           LE CODE GEOMETRIQUE  17 : 'SURFACE B-SPLINE'
            MCN(J+I+NBEFOB) = 17
 270     CONTINUE
      ENDIF
C
C     LE NOMBRE DE TG PAR EF
      MCN( MNFASU + WBTGEF) = NBTGEF
C     LE NOMBRE D'EF A TG
      MCN( MNFASU + WBEFTG) = NBEFTG
C     LE NOMBRE D'EF AVEC POINTEUR
      MCN( MNFASU + WBEFAP) = NBEFAP
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOSU + WNBSOM ) = NBSOSU
C     LE NOMBRE DE TANGENTES CALCULEES
      MCN( MNSOSU + WNBTGS ) = NBTGS
C     LE NOMBRE DE COORDONNEES D'UN SOMMET
      MCN( MNSOSU + WBCOOR ) = 3
C
C     AJOUT DE LA DATE
 9000 CALL ECDATE( MCN(MNFASU) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOSU) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOSU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
      IERR = 0
C
C     SUPPRESSION DES TABLEAUX AUXILIAIRES
 9900 IF( MNXYZL .GT. 0 ) CALL TNMCDS( 'REEL', NBPIEY, MNXYZL )
      IF( MNBX   .GT. 0 ) CALL TNMCDS( 'REEL', KX*(LDEXSB+2),MNBX )
      IF( MNBY   .GT. 0 ) CALL TNMCDS( 'REEL', KY*(LDEYSB+2),MNBY )
      IF( MNAX   .GT. 0 ) CALL TNMCDS( 'REEL', KX * KX, MNAX )
      IF( MNAY   .GT. 0 ) CALL TNMCDS( 'REEL', KY * KY, MNAY )
      IF( MNFACM .GT. 0 ) CALL TNMCDS( 'REEL', MAX(KX,KY), MNFACM )
      IF( MNSPQB .GT. 0 ) CALL TNMCDS( 'REEL', MOSSPL, MNSPQB )
      IF( MNUX   .GT. 0 ) CALL TNMCDS( 'REEL', MOTUX, MNUX )
      IF( MNUY   .GT. 0 ) CALL TNMCDS( 'REEL', MOTUY, MNUY )
      IF( MNTX   .GT. 0 ) CALL TNMCDS( 'REEL', MOTTX, MNTX )
      IF( MNTY   .GT. 0 ) CALL TNMCDS( 'REEL', MOTTY, MNTY )
      IF( MNRX   .GT. 0 ) CALL TNMCDS( 'REEL', LRX+1, MNRX )
      IF( MNRY   .GT. 0 ) CALL TNMCDS( 'REEL', LRY+1, MNRY )
      IF( MNNORX .GT. 0 ) CALL TNMCDS( 'ENTIER', LRX+1, MNNORX )
      IF( MNNORY .GT. 0 ) CALL TNMCDS( 'ENTIER', LRY+1, MNNORY )
      IF( MNPTSB .GT. 0 ) CALL TNMCDS( 'REEL', 3*NBINCS, MNPTSB )
C
C     ERREUR
C     ======
 9999 RETURN
      END
