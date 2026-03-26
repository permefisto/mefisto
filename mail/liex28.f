      SUBROUTINE LIEX28( NTLXLI , LADEFI , RADEFI ,
     %                   NTARLI , MNARLI , NTSOLI , MNSOLI , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UNE LIGNE CIRCULAIRE JOIGNANT DEUX
C -----    SEGMENTS PAR TANGENCE
C
C ENTREES:
C --------
C NTLXLI : NUMERO DU TABLEAU TS DU LEXIQUE DE LA LIGNE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA LIGNE
C          CF ~/TD/D/A_LIGNE__DEFINITION
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA LIGNE
C          CES 2 TABLEAUX LADEFI ET RADEFI ONT MEME ADRESSE A L'APPEL
C
C SORTIES:
C -------
C NTARLI : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ARETES
C MNARLI : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ARETES
C          CF ~/TD/D/A___NSEF
C NTSOLI : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC  PARIS      JUIN 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      REAL              XYZ(3), XYZTG(3), XYZTG1(3), XYZTG2(3),
     %                  XYZ1(3), XYZ2(3), XYZ3(3), XYZ3PT(3,3)
      EQUIVALENCE       (XYZ3PT(1,1),XYZ1(1)),
     %                  (XYZ3PT(1,2),XYZ2(1)),
     %                  (XYZ3PT(1,3),XYZ3(1))
      REAL              XY(2), XY1(2), XY2(2), XY3(2), XYF(2),
     %                  PTG1(2), PTG2(2)
      DOUBLE PRECISION  D2D3(3,3), XYZD(3), DTAILL, DTAIL2, HHH
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE(XYZ)' DES ARETES
C     ===============================================================
C     ICI LA CARTE EST SUPPOSEE ISOTROPE
      NOFOTI = NOFOTIEL()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
C
C     VERIFICATIONS
C     -------------
C     LE NOMBRE D'ARETES ET DE SOMMETS DE LA LIGNE
      NBARLI = LADEFI(WBARLI)
      IF( NOFOTI .GT. 0 ) THEN
C        LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' EXISTE
         IF( NBARLI .LT. 2 ) THEN
            NBARLI = 3
         ENDIF
      ELSE IF( NBARLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE INCORRECT D''ARETES'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
      NBSOLI = NBARLI + 1
      NBTGS  = 0
      NA23   = 0
C
C     LE RAYON
      RAYCOL = RADEFI( WAYCOL )
      IF( RAYCOL .LE. 0. ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAYON=<0'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     LE NUMERO UTILISATEUR DES 3 POINTS DE DEFINITION DU CONGE
C     ---------------------------------------------------------
      DO 5 I=1,3
         IF( I .EQ. 1 ) THEN
            N = LADEFI( WUP1CL )
         ELSE IF( I .EQ. 2 ) THEN
            N = LADEFI( WUP2CL )
         ELSE IF( I .EQ. 3 ) THEN
            N = LADEFI( WUP3CL )
         ENDIF
C        RECUPERATION DES COORDONNEES DU POINT
         CALL LXNLOU( NTPOIN , N , NTLXPO , MN )
         IF( NTLXPO .LE. 0 ) THEN
             NBLGRC(NRERR) = 1
             WRITE(KERR(MXLGER)(1:4),'(I4)') I
             KERR(1) = 'POINT CONGE ' // KERR(MXLGER)(1:4) //' INCONNU'
             CALL LEREUR
             IERR = 6
             GOTO 9999
          ENDIF
          CALL LXTSOU( NTLXPO , 'XYZSOMMET' , NT , MN )
          MN = MN + WYZSOM - 1
          DO 4 J=1,3
             XYZ3PT(J,I) = RMCN( MN + J )
 4        CONTINUE
 5    CONTINUE
C
C     DISTANCE ENTRE CES 3 POINTS
      DO 10 I=1,3
         IF( I.NE. 3 ) THEN
            N = I+1
         ELSE
            N = 1
         ENDIF
         H = DIST2P( XYZ3PT(1,I) , XYZ3PT(1,N) )
         IF( H .LE. 0. ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = '2 DES 3 POINTS DU CONGE SONT CONFONDUS'
            CALL LEREUR
            IERR = 3
            GOTO 9999
         ENDIF
 10   CONTINUE
C
C     PASSAGE AUX COORDONNEES DES 3 POINTS DANS LEUR PLAN
C     ---------------------------------------------------
      CALL DF3D2D( XYZ3PT(1,1), XYZ3PT(1,2), XYZ3PT(1,3),
     %             D2D3, IERR )
      IF( IERR .GT. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = '3 POINTS DU CONGE COLINEAIRES'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C
C     LES 3 POINTS DANS LEUR PLAN
C     ---------------------------
C     LE POINT INITIAL A POUR COORDONNEES 2D : 0. 0.
      XY1(1) = 0.
      XY1(2) = 0.
      CALL CH3D2D( XYZ1, D2D3, XYZ2, XY2 )
      CALL CH3D2D( XYZ1, D2D3, XYZ3, XY3 )
C
C     CALCUL DE L'ANGLE DU CONGE
C     --------------------------
C     IL NE PEUT ETRE NUL CAR LES DISTANCES ENTRE POINTS SONT >0
      PI   = ATAN(1.0) * 4
      X2X3 = XY3(1) - XY2(1)
      X2X1 = XY1(1) - XY2(1)
      Y2Y3 = XY3(2) - XY2(2)
      Y2Y1 = XY1(2) - XY2(2)
      D23  = SQRT( X2X3 * X2X3 + Y2Y3 * Y2Y3 )
      D21  = SQRT( X2X1 * X2X1 + Y2Y1 * Y2Y1 )
      ALPHA = ( X2X3 * X2X1 + Y2Y3 * Y2Y1 ) / (D23 * D21)
      ALPHA = ACOS( ALPHA )
C
C     ANGLE AU CENTRE DU CERCLE = SUPPLEMENTAIRE DE ALPHA ANGLE DU CONGE
      ALPHAC = PI - ALPHA
C
C     CALCUL DES POINTS DE TANGENCE DU CONGE AVEC LES DROITES
C     ET DU CENTRE (XC,YC) DU CERCLE DANS LE PLAN DES 3 POINTS
C     --------------------------------------------------------
      D  = RAYCOL / TAN( ALPHA * 0.5 )
      XC = XY2(1) - D
      IF( XC .LE. 0. .OR. XC .GE. XY2(1) ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAYON TROP GRAND'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
      YC = RAYCOL
      PTG1(1) = XC
      PTG1(2) = 0.0
C     RETOUR AUX COORDONNEES 3D DU POINT PTG1
      CALL CH2D3D( XYZ1, D2D3, PTG1, XYZTG1 )
      PTG2(1) = XY2(1) - D * COS( ALPHA )
      PTG2(2) = D * SIN( ALPHA )
C     RETOUR AUX COORDONNEES 3D DU POINT PTG2
      CALL CH2D3D( XYZ1, D2D3, PTG2, XYZTG2 )
C
C     CALCUL DU NOMBRE D'ARETES POUR LES 2 SEGMENTS ET L'ARC
C     ------------------------------------------------------
C     LA LONGUEUR DES 2 SEGMENTS ET DE L'ARC
      DD = D21 + D23 - D - D + RAYCOL * ALPHAC
      DM = DD / NBARLI
      IF( NOFOTI .LE. 0 ) THEN
C        PAS DE FONCTION TAILLE_IDEALE
         NA12  = NINT( (D21-D) / DM )
         NA23  = NINT( (D23-D) / DM )
         NARC  = NBARLI - NA12 - NA23
         NBTGS = NARC + 1
      ELSE
C        LA FONCTION 'TAILLE_IDEALE' EXISTE
C        CALCUL DU NOMBRE D'ARCS DE CERCLE POUR DECLARER LES TABLEAUX
         NUMARC  = 1
         NBARLI  = 0
         HHH     = D21 - D
C        L'ANGLE AU CENTRE DU CERCLE
         ANGLE   = -PI / 2
         XY(1)   = PTG1(1)
         XY(2)   = PTG1(2)
         XYZD(1) = XYZTG1(1)
         XYZD(2) = XYZTG1(2)
         XYZD(3) = XYZTG1(3)
C        LES 3 PARAMETRES D'APPEL DE LA FONCTION 'TAILLE_IDEALE'
C        AU SOMMET XYZ DU SEGMENT SOMMET1 SOMMET2
 15      CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAILL )
         IF( NCODEV .EQ. 0 ) THEN
C
C           NCODEV  : 0 DTAILL N'EST PAS INITIALISEE EN SORTIE
C                     1 DTAILL   EST     INITIALISEE EN SORTIE
            NBLGRC(NRERR) = 2
            KERR(1) = 'FONCTION TAILLE_IDEALE INCALCULABLE AU POINT'
            KERR(2) = ' '
            WRITE(KERR(2)(1:15), '(E15.7)') XYZD(1)
            WRITE(KERR(2)(18:32),'(E15.7)') XYZD(2)
            WRITE(KERR(2)(35:49),'(E15.7)') XYZD(3)
            CALL LEREUR
            IERR = 5
            GOTO 9999
C
         ELSE
C           ESSAI DE CREER UNE ARETE DE LONGUEUR CETTE TAILLE
            IF( NUMARC .EQ. 1 ) THEN
C
               DTAILL = ABS( DTAILL )
               IF( DTAILL .LT. HHH*0.65 ) THEN
C                 MAILLAGE EN SEGMENTS DROITS DE PTG1 A P1
C                 CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
                  NBARLI = NBARLI + 1
                  IPAS   = 0
C                 LES COORDONNEES 2D DU POINT FINAL DE L'ARETE
 16               XYF(1) = REAL( XY(1) - DTAILL )
                  XYF(2) = XY(2)
C                 RETOUR AUX COORDONNEES 3D
                  CALL CH2D3D( XYZ1, D2D3, XYF, XYZ )
                  IF( IPAS .LE. 4 ) THEN
                     XYZD(1) = XYZ(1)
                     XYZD(2) = XYZ(2)
                     XYZD(3) = XYZ(3)
                     CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                     DTAIL2 = ABS( DTAIL2 )
                     IF( DTAIL2 .LT. DTAILL ) THEN
                        DTAILL = DTAIL2
                        IPAS   = IPAS + 1
                        GOTO 16
                     ENDIF
                  ENDIF
C                 LA LONGUEUR RESTANTE
                  HHH = HHH - DTAILL
                  XY(1) = REAL( XY(1) - DTAILL )
C                 LE SOMMET SUR LE SEGMENT
                  XYZD(1) = XYZ(1)
                  XYZD(2) = XYZ(2)
                  XYZD(3) = XYZ(3)
               ELSE
C                 PASSAGE A L'ARC DE CERCLE DE PTG1 A PTG2
                  NBARLI = NBARLI + 1
                  NA12   = NBARLI
                  NUMARC = 2
                  HHH    = RAYCOL * ALPHAC
                  XY(1)  = PTG1(1)
                  XY(2)  = PTG1(2)
C                 LE PREMIER POINT DE TANGENCE SUR LE CERCLE
                  XYZD(1) = XYZTG1(1)
                  XYZD(2) = XYZTG1(2)
                  XYZD(3) = XYZTG1(3)
               ENDIF
               GOTO 15
C
            ELSE IF( NUMARC .EQ. 2 ) THEN
C              MAILLAGE EN ARCS P3 DE L'ARC PTG1 A PTG2
               IF( DTAILL .LT. HHH*0.65 ) THEN
C                 CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
                  NBARLI = NBARLI + 1
                  IPAS   = 0
                  ANGLE0 = ANGLE
C                 LES COORDONNEES 2D DU POINT FINAL DE L'ARC
 17               ANGLE = REAL( ANGLE0 + DTAILL / RAYCOL )
                  XY(1) = REAL( XC + RAYCOL * COS( ANGLE ) )
                  XY(2) = REAL( YC + RAYCOL * SIN( ANGLE ) )
C                 RETOUR AUX COORDONNEES 3D
                  CALL CH2D3D( XYZ1, D2D3 , XY , XYZ )
                  IF( IPAS .LE. 4 ) THEN
                     XYZD(1) = XYZ(1)
                     XYZD(2) = XYZ(2)
                     XYZD(3) = XYZ(3)
                     CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                     DTAIL2 = ABS( DTAIL2 )
                     IF( DTAIL2 .LT. DTAILL ) THEN
                        DTAILL = DTAIL2
                        IPAS   = IPAS + 1
                        GOTO 17
                     ENDIF
                  ENDIF
C                 LA LONGUEUR RESTANTE
                  HHH = HHH - DTAILL
C                 LE SOMMET SUR LE CERCLE
                  XYZD(1) = XYZ(1)
                  XYZD(2) = XYZ(2)
                  XYZD(3) = XYZ(3)
               ELSE
C                 PASSAGE AU SEGMENT PTG2 P3
                  NBARLI = NBARLI + 1
                  NUMARC = 3
                  HHH    = D23 - D
                  NARC   = NBARLI - NA12
                  XY(1)  = PTG2(1)
                  XY(2)  = PTG2(2)
C                 LE SECOND POINT DE TANGENCE SUR LE CERCLE
                  XYZD(1) = XYZTG2(1)
                  XYZD(2) = XYZTG2(2)
                  XYZD(3) = XYZTG2(3)
                  XDIR = XY3(1) - PTG2(1)
                  YDIR = XY3(2) - PTG2(2)
                  DD   = SQRT( XDIR*XDIR + YDIR*YDIR )
C                 LA TANGENTE UNITAIRE DE PTG2 VERS P3
                  XDIR = XDIR / DD
                  YDIR = YDIR / DD
               ENDIF
               GOTO 15
C
            ELSE IF( NUMARC .EQ. 3 ) THEN
C
               IF( DTAILL .LT. HHH*0.65 ) THEN
C                 MAILLAGE EN SEGMENTS DROITS DE PTG1 A P1
C                 CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
                  NBARLI = NBARLI + 1
                  IPAS   = 0
C                 LES COORDONNEES 2D DU POINT FINAL DE L'ARETE
 18               XYF(1) = REAL( XY(1) + DTAILL * XDIR )
                  XYF(2) = REAL( XY(2) + DTAILL * YDIR )
C                 RETOUR AUX COORDONNEES 3D
                  CALL CH2D3D( XYZ1, D2D3, XYF, XYZ )
                  IF( IPAS .LE. 4 ) THEN
                     XYZD(1) = XYZ(1)
                     XYZD(2) = XYZ(2)
                     XYZD(3) = XYZ(3)
                     CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                     DTAIL2 = ABS( DTAIL2 )
                     IF( DTAIL2 .LT. DTAILL ) THEN
                        DTAILL = DTAIL2
                        IPAS   = IPAS + 1
                        GOTO 18
                     ENDIF
                  ENDIF
C                 LA LONGUEUR RESTANTE
                  HHH = HHH - DTAILL
                  XY(1) = REAL( XY(1) + DTAILL * XDIR )
                  XY(2) = REAL( XY(2) + DTAILL * YDIR )
C                 LE SOMMET SUR LE SEGMENT
                  XYZD(1) = XYZ(1)
                  XYZD(2) = XYZ(2)
                  XYZD(3) = XYZ(3)
                  GOTO 15
               ELSE
C                 FIN DE LA LIGNE
                  NBARLI = NBARLI + 1
                  NA23   = NBARLI - NA12 - NARC
                  NBTGS  = 2 * NARC
               ENDIF
            ENDIF
C
         ENDIF
         IF( NBARLI .LE. 1 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  'NOMBRE INCORRECT D''ARETES'
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
      ENDIF
C
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBSOLI = NBARLI + 1
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     -----------------------------------
      CALL LXTNDC( NTLXLI , 'XYZSOMMET' , 'ENTIER' ,
     %             WYZSOM+3*NBSOLI+3*NBTGS )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOLI + WNBSOM ) = NBSOLI
C
C     LE NOMBRE DE TANGENTES
      MCN( MNSOLI + WNBTGS ) = NBTGS
C
C     LE PREMIER SOMMET DE LA LIGNE EST P1
      MNS  = MNSOLI + WYZSOM
      MNTG = MNS + 3 * NBSOLI
      RMCN(MNS  ) = XYZ1(1)
      RMCN(MNS+1) = XYZ1(2)
      RMCN(MNS+2) = XYZ1(3)
      MNS = MNS + 3
C
      IF( NOFOTI .LE. 0 ) THEN
C
C        PAS DE FONCTION TAILLE_IDEALE => ARETES DE TAILLES EGALES
C        *****************************
         XDIR = ( XYZTG1(1) - XYZ1(1) ) / NA12
         YDIR = ( XYZTG1(2) - XYZ1(2) ) / NA12
         ZDIR = ( XYZTG1(3) - XYZ1(3) ) / NA12
         DO 50 I=1,NA12-1
C           LE SOMMET FINAL DE L'ARETE I DU SEGMENT P1 => PTG1
            RMCN(MNS  ) = XYZ1(1) + I * XDIR
            RMCN(MNS+1) = XYZ1(2) + I * YDIR
            RMCN(MNS+2) = XYZ1(3) + I * ZDIR
            MNS = MNS + 3
 50      CONTINUE
C
C        LE PREMIER SOMMET DE L'ARC DE CERCLE EST LE POINT DE TANGENCE PTG1
         RMCN(MNS  ) = XYZTG1(1)
         RMCN(MNS+1) = XYZTG1(2)
         RMCN(MNS+2) = XYZTG1(3)
C        LA TANGENTE AU POINT DE TANGENCE 1 EST HORIZONTALE EN 2D
         ANGLE  = -PI / 2
         ALFARC = ALPHAC / NARC
C        LE POINT PTG1 + TG(PTG1)
         XY2(1) = PTG1(1) + RAYCOL * ALFARC
         XY2(2) = PTG1(2)
C        RETOUR AUX COORDONNEES 3D
         CALL CH2D3D( XYZ1, D2D3, XY2, XYZ )
         RMCN(MNTG  ) = XYZ(1) - RMCN(MNS  )
         RMCN(MNTG+1) = XYZ(2) - RMCN(MNS+1)
         RMCN(MNTG+2) = XYZ(3) - RMCN(MNS+2)
         MNTG = MNTG + 3
         MNS  = MNS  + 3
C
         XY(1)   = PTG1(1)
         XY(2)   = PTG1(2)
         DO 60 I=1,NARC-1
C           LES COORDONNEES 2D DU POINT FINAL DE L'ARC
            ANGLE = ANGLE + ALFARC
            XY(1) = XC + RAYCOL * COS( ANGLE )
            XY(2) = YC + RAYCOL * SIN( ANGLE )
C           RETOUR AUX COORDONNEES 3D
            CALL CH2D3D( XYZ1, D2D3, XY, RMCN(MNS) )
C           LA TANGENTE AU POINT FINAL DE L'ARC I
            XY2(1) = XY(1) - RAYCOL * ALFARC * SIN( ANGLE )
            XY2(2) = XY(2) + RAYCOL * ALFARC * COS( ANGLE )
C           RETOUR AUX COORDONNEES 3D
            CALL CH2D3D( XYZ1, D2D3, XY2, XYZ )
            RMCN(MNTG  ) = XYZ(1) - RMCN(MNS  )
            RMCN(MNTG+1) = XYZ(2) - RMCN(MNS+1)
            RMCN(MNTG+2) = XYZ(3) - RMCN(MNS+2)
            MNTG = MNTG + 3
            MNS  = MNS  + 3
 60      CONTINUE
C
C        LE DERNIER SOMMET DE L'ARC EST LE POINT DE TANGENCE PTG2
         RMCN(MNS  ) = XYZTG2(1)
         RMCN(MNS+1) = XYZTG2(2)
         RMCN(MNS+2) = XYZTG2(3)
C        LE POINT PTG2 + TG(PTG2)
         ANGLE  = ANGLE + ALFARC
         XY2(1) = PTG2(1) - RAYCOL * ALFARC * SIN( ANGLE )
         XY2(2) = PTG2(2) + RAYCOL * ALFARC * COS( ANGLE )
C        RETOUR AUX COORDONNEES 3D
         CALL CH2D3D( XYZ1, D2D3, XY2, XYZ )
         RMCN(MNTG  ) = XYZ(1) - RMCN(MNS  )
         RMCN(MNTG+1) = XYZ(2) - RMCN(MNS+1)
         RMCN(MNTG+2) = XYZ(3) - RMCN(MNS+2)
         MNTG = MNTG + 3
         MNS  = MNS  + 3
C
         XDIR = ( XYZTG2(1) - XYZ3(1) ) / NA23
         YDIR = ( XYZTG2(2) - XYZ3(2) ) / NA23
         ZDIR = ( XYZTG2(3) - XYZ3(3) ) / NA23
         DO 70 I=NA23-1,1,-1
C           LE SOMMET FINAL DE L'ARETE I DU SEGMENT PTG2 => P3
            RMCN(MNS  ) = XYZ3(1) + I * XDIR
            RMCN(MNS+1) = XYZ3(2) + I * YDIR
            RMCN(MNS+2) = XYZ3(3) + I * ZDIR
            MNS = MNS + 3
 70      CONTINUE
C
      ELSE
C
C        LA FONCTION 'TAILLE_IDEALE' EXISTE
C        **********************************
C        CALCUL DU NOMBRE D'ARCS DE CERCLE POUR DECLARER LES TABLEAUX
         NUMARC  = 1
         N       = NA12
         HHH     = D21 - D
C        L'ANGLE AU CENTRE DU CERCLE
         ANGLE   = -PI / 2
         XY(1)   = PTG1(1)
         XY(2)   = PTG1(2)
C        RETOUR AUX COORDONNEES 3D
         CALL CH2D3D( XYZ1, D2D3, PTG1, XYZ )
         XYZD(1) = XYZ(1)
         XYZD(2) = XYZ(2)
         XYZD(3) = XYZ(3)
C        LES 3 PARAMETRES D'APPEL DE LA FONCTION 'TAILLE_IDEALE'
C        AU SOMMET XYZ DU SEGMENT SOMMET1 SOMMET2
 115     CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAILL )
         DTAILL = ABS( DTAILL )
C        PAS DE TEST SUR NCODEV CAR DEJA FAIT POUR CALCULER NBARLI
C        CREATION D'UNE ARETE DE LONGUEUR CETTE TAILLE
         IF( NUMARC .EQ. 1 ) THEN
C
            IF( DTAILL .LT. HHH*0.65 ) THEN
C              MAILLAGE EN SEGMENTS DROITS DE PTG1 A P1
               IPAS   = 0
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARETE
 116           XYF(1) = REAL( XY(1) - DTAILL )
               XYF(2) = XY(2)
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3, XYF, XYZ )
               IF( IPAS .LE. 4 ) THEN
                  XYZD(1) = XYZ(1)
                  XYZD(2) = XYZ(2)
                  XYZD(3) = XYZ(3)
                  CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                  DTAIL2 = ABS( DTAIL2 )
                  IF( DTAIL2 .LT. DTAILL ) THEN
                     DTAILL = DTAIL2
                     IPAS   = IPAS + 1
                     GOTO 116
                  ENDIF
               ENDIF
C
               N = N - 1
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARETE
               XY(1) = REAL( XY(1) - DTAILL )
C              XY(2) = 0
C              RETOUR AUX COORDONNEES 3D
               MN = MNSOLI + WYZSOM + 3 * N
               CALL CH2D3D( XYZ1, D2D3, XY, RMCN(MN) )
C              LA LONGUEUR RESTANTE
               HHH = HHH - DTAILL
C              LE SOMMET SUR LE SEGMENT
               XYZD(1) = RMCN(MN  )
               XYZD(2) = RMCN(MN+1)
               XYZD(3) = RMCN(MN+2)
            ELSE
C              PASSAGE A L'ARC DE CERCLE DE PTG1 A PTG2
               NUMARC = 2
               HHH    = RAYCOL * ALPHAC
C              LE PREMIER SOMMET DE L'ARC EST LE POINT DE TANGENCE PTG1
               XYZD(1) = XYZTG1(1)
               XYZD(2) = XYZTG1(2)
               XYZD(3) = XYZTG1(3)
               MNS     = MNSOLI + WYZSOM + 3 * NA12
               RMCN(MNS  ) = XYZTG1(1)
               RMCN(MNS+1) = XYZTG1(2)
               RMCN(MNS+2) = XYZTG1(3)
C              LE PREMIER SOMMET DE L'ARC EST LE POINT DE TANGENCE PTG1
C              LA TANGENTE AU POINT DE TANGENCE 1 EST HORIZONTALE EN 2D
               ANGLE = -PI / 2
C              LE POINT PTG1 + TG(PTG1)
               XY2(1) = PTG1(1) + RAYCOL
               XY2(2) = PTG1(2)
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3, XY2, XYZ )
               RMCN(MNTG  ) = XYZ(1) - RMCN(MNS  )
               RMCN(MNTG+1) = XYZ(2) - RMCN(MNS+1)
               RMCN(MNTG+2) = XYZ(3) - RMCN(MNS+2)
               MNS = MNS + 3
            ENDIF
            GOTO 115
C
         ELSE IF( NUMARC .EQ. 2 ) THEN
C           MAILLAGE EN ARCS P3 DE L'ARC PTG1 A PTG2
            IF( DTAILL .LT. HHH*0.65 ) THEN
C              CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
               IPAS   = 0
               ANGLE0 = ANGLE
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARC
 117           ANGLE  = REAL( ANGLE0 + DTAILL / RAYCOL )
               XY(1) = XC + RAYCOL * COS( ANGLE )
               XY(2) = YC + RAYCOL * SIN( ANGLE )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3 , XY , XYZ )
               IF( IPAS .LE. 4 ) THEN
                  XYZD(1) = XYZ(1)
                  XYZD(2) = XYZ(2)
                  XYZD(3) = XYZ(3)
                  CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                  DTAIL2 = ABS( DTAIL2 )
                  IF( DTAIL2 .LT. DTAILL ) THEN
                     DTAILL = DTAIL2
                     IPAS   = IPAS + 1
                     GOTO 117
                  ENDIF
               ENDIF
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARC
               DIFANG = REAL( DTAILL / RAYCOL )
               ANGLE  = ANGLE0 + DIFANG
               XY(1)  = XC + RAYCOL * COS( ANGLE )
               XY(2)  = YC + RAYCOL * SIN( ANGLE )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3, XY, RMCN(MNS) )
C              LA LONGUEUR RESTANTE
               HHH = HHH - DTAILL
C              LE SOMMET SUR LE CERCLE
               XYZD(1) = RMCN(MNS  )
               XYZD(2) = RMCN(MNS+1)
               XYZD(3) = RMCN(MNS+2)
C
C              LA TANGENTE INITIALE EST MULTIPLIEE PAR L'ANGLE DE L'ARC=ARETE
               RMCN(MNTG  ) = RMCN(MNTG  ) * DIFANG
               RMCN(MNTG+1) = RMCN(MNTG+1) * DIFANG
               RMCN(MNTG+2) = RMCN(MNTG+2) * DIFANG
               MNTG = MNTG + 3
C
C              LA TANGENTE AU POINT FINAL DE L'ARC I (SANS L'ANGLE)
               XY2(1) = XY(1) - RAYCOL * SIN( ANGLE )
               XY2(2) = XY(2) + RAYCOL * COS( ANGLE )
C              RETOUR AUX COORDONNEES 3D DE LA TANGENTE
               CALL CH2D3D( XYZ1, D2D3, XY2, XYZTG )
               XYZTG(1) = XYZTG(1) - RMCN(MNS  )
               XYZTG(2) = XYZTG(2) - RMCN(MNS+1)
               XYZTG(3) = XYZTG(3) - RMCN(MNS+2)
C              LA TANGENTE FINALE DE L'ARC
               RMCN(MNTG  ) = - XYZTG(1) * DIFANG
               RMCN(MNTG+1) = - XYZTG(2) * DIFANG
               RMCN(MNTG+2) = - XYZTG(3) * DIFANG
               MNTG = MNTG + 3
C              LA TANGENTE INITIALE DE L'ARC SUIVANT
               RMCN(MNTG  ) = XYZTG(1)
               RMCN(MNTG+1) = XYZTG(2)
               RMCN(MNTG+2) = XYZTG(3)
               MNS = MNS  + 3
            ELSE
C              PASSAGE AU SEGMENT PTG2 P3
               NUMARC = 3
               XY(1)  = PTG2(1)
               XY(2)  = PTG2(2)
               HHH    = D23 - D
C              LA TANGENTE INITIALE DU DERNIER ARC DE CERCLE
C              EST MULTIPLIEE PAR L'ANGLE DE L'ARC=ARETE
               DIFANG = ALPHAC - ( ANGLE + PI / 2 )
               RMCN(MNTG  ) = RMCN(MNTG  ) * DIFANG
               RMCN(MNTG+1) = RMCN(MNTG+1) * DIFANG
               RMCN(MNTG+2) = RMCN(MNTG+2) * DIFANG
               MNTG = MNTG + 3
C              LE DERNIER SOMMET DE L'ARC DE CERCLE EST LE POINT DE TANGENCE PTG
               XYZD(1) = XYZTG2(1)
               XYZD(2) = XYZTG2(2)
               XYZD(3) = XYZTG2(3)
               RMCN(MNS  ) = XYZTG2(1)
               RMCN(MNS+1) = XYZTG2(2)
               RMCN(MNS+2) = XYZTG2(3)
C              LE POINT PTG2 + TG(PTG2)
               ANGLE  = ALPHAC - PI/2
               XY2(1) = PTG2(1) - RAYCOL * SIN( ANGLE )
               XY2(2) = PTG2(2) + RAYCOL * COS( ANGLE )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3, XY2, XYZ )
C              LA DERNIERE TANGENTE DU DERNIER SOMMET DE L'ARC
               RMCN(MNTG  ) = -( XYZ(1) - RMCN(MNS  ) ) * DIFANG
               RMCN(MNTG+1) = -( XYZ(2) - RMCN(MNS+1) ) * DIFANG
               RMCN(MNTG+2) = -( XYZ(3) - RMCN(MNS+2) ) * DIFANG
               MNS = MNS  + 3
C
               XDIR = XY3(1) - PTG2(1)
               YDIR = XY3(2) - PTG2(2)
               DD   = SQRT( XDIR*XDIR + YDIR*YDIR )
C              LA TANGENTE UNITAIRE DE PTG2 VERS P3
               XDIR = XDIR / DD
               YDIR = YDIR / DD
            ENDIF
            GOTO 115
C
         ELSE IF( NUMARC .EQ. 3 ) THEN
C
C           MAILLAGE EN SEGMENTS DROITS DE PTG2 A P3
            IF( DTAILL .LT. HHH*0.65 ) THEN
               IPAS   = 0
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARETE
 118           XYF(1) = REAL( XY(1) + DTAILL * XDIR )
               XYF(2) = REAL( XY(2) + DTAILL * YDIR )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3, XYF, XYZ )
               IF( IPAS .LE. 4 ) THEN
                  XYZD(1) = XYZ(1)
                  XYZD(2) = XYZ(2)
                  XYZD(3) = XYZ(3)
                  CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                  DTAIL2 = ABS( DTAIL2 )
                  IF( DTAIL2 .LT. DTAILL ) THEN
                     DTAILL = DTAIL2
                     IPAS   = IPAS + 1
                     GOTO 118
                  ENDIF
               ENDIF
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARETE
               XY(1) = REAL( XY(1) + DTAILL * XDIR )
               XY(2) = REAL( XY(2) + DTAILL * YDIR )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( XYZ1, D2D3, XY, RMCN(MNS) )
C              LA LONGUEUR RESTANTE
               HHH = HHH - DTAILL
C              LE SOMMET SUR LE SEGMENT DE DROITE PTG2 P3
               XYZD(1) = RMCN(MNS)
               XYZD(2) = RMCN(MNS+1)
               XYZD(3) = RMCN(MNS+2)
               MNS = MNS + 3
               GOTO 115
            ENDIF
         ENDIF
      ENDIF
C
C     LE DERNIER SOMMET DE LA LIGNE EST P3
      MNS = MNSOLI + WYZSOM + 3 * NBSOLI - 3
      RMCN(MNS  ) = XYZ3(1)
      RMCN(MNS+1) = XYZ3(2)
      RMCN(MNS+2) = XYZ3(3)
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOLI + WBCOOR ) = 3
      MCN( MNSOLI + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' LIGNE STRUCTUREE
C     -----------------------------------------------
      CALL LXTNDC( NTLXLI , 'NSEF' , 'ENTIER' ,
     %             WBARSE+1+NBARLI+3*NARC )
      CALL LXTSOU( NTLXLI , 'NSEF' ,  NTARLI  , MNARLI )
C
C     LE TYPE DE L'OBJET : ICI LIGNE
      MCN( MNARLI + WUTYOB ) = 2
C     AJOUT DE LIGNE NON-FERMEE
      MCN( MNARLI + WUTFMA ) = 0
C     AJOUT DU NOMBRE DE SOMMET PAR EF
      MCN( MNARLI + WBSOEF ) = 2
C     AJOUT DU NOMBRE DE TANGENTES PAR EF
      MCN( MNARLI + WBTGEF ) = 2
C     AJOUT DU NOMBRE EF DU PLSV
      MCN ( MNARLI + WBEFOB ) = NBARLI
C     AJOUT DU NOMBRE D'EF A TG
      MCN ( MNARLI + WBEFTG ) = NARC
C     AJOUT DU NOMBRE D'EF POINTES
      MCN ( MNARLI + WBEFAP ) = NBARLI
C     LE TYPE DU MAILLAGE : ICI SEGMENT STRUCTURE
      MCN( MNARLI + WUTYMA ) = 2
C     LE NOMBRE D'ARETES DU SEGMENT STRUCTURE
      MCN( MNARLI + WBARSE ) = NBARLI
C
C     LE POINTEUR SUR CHAQUE EF
      MN  = MNARLI + WBARSE
      MNG = MN + NBARLI
      DO 146 I=1,NA12
C        EF SANS TG
         MCN(MN+I) = 0
 146  CONTINUE
      DO 147 I=1,NARC
C        EF AVEC TG
         MCN(MN+NA12+I) = I
 147  CONTINUE
      DO 148 I=NA12+NARC+1,NBARLI
C        EF SANS TG
         MCN(MN+I) = 0
 148  CONTINUE
C
C     LE CODE GEOMETRIQUE DES EF A TG
      MN = MN + NBARLI
      DO 149 I=1,NARC
C        LE CODE GEOMETRIQUE : CERCLE => 1
         MCN(MN+I) = 1
 149  CONTINUE
C
C     AJOUT DU NUMERO DES TANGENTES DES EF A TG
C     L'ARC DE CERCLE ET SES 2 TANGENTES PAR ARETE
      MN = MN + NARC
      IF( NBTGS .EQ. NARC+1 ) THEN
C        ARCS EGAUX
         MN = MN + 1
         DO 151 I=1,NARC
            MCN(MN  ) =   I
            MCN(MN+1) = -(I+1)
            MN = MN + 2
 151     CONTINUE
      ELSE
C        ARCS NON EGAUX
         DO 152 I=1,2*NARC
            MCN(MN+I) = I
 152     CONTINUE
      ENDIF
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARLI) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARLI + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
      IF( IERR .NE. 0 ) THEN
         CALL LXTSDS( NTLXLI , 'XYZSOMMET' )
         CALL LXTSDS( NTLXLI , 'NSEF' )
         GOTO 9999
      ENDIF
C
C     ARRIVEE ICI EN CAS D'ERREUR
9999  RETURN
      END
