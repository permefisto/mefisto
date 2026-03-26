      SUBROUTINE SUEX48( NTLXVS,  LADEFI,
     %                   NTNSEFS, MNNSEFS, NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPOSER DES POINTS COMME SOMMETS D'UN MAILLAGE DE SURFACE
C ----- LE PLUS PROCHE SOMMET DU POINT HERITE DES COORDONNEES DU POINT
C       ATTENTION, LA QUALITE DES ELEMENTS FINIS SERA MODIFIEE
C
C ENTREES:
C --------
C NTLXVS : NUMERO DU TABLEAU TMS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF ~/TD/D/A_SURFACE__DEFINITION
C
C SORTIES:
C --------
C NTNSEFS: NUMERO      DU TMS 'NSEF' DE LA SURFACE
C MNNSEFS: ADRESSE MCN DU TMS 'NSEF' DE LA SURFACE
C          CF ~/TD/D/A___NSEF
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~/TD/D/A___XYZSOMMET
C IERR   : 0 SI PAS D'ERREUR
C          1 SI SURFACE INITIAL SANS NSEF
C          2 SI SURFACE INITIAL SANS XYZSOMMET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC SAINT PIERRE DU PERRAY Juillet 2012
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      DOUBLE PRECISION  DMIN, D
C
C     LA SURFACE INITIALE
C     ===================
C     LE NOM DE CETTE SURFACE
      NUSUIN = LADEFI( WUSUIN )
C     LE TABLEAU LEXIQUE DE CETTE SURFACE
      CALL LXNLOU( NTSURF, NUSUIN, NTLXSI, MN )
      IF( NTLXSI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE INITIALE INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN INITIAL SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSI, 'NSEF', NTNSEFI, MNNSEFI )
      IF( NTNSEFI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE SANS TMS NSEF'
         ELSE
            KERR(1) = 'SURFACE WITHOUT TMS NSEF'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSI, 'XYZSOMMET', NTXYZI, MNXYZI )
      IF( NTXYZI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE SANS TMS XYZSOMMET'
         ELSE
            KERR(1) = 'SURFACE WITHOUT TMS XYZSOMMET'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      NBSOM  = MCN( MNXYZI + WNBSOM )
      NBTGS  = MCN( MNXYZI + WNBTGS )
      NBCOOR = MCN( MNXYZI + WBCOOR )
C
C     COPIE DU TMS XYZSOMMET
      MOXYZS = WYZSOM + NBCOOR * ( NBSOM + NBTGS )
      CALL LXTNDC( NTLXVS, 'XYZSOMMET', 'MOTS', MOXYZS )
      CALL LXTSOU( NTLXVS, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .LE. 0 ) GOTO 9900
      CALL TRTATA( MCN(MNXYZI), MCN(MNXYZS), MOXYZS )
C
C     LE NOMBRE DE POINTS A IMPOSER COMME SOMMETS
      NBPOSS = LADEFI( WBPOSS )
      DO K = 1, NBPOSS
C
C        LE NOM DU POINT DEVENANT SOMMET
         NUPOSS = LADEFI( WUPOSS - 1 + K )
C
C        LE TABLEAU LEXIQUE DE CE POINT
         CALL LXNLOU( NTPOIN, NUPOSS, NTPOSS, MNPOSS )
         IF( NTPOSS .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'POINT INCONNU'
            ELSE
               KERR(1) = 'UNKNOWN POINT'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9999
         ENDIF
C
C        LE TABLEAU 'XYZSOMMET' DE CE POINT
         CALL LXTSOU( NTPOSS, 'XYZSOMMET', NTXYZP, MNXYZP )
         IF( NTXYZP .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'POINT SANS TMS XYZSOMMET'
            ELSE
               KERR(1) = 'POINT WITHOUT TMS XYZSOMMET'
            ENDIF
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF
C
C        RECHERCHE DU PLUS PROCHE SOMMET DU POINT
         MN = MNXYZP + WYZSOM
         XP = RMCN( MN     )
         YP = RMCN( MN + 1 )
         ZP = RMCN( MN + 2 )
         DMIN = 1D100
         MN = MNXYZS + WYZSOM
         DO N = 1, NBSOM
            D = ( XP - RMCN(MN)   ) ** 2
     %        + ( YP - RMCN(MN+1) ) ** 2
     %        + ( ZP - RMCN(MN+2) ) ** 2
            MN = MN + NBCOOR
            IF( D .LT. DMIN ) THEN
               DMIN = D
               NMIN = N
            ENDIF
         ENDDO
C
C        LES 3 COORDONNEES DU POINT SONT IMPOSEES AU SOMMET NMIN
         MN = MNXYZS + WYZSOM + NBCOOR * (NMIN-1)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19900) NMIN,(RMCN(MN+L),L=0,2),XP,YP,ZP
         ELSE
            WRITE(IMPRIM,29900) NMIN,(RMCN(MN+L),L=0,2),XP,YP,ZP
         ENDIF
         RMCN( MN     ) = XP
         RMCN( MN + 1 ) = YP
         RMCN( MN + 2 ) = ZP
C
      ENDDO
      CALL ECDATE( MCN(MNXYZS) )
19900 FORMAT(' SOMMET',I10,T19,': XYZ INITIAL=',3G15.6/
     %                     T19,'  XYZ FINAL  =',3G15.6)
29900 FORMAT(' VERTEX',I10,T19,': INITIAL XYZ=',3G15.6/
     %                     T19,'  FINAL   XYZ=',3G15.6)
C
C     COPIE DU TMS NSEF
      CALL NBMONSEF( MCN(MNNSEFI), MONSEF )
      CALL LXTNDC( NTLXVS, 'NSEF', 'MOTS', MONSEF )
      CALL LXTSOU( NTLXVS, 'NSEF', NTNSEFS, MNNSEFS )
      IF( NTNSEFS .LE. 0 ) GOTO 9900
      CALL TRTATA( MCN(MNNSEFI), MCN(MNNSEFS), MONSEF )
      CALL ECDATE( MCN(MNNSEFS) )
      GOTO 9999
C
C     PAS ASSEZ DE MEMOIRE
 9900 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'MEMOIRE MCN SATUREE'
      ELSE
         KERR(1) = 'MAIN MEMORY SATURATED'
      ENDIF
      CALL LEREUR
      IERR = 1
C
 9999 RETURN
      END
