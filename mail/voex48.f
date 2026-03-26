      SUBROUTINE VOEX48( NTLXVS,  LADEFI,
     %                   NTNSEFS, MNNSEFS, NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPOSER DES POINTS COMME SOMMETS D'UN MAILLAGE DE VOLUME
C ----- LE PLUS PROCHE SOMMET DU POINT HERITE DES COORDONNEES DU POINT
C       ATTENTION, LA QUALITE DES ELEMENTS FINIS SERA MODIFIEE
C
C ENTREES:
C --------
C NTLXVS : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME
C LADEFI : TABLEAU ENTIER DE DEFINITION DU VOLUME
C          CF ~/TD/D/A_VOLUME__DEFINITION
C
C SORTIES:
C --------
C NTNSEFS: NUMERO      DU TMS 'NSEF' DU VOLUME
C MNNSEFS: ADRESSE MCN DU TMS 'NSEF' DU VOLUME
C          CF ~/TD/D/A___NSEF
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF ~/TD/D/A___XYZSOMMET
C IERR   : 0 SI PAS D'ERREUR
C          1 SI VOLUME INITIAL SANS NSEF
C          2 SI VOLUME INITIAL SANS XYZSOMMET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC et SAINT PIERRE DU PERRAY  MAI 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_volume__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              DMIN, D
      CHARACTER*24      NOMPT
C
C     LE VOLUME INITIAL
C     =================
C     LE NOM DE CE VOLUME
      NUVOIN = LADEFI( WUVOIN )
C     LE TABLEAU LEXIQUE DE CE VOLUME
      CALL LXNLOU( NTVOLU, NUVOIN, NTLXVI, MN )
      IF( NTLXVI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME INITIAL INCONNU'
         ELSE
            KERR(1) = 'UNKNOWN INITIAL VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTSOU( NTLXVI, 'NSEF', NTNSEFI, MNNSEFI )
      IF( NTNSEFI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS TMS NSEF'
         ELSE
            KERR(1) = 'VOLUME WITHOUT TMS NSEF'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CE VOLUME
      CALL LXTSOU( NTLXVI, 'XYZSOMMET', NTXYZI, MNXYZI )
      IF( NTXYZI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS TMS XYZSOMMET'
         ELSE
            KERR(1) = 'VOLUME WITHOUT TMS XYZSOMMET'
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
      NBPOST = LADEFI( WBPOST )
      DO K = 1, NBPOST
C
C        LE NOM DU POINT DEVENANT SOMMET
         NUPOST = LADEFI( WUPOST - 1 + K )
         CALL NMOBNU( 'POINT', NUPOST, NOMPT )
         NCPT = NUDCNB( NOMPT )
C
C        LE TABLEAU LEXIQUE DE CE POINT
         CALL LXNLOU( NTPOIN, NUPOST, NTPOST, MNPOST )
         IF( NTPOST .LE. 0 ) THEN
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
         CALL LXTSOU( NTPOST, 'XYZSOMMET', NTXYZP, MNXYZP )
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
C        RECHERCHE DU PLUS PROCHE SOMMET DU POINT (XP,YP,ZP)
         MN = MNXYZP + WYZSOM
         XP = RMCN( MN     )
         YP = RMCN( MN + 1 )
         ZP = RMCN( MN + 2 )
         print *
         print *,'POINT ',NOMPT(1:NCPT),' XYZ IMPOSEES=',XP,YP,ZP
C
         DMIN = 1E28
         MN   = MNXYZS + WYZSOM
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
