      SUBROUTINE SUEX49( NTLXSU , LADEFI ,
     %                   NTFASU , MNFASU , NTSOFA , MNSOFA , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DESTRUCTURATION D'UNE SURFACE STRUCTUREE EN SURFACE
C -----    NON STRUCTUREE
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF ~td/d/a_surface__definition
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C          CF ~td/d/a___nsef
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC  PARIS        JUILLET 1989
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
C
      IERR   = 0
      NBEFOB = 0
C
C     LA SURFACE A DESTRUCTURER
C     =========================
C     LE NOM DE CETTE SURFACE
      NUSUIN = LADEFI(WUSUIN)
C     LE TABLEAU LEXIQUE DE CETTE SURFACE
      CALL LXNLOU( NTSURF , NUSUIN , NTLXSF , MNLXSU )
      IF( NTLXSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  ' LA SURFACE EST INCONNUE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF , 'NSEF' , NTFAIN , MNFAIN )
      IF( NTFAIN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  ' SURFACE SANS NSEF'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF , 'XYZSOMMET' , NTSOF0 , MNSOF0 )
      IF( NTSOF0 .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  ' SURFACE SANS SOMMETS'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LA SURFACE EST ELLE STRUCTUREE ?
      NUTYSU = MCN( MNFAIN + WUTYMA )
      IF( NUTYSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:6),'(I6)') NUSUIN
         KERR(1) = ' ERREUR SUEX49 ; LA SURFACE' // KERR(MXLGER)(1:6)
     %             // ' N''EST PAS STRUCTUREE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
      IF (NUTYSU.EQ.3) THEN
C        LA SURFACE EST UN TRIANGLE
         NBARTR = MCN ( MNFAIN + WBARTR )
         NBEFOB = NBARTR * NBARTR
      ELSE IF (NUTYSU.EQ.4) THEN
C        LA SURFACE EST UN QUADRILATERE
         NBARXQ = MCN ( MNFAIN + WBARXQ )
         NBARYQ = MCN ( MNFAIN + WBARYQ )
         NBEFOB = NBARXQ * NBARYQ
      ENDIF
C     TEST SUR LE NOMBRE DE SOMMETS
      NBSOFA = MCN ( MNSOF0 + WNBSOM )
      IF( NBSOFA .LE. 1 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)')  NBSOFA
         KERR(1) = ' ERREUR NOMBRE DE SOMMETS' // KERR(MXLGER)(1:4)
     %             // ' AU LIEU DE > 1'
         CALL LEREUR
        IERR = 1
        RETURN
      ENDIF
C
C     CREATION DE LA SURFACE NON STRUCTUREE
C     =====================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     ---------------------------------
      L = WYZSOM + 3 * NBSOFA
      CALL LXTNDC( NTLXSU , 'XYZSOMMET' , 'MOTS' , L )
      CALL LXTSOU( NTLXSU , 'XYZSOMMET' ,  NTSOFA  , MNSOFA )
C     COPIE DU TABLEAU 'XYZSOMMET' STRUCTURE DANS CELUI NON STRUCTURE
      CALL TRTATA( MCN(MNSOF0) , MCN(MNSOFA) , L )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOFA) )
C
C     CONSTRUCTION DU TABLEAU 'NSEF'
C     ------------------------------------
      CALL LXTNDC( NTLXSU , 'NSEF' , 'ENTIER' ,
     S             WUSOEF + NBEFOB * 4 )
      CALL LXTSOU( NTLXSU , 'NSEF' ,  NTFASU  , MNFASU )
C
C     LE TYPE DE L'OBJET : ICI SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     LE TYPE DU MAILLAGE : ICI NON STRUCTURE
      MCN( MNFASU + WUTYMA ) = 0
C     LE NOMBRE DE SOMMETS PAR SOUS-OBJET
      MCN( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE NSEF
      MCN( MNFASU + WBEFOB ) = NBEFOB
C
      IF (NUTYSU.EQ.3) THEN
         CALL SUEXT3(MNFASU,MNFASU+WUSOEF)
      ELSE IF(NUTYSU.EQ.4) THEN
         IA = MNFASU + WUSOEF - 1
         NBX = NBARXQ + 1
         DO NY=1,NBARYQ
            DO NX=1,NBARXQ
               MCN( IA + 1 ) = (NY-1)*NBX+NX
               MCN( IA + 2 ) = (NY-1)*NBX+NX+1
               MCN( IA + 3 ) = NY*NBX+NX+1
               MCN( IA + 4 ) = NY*NBX+NX
               IA            = IA + 4
            ENDDO
         ENDDO
      ENDIF
C
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = 0
C     PAS DE TANGENTES STOCKEES
      MCN( MNFASU + WBTGEF ) = 0
      MCN( MNFASU + WBEFAP ) = 0
      MCN( MNFASU + WBEFTG ) = 0
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNFASU) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
      IERR = 0
C
      RETURN
      END
