        SUBROUTINE SUEX12( NTLXSU , LADEFI ,
     %                     NTFASU , MNFASU , NTSOFA , MNSOFA , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN TRIANGLE AVEC NBINTR
C -----    ARETES PAR COTE
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU QUADRANGLE
C LADEFI : TABLEAU DE DEFINITION DE LA SURFACE PARTITIONNEE
C          CF '~TD/DAT/A_SURFACES__DEFINITION'
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF '~TD/DAT/A___NSEF'
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~TD/DAT/A___SOMMETS'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS   MARS 1989
C.......................................................................
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
C
      IERR = 0
C
C     NOMBRE D'ARETES
C     ===============
      NBINTR = LADEFI(WBINTR)
      IF( NBINTR .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)')NBINTR
         KERR(1) =  ' NOMBRE D''ARETES PAR COTE INCORRECT ='
     %              // KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS
C     ======================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      NB    =  NBINTR + 1
      NBSOM = ( NB * NB + NB ) / 2
      CALL LXTNDC( NTLXSU , 'XYZSOMMET' , 'MOTS' ,  WYZSOM + 3 * NBSOM )
      CALL LXTSOU( NTLXSU , 'XYZSOMMET' ,  NTSOFA , MNSOFA )
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM) = NBSOM
C     CALCUL DES COORDONNEES DES SOMMETS DU TRIANGLE
      MN = MNSOFA + WYZSOM
      DO 20 I=1,NB
         DO 10 J=1,I
            RMCN( MN     ) = I - J
            RMCN( MN + 1 ) = J - 1
            RMCN( MN + 2 ) = 0.
            MN  = MN + 3
 10      CONTINUE
 20   CONTINUE
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + WNBTGS ) = 0
      MCN( MNSOFA + WBCOOR ) = 3
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     GENERATION DES NSEF
C     ==========================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'NSEF' , 'ENTIER', WBARTR+1 )
      CALL LXTSOU( NTLXSU ,'NSEF' ,  NTFASU , MNFASU )
C     MISE A JOUR DU TABLEAU 'NSEF' DE CETTE SURFACE
C     TYPE DE L'OBJET : SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = 0
C     NUMERO DU TYPE DU MAILLAGE : TRIANGLE STRUCTURE
      MCN( MNFASU + WUTYMA ) = 3
C     NBARTR NOMBRE D'ARETES PAR COTE
      MCN( MNFASU + WBARTR ) = NBINTR
C     LE NOMBRE DE SOMMETS PAR FACE
      MCN ( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C0
      MCN ( MNFASU + WBTGEF ) = 0
      MCN ( MNFASU + WBEFAP ) = 0
      MCN ( MNFASU + WBEFTG ) = 0
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNFASU + WBEFOB ) = NBINTR * NBINTR
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      END
