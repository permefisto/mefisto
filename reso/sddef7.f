      SUBROUTINE SDDEF7( NTLXOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IMPOSER DES CONDITIONS AUX LIMITES DE DIRICHLET
C -----    A UN OBJET AUX LIMITES INCLUS DANS L'INTERFACE
C
C ENTREE :
C --------
C NTLXOB : NUMERO DU LEXIQUE DE L'OBJET AUX LIMITES
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  FEVRIER 1990
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
      COMMON  MCN(MOTMCN)
      REAL   RMCN(1)
      EQUIVALENCE ( MCN(1) , RMCN(1) )
C
C     LE TABLEAU 'XYZSOMMET' DE L'OBJET A RETROUVER
      CALL LXTSOU( NTLXOB , 'XYZSOMMET' , NT , MN )
      IF( NT .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : TABLEAU SOMMET NON RETROUVE'
         CALL LEREUR
         CALL LXIM( NTLXOB )
         RETURN
      ENDIF
      NBSOOB = MCN(MN+WNBSOM)
      CALL DIMCOO( NBSOOB , MCN(MN+WYZSOM) , NDIM )
C
C     DECLARATION DU TABLEAU FIXATION
      CALL LXTSOU( NTLXOB , 'FIXATION' , NT , MN )
      IF( NT .GT. 0 ) THEN
C        LE TABLEAU EXISTANT EST DETRUIT
         CALL LXTSDS( NTLXOB , 'FIXATION' )
      ENDIF
      LTFIXA = 1
      NBCOFI = NDIM
      FIXAT  = 0.
      LOCOFI = WUCOFI + NBCOFI * 2
      CALL LXTNDC( NTLXOB , 'FIXATION' , 'MOTS' , LOCOFI )
C     CREATION DU TABLEAU FIXATION
      CALL LXTSOU( NTLXOB , 'FIXATION' , NT , MN )
      MCN(MN+WTFIXA) = LTFIXA
      MCN(MN+WBCOFI) = NBCOFI
      DO 211 N0 = 0 , NBCOFI - 1
         MCN(MN+WUCOFI+N0)         = N0 + 1
         RMCN(MN+WUCOFI+NBCOFI+N0) = FIXAT
 211  CONTINUE
      CALL ECDATE( MCN(MN) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MOREE2 = MOTVAR(6)
      MCN( MN + MOREE2 ) = NONMTD( '~>>>FIXATION' )
C
      RETURN
      END
