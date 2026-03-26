      SUBROUTINE POEX01( NTLXPO , RADEFI ,
     %                   NTSOPO , MNSOPO , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE TABLEAU SOMMET DU POINT ( TYPE DE POINT:1 )
C -----
C ENTREES:
C --------
C NTLXPO : NUMERO DU TABLEAU MS DU LEXIQUE DU  POINT
C LADEFI : LADEFI(0:1)    = DATE   DE CREATION DU TABLEAU 'DEFINITION'
C          LADEFI( 2 )    = NOTD   NUMERO DU TABLEAU DESCRIPTEUR
C          LADEFI(WTYTRP) = NTYTRP NUMERO DU TYPE DE LA TRANSFORMATION DU POINT
C          LADEFI(WUTYPO) = NUTYPO NUMERO DU TYPE DU POINT  ICI 1
C          RADEFI(WOORPO  ) = X DU POINT
C          RADEFI(WOORPO+1) = Y DU POINT
C          RADEFI(WOORPO+2) = Z DU POINT
C RADEFI : TABLEAU LADEFI SOUS FORME DE REELS . MEME ADRESSE A L'APPEL
C
C SORTIES:
C --------
C NTSOPO : NUMERO DU TABLEAU 'XYZSOMMET' DU POINT
C MNSOPO : ADRESSE DANS MCN DE CE TABLEAU MS
C          MCN(MNSOPO:MNSOPO+1) = DATE   DE CREATION (REEL DOUBLE PRECISION)
C          MCN(MNSOPO+2)        = NOTD   NUMERO DU TABLEAU DESCRIPTEUR
C          MCN(MNSOPO+WNBSOM)   = NNBSOM NOMBRE DE SOMMETS DU POINT  (=1)
C          MCN(MNSOPO+WYZSOM:MNSOPO+WYZSOM+2)
C                               = CYZSOM TABLEAU DES 3 COORDONNEES XYZ
C IERR   : 0 SI PAS D'ERREUR
C          1 SI NUMERO DE SOMMET INCORRECT
C          2 MCN TROP PETIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C.............................................................................
      IMPLICIT INTEGER (W)
      include"./incl/a_point__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      REAL              RADEFI(0:*)
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
      CALL LXTNDC( NTLXPO , 'XYZSOMMET' , 'MOTS'  , WYZSOM + 3 )
      CALL LXTSOU( NTLXPO , 'XYZSOMMET' ,  NTSOPO , MNSOPO )
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOPO + WNBSOM ) = 1
C
C     X Y Z
      RMCN( MNSOPO + WYZSOM     ) = RADEFI( WOORPO )
      RMCN( MNSOPO + WYZSOM + 1 ) = RADEFI( WOORPO + 1 )
      RMCN( MNSOPO + WYZSOM + 2 ) = RADEFI( WOORPO + 2 )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOPO) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOPO + WNBTGS ) = 0
      MCN( MNSOPO + WBCOOR ) = 3
      MCN( MNSOPO + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
      IERR = 0
      END
