      SUBROUTINE DI2PUT( KNOMP1, KNOMP2, DISTAN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA DISTANCE ENTRE LES POINTS KNOMP1 ET KNOMP2
C -----
C
C ENTREES:
C --------
C KNOMP1 : NOM DU PREMIER POINT
C KNOMP2 : NOM DU SECOND  POINT
C
C SORTIES:
C --------
C DISTAN : DISTANCE ENTRE LES 2 POINTS
C          -1D0 SI L'UN DES 2 POINTS EST INCONNU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     FEVRIER 1994
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
      COMMON            RMCN(MOTMCN)
      REAL              XYZ1(3),XYZ2(3)
      DOUBLE PRECISION  DISTAN
      CHARACTER*(*)     KNOMP1,KNOMP2
C
C     COORDONNEES DU POINT INITIAL
      IF( NTPOIN .LE. 0 ) GOTO 9999
      CALL LXLXOU( NTPOIN , KNOMP1, NTPOI , MN )
      IF( NTPOI .LE. 0 ) GOTO 9999
      CALL LXTSOU( NTPOI  , 'XYZSOMMET' ,  NTSOM , MNSOM )
      IF( NTSOM .LE. 0 ) GOTO 9999
      MN = MNSOM + WYZSOM
      XYZ1( 1 ) = RMCN( MN )
      XYZ1( 2 ) = RMCN( MN + 1 )
      XYZ1( 3 ) = RMCN( MN + 2 )
C
C     COORDONNEES DU POINT FINAL
      CALL LXLXOU( NTPOIN , KNOMP2, NTPOI , MN )
      IF( NTPOI .LE. 0 ) GOTO 9999
      CALL LXTSOU( NTPOI  , 'XYZSOMMET' ,  NTSOM , MNSOM )
      IF( NTSOM .LE. 0 ) GOTO 9999
      MN = MNSOM + WYZSOM
      XYZ2( 1 ) = RMCN( MN )
      XYZ2( 2 ) = RMCN( MN + 1 )
      XYZ2( 3 ) = RMCN( MN + 2 )
C
C     DISTANCE ENTRE CES 2 POINTS
      DISTAN = DIST2P( XYZ1 , XYZ2 )
      RETURN
C
C     ERREUR
 9999 DISTAN = -1D0
      END
