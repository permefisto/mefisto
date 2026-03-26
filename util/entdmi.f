      SUBROUTINE ENTDMI( NBP , PXYZD , DISMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCUL DE LA DISTANCE MINIMALE ENTRE 2 POINTS
C -----     DES NBP POINTS DEFINIS PAR PXYZD
C
C ENTREES :
C ---------
C NBP     : NOMBRE DE POINTS
C PXYZD   : 3 COORDONNEES ET DISTANCE SOUHAITEE DES NBP POINTS
C
C SORTIE  :
C ---------
C DISMIN  : DISTANCE MINMALE ENTRE 2 POINTS DU NUAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      REAL              PXYZD(4,NBP)
C
C     INITIALISATION
      DISMIN = ( PXYZD(1,2) - PXYZD(1,1) ) ** 2 +
     %         ( PXYZD(2,2) - PXYZD(2,1) ) ** 2 +
     %         ( PXYZD(3,2) - PXYZD(3,1) ) ** 2
C
      DO 20 J=1,NBP-1
         X = PXYZD(1,J)
         Y = PXYZD(2,J)
         Z = PXYZD(3,J)
         DO 10 I=J+1,NBP
            D = ( PXYZD(1,I) - X ) ** 2 +
     %          ( PXYZD(2,I) - Y ) ** 2 +
     %          ( PXYZD(3,I) - Z ) ** 2
            IF( D .LT. DISMIN ) THEN
C
               IF( D .LE. 0. ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1)='ENTDMI: 2 POINTS'
                  CALL LEREUR
                  WRITE(IMPRIM,*) '2 POINTS IDENTIQUES ',J,I,
     %                            ' IDENTIQUES',X,Y,Z
               ENDIF
C
               DISMIN = D
            ENDIF
 10      CONTINUE
 20   CONTINUE
C
      DISMIN = SQRT( DISMIN )
      END
