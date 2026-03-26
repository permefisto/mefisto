       SUBROUTINE TAILIDEAR( NOFOT, XYZR, NCODEV, TAID )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER LA VALEUR>0 DE LA TAILLE IDEALE DES ARETES AUTOUR DE XYZR
C -----  SOIT A PARTIR DE DARETE LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE
C        SOIT A PARTIR DE LA FONCTION DONNEE PAR L'UTILISATEUR
C                            TAILLE_IDEALE(X,Y,Z) ou EDGE_LENGTH(X,Y,Z)

C ENTREES:
C --------
C DARETE : DANS ./incl/darete.inc   (OPTION 0: DU MAILLAGE)
C NOFOT  : INACTIVE le 27/05/2020 AU PROFIT de NOFOTI de ./incl/darete.inc
C          NUMERO DE LA FONCTION UTILISATEUR
C          TAILLE_IDEALE(X,Y,Z) ou EDGE_LENGTH(X,Y,Z)
C XYZR   : XYZ Du POINT de CALCUL de la TAILLE SOUHAITEE DES ARETES ISSUES

C SORTIES:
C --------
C NCODEV : 0 TAID N'EST PAS INITIALISEE EN SORTIE et TAID=0D0
C          1 TAID   EST     INITIALISEE EN SORTIE
C TAID   : SI NCODEV=1 ALORS la TAILLE_IDEALE>0 au POINT XYZR
C                      SINON 0D0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET LJLL UPMC PARIS & Saint PIERRE DU PERRAY   Juin 2017
C2345X7..............................................................012
       include"./incl/darete.inc"
       REAL              XYZR(3)
       DOUBLE PRECISION  XYZD(3), TAID
       INTRINSIC         DBLE

C      PASSAGE DE XYZR REAL en DOUBLE PRECISION XYZD
       XYZD( 1 ) = DBLE( XYZR( 1 ) )
       XYZD( 2 ) = DBLE( XYZR( 2 ) )
       XYZD( 3 ) = DBLE( XYZR( 3 ) )

C      APPEL DE LA VERSION DOUBLE PRECISION
       CALL TAILIDEA( NOFOTI, XYZD,  NCODEV, TAID )

       RETURN
       END
