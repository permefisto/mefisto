       SUBROUTINE TAILIDEA( NOFOT, XYZD, NCODEV, TAID )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCULER LA VALEUR>0 DE LA TAILLE IDEALE DES ARETES AUTOUR DE XYZD
C -----  SOIT A PARTIR DE DARETE LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE
C        SOIT A PARTIR DE LA FONCTION DONNEE PAR L'UTILISATEUR
C                            TAILLE_IDEALE(X,Y,Z) ou EDGE_LENGTH(X,Y,Z)

C ENTREES:
C --------
C DARETE : DANS ./incl/darete.inc   (OPTION 0: DU MAILLAGE)
C NOFOT  : INACTIVE le 27/05/2020 AU PROFIT de NOFOTI de ./incl/darete.inc
C          NUMERO DE LA FONCTION UTILISATEUR
C          TAILLE_IDEALE(X,Y,Z) ou EDGE_LENGTH(X,Y,Z)
C XYZD   : XYZ Du POINT de CALCUL de la TAILLE SOUHAITEE DES ARETES ISSUES

C SORTIES:
C --------
C NCODEV : =0 TAID N'EST PAS INITIALISEE EN SORTIE et TAID=0D0
C          =1 TAID   EST     INITIALISEE EN SORTIE
C TAID   : SI NCODEV=1 ALORS la TAILLE_IDEALE>0 au POINT XYZD
C                      SINON DARETE la TAILLE par DEFAUT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET LJLL UPMC PARIS & Saint PIERRE DU PERRAY   Juin 2017
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/darete.inc"
      DOUBLE PRECISION  XYZD(3), TAID

      IF( NOFOTI .GT. 0 ) THEN

C         CALCUL DE TAILLE_IDEALE(X,Y,Z) AU SOMMET XYZD
          CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, TAID )

       ELSE

C         LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE
          TAID   = DARETE
          NCODEV = 1

       ENDIF

       IF( TAID .LE. 0D0 ) THEN

          IF( DARETE .GT. 0D0 ) THEN
C            LA VALEUR PAR DEFAUT DES ARETES DU MAILLAGE
             TAID   = DARETE
             NCODEV = 1
          ELSE
C            IMPOSSIBLE D'AVOIR UNE VALEUR RAISONNABLE
C            CAR DARETE VALEUR PAR DEFAUT EST INCORRECTE
             TAID   = 0D0
             NCODEV = 0
             IF( LANGAG .EQ. 0 ) THEN
             PRINT*,'INITIALISER LA TAILLE D''ARETE PAR DEFAUT DARETE=',
     %               DARETE,' par l''OPTION 0; du MENU Debut'
              ELSE
              PRINT*,'GIVE a VALUE to DEFAULT EDGE_LENGTH by the OPTION 
     %0; of Debut MENU'
              ENDIF
          ENDIF

       ENDIF

       RETURN
       END
