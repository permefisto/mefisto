      SUBROUTINE REVITANG( NBNOEF, XYZEF,
     %                     NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                     Omega )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECUPERER LA VITESSE ANGULAIRE DE ROTATION DU FLUIDE
C -----    A L'INSTANT TEMPS STOCKE DANS LA VARIABLE TEMPS
C          DE ./incl/ctemps.inc
C
C ENTREES:
C --------
C NBNOEF : NOMBRE DE SOMMETS D'UN EF
C XYZEF  : 3 COORDONNEES DES NBNOEF SOMMETS DE L'EF
C NOOBVC : NUMERO DE VOLUME DU TETRAEDRE
C NUMIVO : NUMERO MINIMAL DES VOLUMES DE L'OBJET
C NUMAVO : NUMERO MAXIMAL DES VOLUMES DE L'OBJET
C LTDEVO : ADRESSE MCN DES TABLEAUX DES DONNEES DU FLUIDE
C
C SORTIES:
C --------
C Omega  : VITESSE ANGULAIRE DE LA ROTATION OMEGA DE L'OBJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Mai 2011
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
      INTEGER           LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
      REAL              XYZEF(NBNOEF,3)
      DOUBLE PRECISION  Omega(3), XYZD(1:3), S
C
C     BARYCENTRE DE L'EF
      DO K=1,3
         S = 0D0
         DO N=1,NBNOEF
            S = S + XYZEF(N,K)
         ENDDO
         XYZD(K)= S / NBNOEF
      ENDDO
C
C     VECTEUR(3) DE VITESSE ANGULAIRE AU BARYCENTRE DE L'EF
      MNVIAN = LTDEVO(LPVIAN,NOOBVC)
      IF( MNVIAN .GT. 0 ) THEN
         CALL REVIAN( 4, NOOBVC, XYZD(1), XYZD(2), XYZD(3),
     %                MNVIAN, Omega )
      ELSE
         Omega(1) = 0D0
         Omega(2) = 0D0
         Omega(3) = 0D0
      ENDIF
C
      RETURN
      END
