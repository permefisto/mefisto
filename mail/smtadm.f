      SUBROUTINE SMTADM( NS1, NS2, NBSTIS, NOSTIS, PTXYZD,
     %                   NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SAVOIR SI L'ARETE DE SOMMETS NS1 NS2 DANS PTXYZD
C -----    EST ASSEZ ELOIGNEE DES NBSTIS POINTS ISOLES
C
C ENTREES:
C --------
C NS1    : LE NUMERO PTXYZD DU PREMIER SOMMET DE L'ARETE
C NS2    : LE NUMERO PTXYZD DU SECOND  SOMMET DE L'ARETE
C NBSTIS : NOMBRE DE POINTS ISOLES
C NOSTIS : NUMERO PTXYZD DES POINTS ISOLES
C PTXYZD : TABLEAU DES COORDONNEES DES SOMMETS
C
C SORTIE :
C --------
C NONOUI : 1 SI L'ARETE EST ELOIGNEE DES POINTS ISOLES
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1993
C2345X7..............................................................012
      DOUBLE PRECISION  PTXYZD(1:4,1:*),K,D
      INTEGER           NOSTIS(NBSTIS)
C
      DO 100 I=1,NBSTIS
         NS = NOSTIS(I)
         IF( NS .EQ. NS1 .OR. NS .EQ. NS2 ) GOTO 100
C
C        NS N'EST PAS UN POINT ISOLE
C        M LE POINT PROJETE DE NS SUR NS1 NS2
C        CALCUL DE K TEL QUE NS1 M = K NS1 NS2
C        ET                  NS1 M ORTHOGONAL A NS M
         D =  (PTXYZD(1,NS2)-PTXYZD(1,NS1))**2
     %     +  (PTXYZD(2,NS2)-PTXYZD(2,NS1))**2
     %     +  (PTXYZD(3,NS2)-PTXYZD(3,NS1))**2
         K =((PTXYZD(1,NS)-PTXYZD(1,NS1))*(PTXYZD(1,NS2)-PTXYZD(1,NS1))
     %     + (PTXYZD(2,NS)-PTXYZD(2,NS1))*(PTXYZD(2,NS2)-PTXYZD(2,NS1))
     %     + (PTXYZD(3,NS)-PTXYZD(3,NS1))*(PTXYZD(3,NS2)-PTXYZD(3,NS1)))
     %     / D
C
         IF( 0.001D0 .LE. K .AND. K .LE. 1.001D0 ) THEN
C           M LA PROJECTION DE NS EST ENTRE NS1 NS2 OU TRES PROCHE
C           DE NS1 OU DE NS2
C           CALCUL DE LA DISTANCE M NS
            K =  (PTXYZD(1,NS)-PTXYZD(1,NS1))**2
     %        +  (PTXYZD(2,NS)-PTXYZD(2,NS1))**2
     %        +  (PTXYZD(3,NS)-PTXYZD(3,NS1))**2
     %        -  K * K * D
            IF( D .GE. 400D0*K ) THEN
C              NS1 NS2 EST 20 FOIS PLUS GRAND QUE NS M
               NONOUI = 0
               RETURN
            ENDIF
          ENDIF
 100  CONTINUE
C
C     PAS DE POINTS ISOLES TROP PROCHES DE L'ARETE NS1 NS2
      NONOUI = 1
      END
