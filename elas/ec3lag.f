      SUBROUTINE EC3LAG( NUELEM,NBELEM,NBPOLY,NUNDEL,NPI,NDSM,NDSMCI,
     %                  POLY,DP,F,
     %                  NOOBVC,NUMIVO,NUMAVO,LTDEVO,
     %                  MNTEMP,NTDLT,TEMPER,NTDLE,DEPLAC,
     %                  ELAS,STRELT,CONINI,
     %                  COOPTC,STRESS)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU TENSEUR DES CONTRAINTES AUX POINTS D'INTEGRATION
C -----    POUR LES ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2 EN 3D
C          SAUF TETR 3P1D
C
C ENTREES:
C --------
C NUELEM : NUMERO DE L'ELEMENT FINI POUR CE TYPE D'ELEMENT
C NBELEM : NOMBRE TOTAL D'ELEMENTS FINIS DE CE TYPE
C NBPOLY : NOMBRE DE POLYNOMES
C NUNDEL : NUMERO DES NBPOLY NOEUDS DES ELEMENTS FINIS
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DE L ELEMENT
C NDSM   : NOMBRE DE CAS DE CHARGE OU SECONDS MEMBRES DU SYSTEME
C NDSMCI : NOMBRE DE JEUX DE CONTRAINTES INITIALES = 0 OU 1 OU NDSM
C POLY   : POLY(I,L) = PI(XL,YL)
C DP     : GRADIENT DES FONCTIONS DE BASE AUX POINTS D INTEGRATION
C F      : X Y Z DES NPI POINTS D INTEGRATION
C
C NOOBVC : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMIVO : NUMERO MINIMAL DES OBJETS SURFACES
C NUMAVO : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES VOLUMES DE L'OBJET
C
C MNTEMP : >0 INDICATEUR DE CONTRAINTES THERMIQUES
C          =0 SI PAS DE  DE CONTRAINTES THERMIQUES
C NTDLT  : NOMBRE TOTAL DE DL EN TEMPERATURE DU MAILLAGE
C TEMPER : TEMPERATURES (NTDLT,NDSM)
C NTDLE  : NOMBRE TOTAL DE DL EN DEPLACEMENTS DU MAILLAGE
C DEPLAC : DEPLACEMENTS (NTDLE,NDSM)
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C ELAS   : TENSEUR DE L ELASTICITE SYMETRIQUE  CALCULE DANS CE SP
C STRELT : TABLEAU AUXILIAIRE (6,3*NBPOLY)
C CONINI : TABLEAU AUXILIAIRE (NDSM,6)
C
C SORTIES:
C --------
C COOPTC : LES COORDONNEES DES POINTS DE CALCUL DES CONTRAINTES
C STRESS : LES CONTRAINTES AUX POINTS ET POUR CHAQUE CAS DE CHARGE
C          STRESS ( 6 , NPI , NDSM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1998
C23456...............................................................012
      include"./incl/donela.inc"
      DOUBLE PRECISION POLY(NBPOLY,NPI),DP(3,NBPOLY,NPI),
     %                 F(NPI,3),STRELT(6,3*NBPOLY),CONINI(NDSM,6),
     %                 TEMPER(NTDLT,NDSM),DEPLAC(NTDLE,NDSM),
     %                 DILATA,STRESS(6,NPI,NDSM)
      REAL             COOPTC(NBELEM,NPI,3)
      INTEGER          NUNDEL(NBELEM,NBPOLY)
      INTEGER          LTDEVO(1:MXDOEL,NUMIVO:NUMAVO)
      DOUBLE PRECISION ELAS(21),ED(6,3),ED3(6),D,TEMPIN
      EQUIVALENCE      (ED,ED3)
      INTRINSIC        REAL
C
C     INITIALISATION A ZERO DES CONTRAINTES
C     -------------------------------------
      CALL AZEROD( 6*NPI*NDSM, STRESS )
C
C     BOUCLE SUR LES POINTS D'INTEGRATION
C     ===================================
      DO 100 L=1,NPI
C
C        LE TABLEAU COOPTC COORDONNEES DES POINTS DE CALCUL
C        --------------------------------------------------
         COOPTC(NUELEM,L,1) = REAL( F(L,1) )
         COOPTC(NUELEM,L,2) = REAL( F(L,2) )
         COOPTC(NUELEM,L,3) = REAL( F(L,3) )
C
C        LE TENSEUR DE L'ELASTICITE EN CE POINT D'INTEGRATION
C        ----------------------------------------------------
         CALL REELAS( 4,NOOBVC,3,F(L,1),F(L,2),F(L,3),
     %                LTDEVO(LPYOUN,NOOBVC),
     %                ELAS )
C
C        CALCUL DES CONTRAINTES ELASTIQUES
C        ---------------------------------
         ED(1,1)=ELAS( 1)
         ED(2,1)=ELAS( 2)
         ED(3,1)=ELAS( 4)
         ED(4,1)=ELAS( 7)
         ED(5,1)=ELAS(11)
         ED(6,1)=ELAS(16)
C
         ED(1,2)=ELAS( 7)
         ED(2,2)=ELAS( 8)
         ED(3,2)=ELAS( 9)
         ED(4,2)=ELAS(10)
         ED(5,2)=ELAS(14)
         ED(6,2)=ELAS(19)
C
         ED(1,3)=ELAS(16)
         ED(2,3)=ELAS(17)
         ED(3,3)=ELAS(18)
         ED(4,3)=ELAS(19)
         ED(5,3)=ELAS(20)
         ED(6,3)=ELAS(21)
C
C        CALCUL DU TABLEAU (ED1)*(DP)
         CALL AB0D(6,3,NBPOLY,ED,DP(1,1,L),STRELT(1,1))
C
         ED(1,1)=ELAS( 7)
         ED(2,1)=ELAS( 8)
         ED(3,1)=ELAS( 9)
         ED(4,1)=ELAS(10)
         ED(5,1)=ELAS(14)
         ED(6,1)=ELAS(19)
C
         ED(1,2)=ELAS( 2)
         ED(2,2)=ELAS( 3)
         ED(3,2)=ELAS( 5)
         ED(4,2)=ELAS( 8)
         ED(5,2)=ELAS(12)
         ED(6,2)=ELAS(17)
C
         ED(1,3)=ELAS(11)
         ED(2,3)=ELAS(12)
         ED(3,3)=ELAS(13)
         ED(4,3)=ELAS(14)
         ED(5,3)=ELAS(15)
         ED(6,3)=ELAS(20)
C
C        CALCUL DU TABLEAU (ED2)*(DP)
         CALL AB0D(6,3,NBPOLY,ED,DP(1,1,L),STRELT(1,1+NBPOLY))
C
         ED(1,1)=ELAS(16)
         ED(2,1)=ELAS(17)
         ED(3,1)=ELAS(18)
         ED(4,1)=ELAS(19)
         ED(5,1)=ELAS(20)
         ED(6,1)=ELAS(21)
C
         ED(1,2)=ELAS(11)
         ED(2,2)=ELAS(12)
         ED(3,2)=ELAS(13)
         ED(4,2)=ELAS(14)
         ED(5,2)=ELAS(15)
         ED(6,2)=ELAS(20)
C
         ED(1,3)=ELAS( 4)
         ED(2,3)=ELAS( 5)
         ED(3,3)=ELAS( 6)
         ED(4,3)=ELAS( 9)
         ED(5,3)=ELAS(13)
         ED(6,3)=ELAS(18)
C
C        CALCUL DU TABLEAU (ED3)*(DP)
         CALL AB0D(6,3,NBPOLY,ED,DP(1,1,L),STRELT(1,1+2*NBPOLY))
C
C        CONTRIBUTION DES DEPLACEMENTS AU TENSEUR DES CONTRAINTES
         DO 9 K=1,6
            ID = 0
            DO 8 I=1,3
               DO 7 J=1,NBPOLY
                  NDL = 3 * NUNDEL( NUELEM , J ) - 3
                  DO 6 M=1,NDSM
                     STRESS(K,L,M) = STRESS(K,L,M)  +
     %                               STRELT(K,ID+J) * DEPLAC(I+NDL,M)
    6             CONTINUE
    7          CONTINUE
               ID = ID + NBPOLY
    8       CONTINUE
    9    CONTINUE
C
C        CALCUL DES CONTRAINTES THERMIQUES
C        ---------------------------------
         IF ( MNTEMP .GT. 0 ) THEN
            IF( LTDEVO(LPDILA,NOOBVC) .GT. 0 ) THEN
C              CALCUL DE -(E)*(DILATA(NDSDE))=(ED3)
C              LE COEFFICIENT DE DILATATION THERMIQUE
               CALL REDILA(4,NOOBVC,F(L,1),F(L,2),F(L,3),
     %                     LTDEVO(LPDILA,NOOBVC),
     %                     DILATA , TEMPIN )
C
C              CALCUL DE ED3 = - (E) (DILATA)
               D      = - DILATA
               ED3(1) = D * ( ELAS( 1) + ELAS( 2) + ELAS( 4) )
               ED3(2) = D * ( ELAS( 2) + ELAS( 3) + ELAS( 5) )
               ED3(3) = D * ( ELAS( 4) + ELAS( 5) + ELAS( 6) )
               ED3(4) = D * ( ELAS( 7) + ELAS( 8) + ELAS( 9) )
               ED3(5) = D * ( ELAS(11) + ELAS(12) + ELAS(13) )
               ED3(6) = D * ( ELAS(16) + ELAS(17) + ELAS(18) )
C
C              CALCUL DE STRELT=(ED3)*(P)
               CALL AB0D(6,1,NBPOLY,ED3,POLY(1,L),STRELT(1,1))
               DO 19 K=1,6
                  DO 18 J=1,NBPOLY
                     DO 17 M=1,NDSM
                        STRESS(K,L,M) = STRESS(K,L,M) +
     %                   STRELT(K,J)*(TEMPER(NUNDEL(NUELEM,J),M)-TEMPIN)
   17                CONTINUE
   18             CONTINUE
   19          CONTINUE
            ENDIF
         ENDIF
C
C        CALCUL DES CONTRAINTES INITIALES
C        --------------------------------
         IF( NDSMCI .GT. 0 .AND. LTDEVO(LPCOIN,NOOBVC) .GT. 0 ) THEN
C           LA VALEUR DES CONTRAINTES INITIALES EN CE POINT
            CALL RECOIN( 4,NOOBVC,6,F(L,1),F(L,2),F(L,3),
     %                   LTDEVO(LPCOIN,NOOBVC), CONINI )
            DO 29 M=1,NDSM
               IF( NDSMCI .EQ. 1 ) THEN
                  MM = 1
               ELSE
                  MM = M
               ENDIF
               DO 27 K=1,6
                  STRESS(K,L,M) = STRESS(K,L,M) + CONINI(MM,K)
  27           CONTINUE
  29        CONTINUE
         ENDIF
 100  CONTINUE
      END
