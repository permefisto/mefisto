      SUBROUTINE EC2LAG( NUELEM,NBELEM,NBPOLY,NUNDEL,NPI,NDSM,NDSMCI,
     %                   POLY,DP,F1,F2,
     %                   NOOBSF,NUMISU,NUMASU,LTDESU,
     %                   MNTEMP,NTDLT,TEMPER,NTDLE,DEPLAC,
     %                   ELAS,STRELT,CONINI,
     %                   COOPTC,STRESS )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU TENSEUR DES CONTRAINTES AUX POINTS D'INTEGRATION
C -----    POUR LES ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2
C
C ENTREES :
C ---------
C NUELEM : NUMERO DE L'ELEMENT FINI POUR CE TYPE D'ELEMENT
C NBELEM : NOMBRE TOTAL D'ELEMENTS FINIS DE CE TYPE
C NBPOLY : NOMBRE DE POLYNOMES
C NUNDEL : NUMERO DES NBPOLY NOEUDS DES ELEMENTS FINIS
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DE L ELEMENT
C NDSM   : NOMBRE DE CAS DE CHARGE OU SECONDS MEMBRES DU SYSTEME
C NDSMCI : NOMBRE DE JEUX DE CONTRAINTES INITIALES = 0 OU 1 OU NDSM
C POLY   : POLY(I,L) = PI(XL,YL)
C DP     : GRADIENT DES FONCTIONS DE BASE AUX POINTS D INTEGRATION
C F1,F2  : X Y DES NPI POINTS D INTEGRATION
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES OBJETS SURFACES
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
C STRELT : TABLEAU AUXILIAIRE (3,2*NBPOLY)
C CONINI : TABLEAU AUXILIAIRE (NDSM,3)
C
C PARAMETRE RESULTAT :
C --------------------
C COOPTC : LES COORDONNEES DES POINTS DE CALCUL DES CONTRAINTES
C STRESS : LES CONTRAINTES AUX POINTS ET POUR CHAQUE CAS DE CHARGE
C          STRESS ( 3 , NPI , NDSM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456...............................................................012
      include"./incl/donela.inc"
      DOUBLE PRECISION POLY(NBPOLY,NPI),DP(2,NBPOLY,NPI),
     %                 F1(NPI),F2(NPI),
     %                 STRELT(3,2*NBPOLY),CONINI(NDSM,3),
     %                 TEMPER(NTDLT,NDSM),DEPLAC(NTDLE,NDSM),
     %                 DILATA,STRESS(3,NPI,NDSM)
      REAL             COOPTC(NBELEM,NPI,2)
      INTEGER          NUNDEL(NBELEM,NBPOLY)
      INTEGER          LTDESU(1:MXDOEL,NUMISU:NUMASU)
      DOUBLE PRECISION ELAS(6),ED1(3,2),ED2(3,2),ED3(3),D,TEMPIN
      INTRINSIC        REAL
C
C     INITIALISATION A ZERO DES CONTRAINTES
C     -------------------------------------
      CALL AZEROD( 3*NPI*NDSM , STRESS )
C
C     BOUCLE SUR LES POINTS D'INTEGRATION
C     ===================================
      DO 100 L=1,NPI
C
C        LE TABLEAU COOPTC COORDONNEES DES POINTS DE CALCUL
C        --------------------------------------------------
         COOPTC(NUELEM,L,1) = REAL( F1(L) )
         COOPTC(NUELEM,L,2) = REAL( F2(L) )
C
C        LE TENSEUR DE L'ELASTICITE EN CE POINT D'INTEGRATION
C        ----------------------------------------------------
         CALL REELAS( 3,NOOBSF,2,F1(L),F2(L),0D0,
     %                LTDESU(LPYOUN,NOOBSF),
     %                ELAS )
C
C        CALCUL DES CONTRAINTES ELASTIQUES
C        ---------------------------------
         ED1(1,1)=ELAS(1)
         ED1(1,2)=ELAS(4)
         ED1(2,1)=ELAS(2)
         ED1(2,2)=ELAS(5)
         ED1(3,1)=ELAS(4)
         ED1(3,2)=ELAS(6)
C
         ED2(1,1)=ELAS(4)
         ED2(1,2)=ELAS(2)
         ED2(2,1)=ELAS(5)
         ED2(2,2)=ELAS(3)
         ED2(3,1)=ELAS(6)
         ED2(3,2)=ELAS(5)
C
C        CALCUL DU TABLEAU (ED1)*(DP)
         CALL AB0D(3,2,NBPOLY,ED1,DP(1,1,L),STRELT(1,1))
C
C        CALCUL DU TABLEAU (ED2)*(DP)
         CALL AB0D(3,2,NBPOLY,ED2,DP(1,1,L),STRELT(1,NBPOLY+1))
C
C        CONTRIBUTION DES DEPLACEMENTS AU TENSEUR DES CONTRAINTES
         DO 9 K=1,3
            ID = 0
            DO 8 I=1,2
               DO 7 J=1,NBPOLY
                  NDL = 2 * NUNDEL( NUELEM , J ) - 2
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
            IF( LTDESU(LPDILA,NOOBSF) .GT. 0 ) THEN
C              CALCUL DE -(E)*(DILATA(NDSDE))=(ED3)
C              LE COEFFICIENT DE DILATATION THERMIQUE
               CALL REDILA(3,NOOBSF,F1(L),F2(L),0D0,
     %                     LTDESU(LPDILA,NOOBSF),
     %                     DILATA , TEMPIN )
C              CALCUL DE ED1 = - (E) (DILATA)
               D      = - DILATA
               ED3(1) = D * ( ELAS(1)+ELAS(2)+ELAS(4) )
               ED3(2) = D * ( ELAS(2)+ELAS(3)+ELAS(5) )
               ED3(3) = D * ( ELAS(4)+ELAS(5)+ELAS(6) )
C              CALCUL DE STRELT=(ED3)*(P)
               CALL AB0D(3,1,NBPOLY,ED3,POLY(1,L),STRELT(1,1))
               DO 19 K=1,3
                  DO 18 J=1,NBPOLY
                     DO 17 M=1,NDSM
                        STRESS(K,L,M) = STRESS(K,L,M) +
     %                   STRELT(K,J)*(TEMPER(NUNDEL(NUELEM,J),M)-TEMPIN)
  17                 CONTINUE
  18              CONTINUE
  19           CONTINUE
            ENDIF
         ENDIF
C
C        CALCUL DES CONTRAINTES INITIALES
C        --------------------------------
         IF( NDSMCI .GT. 0 .AND. LTDESU(LPCOIN,NOOBSF) .GT. 0 ) THEN
C           LA VALEUR DES CONTRAINTES INITIALES EN CE POINT
            CALL RECOIN( 3,NOOBSF,3,F1(L),F2(L),0.D0,
     %                   LTDESU(LPCOIN,NOOBSF), CONINI )
            DO 29 M=1,NDSM
               IF( NDSMCI .EQ. 1 ) THEN
                  MM = 1
               ELSE
                  MM = M
               ENDIF
               DO 27 K=1,3
                  STRESS(K,L,M) = STRESS(K,L,M) + CONINI(MM,K)
  27           CONTINUE
  29        CONTINUE
         ENDIF
 100  CONTINUE
      END
