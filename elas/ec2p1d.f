      SUBROUTINE EC2P1D( X,      NDSM,   NDSMCI, NUELEM, NBELEM, NUNDEL,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   MNTEMP, NTDLT,  TEMPER,
     %                   NTDLE,  DEPLAC,
     %                   ELAS,   CONINI,
     %                   COOPTC, STRESS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU TENSEUR DES CONTRAINTES AUX POINTS D'INTEGRATION
C -----    POUR LES ELEMENTS FINIS LAGRANGE DE DEGRE 1 OU 2
C
C ENTREES :
C ---------
C X      : 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NDSMCI : NOMBRE DE JEUX DE CONTRAINTES INITIALES = 0 OU 1 OU NDSM
C NUELEM : NUMERO DE L'ELEMENT FINI
C NBELEM : NOMBRE D'ELEMENTS FINIS
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES OBJETS SURFACES
C
C MNTEMP : >0 CALCUL DEMANDE DES CONTRAINTES THERMIQUES,=<0 SINON
C NTDLT  : NOMBRE DE NOEUDS THERMIQUES OU DL THERMIQUES DU MAILLAGE
C TEMPER : TEMPERATURES AUX NOEUDS DU MAILLAGE TEMPER(NTDLT,NDSM)
C NTDLE  : NOMBRE TOTAL DE DL EN DEPLACEMENTS DU MAILLAGE
C DEPLAC : DEPLACEMENTS (NTDLE,NDSM)
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C ELAS   : TENSEUR DE L ELASTICITE SYMETRIQUE  CALCULE DANS CE SP
C CONINI : TABLEAU AUXILIAIRE (NDSM,3)
C
C SORTIES:
C --------
C COOPTC : LES COORDONNEES DU BARYCENTRE DU TRIANGLE POINT DE CALCUL DES CONTRAI
C STRESS : LES CONTRAINTES AU BARYCENTRE ET POUR CHAQUE CAS DE CHARGE
C          STRESS( 3, NDSM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1998
C23456...............................................................012
      include"./incl/donela.inc"
      DOUBLE PRECISION  TEMPER(NTDLT,NDSM),
     %                  DEPLAC(NTDLE,NDSM),
     %                  ELAS(6),
     %                  CONINI(NDSM,3),
     %                  STRESS(3,NDSM)
      REAL              X(3,2)
      REAL              COOPTC(NBELEM,2)
      INTEGER           NUNDEL(NBELEM,3)
      INTEGER           LTDESU(1:MXDOEL,NUMISU:NUMASU)
      DOUBLE PRECISION  DILATA,TEMPIN
C
      DOUBLE PRECISION  D, DELTA, XD, YD
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32
      INTRINSIC         REAL
C
C     INITIALISATION A ZERO DES CONTRAINTES
C     -------------------------------------
      CALL AZEROD( 3*NDSM, STRESS )
C
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)
C
      DELTA = ABS( X21 * Y31 - X31 * Y21 )
C
C     LE NUMERO DES 6 DEGRES DE LIBERTE
      NU21 = NUNDEL(NUELEM,1) * 2
      NU22 = NUNDEL(NUELEM,2) * 2
      NU23 = NUNDEL(NUELEM,3) * 2
      NU11 = NU21 - 1
      NU12 = NU22 - 1
      NU13 = NU23 - 1
C
C     LE TABLEAU COOPTC COORDONNEES DES POINTS DE CALCUL
C     --------------------------------------------------
      XD = (X(1,1)+X(2,1)+X(3,1)) / 3D0
      YD = (X(1,2)+X(2,2)+X(3,2)) / 3D0
      COOPTC(NUELEM,1) = REAL( XD )
      COOPTC(NUELEM,2) = REAL( YD )
C
C     CALCUL DES CONTRAINTES ELASTIQUES
C     ---------------------------------
C     LE TENSEUR DE L'ELASTICITE EN CE BARYCENTRE
      CALL REELAS( 3,NOOBSF,2,XD,YD,0D0,LTDESU(LPYOUN,NOOBSF),
     %             ELAS )
      DO 10 N=1,NDSM
         STRESS(1,N)=((-ELAS(1) * Y32 +ELAS(4) * X32 ) * DEPLAC(NU11,N)
     %               +( ELAS(1) * Y31 -ELAS(4) * X31 ) * DEPLAC(NU12,N)
     %               +(-ELAS(1) * Y21 +ELAS(4) * X21 ) * DEPLAC(NU13,N)
     %               +(-ELAS(4) * Y32 +ELAS(2) * X32 ) * DEPLAC(NU21,N)
     %               +( ELAS(4) * Y31 -ELAS(2) * X31 ) * DEPLAC(NU22,N)
     %               +(-ELAS(4) * Y21 +ELAS(2) * X21 ) * DEPLAC(NU23,N))
     %               / DELTA
         STRESS(2,N)=((-ELAS(2) * Y32 +ELAS(5) * X32 ) * DEPLAC(NU11,N)
     %               +( ELAS(2) * Y31 -ELAS(5) * X31 ) * DEPLAC(NU12,N)
     %               +(-ELAS(2) * Y21 +ELAS(5) * X21 ) * DEPLAC(NU13,N)
     %               +(-ELAS(5) * Y32 +ELAS(3) * X32 ) * DEPLAC(NU21,N)
     %               +( ELAS(5) * Y31 -ELAS(3) * X31 ) * DEPLAC(NU22,N)
     %               +(-ELAS(5) * Y21 +ELAS(3) * X21 ) * DEPLAC(NU23,N))
     %               / DELTA
         STRESS(3,N)=((-ELAS(4) * Y32 +ELAS(6) * X32 ) * DEPLAC(NU11,N)
     %               +( ELAS(4) * Y31 -ELAS(6) * X31 ) * DEPLAC(NU12,N)
     %               +(-ELAS(4) * Y21 +ELAS(6) * X21 ) * DEPLAC(NU13,N)
     %               +(-ELAS(6) * Y32 +ELAS(5) * X32 ) * DEPLAC(NU21,N)
     %               +( ELAS(6) * Y31 -ELAS(5) * X31 ) * DEPLAC(NU22,N)
     %               +(-ELAS(6) * Y21 +ELAS(5) * X21 ) * DEPLAC(NU23,N))
     %               / DELTA
 10   CONTINUE
C
C     CALCUL DES CONTRAINTES THERMIQUES
C     ---------------------------------
      IF( MNTEMP .GT. 0 ) THEN
         IF( LTDESU(LPDILA,NOOBSF) .GT. 0 ) THEN
C           CALCUL DE -(E)*(DILATA(NDSDE))
C           LE COEFFICIENT DE DILATATION THERMIQUE
            CALL REDILA( 3,NOOBSF,XD,YD,0D0,LTDESU(LPDILA,NOOBSF),
     %                   DILATA , TEMPIN )
            DO 20 N=1,NDSM
               D = ( ( TEMPER(NUNDEL(NUELEM,1),N)
     %               + TEMPER(NUNDEL(NUELEM,2),N)
     %               + TEMPER(NUNDEL(NUELEM,2),N) ) / 3D0 - TEMPIN )
     %             * DILATA
               STRESS(1,N) = STRESS(1,N) - (ELAS(1)+ELAS(2)) * D
               STRESS(2,N) = STRESS(2,N) - (ELAS(2)+ELAS(3)) * D
               STRESS(3,N) = STRESS(3,N) - (ELAS(4)+ELAS(5)) * D
 20         CONTINUE
         ENDIF
      ENDIF
C
C     CALCUL DES CONTRAINTES INITIALES
C     --------------------------------
      IF( NDSMCI .GT. 0 .AND. LTDESU(LPCOIN,NOOBSF) .GT. 0 ) THEN
C        LA VALEUR DES CONTRAINTES INITIALES EN CE POINT
         CALL RECOIN( 3,NOOBSF,3,XD,YD,0.D0,
     %                LTDESU(LPCOIN,NOOBSF), CONINI )
         DO 40 N=1,NDSM
            IF( NDSMCI .EQ. 1 ) THEN
               MM = 1
            ELSE
               MM = N
            ENDIF
            DO 30 K=1,3
               STRESS(K,N) = STRESS(K,N) + CONINI(MM,K)
 30         CONTINUE
 40      CONTINUE
      ENDIF
C
      RETURN
      END
