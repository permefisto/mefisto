      PROGRAM PPPOBA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C BUT :  GENERATION DU FICHIER NFPOBA DES VALEURS
C -----  DES POIDS ET COORDONNEES DES POINTS DES FORMULES
C        D INTEGRATION NUMERIQUE SUR LES ELEMENTS UNITE
C        DES VALEURS DES POLYNOMES DE BASE
C        ET DE LEURS DERIVEES AUX POINTS D INTEGRATION NUMERIQUE
C        DES ELEMENTS FINIS DE TABLEAUX ENCOMBRANTS
C        *************************************************************
C        *  CE PROGRAMME DOIT ETRE EXECUTE AVANT TOUT APPEL A UN TEL *
C        *  ELEMENT FINI SOUS PEINE D ERREUR.......                  *
C        *************************************************************
C
C        GENERATION of Omegal, P(bl), DP(bl) of every finite element
C        for quadrature formula on the reference finite element
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1989
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1990
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C AJOUTS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2005
C AJOUTS : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris    JUIN 2007
C AJOUTS : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris OCTOBRE 2007
C AJOUTS : ALAIN PERRONNET Saint Pierre du Perray & LJLL   DECEMBRE 2009
C AJOUTS : ALAIN PERRONNET Saint Pierre du Perray & LJLL   DECEMBRE 2012
C23456---------------------------------------------------------------012
      include"./incl/homdir.inc"
      include"./incl/motmcg.inc"
      include"./incl/pppoba.inc"
      include"./incl/msvaau.inc"
      include"./incl/ppmck.inc"
      COMMON /        / MCN( MOTMCN )
      REAL              RMCN(MOTMCN )
      DOUBLE PRECISION  DMCN(MOTMCN/2)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
C
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
C
      CHARACTER*160     KNOM
      CHARACTER*4       KPPP
      INTEGER          IKPPP
      EQUIVALENCE      (KPPP,IKPPP)
      DOUBLE PRECISION  PP, V, TABLO(10000)
C
C     LES POIDS ET COORDONNEES DES FORMULES D INTEGRATION
C     ===================================================
C     PO COMME POINT
C     SE COMME SEGMENT
C     TR COMME TRIANGLE
C     QU COMME QUADRANGLE
C     TE COMME TETRAEDRE
C     PE COMME PENTAEDRE
C     HE COMME HEXAEDRE
C     PY COMME PYRAMIDE
C     P3 COMME POLYNOME DE DEGRE 3 EN (X,Y,Z) OU (X,Y)
C     R3 COMME POLYNOME DE DEGRE 3 EN (X,Y),Z
C     Q3 COMME POLYNOME DE DEGRE 3 EN X ET Y ET Z ...
      DOUBLE PRECISION RAC,D,DD,S
      DOUBLE PRECISION POSEP1( 2), COSEP1( 2 )
      DOUBLE PRECISION POSEP3( 2), COSEP3( 2 )
      DOUBLE PRECISION POSEP5( 3), COSEP5( 3 )
      DOUBLE PRECISION POSEP7( 4), COSEP7( 4 )
      DOUBLE PRECISION POSEP9( 5), COSEP9( 5 )

      DOUBLE PRECISION POTRP1( 3), COTRP1( 2,  3 )
      DOUBLE PRECISION POTRP2( 3), COTRP2( 2,  3 )
      DOUBLE PRECISION POTRP5( 7), COTRP5( 2,  7 )

      DOUBLE PRECISION POQUQ1( 4), COQUQ1( 2,  4 )
      DOUBLE PRECISION POQUQ3( 4), COQUQ3( 2,  4 )
      DOUBLE PRECISION POQUQ5( 9), COQUQ5( 2,  9 )
      DOUBLE PRECISION POQUQ7(16), COQUQ7( 2, 16 )
      DOUBLE PRECISION POQUQ9(25), COQUQ9( 2, 25 )

      DOUBLE PRECISION POTEP1( 4), COTEP1( 3,  4 )
      DOUBLE PRECISION POTEP5(15), COTEP5( 3, 15 )

      DOUBLE PRECISION POPER1( 6), COPER1( 3,  6 )
      DOUBLE PRECISION POPER2( 6), COPER2( 3,  6 )
      DOUBLE PRECISION POPER5(21), COPER5( 3, 21 )

      DOUBLE PRECISION POHEQ1(  8), COHEQ1( 3,   8 )
      DOUBLE PRECISION POHEQ3(  8), COHEQ3( 3,   8 )
      DOUBLE PRECISION POHEQ5( 27), COHEQ5( 3,  27 )
      DOUBLE PRECISION POHEQ7( 64), COHEQ7( 3,  64 )
      DOUBLE PRECISION POHEQ9(125), COHEQ9( 3, 125 )

      DOUBLE PRECISION POPYQ1( 5), COPYQ1( 3,  5 )
      DOUBLE PRECISION POPYQ3( 8), COPYQ3( 3,  8 )
      DOUBLE PRECISION POPYQ5(27), COPYQ5( 3, 27 )

      DOUBLE PRECISION PO6CQ1(64), CO6CQ1( 6, 64 )
C
C     LES COEFFICIENTS DES POLYNOMES DE BASE DES ELEMENTS
C     ===================================================
C     POLYNOME SEP1 LAGRANGE P1 SUR LE SEGMENT UNITE
      DOUBLE PRECISION SEP1(2,2)
C
C     POLYNOME SEP2 LAGRANGE P2 SUR LE SEGMENT UNITE
      DOUBLE PRECISION SEP2(3,3)
C
C     POLYNOME TRP1 LAGRANGE P1 SUR LE TRIANGLE UNITE
      DOUBLE PRECISION TRP1(2,2,3)
C
C     POLYNOME TRP2 LAGRANGE P2 SUR LE TRIANGLE UNITE
      DOUBLE PRECISION TRP2(3,3,6)
C
C     POLYNOME TRP1Bulle BREZZI-FORTIN SUR LE TRIANGLE UNITE
      DOUBLE PRECISION TRP1Bulle(4,4,4)
C
C     POLYNOME QUQ1 LAGRANGE Q1 SUR LE QUADRANGLE UNITE
      DOUBLE PRECISION QUQ1(2,2,4)
C
C     POLYNOME QUQ2 LAGRANGE Q2 SUR LE QUADRANGLE UNITE
      DOUBLE PRECISION QUQ2(3,3,8)
C
C     POLYNOME TEP1 LAGRANGE P1 SUR LE TETRAEDRE UNITE
      DOUBLE PRECISION TEP1(2,2,2,4)
C
C     POLYNOME TEP2 LAGRANGE P2 SUR LE TETRAEDRE UNITE
      DOUBLE PRECISION TEP2(3,3,3,10)
C
C     POLYNOME TEP1Bulle BREZZI-FORTIN SUR LE TETRAEDRE UNITE
      DOUBLE PRECISION TEP1Bulle(5,5,5, 5)
C
C     TABLEAU AUXILIAIRE (LES DERIVEES PREMIERES ET SECONDES D UN POLYNOME)
      DOUBLE PRECISION DPXYZ(6*125), DDPXYZ(21*125)
C
cccC     TABLEAU DES INTEGRALES DES DERIVEES DES POLYNOMES SUR TETRAEDRE P2
ccc      DOUBLE PRECISION DPDP(3,10,3,10)
C
C     TABLEAUX DES VALEURS DES POLYNOMES (POLYPI) DES DERIVEES (DPOLYP)
C     AUX POINTS D INTEGRATION NUMERIQUE (MAXIMUM POUR 6-CUBE 6Q1C)
C     =================================================================
      PARAMETER        (MXDIM=6,  MXDIM2=MXDIM*(MXDIM+1)/2,
     %                  MXPOL=64, MXPTI=125)
      DOUBLE PRECISION  POLYPI(       MXPOL, MXPTI),
     %                  DPOLYP(MXDIM, MXPOL, MXPTI),
     %                  DDPOLY(MXDIM2,MXPOL, MXPTI)
C     SI IPOIDS = 1     POIDS(NPI) = POIDS DES POINTS
C     SI IPOLY  = 1     POLY(NBPOLY,NPI) =   PI EN (XL,YL)
C                       POLY(  I   , L ) =
C     SI IDPOLY = 1     DPOLYP(NDIM,NBPOLY,NPI)
C                       DPOLYP( I  , J    , L ) =DPJ/DXI (XL,YL)
C     SI IDDPOL = 1     DDPOLY(NDIM2,NBPOLY,NPI)=DDPJ/DXXI (XL,YL)
C       pour NDIM =3    DDPOLY( 1, J, L )=DDPJ/DDXX (XL,YL)
C            NDIM2=6    DDPOLY( 2, J, L )=DDPJ/DDYX (XL,YL)
C                       DDPOLY( 3, J, L )=DDPJ/DDYY (XL,YL)
C                       DDPOLY( 4, J, L )=DDPJ/DDZX (XL,YL,ZL)
C                       DDPOLY( 5, J, L )=DDPJ/DDZY (XL,YL,ZL)
C                       DDPOLY( 6, J, L )=DDPJ/DDZZ (XL,YL,ZL) ...
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      INTEGER           ITA(4,127)
C
      COMMON / TRAVA1 / NTITRE(20),NDATE(2),NOMCRE(6)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON /LIBINT/ MPFIC(16),MPLIG(16),IPILE,ILIGNE,ISOUFF,IRETAP
C
      COMMON / ELFINI / NELFI(512)
      EQUIVALENCE(NELFI(1),NBTAST),(NELFI(2),NODEPA),(NELFI(3),ITA(1,1))

      double precision psebl(3,3),p
C
      DATA    SEP1 / 1.D0 ,  -1.D0 ,
     %               0.D0 ,   1.D0 /
C
      DATA    SEP2 / 1.D0 ,  -3.D0 ,  2D0 ,
     %               0.D0 ,  -1.D0 ,  2D0 ,
     %               0.D0 ,   4.D0 , -4D0 /
C
      DATA    TRP1 /
     %        1D0 , -1D0 , -1D0 ,  0D0 ,
     %        0D0 ,  1D0 ,  0D0 ,  0D0 ,
     %        0D0 ,  0D0 ,  1D0 ,  0D0 /
C
      DATA    TRP2 /
     %        1D0 , -3D0 , 2D0 , -3D0 , 4D0 , 0D0 , 2D0 , 0D0 , 0D0 ,
     %        0D0 , -1D0 , 2D0 ,  0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 ,
     %        0D0 ,  0D0 , 0D0 ,-1.D0 , 0D0 , 0D0 , 2D0 , 0D0 , 0D0 ,
     %        0D0 ,  4D0 ,-4D0 ,  0D0 ,-4D0 , 0D0 , 0D0 , 0D0 , 0D0 ,
     %        0D0 ,  0D0 , 0D0 ,  0D0 , 4D0 , 0D0 , 0D0 , 0D0 , 0D0 ,
     %        0D0 ,  0D0 , 0D0 ,  4D0 ,-4D0 , 0D0 ,-4D0 , 0D0 , 0D0/
C
C     LES 4 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TRIANGLE REFERENCE
      DATA    TRP1Bulle /
     %        1D0,-1D0,0D0,0D0,  -1D0,-9D0, 9D0,0D0,  0D0,9D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,
     %        0D0, 1D0,0D0,0D0,   0D0,-9D0, 9D0,0D0,  0D0,9D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,   1D0,-9D0, 9D0,0D0,  0D0,9D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,   0D0,27D0,-27D0,0D0, 0D0,-27D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0 /
C
      DATA    QUQ1 /
     %        1D0 , -1D0 , -1D0 ,  1D0 ,
     %        0D0 ,  1D0 ,  0D0 , -1D0 ,
     %        0D0 ,  0D0 ,  0D0 ,  1D0 ,
     %        0D0 ,  0D0 ,  1D0 , -1D0 /
C
      DATA    QUQ2 /
     %      1D0 , -3D0 ,  2D0 , -3D0 ,  5D0 , -2D0 ,  2D0 , -2D0 , 0D0 ,
     %      0D0 , -1D0 ,  2D0 ,  0D0 , -1D0 , -2D0 ,  0D0 ,  2D0 , 0D0 ,
     %      0D0 ,  0D0 ,  0D0 ,  0D0 , -3D0 ,  2D0 ,  0D0 ,  2D0 , 0D0 ,
     %      0D0 ,  0D0 ,  0D0 , -1D0 , -1D0 ,  2D0 ,  2D0 , -2D0 , 0D0 ,
     %      0D0 ,  4D0 , -4D0 ,  0D0 , -4D0 ,  4D0 ,  0D0 ,  0D0 , 0D0 ,
     %      0D0 ,  0D0 ,  0D0 ,  0D0 ,  4D0 ,  0D0 ,  0D0 , -4D0 , 0D0 ,
     %      0D0 ,  0D0 ,  0D0 ,  0D0 ,  4D0 , -4D0 ,  0D0 ,  0D0 , 0D0 ,
     %      0D0 ,  0D0 ,  0D0 ,  4D0 , -4D0 ,  0D0 , -4D0 ,  4D0 , 0D0 /
C     QUADRANGLE Serendipity Valeur imposee au barycentre
C     Cf Cours Thomas Raviart 1972  page X-22
C     P(Barycentre) = - SOM 1 a 4 P(Sommet) / 4 + Som 5 a 8 P(Milieu) / 2
C
      DATA    TEP1 /
     %        1D0, -1D0, -1D0,  0D0,  -1D0, 0D0, 0D0, 0D0,
     %        0D0,  1D0,  0D0,  0D0,   0D0, 0D0, 0D0, 0D0,
     %        0D0,  0D0,  1D0,  0D0,   0D0, 0D0, 0D0, 0D0,
     %        0D0,  0D0,  0D0,  0D0,   1D0, 0D0, 0D0, 0D0 /
C
      DATA    TEP2 /
     %        1D0, -3D0, 2D0,  -3D0, 4D0, 0D0,  2D0, 0D0, 0D0,
     %       -3D0,  4D0, 0D0,   4D0, 0D0, 0D0,  0D0, 0D0, 0D0, 
     %        2D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0, -1D0, 2D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,  -1D0, 0D0, 0D0,  2D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %       -1D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        2D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  4D0,-4D0,   0D0,-4D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0, -4D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,   0D0, 4D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,   4D0,-4D0, 0D0, -4D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,  -4D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        4D0, -4D0, 0D0,  -4D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %       -4D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  4D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,

     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   4D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     %        0D0,  0D0, 0D0,   0D0, 0D0, 0D0,  0D0, 0D0, 0D0 /
C
C     LES 5 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TETRAEDRE REFERENCE
      DATA TEP1Bulle/
     %        1D0,   -1D0,  3*0D0,   -1D0, 19*0D0,
     %       -1D0,  5*0D0,  -64D0,   64D0,  3*0D0,  64D0,  13*0D0,
     %      6*0D0,   64D0, 68*0D0,

     %        0D0,    1D0, 23*0D0,
     %      6*0D0,  -64D0,   64D0,  3*0D0,   64D0,  13*0D0,
     %      6*0D0,   64D0, 68*0D0,

     %      5*0D0,    1D0,  19*0D0,
     %      6*0D0,  -64D0,    64D0, 3*0D0,   64D0,  13*0D0,
     %      6*0D0,   64D0,  68*0D0,

     %     25*0D0,
     %        1D0,  5*0D0,  -64D0,   64D0,  3*0D0,  64D0,  13*0D0,
     %      6*0D0,   64D0, 68*0D0,

     %     31*0D0,  256D0,  -256D0, 3*0D0, -256D0, 13*0D0,
     %      6*0D0, -256D0,  68*0D0 /
C
 9999 FORMAT(/' LE CONTENU DU COMMON / ELFINI /'/1X,130(1H=)/
     %' NBTAST=',I3,' NODEPA=',I5//5(1X,A4,I7,2I5,4X))
   75 FORMAT(2D25.17)
10100 FORMAT(/' TABLEAU ',A4,' ADRESSE DANS M ',I8,3X,' NBRE DE MOTS ',
     %I7,3X,' NO DE SA 1-ERE PAGE ',I7/1X,130(1H=)/(10G13.5))
10150 FORMAT(/' TABLEAU ',A4,' ADRESSE DANS M ',I9,3X,' NBRE DE MOTS ',
     %I7,3X,' NO DE SA 1-ERE PAGE ',I6/1X,130(1H=)/
     %' NDIM =',I5,' DEGRE + 1 =',I5,' P(0) OU Q(1) OU R(2) =',I5,
     %' NBRE POLYNOMES =',I5/(10G13.5))
19997 FORMAT(' MAIN : FIN NORMALE DE POBA'/)
C
C     ==================================================================
C     CALCUL D'INTEGRALES DE POLYNOMES ET DERIVEES et AFFICHAGE SIMPLE
C     RECUPERATION POUR FICHIER INCLUDE DANS LES SUBROUTINES DE CALCUL
C     ==================================================================
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de DP1Bulle dx dy
      CALL PN2DPN( 3, 4, TRP1Bulle, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P1 DP1Bulle dx dy
      CALL PN2PMDPN( 1, 3, TRP1, 3, 4, TRP1Bulle, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de DP2 dx dy
      CALL PN2DPN( 2, 6, TRP2, TABLO )
C
C     TABLEAU DDP2/dxi dxj sur le TRIANGLE UNITE
      CALL DDP2TR( TRP2, TABLO )
C
C     TABLEAU Integrale P1i P2j P1k sur le TRIANGLE UNITE
      CALL P1P2P1TR( TRP1, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P2 P1 dx dy
      CALL PN2PMPN( 2, 6, TRP2, 1, 3, TRP1, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P2 P2 dx dy
      CALL PN2PMPN( 2, 6, TRP2, 2, 6, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P2 DP2 dx dy
      CALL PN2PMDPN( 2, 6, TRP2, 2, 6, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de DP2 DP2 dx dy
      CALL PN2DPMDPM( 2, 6, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE de P2 P2 DP2 dx dy
      CALL PN2PNPNDPN( 2, 6, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P1 DP2 dx dy
      CALL PN2PMDPN( 1, 3, TRP1, 2, 6, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P1 DP2 DP2 dx dy
      CALL PN2PMDPNDPN( 1, 3, TRP1, 2, 6, TRP2, TABLO )
C
C     TABLEAU Integrale sur le TRIANGLE UNITE de P1 P2 DDP2 dx dy
      CALL PN2PMPNDDPN( 1, 3, TRP1, 2, 6, TRP2, TABLO )
C
C     -----------------------------------------------------------------
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de DP1Bulle dx dy dz
      CALL PN3DPN( 4, 5, TEP1Bulle, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de P1 DP1Bulle dx dy dz
      CALL PN3PMDPN( 1, 4, TEP1, 4, 5, TEP1Bulle, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de DP2 dx dy dz
      CALL PN3DPN( 2, 10, TEP2, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de P2 P1 dx dy dz
      CALL PN3PMPN( 2, 10, TEP2, 1, 4, TEP1, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de P2 P2 dx dy dz
      CALL PN3PMPN( 2, 10, TEP2, 2, 10, TEP2, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de P2 DP2 dx dy dz
      CALL PN3PMDPN( 2, 10, TEP2, 2, 10, TEP2, TABLO )

C     TABLEAU Integrale sur le TETRAEDRE UNITE de DP2 DP2 dx dy dz
      CALL PN3DPMDPM( 2, 10, TEP2, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE UNITE de P1 DP2 dx dy dz
      CALL PN3PMDPN( 1, 4, TEP1, 2, 10, TEP2, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE de P2 P2 DP2 dx dy
      CALL PN3PNPNDPN( 2, 10, TEP2, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE de P1 DP2 DP2 dx dy
      CALL PN3PMDPNDPN( 1, 4, TEP1, 2, 10, TEP2, TABLO )
C
C     TABLEAU Integrale sur le TETRAEDRE de P1 P2 DDP2 dx dy
      CALL PN3PMPNDDPN( 1, 4, TEP1, 2, 10, TEP2, TABLO )
C
C
C     ==================================================================
C     INITIALISATION DES TABLEAUX ET DU FICHIER POBA
C     ==================================================================
      KPPP = '....'
C
C     LE FORMAT LIBRE , LES NUMEROS DES UNITES
      CALL INITIA( .FALSE. )
C
C     INITIALISATION DE VARIABLES DE LA MS
      CALL MSINIT
C
C     INITIALISATION DES TABLES DES TABLEAUX MC
C
C     ATTENTION : MOTSMG ET MOTSMK DOIVENT CORRESPONDRE A LA DECLARATION
C     =========== DU NOMBRE DE MOTS DES TABLEAUX MCG ET MCK DES
C                 COMMON / MSMCG / ET / MSMCK /
      MOTSMG = MOTMCG
      MOTSMN = MOTMCN
      MOTSMK = MOTMCK
C
C     INITIALISATION DU REPERTOIRE DES ZONES LIBRES DANS MCG
C     ======================================================
C     LE NOMBRE MAXIMAL DE ZONES LIBRES PERMISES ENTRE LES TABLEAUX
C     STOCKES DANS MCG
      M2ZLMG = 16
      M2ZLMK = 8
      M2ZLMN = 128
C
C     LE NOMBRE DE MOTS DU REPERTOIRE REZLMN REZLMG ET REZLMK
      MOZLMG = 3 * M2ZLMG + 3
      MOZLMK = 3 * M2ZLMK + 3
      MOZLMN = 3 * M2ZLMN + 3
C     Y A-T-IL ASSEZ DE PLACE EN MCG ?
      I = MOZLMN + MOZLMG + MOZLMK
      IF( I .GT. MOTSMG ) THEN
C        NON:ERREUR
         WRITE(IMPRIM,10006) MOTSMG,I
10006 FORMAT(' PPPOBA:ERREUR LE NOMBRE DE MOTS DE MCG=',I12,
     %' DOIT ETRE SUPERIEUR A ',I12/)
         CALL ARRET( 100 )
      ENDIF
C
C     LES TABLEAUX REZLMN,G ET K SONT AU BOUT DU SUPER-TABLEAU MCG
      MGZLMG = MOTSMG - MOZLMG - MOZLMN - MOZLMK + 1
      CALL TAMCIN( MCG,MOTSMG,MOTSMG-MOZLMN-MOZLMG-MOZLMK,M2ZLMG,MGZLMG)
      MGZLMK = MGZLMG + MOZLMG
      CALL TAMCIN( MCG , MOTSMG , MOTSMK , M2ZLMK , MGZLMK )
      MGZLMN = MGZLMK + MOZLMK
      CALL TAMCIN( MCG , MOTSMG , MOTSMN , M2ZLMN , MGZLMN )
C
C     LES CARACTERISTIQUES DU FICHIER NFPOBA EN ACCES DIRECT
C     ------------------------------------------------------
      CALL TRUNIT( NFPOBA )
CCC      NBPAGE = 130
      MOPAGE = 512
C
C     DEFINITION ET OUVERTURE DE NFPOBA
C     ---------------------------------
      KNOM = HOMDIR // '/pp/pxyz'
      WRITE(IMPRIM,*) 'OPEN of ',KNOM 
      OPEN( UNIT=NFPOBA, ERR=9900, STATUS='NEW',
     %      FILE=KNOM  , ACCESS='DIRECT',
     %      FORM='UNFORMATTED', IOSTAT=IOERR,
     %      RECL= MOPAGE * NBCHMO )
C
C     INITIALISATION DU COMMON/ELFINI/ ET PROTECTION DANS PAGE 1
C     ----------------------------------------------------------
      DO 2 I=1,MOPAGE
         NELFI(I) = 0
    2 ENDDO
      L = (MOPAGE - 2) / 4
      DO 3 I=1,L
         ITA(1,I) = IKPPP
    3 ENDDO
C     ECRITURE SUR LE FICHIER POBA DE LA PAGE 1
      NOPAGE = 1
      CALL FIPAEC( NFPOBA , NOPAGE , MOPAGE , NELFI )
c
C     LES DONNEES SONT LUES SUR LE FICHIER ./td/p/poba
C     ------------------------------------------------
      CALL TRUNIT( NFDONN )
      KNOM = HOMDIR // '/td/p/poba'
      WRITE(IMPRIM,*) 'OPEN of ',KNOM 
      OPEN( UNIT=NFDONN , ERR=9900 , STATUS='OLD' ,
     %      FILE=KNOM   , ACCESS='SEQUENTIAL' ,
     %      FORM='FORMATTED' , IOSTAT=IOERR )
C
C     LE FICHIER DES DONNEES EST NFDONN
C     MISE A JOUR DANS LE FORMAT LIBRE
      IPILE  = IPILE + 1
      MPFIC(IPILE) = NFDONN
      MPLIG(IPILE) = 0
      LECTEU = NFDONN
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     +++++            LES FORMULES D INTEGRATION NUMERIQUE        +++++
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     AUX EXTREMITES DE L'INTERVALLE [0,1]
C     A NPISE1=2 POINTS SUR LE SEGMENT UNITE EXACTE POUR P1
C     ==================================================================
      POSEP1(1) = 0.5D0
      POSEP1(2) = 0.5D0
      COSEP1(1) = 0.D0
      COSEP1(2) = 1.D0
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     DE GAUSS-LEGENDRE
C     A NPISE3=2 POINTS SUR LE SEGMENT UNITE EXACTE POUR P3
C     ==================================================================
      POSEP3(1) = 0.5D0
      POSEP3(2) = 0.5D0

      RAC       = 0.5D0 / DSQRT(3.D0)
      COSEP3(1) = 0.5D0 - RAC
      COSEP3(2) = 0.5D0 + RAC
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     DE GAUSS-LEGENDRE
C     A NPISE5=3 POINTS SUR LE SEGMENT UNITE EXACTE POUR P5
C     ==================================================================
      POSEP5(1) = 5.0D0 / 18.D0
      POSEP5(2) = 8.0D0 / 18.D0
      POSEP5(3) = POSEP5(1)

      RAC       = 0.5D0 * DSQRT(3.D0 / 5.D0)
      COSEP5(1) = 0.5D0 - RAC
      COSEP5(2) = 0.5D0
      COSEP5(3) = 0.5D0 + RAC
c
cccc     Pol no(Bl) 1:S1 2:S2 3:Milieu (Calcul pour f2sp2p1.f)
ccc      do l=1,3
ccc         p = COSEP5(l)
ccc         PSEBL( 1, L ) = ( 1D0 - P ) * ( 1D0 - 2D0 * P )
ccc         PSEBL( 2, L ) =  P * ( 2D0 * P - 1D0 )
ccc         PSEBL( 3, L ) =  4D0 * P * ( 1D0 - P )
ccc      enddo
ccc      print 10101, posep5
ccc      print 10101, cosep5
ccc      print 10101, psebl
ccc10101 format('     %',D24.17,', ',D24.17,',')
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     DE GAUSS-LEGENDRE
C     A NPISE7=4 POINTS SUR LE SEGMENT UNITE EXACTE POUR P7
C     ==================================================================
      POSEP7(1) = ( 18D0 + SQRT(30D0) ) / 72D0
      POSEP7(2) = POSEP7(1)
      POSEP7(3) = ( 18D0 - SQRT(30D0) ) / 72D0
      POSEP7(4) = POSEP7(3)

      COSEP7(1) = (1D0 - SQRT( (3D0 - 2D0 * SQRT(6D0/5D0) ) /7D0 ))/ 2D0
      COSEP7(2) = (1D0 + SQRT( (3D0 - 2D0 * SQRT(6D0/5D0) ) /7D0 ))/ 2D0
      COSEP7(3) = (1D0 - SQRT( (3D0 + 2D0 * SQRT(6D0/5D0) ) /7D0 ))/ 2D0
      COSEP7(4) = (1D0 + SQRT( (3D0 + 2D0 * SQRT(6D0/5D0) ) /7D0 ))/ 2D0
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     DE GAUSS-LEGENDRE
C     A NPISE9=5 POINTS SUR LE SEGMENT UNITE EXACTE POUR P9
C     ==================================================================
      POSEP9(1) = ( 322D0 + 13 * SQRT(70D0) ) / 1800D0
      POSEP9(2) = POSEP9(1)
      POSEP9(3) = 128D0 / 450D0
      POSEP9(4) = ( 322D0 - 13 * SQRT(70D0) ) / 1800D0
      POSEP9(5) = POSEP9(4)

      COSEP9(1) =( 1D0 - SQRT( 5D0 - 2D0 * SQRT(10D0/7D0) ) /3D0 ) / 2D0
      COSEP9(2) =( 1D0 + SQRT( 5D0 - 2D0 * SQRT(10D0/7D0) ) /3D0 ) / 2D0
      COSEP9(3) =  0.5D0
      COSEP9(4) =( 1D0 - SQRT( 5D0 + 2D0 * SQRT(10D0/7D0) ) /3D0 ) / 2D0
      COSEP9(5) =( 1D0 + SQRT( 5D0 + 2D0 * SQRT(10D0/7D0) ) /3D0 ) / 2D0
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=3 POINTS SUR LE TRIANGLE UNITE EXACTE POUR P1
C     ==================================================================
C     LES COORDONNEES
C     ---------------
      COTRP1(1,1)  = 0.D0
      COTRP1(2,1)  = 0.D0
      COTRP1(1,2)  = 1.D0
      COTRP1(2,2)  = 0.D0
      COTRP1(1,3)  = 0.D0
      COTRP1(2,3)  = 1.D0
C
C     LES POIDS
C     ---------
      RAC       = 1.D0 / 6.D0
      POTRP1(1) = RAC
      POTRP1(2) = RAC
      POTRP1(3) = RAC
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=3 POINTS SUR LE TRIANGLE UNITE EXACTE POUR P2
C     ==================================================================
C     LES COORDONNEES
C     ---------------
      COTRP2(1,1)  = 1.D0 / 6.D0
      COTRP2(2,1)  = COTRP2(1,1)
      COTRP2(1,2)  = 2.D0 / 3.D0
      COTRP2(2,2)  = COTRP2(1,1)
      COTRP2(1,3)  = COTRP2(1,1)
      COTRP2(2,3)  = COTRP2(1,2)
C
C     LES POIDS
C     ---------
      POTRP2(1)  = 1.D0 / 6.D0
      POTRP2(2)  = POTRP2(1)
      POTRP2(3)  = POTRP2(1)
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPITR5=7 POINTS SUR LE TRIANGLE UNITE EXACTE POUR P5
C     ==================================================================
      RAC = DSQRT(15.D0)
C
C     LES COORDONNEES
C     ---------------
      COTRP5(1,1)  = (6.D0 -        RAC) / 21.D0
      COTRP5(2,1)  = COTRP5(1,1)
      COTRP5(1,2)  = (9.D0 + 2.D0 * RAC) / 21.D0
      COTRP5(2,2)  = COTRP5(1,1)
      COTRP5(1,3)  = COTRP5(1,1)
      COTRP5(2,3)  = COTRP5(1,2)
      COTRP5(1,4)  = (6.D0 +        RAC) / 21.D0
      COTRP5(2,4)  = (9.D0 - 2.D0 * RAC) / 21.D0
      COTRP5(1,5)  = COTRP5(1,4)
      COTRP5(2,5)  = COTRP5(1,4)
      COTRP5(1,6)  = COTRP5(2,4)
      COTRP5(2,6)  = COTRP5(1,4)
      COTRP5(1,7)  = 1.D0 / 3.D0
      COTRP5(2,7)  = COTRP5(1,7)
ccc      print 10101, cotrp5
ccc10101 format('     %',D24.17,', ',D24.17,',')
C
C     LES POIDS
C     ---------
      POTRP5(1)  = (155.D0 - RAC) / 2400.D0
      POTRP5(2)  = POTRP5(1)
      POTRP5(3)  = POTRP5(1)
      POTRP5(4)  = (155.D0 + RAC) / 2400.D0
      POTRP5(5)  = POTRP5(4)
      POTRP5(6)  = POTRP5(4)
      POTRP5(7)  = 9.D0 / 80.D0
ccc      print 10101, potrp5
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=4 POINTS SUR LE QUADRANGLE UNITE , EXACTE POUR Q1
C     ==================================================================
C     LES COORDONNEES COQUQ1 ET LES POIDS POQUQ1
C     ------------------------------------------
      COQUQ1(1,1) = 0.D0
      COQUQ1(2,1) = 0.D0
      COQUQ1(1,2) = 1.D0
      COQUQ1(2,2) = 0.D0
      COQUQ1(1,3) = 0.D0
      COQUQ1(2,3) = 1.D0
      COQUQ1(1,4) = 1.D0
      COQUQ1(2,4) = 1.D0
      POQUQ1( 1 ) = 0.25D0
      POQUQ1( 2 ) = 0.25D0
      POQUQ1( 3 ) = 0.25D0
      POQUQ1( 4 ) = 0.25D0
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=4 POINTS SUR LE QUADRANGLE UNITE , EXACTE POUR Q3
C     ==================================================================
C     LES COORDONNEES COQUQ3 ET LES POIDS POQUQ3
C     ------------------------------------------
      L = 0
      DO 100 I=1,2
         DO 90 J=1,2
            L = L + 1
            POQUQ3( L )  = POSEP3(I) * POSEP3(J)
            COQUQ3(1,L)  = COSEP3(J)
            COQUQ3(2,L)  = COSEP3(I)
 90      ENDDO
 100  ENDDO
C
C     PERMUTATION DES POINTS 3 ET 4
      DO 105 I=1,2
         D           = COQUQ3(I,3)
         COQUQ3(I,3) = COQUQ3(I,4)
         COQUQ3(I,4) = D
 105  ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=9 POINTS SUR LE QUADRANGLE UNITE EXACTE POUR Q5
C     ==================================================================
C     LES COORDONNEES COQUQ5 ET LES POIDS POQUQ5
C     ------------------------------------------
      L = 0
      DO 120 I=1,3
         DO 110 J=1,3
            L = L + 1
            POQUQ5( L ) = POSEP5(I) * POSEP5(J)
            COQUQ5(1,L) = COSEP5(J)
            COQUQ5(2,L) = COSEP5(I)
 110     ENDDO
 120  ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=16 POINTS SUR LE QUADRANGLE UNITE EXACTE POUR Q7
C     ==================================================================
C     LES COORDONNEES COQUQ7 ET LES POIDS POQUQ7
C     ------------------------------------------
      L = 0
      DO 140 I=1,4
         DO 130 J=1,4
            L = L + 1
            POQUQ7( L ) = POSEP7(I) * POSEP7(J)
            COQUQ7(1,L) = COSEP7(J)
            COQUQ7(2,L) = COSEP7(I)
 130     ENDDO
 140  ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=9 POINTS SUR LE QUADRANGLE UNITE EXACTE POUR Q9
C     ==================================================================
C     LES COORDONNEES COQUQ9 ET LES POIDS POQUQ9
C     ------------------------------------------
      L = 0
      DO 160 I=1,5
         DO 150 J=1,5
            L = L + 1
            POQUQ9( L ) = POSEP9(I) * POSEP9(J)
            COQUQ9(1,L) = COSEP9(J)
            COQUQ9(2,L) = COSEP9(I)
 150     ENDDO
 160  ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=4 POINTS SUR LE TETRAEDRE UNITE EXACTE POUR P1
C     ==================================================================
      POTEP1(1)   = 1D0/24D0
      POTEP1(2)   = 1D0/24D0
      POTEP1(3)   = 1D0/24D0
      POTEP1(4)   = 1D0/24D0
C
      COTEP1(1,1) = 0.D0
      COTEP1(2,1) = 0.D0
      COTEP1(3,1) = 0.D0
C
      COTEP1(1,2) = 1.D0
      COTEP1(2,2) = 0.D0
      COTEP1(3,2) = 0.D0
C
      COTEP1(1,3) = 0.D0
      COTEP1(2,3) = 1.D0
      COTEP1(3,3) = 0.D0
C
      COTEP1(1,4) = 0.D0
      COTEP1(2,4) = 0.D0
      COTEP1(3,4) = 1.D0
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=15 POINTS SUR LE TETRAEDRE UNITE EXACTE POUR P5
C     ==================================================================
      POTEP5(1)   = 8D0/405D0
      COTEP5(1,1) = 0.25D0
      COTEP5(2,1) = 0.25D0
      COTEP5(3,1) = 0.25D0
C
      RAC         = DSQRT(15D0)
      D           = (7D0 - RAC)/34D0
      DD          = (13D0 + 3D0 * RAC)/34D0
      S           = (2665D0 + 14D0 * RAC)/226800D0
      J           = 0
C
  295 DO 300 I=1,4
         K             = J + 2
         POTEP5(K)     = S
         COTEP5(1,K)   = D
         COTEP5(2,K)   = D
         COTEP5(3,K)   = D
         K             = J + 3
         POTEP5(K)     = S
         COTEP5(1,K)   = DD
         COTEP5(2,K)   = D
         COTEP5(3,K)   = D
         K             = J + 4
         POTEP5(K)     = S
         COTEP5(1,K)   = D
         COTEP5(2,K)   = DD
         COTEP5(3,K)   = D
         K             = J + 5
         POTEP5(K)     = S
         COTEP5(1,K)   = D
         COTEP5(2,K)   = D
         COTEP5(3,K)   = DD
  300 ENDDO
C
      IF( J .EQ. 4 ) GOTO 303
      D  = (7D0 + RAC) / 34D0
      DD = (13D0 - 3D0 * RAC) / 34D0
      S  = (2665D0 - 14D0 * RAC) / 226800D0
      J  = 4
      GOTO 295
C
  303 D  = (10D0 - 2D0 * RAC) / 40D0
      DD = (10D0 + 2D0 * RAC) / 40D0
      S  = 5D0 / 567D0
C
      DO 310 I=10,15
         POTEP5(I) = S
  310 ENDDO
      COTEP5(1,10) = D
      COTEP5(2,10) = DD
      COTEP5(3,10) = DD
      COTEP5(1,11) = D
      COTEP5(2,11) = D
      COTEP5(3,11) = DD
      COTEP5(1,12) = DD
      COTEP5(2,12) = D
      COTEP5(3,12) = D
      COTEP5(1,13) = DD
      COTEP5(2,13) = DD
      COTEP5(3,13) = D
      COTEP5(1,14) = DD
      COTEP5(2,14) = D
      COTEP5(3,14) = DD
      COTEP5(1,15) = D
      COTEP5(2,15) = DD
      COTEP5(3,15) = D
C
ccc      print 10101, potep5
ccc      print *
ccc      print 10101, cotep5
ccc10101 format('     %',D24.17,', ',D24.17,',')
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=6  POINTS SUR LE PENTAEDRE UNITE EXACTE POUR R1
C     ==================================================================
      COPER1(1,1) = 0D0
      COPER1(2,1) = 0D0
      COPER1(3,1) = 0D0
C
      COPER1(1,2) = 1D0
      COPER1(2,2) = 0D0
      COPER1(3,2) = 0D0
C
      COPER1(1,3) = 0D0
      COPER1(2,3) = 1D0
      COPER1(3,3) = 0D0
C
      COPER1(1,4) = 0D0
      COPER1(2,4) = 0D0
      COPER1(3,4) = 1D0
C
      COPER1(1,5) = 1D0
      COPER1(2,5) = 0D0
      COPER1(3,5) = 1D0
C
      COPER1(1,6) = 0D0
      COPER1(2,6) = 1D0
      COPER1(3,6) = 1D0
C
      DO 315 I=1,6
         POPER1(I) = 1D0 / 6D0
  315 ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI= 6 POINTS SUR LE PENTAEDRE UNITE EXACTE POUR R2
C     ==================================================================
      L = 0
      DO 312 J=1,2
         DO 311 I=1,3
            L = L + 1
            POPER2( L ) = POTRP2(I) * POSEP3(J)
            COPER2(1,L) = COTRP2(1,I)
            COPER2(2,L) = COTRP2(2,I)
            COPER2(3,L) = COSEP3( J )
  311    ENDDO
  312 ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=21 POINTS SUR LE PENTAEDRE UNITE EXACTE POUR R5
C     ==================================================================
      L = 0
      DO 330 J=1,3
         DO 320 I=1,7
            L = L + 1
            POPER5(L)   = POTRP5(I) * POSEP5(J)
            COPER5(1,L) = COTRP5(1,I)
            COPER5(2,L) = COTRP5(2,I)
            COPER5(3,L) = COSEP5( J )
  320    ENDDO
  330 ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=8  POINTS SUR L HEXAEDRE UNITE EXACTE POUR Q1
C     ==================================================================
      COHEQ1(1,1) = 0D0
      COHEQ1(2,1) = 0D0
      COHEQ1(3,1) = 0D0
C
      COHEQ1(1,2) = 1D0
      COHEQ1(2,2) = 0D0
      COHEQ1(3,2) = 0D0
C
      COHEQ1(1,3) = 0D0
      COHEQ1(2,3) = 1D0
      COHEQ1(3,3) = 0D0
C
      COHEQ1(1,4) = 1D0
      COHEQ1(2,4) = 1D0
      COHEQ1(3,4) = 0D0
C
      COHEQ1(1,5) = 0D0
      COHEQ1(2,5) = 0D0
      COHEQ1(3,5) = 1D0
C
      COHEQ1(1,6) = 1D0
      COHEQ1(2,6) = 0D0
      COHEQ1(3,6) = 1D0
C
      COHEQ1(1,7) = 0D0
      COHEQ1(2,7) = 1D0
      COHEQ1(3,7) = 1D0
C
      COHEQ1(1,8) = 1D0
      COHEQ1(2,8) = 1D0
      COHEQ1(3,8) = 1D0
C
      DO 335 I=1,8
         POHEQ1(I) = 0.125D0
  335 ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI= 8 POINTS SUR L HEXAEDRE UNITE EXACTE POUR Q3
C     ==================================================================
      L = 0
      DO 333 K=1,2
         DO 332 J=1,2
            DO 331 I=1,2
               L = L + 1
               POHEQ3( L ) = POSEP3(I) * POSEP3(J) * POSEP3(K)
               COHEQ3(1,L) = COSEP3(I)
               COHEQ3(2,L) = COSEP3(J)
               COHEQ3(3,L) = COSEP3(K)
  331       ENDDO
  332    ENDDO
  333 ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=27 POINTS SUR L HEXAEDRE UNITE EXACTE POUR Q5
C     ==================================================================
      L = 0
      DO 360 K=1,3
         DO 350 J=1,3
            DO 340 I=1,3
               L = L + 1
               POHEQ5( L ) = POSEP5(I) * POSEP5(J) * POSEP5(K)
               COHEQ5(1,L) = COSEP5(I)
               COHEQ5(2,L) = COSEP5(J)
               COHEQ5(3,L) = COSEP5(K)
  340       ENDDO
  350    ENDDO
  360 ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=64 POINTS SUR L HEXAEDRE UNITE EXACTE POUR Q7
C     ==================================================================
      L = 0
      DO 365 K=1,4
         DO 363 J=1,4
            DO 361 I=1,4
               L = L + 1
               POHEQ7( L ) = POSEP7(I) * POSEP7(J) * POSEP7(K)
               COHEQ7(1,L) = COSEP7(I)
               COHEQ7(2,L) = COSEP7(J)
               COHEQ7(3,L) = COSEP7(K)
 361        ENDDO
 363     ENDDO
 365  ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=125 POINTS SUR L HEXAEDRE UNITE EXACTE POUR Q9
C     ==================================================================
      L = 0
      DO 369 K=1,5
         DO 367 J=1,5
            DO 366 I=1,5
               L = L + 1
               POHEQ9( L ) = POSEP9(I) * POSEP9(J) * POSEP9(K)
               COHEQ9(1,L) = COSEP9(I)
               COHEQ9(2,L) = COSEP9(J)
               COHEQ9(3,L) = COSEP9(K)
 366        ENDDO
 367     ENDDO
 369  ENDDO
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=5  POINTS SUR LA PYRAMIDE UNITE EXACTE POUR Q1xP1
C     ==================================================================
      COPYQ1(1,1) = 0D0
      COPYQ1(2,1) = 0D0
      COPYQ1(3,1) = 0D0
C
      COPYQ1(1,2) = 1D0
      COPYQ1(2,2) = 0D0
      COPYQ1(3,2) = 0D0
C
      COPYQ1(1,3) = 1D0
      COPYQ1(2,3) = 1D0
      COPYQ1(3,3) = 0D0
C
      COPYQ1(1,4) = 0D0
      COPYQ1(2,4) = 1D0
      COPYQ1(3,4) = 0D0
C
      COPYQ1(1,5) = 0D0
      COPYQ1(2,5) = 0D0
      COPYQ1(3,5) = 1D0
C
      DO 370 I=1,5
         POPYQ1(I) = 1D0 / 15D0
 370  ENDDO
C
      CALL VALXYZ( 5, POPYQ1, 3, COPYQ1, V )
ccc      print *,'Integrale x2 y2 z2 sur Pyramide PYQ1=',V
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=8 POINTS SUR LA PYRAMIDE UNITE EXACTE POUR Q3xP3
C     ==================================================================
      L = 0
      DO 373 K=1,2
         DO 372 J=1,2
            DO 371 I=1,2
               L = L + 1
               POPYQ3( L ) = POSEP3(I) * POSEP3(J) * POSEP3(K)
               COPYQ3(1,L) = COSEP3(I) * COSEP3(K)
               COPYQ3(2,L) = COSEP3(J) * COSEP3(K)
               COPYQ3(3,L) = COSEP3(K)
  371       ENDDO
  372    ENDDO
  373 ENDDO
C
      CALL VALXYZ( 8, POPYQ3, 3, COPYQ3, V )
ccc      print *,'Integrale x2 y2 z2 sur Pyramide PYQ3=',V
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=27 POINTS SUR LA PYRAMIDE UNITE EXACTE POUR Q5xP5
C     ==================================================================
      L = 0
      DO 390 K=1,3
         DO 385 J=1,3
            DO 380 I=1,3
               L = L + 1
               POPYQ5( L ) = POSEP5(I) * POSEP5(J) * POSEP5(K)
               COPYQ5(1,L) = COSEP5(I) * COSEP5(K)
               COPYQ5(2,L) = COSEP5(J) * COSEP5(K)
               COPYQ5(3,L) = COSEP5(K)
 380        ENDDO
 385     ENDDO
 390  ENDDO
C
      CALL VALXYZ( 27, POPYQ5, 3, COPYQ5, V )
ccc      print *,'Integrale x2 y2 z2 sur Pyramide PYQ5=',V
C
C     ==================================================================
C     CALCUL DES POIDS ET COORDONNEES D UNE FORMULE D INTEGRATION
C     A NPI=64  POINTS SUR LE 6-CUBE XYZUVW OBTENUE PAR PRODUIT TENSORIEL DE
C     CELLE SUR L HEXAEDRE UNITE XYZ EXACTE POUR Q1 A 8 POINTS=SOMMETS ET
C     CELLE SUR L HEXAEDRE UNITE UVW EXACTE POUR Q3 A 8 POINTS DE GAUSS
C     ==================================================================
      LL = 0
      DO 436 N=1,2
         DO 435 M=1,2
            DO 434 L=1,2
C              LE POIDS POUR LA FORMULE EXACTE POUR Q3(U,V,W)
               PP = POSEP3(L) * POSEP3(M) * POSEP3(N)
CCCC              LE POIDS POUR LA FORMULE EXACTE POUR Q1(U,V,W)
CCC               PP = POSEP1(L) * POSEP1(M) * POSEP1(N)
               DO 433 K=1,2
                  DO 432 J=1,2
                     DO 431 I=1,2
                        LL = LL + 1
C
C                       LE POIDS POUR LA FORMULE EXACTE POUR Q1(X,Y,Z)
                        PO6CQ1( LL ) = PP*POSEP1(I)*POSEP1(J)*POSEP1(K)
CCCC                       LE POIDS POUR LA FORMULE EXACTE POUR Q3(X,Y,Z)
CCC                        PO6CQ1( LL ) = PP*POSEP3(I)*POSEP3(J)*POSEP3(K)
C
C                       LES 6 COORDONNEES
                        CO6CQ1(1,LL) = COSEP1(I)
                        CO6CQ1(2,LL) = COSEP1(J)
                        CO6CQ1(3,LL) = COSEP1(K)
                        CO6CQ1(4,LL) = COSEP3(L)
                        CO6CQ1(5,LL) = COSEP3(M)
                        CO6CQ1(6,LL) = COSEP3(N)
CCC
CCC                        CO6CQ1(1,LL) = COSEP3(I)
CCC                        CO6CQ1(2,LL) = COSEP3(J)
CCC                        CO6CQ1(3,LL) = COSEP3(K)
CCC                        CO6CQ1(4,LL) = COSEP3(L)
CCC                        CO6CQ1(5,LL) = COSEP3(M)
CCC                        CO6CQ1(6,LL) = COSEP3(N)
C
 431                 ENDDO
 432              ENDDO
 433           ENDDO
 434        ENDDO
 435     ENDDO
 436  ENDDO
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     +++++            LA GENERATION DES TABLEAUX NUMERIQUES       +++++
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C
C     ==================================================================
C     CALCUL DES VALEURS DES DERIVEES DES  POLYNOMES EN CERTAINS POINTS
C     et AFFICHAGE SIMPLE
C     RECUPERATION POUR FICHIER INCLUDE DANS LES SUBROUTINES DE CALCUL
C     ==================================================================
C     TABLEAU DP2(3 Sommets) sur le TRIANGLE UNITE
      CALL PN2DPPT( 2,  6, TRP2, 3, COTRP1, TABLO )
C
C     TABLEAU DP2(4 Sommets) sur le TETRAEDRE UNITE
      CALL PN3DPPT( 2, 10, TEP2, 4, COTEP1, TABLO )
C
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     SEGMENT  P1-LAGRANGE .LES COEFFICIENTS DES POLYNOMES
C     --------------------- PRECEDES DU NOMBRE DE VARIABLES
C                                    DU DEGRE+1 DES POLYNOMES
C                                    DU TYPE DE POLYNOMES (P)
C                                    DU NOMBRE DE POLYNOMES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL COPOLY('SEP1',IA,L,1,1,0,2,SEP1)
      CALL TRPOBA('SEP1',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION DES CARACTERISTIQUES ET VALEURS DES POLYNOMES
      IA1 = (IA + 3) / 2
      LL  = (L  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     SEGMENT  P1-LAGRANGE .LES POLYNOMES
C    ---------------------  SONT CALCULES AUX 2 POINTS DE LA
C                           FORMULE D INTEGRATION NUMERIQUE DE
C                           GAUSS-LEGENDRE EXACTE POUR LES POLYNOMES P3
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('1P13',NFPOBA,MOPAGE,NOPAGE,
     %            1,1,1,1,1,
     %            1,2,2,POSEP3,COSEP3,SEP1,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERSOM('2F11',2,2,1,SEP1,DPXYZ,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     SEGMENT  P2-LAGRANGE .LES COEFFICIENTS DES POLYNOMES
C     --------------------- PRECEDES DU NOMBRE DE VARIABLES
C                                    DU DEGRE+1 DES POLYNOMES
C                                    DU TYPE DE POLYNOMES (P)
C                                    DU NOMBRE DE POLYNOMES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL COPOLY('SEP2',IA,L,1,2,0,3,SEP2)
      CALL TRPOBA('SEP2',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION DES CARACTERISTIQUES ET VALEURS DES POLYNOMES
      IA1 = (IA + 3) / 2
      LL  = (L  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     SEGMENT  P2-LAGRANGE .LES POLYNOMES
C     --------------------- SONT CALCULES AUX 3 POINTS DE LA
C                           FORMULE D INTEGRATION NUMERIQUE DE
C                           GAUSS-LEGENDRE EXACTE POUR LES POLYNOMES P5
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('1P25',NFPOBA,MOPAGE,NOPAGE,
     %            1,1,1,1,1,
     %            2,3,3,POSEP5,COSEP5,SEP2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERSOM('2F21',2,3,2,SEP2,DPXYZ,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     SEGMENT  P2-LAGRANGE .LES POLYNOMES
C     --------------------- SONT CALCULES AUX 4 POINTS DE LA
C                           FORMULE D INTEGRATION NUMERIQUE DE
C                           GAUSS-LEGENDRE EXACTE POUR LES POLYNOMES P7
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('1P27',NFPOBA,MOPAGE,NOPAGE,
     %            1,1,1,1,1,
     %            2,3,4,POSEP7,COSEP7,SEP2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     SEGMENT  P2-LAGRANGE .LES POLYNOMES
C     --------------------- SONT CALCULES AUX 5 POINTS DE LA
C                           FORMULE D INTEGRATION NUMERIQUE DE
C                           GAUSS-LEGENDRE EXACTE POUR LES POLYNOMES P9
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('1P29',NFPOBA,MOPAGE,NOPAGE,
     %            1,1,1,1,1,
     %            2,3,5,POSEP9,COSEP9,SEP2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE P1-LAGRANGE .LES COEFFICIENTS DES POLYNOMES
C     --------------------- PRECEDES DU NOMBRE DE VARIABLES
C                                    DU DEGRE+1 DES POLYNOMES
C                                    DU TYPE DE POLYNOMES (P OU Q OU R)
C                                    DU NOMBRE DE POLYNOMES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL COPOLY('TRP1',IA,L,2,1,0,3,TRP1)
      CALL TRPOBA('TRP1',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION DES CARACTERISTIQUES ET VALEURS DES POLYNOMES
      IA1 = (IA + 3) / 2
      LL  = (L  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE P2-LAGRANGE .LES COEFFICIENTS DES POLYNOMES
C     --------------------- PRECEDES DU NOMBRE DE VARIABLES
C                                    DU DEGRE+1 DES POLYNOMES
C                                    DU TYPE DE POLYNOMES (P OU Q OU R)
C                                    DU NOMBRE DE POLYNOMES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CCC   ECRASEMENT INCONNU SUR ISIS2.EDCSM et PAS SUR LES AUTRES!...
      TRP2(1,1,1)=1D0
CCC   CORRIGE CET ECRASEMENT...
      CALL COPOLY('TRP2',IA,L,2,2,0,6,TRP2)
      CALL TRPOBA('TRP2',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION DES CARACTERISTIQUES ET VALEURS DES POLYNOMES
      IA1 = (IA + 3) / 2
      LL  = (L  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE Q1-LAGRANGE .LES COEFFICIENTS DES POLYNOMES
C     ----------------------- PRECEDES DU NOMBRE DE VARIABLES
C                                    DU DEGRE+1 DES POLYNOMES
C                                    DU TYPE DE POLYNOMES (P OU Q OU R)
C                                    DU NOMBRE DE POLYNOMES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL COPOLY('QUQ1',IA,L,2,1,1,4,QUQ1)
      CALL TRPOBA('QUQ1',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION DES CARACTERISTIQUES ET VALEURS DES POLYNOMES
      IA1 = (IA + 3) / 2
      LL  = (L  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE Q2-LAGRANGE .LES COEFFICIENTS DES POLYNOMES
C     ----------------------- PRECEDES DU NOMBRE DE VARIABLES
C                                    DU DEGRE+1 DES POLYNOMES
C                                    DU TYPE DE POLYNOMES (P OU Q OU R)
C                                    DU NOMBRE DE POLYNOMES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL COPOLY('QUQ2',IA,L,2,2,1,8,QUQ2)
      CALL TRPOBA('QUQ2',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION DES CARACTERISTIQUES ET VALEURS DES POLYNOMES
      IA1 = (IA + 3) / 2
      LL  = (L  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE P1-LAGRANGE .LES POLYNOMES ET LEURS DERIVEES
C     --------------------- SONT CALCULES AUX 3 POINTS D UNE
C                           FORMULE D INTEGRATION NUMERIQUE
C                           EXACTE POUR LES POLYNOMES P2
C                           DERIVEES DES POLYNOMES AUX POINTS D
C                           INTEGRATION DES 3 COTES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('2P12',NFPOBA,MOPAGE,NOPAGE,
     %            2,1,1,1,1,
     %            1,3,3,POTRP2,COTRP2,TRP1,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERARE('3F13',3,3,1,TRP1,DPXYZ,2,COSEP3,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE P2-LAGRANGE .LES POLYNOMES ET LEURS DERIVEES
C     --------------------- SONT CALCULES AUX 7 POINTS D UNE
C                           FORMULE D INTEGRATION NUMERIQUE
C                           EXACTE POUR LES POLYNOMES P5
C                           DERIVEES DES POLYNOMES AUX POINTS D
C                           INTEGRATION DES 3 COTES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('2P25',NFPOBA,MOPAGE,NOPAGE,
     %            2,1,1,1,1,
     %            2,6,7,POTRP5,COTRP5,TRP2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERARE('3F25',3,6,2,TRP2,DPXYZ,3,COSEP5,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE Q1-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 4 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES Q3
C                             DERIVEES DES POLYNOMES AUX POINTS D
C                             INTEGRATION DES 4 COTES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('2Q13',NFPOBA,MOPAGE,NOPAGE,
     %            2,1,1,1,1,
     %            1,4,4,POQUQ3,COQUQ3,QUQ1,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERARE('4F13',4,4,1,QUQ1,DPXYZ,2,COSEP3,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE Q2LAGRANGE LES POLYNOMES ET LEURS DERIVEES
C     --------------------- SONT CALCULES AUX 9 POINTS D UNE
C                           FORMULE D INTEGRATION NUMERIQUE
C                           EXACTE POUR LES POLYNOMES Q5
C                           DERIVEES DES POLYNOMES AUX POINTS D
C                           INTEGRATION DES 4 COTES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('2Q25',NFPOBA,MOPAGE,NOPAGE,
     %            2,1,1,1,1,
     %            2,8,9,POQUQ5,COQUQ5,QUQ2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERARE('4F25',4,8,2,QUQ2,DPXYZ,3,COSEP5,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE Q2LAGRANGE LES POLYNOMES ET LEURS DERIVEES
C     --------------------- SONT CALCULES AUX 16 POINTS D UNE
C                           FORMULE D INTEGRATION NUMERIQUE
C                           EXACTE POUR LES POLYNOMES Q7
C                           DERIVEES DES POLYNOMES AUX POINTS D
C                           INTEGRATION DES 4 COTES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('2Q27',NFPOBA,MOPAGE,NOPAGE,
     %            2,1,1,1,1,
     %            2,8,16,POQUQ7,COQUQ7,QUQ2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERARE('4F27',4,8,2,QUQ2,DPXYZ,4,COSEP7,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE Q2LAGRANGE LES POLYNOMES ET LEURS DERIVEES
C     --------------------- SONT CALCULES AUX 25 POINTS D UNE
C                           FORMULE D INTEGRATION NUMERIQUE
C                           EXACTE POUR LES POLYNOMES Q9
C                           DERIVEES DES POLYNOMES AUX POINTS D
C                           INTEGRATION DES 4 COTES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('2Q29',NFPOBA,MOPAGE,NOPAGE,
     %            2,1,1,1,1,
     %            2,8,25,POQUQ9,COQUQ9,QUQ2,DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERARE('4F29',4,8,2,QUQ2,DPXYZ,5,COSEP9,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE MT21 MIXTE   FLUX P1 PAR ARETE
C     --------------------- TEMPERATURE P1 SUR LE TRIANGLE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      L  = 252
      IA =   0
      CALL TNMCDC( 'REEL2' , L , IA )
C
      CALL TMXT21( DMCN((IA+1)/2),     DMCN((IA+1+72)/2) ,
     %             DMCN((IA+1+144)/2), DMCN((IA+1+216)/2) )
C
      CALL TRPOBA('MT21',IA,L,NFPOBA,MOPAGE,NOPAGE)
      IA1 = (IA - 1) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L ,IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE MQ10 MIXTE FLUX P0 PAR ARETE
C     --------------------- TEMPERATURE Q0 SUR LE QUADRANGLE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      L      =  80
      IA     =   0
      CALL TNMCDC( 'REEL2' , L ,IA )
C
      CALL TMXQ10( DMCN((IA+1)/2) )
C
      CALL TRPOBA('MQ10',IA,L,NFPOBA,MOPAGE,NOPAGE)
      IA1 = (IA - 1) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     QUADRANGLE MQ21 MIXTE FLUX P1 PAR ARETE
C     --------------------- TEMPERATURE Q1 SUR LE QUADRANGLE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      L      = 708
      IA     =   0
      CALL TNMCDC( 'REEL2' , L , IA )
C
      CALL TMXQ21( DMCN((IA+1)/2) ,     DMCN((IA+1+216)/2) ,
     %             DMCN((IA+1+432)/2) , DMCN((IA+1+528)/2) )
C
      CALL TRPOBA('MQ21',IA,L,NFPOBA,MOPAGE,NOPAGE)
      IA1 = (IA - 1) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TETRAEDRE  M3T2 MIXTE FLUX P1 PAR FACE
C     --------------------- TEMPERATURE P1 DANS LE TETRAEDRE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      L  = 1560
      IA =   0
      CALL TNMCDC( 'REEL2' , L , IA )
C
      CALL TMX3T2( DMCN((IA+1)/2) ,     DMCN((IA+1+240)/2) ,
     %             DMCN((IA+1+480)/2) , DMCN((IA+1+720)/2) ,
     %             DMCN((IA+1+960)/2) , DMCN((IA+1+1200)/2) ,
     %             DMCN((IA+1+1440)/2) ,
     %             DMCN((IA+1+1536)/2))
C
      CALL TRPOBA('M3T2',IA,L,NFPOBA,MOPAGE,NOPAGE)
      IA1 = (IA - 1) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     HEXAEDRE  M3H1  MIXTE FLUX P0 PAR FACE
C     --------------------- TEMPERATURE P0 DANS L HEXAEDRE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      L  = 1584
      IA =   0
      CALL TNMCDC( 'REEL2' , L , IA )
C
      CALL TMX3H1( DMCN((IA+1)/2) ,     DMCN((IA+1+336)/2) ,
     %             DMCN((IA+1+496)/2) , DMCN((IA+1+656)/2) ,
     %             DMCN((IA+1+816)/2) , DMCN((IA+1+1072)/2),
     %             DMCN((IA+1+ 1328)/2))
C
      CALL TRPOBA('M3H1',IA,L,NFPOBA,MOPAGE,NOPAGE)
      IA1 = (IA - 1) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     HEXAEDRE M3H2  MIXTE FLUX Q1 PAR FACE
C     --------------------- TEMPERATURE Q1 DANS L HEXAEDRE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     NOMBRE DE MOTS DU TABLEAU 'M3H2'
      L  = 12090
C     ADRESSAGE AUTOMATIQUE DU TABLEAU
      IA = 0
C     ADRESSAGE EFFECTIF DANS LE SUPER TABLEAU M
      CALL TNMCDC( 'REEL2' , L , IA )
C
C     CALCUL DES VALEURS DU TABLEAU 'M3H2' COMPOSE DE PLUSIEURS
C     SOUS-TABLEAUX (CF L ADRESSAGE A TRAVERS M)
C
      CALL TMX3H2( DMCN((IA+1)/2)  ,      DMCN((IA+1+1134)/2) ,
     %             DMCN((IA+1+1674)/2) ,  DMCN((IA+1+2214)/2),
     %             DMCN((IA+1+2754)/2) ,  DMCN((IA+1+3618)/2) ,
     %             DMCN((IA+1+4482)/2) ,  DMCN((IA+1+5346)/2) ,
     %             DMCN((IA+1+7290)/2) ,  DMCN((IA+1+9234)/2) ,
     %             DMCN((IA+1+11178)/2) , DMCN((IA+1+11562)/2),
     %             DMCN((IA+1+11754)/2) )
C
C     INITIATIONS OBLIGATOIRES POUR LE TRANSFERT SUR LE FICHIER POBA
C     NO DU TABLEAU SUR POBA ET DANS LE COMMON /ELFINI/
      CALL TRPOBA('M3H2',IA,L,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IA - 1) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
C
C     DESTRUCTION DU TABLEAU 'M3H2' DU SUPER-TABLEAU M
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TETRAEDRE  P1-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 4 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES P1
C                             DERIVEES DES POLYNOMES AUX 3 POINTS
C                             D INTEGRATION DES 4 FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'TEP1' DANS LE SUPER-TABLEAU M => 68 MOTS
C     4 POLYNOMES DE 2**3=8 COEFFICIENTS DOUBLE PRECISION CHACUN
      LTEP1  = 4 + 2 * 4 * 8
      IATEP1 = 0
      CALL TNMCDC( 'REEL2' , LTEP1 , IATEP1 )
C
C     LECTURE DES  8 MONOMES DES  4 POLYNOMES DU TETRAEDRE UNITE P1
C     NOMBRE DE VARIABLES ICI 3 : X , Y , Z
      MCN(IATEP1  ) = 3
C     DEGRE DES POLYNOMES + 1
      MCN(IATEP1+1) = 2
C     NPOUQ : 0 => POLYNOME EN (X,Y,Z)
C             1 => POLYNOME EN (X) , (Y) , (Z)
C             2 => POLYNOME EN (X,Y) , (Z)
      MCN(IATEP1+2) = 0
C     NOMBRE DE POLYNOMES
      MCN(IATEP1+3) = 4
      NCVALS = 0
C
C     LECTURE SUR LE FICHIER $MEFISTO/td/p/poba DES 8 COEFFICIENTS
C     DES 4 POLYNOMES DE BASE DU TETR 3P1C
      L = (LTEP1-4)/2
      CALL LIRTRD( NCVALS, L, DMCN( (IATEP1+1+4)/2 ) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('TEP1',IATEP1,LTEP1,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IATEP1 + 3) / 2
      LL  = (LTEP1  - 4) / 2
      WRITE (IMPRIM,10150)(ITA(I,NBTAST),I=1,4),(MCN(IATEP1-1+I),I=1,4),
     %                    (DMCN(IA1+I),I=1,LL)
C
C      CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C      AUX POINTS D INTEGRATION NUMERIQUE
C
      CALL SUPOBA('3P11',NFPOBA,MOPAGE,NOPAGE,
     %            3,1,1,1,1,
     %            1,4,4,POTEP1,COTEP1,DMCN((IATEP1+1+4)/2),DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('5F12',5,4,1,DMCN((1+IATEP1+4)/2),DPXYZ,
     %            3,COTRP2,3,COTRP2,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TETRAEDRE  P2-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 15 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES P5
C                             DERIVEES DES POLYNOMES AUX 7 POINTS
C                             D INTEGRATION DES 4 FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'TEP2' DANS LE SUPER-TABLEAU M => 544
C     10 POLYNOMES DE 3**3=27 COEFFICIENTS DOUBLE PRECISION CHACUN
      LTEP2  = 4 + 2 * 10 * 27
      IATEP2 = 0
      CALL TNMCDC( 'REEL2' , LTEP2 , IATEP2 )
C
C     LECTURE DES 27 MONOMES DES 10 POLYNOMES DU TETRAEDRE UNITE P2
      MCN(IATEP2  ) = 3
      MCN(IATEP2+1) = 3
      MCN(IATEP2+2) = 0
      MCN(IATEP2+3) = 10
      NCVALS = 0
      L = (LTEP2-4)/2
      CALL LIRTRD( NCVALS, L,  DMCN( (IATEP2+1+4)/2 ) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('TEP2',IATEP2,LTEP2,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IATEP2 + 3) / 2
      L   = (LTEP2  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IATEP2-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,L)
C
C     CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C     AUX POINTS D INTEGRATION NUMERIQUE
      CALL SUPOBA('3P25',NFPOBA,MOPAGE,NOPAGE,
     %            3,1,1,1,1,
     %            2,10,15,POTEP5,COTEP5,DMCN((IATEP2+1+4)/2),
     %            DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
cccC     CALCUL DE L'INTEGRALE DP2 DP2 dx dy dz sur LE TETRAEDRE
cccC     POUR L'EF de TAYLOR-HOOD en FLUIDE
ccc      CALL INTDPDP( 15, POTEP5, 3, 10, DPOLYP, DPDP )
C
      CALL DERFAC('5F25',5,10,2,DMCN((1+IATEP2+4)/2),DPXYZ,
     %            7,COTRP5,7,COTRP5,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     PYRAMIDE Q1P1-LAGRANGE  LES MONOMES DES POLYNOMES
C     ----------------------  DE BASE DE LA PYRAMIDE UNITE SONT
C                             STOCKES SUR LE FICHIER POBA
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'PYR1' DANS LE SUPER-TABLEAU M
C     5 POLYNOMES DE 2**3=8 COEFFICIENTS DOUBLE PRECISION CHACUN
      LPYR1  = 4 + 5 * 8 * 2
      IAPYR1 = 0
      CALL TNMCDC( 'REEL2' , LPYR1 , IAPYR1 )
C
C     LECTURE DES 8 MONOMES DES 5 POLYNOMES DE LA PYRAMIDE UNITE Q1P1
      MCN(IAPYR1  ) = 3
      MCN(IAPYR1+1) = 2
      MCN(IAPYR1+2) = 2
      MCN(IAPYR1+3) = 5
      NCVALS        = 0
      L = (LPYR1-4)/2
      CALL LIRTRD( NCVALS, L, DMCN( (IAPYR1+1+4)/2 ) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA( 'PYR1', IAPYR1, LPYR1, NFPOBA, MOPAGE, NOPAGE )
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IAPYR1 + 3) / 2
      LL  = (LPYR1  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IAPYR1-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     PYRAMIDE Q1P1-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX  8 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES Q1P1
C                             DERIVEES DES POLYNOMES AUX 5 SOMMETS DES
C                             5 FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('3PY1', NFPOBA, MOPAGE, NOPAGE,
     %            3,1,1,1,1,
     %            1,5,8,POPYQ3,COPYQ3,DMCN((IAPYR1+1+4)/2),DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('FPY1',9, 5,1,DMCN((1+IAPYR1+4)/2),DPXYZ,
     %            3,COTRP2,4,COQUQ3,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     PYRAMIDE Q2P2-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 27 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES Q5P5
C                             DERIVEES DES POLYNOMES AUX POINTS
C                             D INTEGRATION DES FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'PYR2' DANS LE SUPER-TABLEAU M
C     13 POLYNOMES DE 3**3=27 COEFFICIENTS DOUBLE PRECISION CHACUN
      IAPYR2 = 0
      LPYR2  = 4 + 13 * 27 * 2
      CALL TNMCDC( 'REEL2', LPYR2, IAPYR2 )
C
C     LECTURE DES 27 MONOMES DES 13 POLYNOMES DU PYRAMIDE UNITE R2
      MCN(IAPYR2  ) = 3
      MCN(IAPYR2+1) = 3
      MCN(IAPYR2+2) = 2
      MCN(IAPYR2+3) = 13
      NCVALS        = 0
      L = (LPYR2-4)/2
      CALL LIRTRD( NCVALS, L, DMCN((IAPYR2+1+4)/2) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('PYR2',IAPYR2,LPYR2,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IAPYR2 + 3) / 2
      L   = 13 * 27
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IAPYR2-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,L)
C
C     CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C     AUX POINTS D INTEGRATION NUMERIQUE
      CALL SUPOBA('3PY2',NFPOBA,MOPAGE,NOPAGE,
     %            3,1,1,1,1,
     %            2,13,27,POPYQ5,COPYQ5,DMCN((IAPYR2+1+4)/2),
     %            DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('FPY2',9,13,2,DMCN((1+IAPYR2+4)/2),DPXYZ,
     %            7,COTRP5,9,COQUQ5,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     PENTAEDRE  R1-LAGRANGE  LES MONOMES DES POLYNOMES
C     ----------------------  DE BASE DU PENTAEDRE UNITE SONT
C                             STOCKES SUR LE FICHIER POBA
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'PER1' DANS LE SUPER-TABLEAU M
C     6 POLYNOMES DE 2**3=8 COEFFICIENTS DOUBLE PRECISION CHACUN
      LPER1  = 4 + 2 * 6 * 8
      IAPER1 = 0
      CALL TNMCDC( 'REEL2' , LPER1 , IAPER1 )
C
C     LECTURE DES  8 MONOMES DES  6 POLYNOMES DU PENTAEDRE UNITE R1
      MCN(IAPER1  ) = 3
      MCN(IAPER1+1) = 2
      MCN(IAPER1+2) = 2
      MCN(IAPER1+3) = 6
      NCVALS        = 0
      L = (LPER1-4)/2
      CALL LIRTRD( NCVALS, L, DMCN((IAPER1+1+4)/2) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('PER1',IAPER1,LPER1,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IAPER1 + 3) / 2
      LL  = (LPER1  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IAPER1-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     PENTAEDRE  R1-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX  6 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES R1
C                             DERIVEES DES POLYNOMES AUX 4 SOMMETS DES
C                             5 FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('3R12',NFPOBA,MOPAGE,NOPAGE,
     %            3,1,1,1,1,
     %            1,6,6,POPER2,COPER2,DMCN((IAPER1+1+4)/2),DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('6F12',6, 6,1,DMCN((1+IAPER1+4)/2),DPXYZ,
     %            3,COTRP2,4,COQUQ3,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     PENTAEDRE  R2-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 21 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES R5
C                             DERIVEES DES POLYNOMES AUX POINTS
C                             D INTEGRATION DES FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'PER2' DANS LE SUPER-TABLEAU M
C     18 POLYNOMES DE 3**3=27 COEFFICIENTS DOUBLE PRECISION CHACUN
      IAPER2 = 0
      LPER2  = 4 + 27 * 18 * 2
      CALL TNMCDC( 'REEL2', LPER2, IAPER2 )
C
C     LECTURE DES 27 MONOMES DES 18 POLYNOMES DU PENTAEDRE UNITE R2
C
C     ATTENTION : IL S AGIT ICI DE L INTERPOLATION R2 COMPLETE
C     ----------- ELLE EST ENSUITE REDUITE PAR COMBINAISONS LINEAIRES
C                 C-A-D PAR SUPPRESSION DES POINTS BARYCENTRE DES CARRES
C
      MCN(IAPER2  ) = 3
      MCN(IAPER2+1) = 3
      MCN(IAPER2+2) = 2
      MCN(IAPER2+3) = 18
      NCVALS        = 0
      L = (LPER2-4)/2
      CALL LIRTRD( NCVALS, L, DMCN((IAPER2+1+4)/2) )
C
C     SUPPRESSION DU BARYCENTRE DES 3 FACES CARRE PAR COMBINAISON
C     LINEAIRE DES POLYNOMES
      CALL CL3R2C( IAPER2 )
C
C     MISE A JOUR DU NOMBRE DE POLYNOMES DE L INTERPOLATION
      MCN(IAPER2+3)=15
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('PER2',IAPER2,LPER2-162,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IAPER2 + 3) / 2
      L   = 15 * 27
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IAPER2-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,L)
C
C      CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C      AUX POINTS D INTEGRATION NUMERIQUE
      CALL SUPOBA('3R25',NFPOBA,MOPAGE,NOPAGE,
     %            3,1,1,1,1,
     %            2,15,21,POPER5,COPER5,DMCN((IAPER2+1+4)/2),
     %            DPXYZ,DDPXYZ,
     %            POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('6F25',6,15,2,DMCN((1+IAPER2+4)/2),DPXYZ,
     %            7,COTRP5,9,COQUQ5,
     %            NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     HEXAEDRE  Q1-LAGRANGE   LES MONOMES DES POLYNOMES
C     ----------------------  DE BASE DE L'HEXAEDRE UNITE SONT
C                             STOCKES SUR LE FICHIER POBA
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'HEQ1' DANS LE SUPER-TABLEAU M => 132 MOTS
C     8 POLYNOMES DE 2**3=8 COEFFICIENTS DOUBLE PRECISION CHACUN
      LHEQ1  = 4 + 2 * 8 * 8
      IAHEQ1 = 0
      CALL TNMCDC( 'REEL2' , LHEQ1 , IAHEQ1 )
C
C     LECTURE DES  8 MONOMES DES  8 POLYNOMES DE L HEXAEDRE UNITE Q1
      MCN(IAHEQ1  ) = 3
      MCN(IAHEQ1+1) = 2
      MCN(IAHEQ1+2) = 1
      MCN(IAHEQ1+3) = 8
      NCVALS = 0
C
C     LECTURE SUR LE FICHIER $MEFISTO/td/p/poba DES 68 COEFFICIENTS
C     DES 8 POLYNOMES DE BASE DU HEXA 3Q1C
      L = (LHEQ1-4)/2
      CALL LIRTRD( NCVALS, L, DMCN((IAHEQ1+1+4)/2) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('HEQ1',IAHEQ1,LHEQ1,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IAHEQ1 + 3) / 2
      LL  = (LHEQ1  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IAHEQ1-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     HEXAEDRE  Q1-LAGRANGE   LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 8 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES Q3
C                             DERIVEES DES POLYNOMES AUX POINTS
C                             D INTEGRATION DES FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA('3Q13',NFPOBA,MOPAGE,NOPAGE,
     %             3,1,1,1,1,
     %            1,8,8,POHEQ3,COHEQ3,DMCN((IAHEQ1+1+4)/2),DPXYZ,DDPXYZ,
     %             POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('7F13',7, 8,1,DMCN((1+IAHEQ1+4)/2),DPXYZ,
     %             4,COQUQ3,4,COQUQ3,
     %             NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     HEXAEDRE  Q2-LAGRANGE   LES POLYNOMES ET LEURS DERIVEES
C     ----------------------  SONT CALCULES AUX 27 POINTS D UNE
C                             FORMULE D INTEGRATION NUMERIQUE
C                             EXACTE POUR LES POLYNOMES Q5
C                             DERIVEES DES POLYNOMES AUX POINTS
C                             D INTEGRATION DES FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'HEQ2' DANS LE SUPER-TABLEAU M => 1462 MOTS
C     27 POLYNOMES DE 3**3=27 COEFFICIENTS DOUBLE PRECISION CHACUN
      LHEQ2  = 4 + 2 * 27 * 27
      IAHEQ2 = 0
      CALL TNMCDC( 'REEL2', LHEQ2, IAHEQ2 )
C
C     LECTURE DES 27 MONOMES DES 27 POLYNOMES DE L HEXAEDRE UNITE Q2
C
C     ATTENTION : IL S AGIT ICI DE L INTERPOLATION Q2 COMPLETE
C     ----------- ELLE EST ENSUITE REDUITE PAR COMBINAISONS LINEAIRES
C                 C-A-D PAR SUPPRESSION DES POINTS BARYCENTRE DES CARRES
      MCN(IAHEQ2  ) = 3
      MCN(IAHEQ2+1) = 3
      MCN(IAHEQ2+2) = 1
      MCN(IAHEQ2+3) = 27
      NCVALS        = 0
C
C     LECTURE SUR LE FICHIER $MEFISTO/td/p/poba DES 27 COEFFICIENTS
C     DES 27 POLYNOMES DE BASE DU HEXA 3Q2C COMPLET
      L = (LHEQ2-4)/2
      CALL LIRTRD( NCVALS, L, DMCN((IAHEQ2+1+4)/2) )
C
C     SUPPRESSION DU BARYCENTRE DES 6 FACES CARRE PAR COMBINAISON
C     LINEAIRE DES POLYNOMES
      CALL CL3Q2C( IAHEQ2 )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('HEQ2',IAHEQ2,LHEQ2-378,NFPOBA,MOPAGE,NOPAGE)
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IAHEQ2 + 3) / 2
      LL  = 20 * 27
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IAHEQ2-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
C
C     CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C     AUX 27 POINTS D INTEGRATION NUMERIQUE EXACTE Q5
      CALL SUPOBA('3Q25',NFPOBA,MOPAGE,NOPAGE,
     %             3,1,1,1,1,
     %             2,20,27,POHEQ5,COHEQ5,DMCN((IAHEQ2+1+4)/2),
     %             DPXYZ,DDPXYZ,
     %             POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('7F25',7,20,2,DMCN((1+IAHEQ2+4)/2),DPXYZ,
     %             9,COQUQ5,9,COQUQ5,
     %             NFPOBA,MOPAGE,NOPAGE)
C
C     CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C     AUX 64 POINTS D INTEGRATION NUMERIQUE EXACTE Q7
      CALL SUPOBA('3Q27',NFPOBA,MOPAGE,NOPAGE,
     %             3,1,1,1,1,
     %             2,20,64,POHEQ7,COHEQ7,DMCN((IAHEQ2+1+4)/2),
     %             DPXYZ,DDPXYZ,
     %             POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('7F27',7,20,2,DMCN((1+IAHEQ2+4)/2),DPXYZ,
     %             16,COQUQ7,16,COQUQ7,
     %             NFPOBA,MOPAGE,NOPAGE)
C
C     CALCUL DES VALEURS DES POLYNOMES ET DE SES DERIVEES
C     AUX 125 POINTS D INTEGRATION NUMERIQUE EXACTE Q9
      CALL SUPOBA('3Q29',NFPOBA,MOPAGE,NOPAGE,
     %             3,1,1,1,1,
     %             2,20,125,POHEQ9,COHEQ9,DMCN((IAHEQ2+1+4)/2),
     %             DPXYZ,DDPXYZ,
     %             POLYPI,DPOLYP,DDPOLY)
C
      CALL DERFAC('7F29',7,20,2,DMCN((1+IAHEQ2+4)/2),DPXYZ,
     %             25,COQUQ9,25,COQUQ9,
     %             NFPOBA,MOPAGE,NOPAGE)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE  HYBRIDE DUAL : DEFINITION D'UNE MATRICE CONSTANTE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU HD06 DANS LE SUPER TABLEAU M
      L  = 288
      IA = 0
      CALL TNMCDC( 'REEL2' , L , IA )
C     GENERATION DU TABLEAU HD06
      CALL TAHD06( DMCN((IA+1)/2) )
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('HD06',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION SUR LE LISTING
      IA1 = ( IA - 1 ) / 2
      LL  = L / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),(DMCN(IA1+1),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     TRIANGLE  EQUILIBRE : DEFINITION D'UNE MATRICE CONSTANTE
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU EQ06 DANS LE SUPER TABLEAU M
      L  = 252
      IA = 0
      CALL TNMCDC( 'REEL2' , L , IA )
C     GENERATION DU TABLEAU EQ06
      CALL TAEQ06( DMCN((IA+1)/2) )
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA('EQ06',IA,L,NFPOBA,MOPAGE,NOPAGE)
C     IMPRESSION SUR LE LISTING
      IA1= (IA-1)/2
      LL = L/2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),(DMCN(IA1+1),I=1,LL)
      CALL TNMCDS( 'REEL2' , L , IA )
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     6-CUBE  Q1-LAGRANGE  LES MONOMES DES POLYNOMES
C     -------------------  DE BASE DU 6-CUBE UNITE SONT
C                          STOCKES SUR LE FICHIER POBA
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     ADRESSAGE DU TABLEAU 'HEQ1' DANS LE SUPER-TABLEAU M
      L6CQ1  = 4 + 2 * 64 * 64
      IA6CQ1 = 0
      CALL TNMCDC( 'REEL2', L6CQ1, IA6CQ1 )
C
C     CALCUL DES 64 MONOMES DES 64 POLYNOMES DU 6-CUBE UNITE Q1
C     DIMENSION
      MCN(IA6CQ1  ) = 6
C     DEGRE + 1
      MCN(IA6CQ1+1) = 2
C     P=>0 Q=>1 R=>2   ICI Q=>1
      MCN(IA6CQ1+2) = 1
C     NOMBRE DE POLYNOMES 2**6=64
      MCN(IA6CQ1+3) = 64
C
C     CONSTRUCTION DIRECTE DES 64 COEFFICIENTS DES 64 POLYNOMES DU 6CUBE 6Q1C
      CALL CF6Q1C( MCN(IA6CQ1+4) )
C
C     TRANSFERT SUR LE FICHIER POBA
      CALL TRPOBA( '6CQ1', IA6CQ1, L6CQ1, NFPOBA, MOPAGE, NOPAGE )
C
C     L IMPRESSION SUR LE LISTING DE SES CARACTERISTIQUES ET VALEURS
      IA1 = (IA6CQ1 + 3) / 2
      LL  = (L6CQ1  - 4) / 2
      WRITE (IMPRIM,10150) (ITA(I,NBTAST),I=1,4),
     %                     (MCN(IA6CQ1-1+I),I=1,4),
     %                     (DMCN(IA1+I),I=1,LL)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     6-CUBE  Q1-LAGRANGE  LES POLYNOMES ET LEURS DERIVEES
C     -------------------  SONT CALCULES AUX 64 POINTS D UNE
C                          FORMULE D INTEGRATION NUMERIQUE
C                          EXACTE POUR LES POLYNOMES Q1
C               ATTENTION: ICI AUCUNE DERIVEE DES POLYNOMES AUX POINTS
C                          D INTEGRATION DES FACES
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL SUPOBA( '6Q13', NFPOBA, MOPAGE, NOPAGE,
     %              6,1,1,1,0,
     %              1,64,64, PO6CQ1,CO6CQ1,
     %              DMCN((IA6CQ1+1+4)/2),
     %              DPXYZ,  DDPXYZ,
     %              POLYPI, DPOLYP, DDPOLY)
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     LE TABLEAU CONTIENT : NDIM  = DIMENSION DE L ESPACE DE L EF
C     SUR LE FICHIER PXYZ   NBPOLY= NBRE DE POLYNOMES DE BASE
C                           NPI   = NBRE DE POINTS D'INTEGRATION
C                           IPOIDS= INDICATEUR DE PRESENCE DES POIDS
C                           IPOLY = INDICATEUR DE PRESENCE DES POLYNOMES
C                           IDPOLY= INDICATEUR DE PRESENCE DES DERIVEES 1
C                           IDDPOL= INDICATEUR DE PRESENCE DES DERIVEES 2
C     SI IPOIDS = 1         POIDS(NPI) = POIDS DES POINTS
C     SI IPOLY  = 1         POLY(NBPOLY,NPI) =   PI EN (XL,YL)
C                           POLY(  I   , L ) =
C     SI IDPOLY = 1         DPOLYP(NDIM,NBPOLY,NPI)=DPJ/DXI (XL,YL)
C                           DPOLYP( I  , J    , L )=
C     SI IDDPOL = 1         DDPOLY(NDIM2,NBPOLY,NPI)=DDPJ/DXXI (XL,YL)
C                           DDPOLY( 1, J, L )=DDPJ/DDXX (XL,YL)
C                           DDPOLY( 2, J, L )=DDPJ/DDYX (XL,YL)
C                           DDPOLY( 3, J, L )=DDPJ/DDYY (XL,YL)
C                           DDPOLY( 4, J, L )=DDPJ/DDZX (XL,YL,ZL)
C                           DDPOLY( 5, J, L )=DDPJ/DDZY (XL,YL,ZL)
C                           DDPOLY( 6, J, L )=DDPJ/DDZZ (XL,YL,ZL)
C     AUCUNE DERIVEE SECONDE EN DIMENSION 6 (NON PROGRAMME)
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C     *****************************************************************
C     *  INSERER ICI LES TABLEAUX DES NOUVEAUX TYPES D'ELEMENTS FINIS *
C     *****************************************************************
      NODEPA = NOPAGE + 1
C
C     MISE A ZERO DE L ADRESSE DANS M DES TABLEAUX
C     MISE A JOUR DE LA 1-ERE PAGE
C
      DO 9000 I=1,NBTAST
            ITA(2,I) = 0
 9000 ENDDO
C
      I = 1
      CALL FIPAEC( NFPOBA , I , MOPAGE , NELFI )
      WRITE (IMPRIM,9999) NELFI
C
C     FERMETURE DU FICHIER NFPOBA
C     ---------------------------
      CLOSE( NFPOBA )
      WRITE(IMPRIM,19997)
      STOP
C
C     ERREUR A L'OUVERTURE
9900  WRITE(IMPRIM,*) ' ERREUR A L''OUVERTURE DU FICHIER',KNOM
      STOP
      END


      SUBROUTINE CL3Q2C( IAHEQ2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFECTUER LES COMBINAISONS LINEAIRES ENTRE LES POLYNOMES DE
C ----- BASE DE L ELEMENT HEXAEDRE Q2-LAGRANGE ISOPARAMETRIQUE COMPLET
C       AFIN DE SUPPRIMER LE BARYCENTRE DES 6 FACES QUADRANGULAIRES
C       ET LE BARYCENTRE DU CUBE
C
C PARAMETRE D ENTREE :
C --------------------
C IAHEQ2 : ADRESSE DANS M DU TABLEAU DES COEFFICIENTS DES POLYNOMES
C          DE BASE DE L HEXAEDRE Q2-LAGRANGE COMPLET
C          EN SORTIE : ADRESSE DES POLYNOMES DE Q2-INCOMPLET (SERENDIP)
C
C PARAMETRE D ENTREE ET RESULTAT :
C --------------------------------
C MCN    : SUPER-TABLEAU DE TRAVAIL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA   AOUT 1981
C.......................................................................
      include"./incl/pppoba.inc"
      COMMON MCN(MOTMCN)
C
C     SUPPRESSION DU BARYCENTRE DES FACES CARRE PAR COMBINAISON
C     LINEAIRE DES POLYNOMES
      IA  = IAHEQ2 - 50
      IA1 = IAHEQ2 +  4
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),-.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+21*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+23*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+25*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+26*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
C
      IA1 = IA1 + 54
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+22*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.50D0,MCN(IA+24*54),MCN(IA1))
      CALL CL2VED(27,1D0,MCN(IA1),+.25D0,MCN(IA+27*54),MCN(IA1))
      RETURN
      END


      SUBROUTINE CL3R2C( IAPER2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFECTUER LES COMBINAISONS LINEAIRES ENTRE LES POLYNOMES DE
C ----- BASE DE L ELEMENT PENTAEDRE R2-LAGRANGE ISOPARAMETRIQUE COMPLET
C       AFIN DE SUPPRIMER LE BARYCENTRE DES 3 FACES QUQDRANGULAIRES
C
C PARAMETRE D ENTREE :
C --------------------
C IAPER2 : ADRESSE DANS M DU TABLEAU DES COEFFICIENTS DES POLYNOMES
C          DE BASE DE L HEXAEDRE Q2-LAGRANGE COMPLET
C          EN SORTIE : ADRESSE DES POLYNOMES DE Q2-INCOMPLET (SERENDIP)
C
C PARAMETRE D ENTREE ET RESULTAT :
C --------------------------------
C MCN    : SUPER-TABLEAU DE TRAVAIL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA   AOUT 1981
C.......................................................................
       include"./incl/pppoba.inc"
      COMMON     MCN( MOTMCN )
C
C     ADRESSAGE D UN TABLEAU AUXILIAIRE DE 3 POLYNOMES
      IAAUX = 0
      LAUX  = 162
      CALL TNMCDC( 'REEL2' , LAUX , IAAUX )
C
C     SUPPRESSION DU BARYCENTRE DES FACES CARRE PAR COMBINAISON
C     LINEAIRE DES POLYNOMES
      IA    = IAPER2 - 50
      CALL CL2VED(27,1D0,MCN(IA+16*54),1D0,MCN(IA+17*54),MCN(IAAUX))
      CALL CL2VED(27,1D0,MCN(IA+17*54),1D0,MCN(IA+18*54),MCN(IAAUX+54))
      CALL CL2VED(27,1D0,MCN(IA+18*54),1D0,MCN(IA+16*54),MCN(IAAUX+108))
C
      CALL CL2VED(27,1D0,MCN(IA+   54),-.25D0,MCN(IAAUX    ),
     %                   MCN(IA+   54))
      CALL CL2VED(27,1D0,MCN(IA+ 2*54),-.25D0,MCN(IAAUX+ 54),
     %                   MCN(IA+ 2*54))
      CALL CL2VED(27,1D0,MCN(IA+ 3*54),-.25D0,MCN(IAAUX+108),
     %                   MCN(IA+ 3*54))
      CALL CL2VED(27,1D0,MCN(IA+ 4*54),-.25D0,MCN(IAAUX    ),
     %                   MCN(IA+ 4*54))
      CALL CL2VED(27,1D0,MCN(IA+ 5*54),-.25D0,MCN(IAAUX+ 54),
     %                   MCN(IA+ 5*54))
      CALL CL2VED(27,1D0,MCN(IA+ 6*54),-.25D0,MCN(IAAUX+108),
     %                   MCN(IA+ 6*54))
      CALL CL2VED(27,1D0,MCN(IA+ 7*54),+0.5D0,MCN(IA+ 17*54),
     %                   MCN(IA+ 7*54))
      CALL CL2VED(27,1D0,MCN(IA+ 8*54),+0.5D0,MCN(IA+ 18*54),
     %                   MCN(IA+ 8*54))
      CALL CL2VED(27,1D0,MCN(IA+ 9*54),+0.5D0,MCN(IA+ 16*54),
     %                   MCN(IA+ 9*54))
      CALL CL2VED(27,1D0,MCN(IA+10*54),+0.5D0,MCN(IAAUX),MCN(IA+10*54))
      CALL CL2VED(27,1D0,MCN(IA+11*54),+0.5D0,MCN(IAAUX+ 54),
     %                   MCN(IA+11*54))
      CALL CL2VED(27,1D0,MCN(IA+12*54),+0.5D0,MCN(IAAUX+108),
     %                   MCN(IA+12*54))
      CALL CL2VED(27,1D0,MCN(IA+13*54),+0.5D0,MCN(IA+ 17*54),
     %                   MCN(IA+13*54))
      CALL CL2VED(27,1D0,MCN(IA+14*54),+0.5D0,MCN(IA+ 18*54),
     %                   MCN(IA+14*54))
      CALL CL2VED(27,1D0,MCN(IA+15*54),+0.5D0,MCN(IA+ 16*54),
     %                   MCN(IA+15*54))
C
      CALL TNMCDS( 'REEL2' , LAUX , IAAUX )
      END



      SUBROUTINE DESOM1( NCOGEL,NBPOLY,K,N1,POLY,DP,
     %                   PSOMT,DPSOMT,NVADPF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL ET STOCKAGE DU TABLEAU DES VALEURS DES
C ----- POLYNOMES AUX SOMMETS DE L EF REFERENCE
C       DERIVEE DES POLYNOMES AUX SOMMETS DE L'EF DE REFERENCE
C
C PARAMETRES D ENTREE :
C ---------------------
C NCOGEL : CODE DE GEOMETRIE DE L ELEMENT (3 TRIANGLE,...)
C NBPOLY : NOMBRE DE POLYNOMES DE L'EF
C K      : NO DU SOMMET TRAITE ( 1 OU 2 )
C NPIS   : NOMBRE DE POINTS D INTEGRATION AU SOMMET => 1
C N1     : DEGRE DES POLYNOMES + 1
C POLY   : POLY(I,K) = COEFFICIENT DE X**(I-1) DU K-EME POLYNOME
C DP     : TABLEAU AUXILIAIRE(N1) REEL DOUBLE PRECISION
C
C PARAMETRE RESULTAT :
C --------------------
C PSOMT  : VALEURS DES POLYNOMES  AUX 2 SOMMETS DE L'EF DE REFERENCE
C DPSOMT : DERIVEES DES POLYNOMES AUX 2 SOMMETS DE L'EF DE REFERENCE
C NVADPF : NOMBRE DE VARIABLES DE PSOMT ET DPSOMT POUR CE SOMMET K
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     JUIN 2009
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PSOMT(NBPOLY),
     %                  DPSOMT(NBPOLY),
     %                  S,
     %                  POLY(N1,NBPOLY),
     %                  DP(*),
     %                  X
C
C     POUR LE SOMMET K
      IF( K .EQ. 1 ) THEN
         X = 0D0
      ELSE
         X = 1D0
      ENDIF
C
C     BOUCLE SUR LES POLYNOMES
C     ========================
      DO 10 J=1,NBPOLY
C
C        VALEUR AU SOMMET K
         CALL PN1DVA( N1, POLY(1,J), X, PSOMT(J) )
C
C        VALEUR DE LA DERIVEE AU SOMMET K
         CALL PN1DDE( N1, POLY(1,J), DP )
         CALL PN1DVA( N1, DP, X, DPSOMT(J) )
C
 10   ENDDO
C
      NVADPF = 2 * NBPOLY
      RETURN
      END


      SUBROUTINE DERSOM( NOMTAB, NCOGEL, NBPOLY, NDEGRE, POLY, DP,
     %                   NFPOBA, MOPAGE, NOPAGE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL ET STOCKAGE DU TABLEAU DE LA DERIVEE DES POLYNOMES AUX
C -----  2 SOMMETS DU SEGMENT UNITE  EF DE REFERENCE 1D
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMTAB : 4 CARACTERES DU NOM DU TABLEAU A STOCKER
C NCOGEL : CODE DE GEOMETRIE DE L ELEMENT (2:SEGMENT)
C NBPOLY : NOMBRE DE POLYNOMES
C NDEGRE : DEGRE DES POLYNOMES
C POLY   : POLY(I,K) = COEFFICIENT DE X**(I-1) DU K-EME POLYNOME
C DP     : TABLEAU AUXILIAIRE(NDEGRE+1) REEL DOUBLE PRECISION
C NPIS   : NOMBRE DE POINTS D INTEGRATION SUR LE SOMMET => 1
C COSEDE : COORDONNEES DES POINTS D INTEGRATION DES SOMMETS
C NFPOBA : FICHIER DE STOCKAGE DES TABLEAUX
C MOPAGE : NOMBRE DE MOTS DE CHAQUE PAGE D ACCES DIRECT DU FICHIER POBA
C
C PARAMETRE MODIFIE :
C -------------------
C NOPAGE : NO DE LA PAGE QUI PRECEDE CELLE DE STOCKAGE DU TABLEAU
C          EN SORTIE NO DE LA DERNIERE PAGE SUPPORT DU TABLEAU SUR POBA
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     JUIN 2009
C23456---------------------------------------------------------------012
      include"./incl/pppoba.inc"
      COMMON             MCN( MOTMCN )
      DOUBLE  PRECISION  DMCN(1)
      EQUIVALENCE       (MCN(1),DMCN(1))
      CHARACTER*(*)      NOMTAB
C
      DOUBLE  PRECISION POLY(1), DP(1)
      COMMON / ELFINI / NBTAST,NODEPA,ITA(4,127),NZ(2)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
10000 FORMAT(/' TABLEAU ',A4,' ADRESSE DANS M ',I8,3X,' NBRE DE MOTS ',
     %I7,3X,' NO DE SA 1-ERE PAGE ',I7/1X,130(1H=)/
     %' ADRESSE DE LA 1-ERE VARIABLE DE CHAQUE SOMMET =',6I5/
     %(10G10.3))
C
C     INITIALISATIONS
C     ===============
      N1 = NDEGRE + 1
C
C     ADRESSAGE DU TABLEAU NOMTAB
C     ===========================
      IAT = 0
      LTT = 6 + 12 * NBPOLY
      CALL TNMCDC( 'REEL2', LTT, IAT )
C
C     NOMBRE DE SOMMETS DE L ELEMENT FINI
C     ===================================
      NBS = 2
C
C     BOUCLE SUR LES SOMMETS DE L ELEMENT FINI
C     ========================================
      IA0 = (IAT + 6) / 2
      IA  = 1
      IAD = 1 + NBPOLY
      DO 100 K=1,NBS
C          CALCUL DE DP AU SOMMET K
           CALL DESOM1( NCOGEL, NBPOLY, K, N1,
     %                  POLY, DP, DMCN(IA0+IA), DMCN(IA0+IAD), NVADPF )
C          MISE A JOUR DE IA
           MCN(IAT-1+K) = IA
           IA  = IA  + NVADPF
           IAD = IAD + NVADPF
  100 ENDDO
C
C     MISE A JOUR DU TABLEAU POINTEUR DE NBS+1 A 6
C     ============================================
      DO 110 K = NBS, 5
           MCN( IAT + K ) = 0
  110 ENDDO
C
C     TRANSFERT SUR LE FICHIER
C     ========================
      IA = IA - 1
      K  = 6 + IA * 2
      CALL TRPOBA( NOMTAB, IAT, K, NFPOBA, MOPAGE, NOPAGE )
C
C     IMPRESSION DES VALEURS STOCKEES
C     ===============================
      WRITE( IMPRIM, 10000 ) (ITA(K,NBTAST),K=1,4),
     %                       (MCN(IAT-1+K),K=1,6),
     %                       (DMCN(IA0+K),K=1,IA)
C
C     DESTRUCTION DU TABLEAU
C     ======================
      CALL TNMCDS( 'REEL2', LTT, IAT )
      RETURN
      END


      SUBROUTINE DERAR1( NCOGEL,NBPOLY,K,NPIA,COSEDE,N1,N2,POLY,DP,
     %                   PARET,DPARET,NVADPF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL ET STOCKAGE DU TABLEAU DES VALEURS DES
C ----- POLYNOMES AUX POINTS D INTEGRATION D UNE ARETE DE L EF REFERENCE
C       DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION D UNE ARETE DE
C       L ELEMENT FINI DE REFERENCE
C
C PARAMETRES D ENTREE :
C ---------------------
C NCOGEL : CODE DE GEOMETRIE DE L ELEMENT (3 TRIANGLE,...)
C NBPOLY : NOMBRE DE POLYNOMES
C K      : NO DE L ARETE TRAITEE
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR L ARETE
C COSEDE : COORDONNEES DES POINTS D INTEGRATION DE L ARETE
C N1     : DEGRE DES POLYNOMES + 1
C N2     : N1 ** 2
C POLY   : POLY(I,K) = COEFFICIENT DE X**(I-1) DU K-EME
C          POLYNOME.LE 1-ER INDICE CONDENSE LES 2 CI-DESSUS
C DP     : TABLEAU AUXILIAIRE(NDEGRE+1,NDEGRE+1) REEL DOUBLE PRECISION
C
C PARAMETRE RESULTAT :
C --------------------
C PARET  : VALEURS DES POLYNOMES AUX POINTS D INTEGRATION DE LA
C          K-EME ARETE DE L'EF
C DPARET : DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION DE LA
C          K-EME ARETE DE L'EF
C NVADPF : NOMBRE DE VARIABLES DE PARET ET DPARET POUR CETTE ARETE K
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  PARET(NBPOLY,NPIA),
     %                  DPARET(2,NBPOLY,NPIA),
     %                  COSEDE(NPIA),S,
     %                  POLY(N2,NBPOLY),
     %                  DP(*),
     %                  X(2),XX(2),ALFA(2,2)
C
      DO 100 L=1,NPIA
C
C          RECHERCHE DES COEFFICIENTS DE PARAMETRAGE DE G
C          PAR RAPPORT A F (G = F RESTREINTE A L ARETE)
C          ==============================================
           CALL GFACE(2,NCOGEL,K,ALFA)
C
C          CALCUL DES COORDONNEES DES POINTS SUR L ARETE K
C          ===============================================
           XX(1) = 1D0
           XX(2) = COSEDE(L)
           DO 20 I=1,2
                S = 0D0
                DO 10 J=1,2
                     S = S + ALFA(I,J) * XX(J)
   10           ENDDO
                X(I) = S
   20      ENDDO
C
C          BOUCLE SUR LES POLYNOMES
C          ========================
           DO 40 J=1,NBPOLY
C
C               VALEUR AUX POINTS
                CALL PN2DVA( N1, POLY(1,J), X(1), X(2), PARET(J,L) )
C               BOUCLE SUR LA DERIVATION
                DO 30 I=1,2
                     CALL PN2DDE( I, N1, POLY(1,J), DP )
                     CALL PN2DVA( N1, DP, X(1), X(2), DPARET(I,J,L) )
   30           ENDDO
   40      ENDDO
  100 ENDDO
C
      NVADPF = 3 * NBPOLY * NPIA
      END


      SUBROUTINE DERARE( NOMTAB, NCOGEL, NBPOLY, NDEGRE, POLY, DP,
     %                   NPIA,   COSEDE, NFPOBA, MOPAGE, NOPAGE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL ET STOCKAGE DU TABLEAU DES DERIVEES DES POLYNOMES AUX
C ----- POINTS D INTEGRATION DE CHAQUE ARETE DE L ELEMENT DE REFERENCE
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMTAB : 4 CARACTERES DU NOM DU TABLEAU A STOCKER
C NCOGEL : CODE DE GEOMETRIE DE L ELEMENT (3 TRIANGLE,...)
C NBPOLY : NOMBRE DE POLYNOMES
C NDEGRE : DEGRE DES POLYNOMES
C POLY   : POLY(I,K) = COEFFICIENT DE X**(I-1) DU K-EME POLYNOME
C DP     : TABLEAU AUXILIAIRE(NDEGRE+1,NDEGRE+1) REEL DOUBLE PRECISION
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR LE SEGMENT UNITE
C COSEDE : COORDONNEES DES POINTS D INTEGRATION DU SEGMENT UNITE
C NFPOBA : FICHIER DE STOCKAGE DES TABLEAUX
C MOPAGE : NOMBRE DE MOTS DE CHAQUE PAGE D ACCES DIRECT DU FICHIER POBA
C
C PARAMETRE MODIFIE :
C -------------------
C NOPAGE : NO DE LA PAGE QUI PRECEDE CELLE DE STOCKAGE DU TABLEAU
C          EN SORTIE NO DE LA DERNIERE PAGE SUPPORT DU TABLEAU SUR POBA
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C ......................................................................
      include"./incl/pppoba.inc"
      COMMON             MCN( MOTMCN )
      DOUBLE  PRECISION  DMCN(1)
      EQUIVALENCE       (MCN(1),DMCN(1))
      CHARACTER*(*)      NOMTAB
C
      DOUBLE  PRECISION POLY(1),DP(1),COSEDE(1)
      COMMON / ELFINI / NBTAST,NODEPA,ITA(4,127),NZ(2)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
10000 FORMAT(/' TABLEAU ',A4,' ADRESSE DANS M ',I8,3X,' NBRE DE MOTS ',
     %I7,3X,' NO DE SA 1-ERE PAGE ',I7/1X,130(1H=)/
     %' ADRESSE DE LA 1-ERE VARIABLE DE CHAQUE ARETE =',6I5/
     %(10G10.3))
C
C     INITIALISATIONS
C     ===============
      N1 = NDEGRE + 1
      N2 = N1 ** 2
C
C     ADRESSAGE DU TABLEAU NOMTAB
C     ===========================
      IAT = 0
      LTT = 6 + 12 * NBPOLY * NPIA
      CALL TNMCDC( 'REEL2' , LTT , IAT )
C
C     RECHERCHE DU NOMBRE DE ARETES DE L ELEMENT FINI
C     ===============================================
      NBARE = NBARET( NCOGEL )
C
C     BOUCLE SUR LES ARETES DE L ELEMENT FINI
C     =======================================
      IA0 = (IAT + 6) / 2
      IA  = 1
      IAD = 1 + NBPOLY * NPIA
      DO 100 K=1,NBARE
C          CALCUL DE DP AUX POINTS D INTEGRATION DE LA K-EME ARETE
           CALL DERAR1(NCOGEL,NBPOLY,K,NPIA,COSEDE,N1,N2,
     %                 POLY,DP,DMCN(IA0+IA),DMCN(IA0+IAD),NVADPF)
C          MISE A JOUR DE IA
           MCN(IAT-1+K) = IA
           IA  = IA  + NVADPF
           IAD = IAD + NVADPF
  100 ENDDO
C
C     MISE A JOUR DU TABLEAU POINTEUR DE NBARE+1 A 6
C     ==============================================
      DO 110 K=NBARE , 5
           MCN( IAT + K ) = 0
  110 ENDDO
C
C     TRANSFERT SUR LE FICHIER
C     ========================
      IA = IA - 1
      K  = 6 + IA * 2
      CALL TRPOBA(NOMTAB,IAT,K,NFPOBA,MOPAGE,NOPAGE)
C
C     IMPRESSION DES VALEURS STOCKEES
C     ===============================
      WRITE( IMPRIM , 10000 ) (ITA(K,NBTAST),K=1,4),
     %                        (MCN(IAT-1+K),K=1,6),
     %                        (DMCN(IA0+K),K=1,IA)
C
C     DESTRUCTION DU TABLEAU
C     ======================
      CALL TNMCDS( 'REEL2' , LTT , IAT )
      RETURN
      END


      SUBROUTINE DERFA1(NCOGEL,NBPOLY,K,NPIF,COFADE,N1,N3,POLY,DP,
     %                  PFACE,DPFACE,NVADPF)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL ET STOCKAGE DU TABLEAU DES VALEURS DES
C ----- POLYNOMES AUX POINTS D INTEGRATION D UNE FACE DE L'EF
C       DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION D UNE FACE DE
C       L'ELEMENT FINI DE REFERENCE
C
C PARAMETRES D ENTREE :
C ---------------------
C NCOGEL : CODE DE GEOMETRIE DE L ELEMENT (3 TRIANGLE,...)
C NBPOLY : NOMBRE DE POLYNOMES
C K      : NO DE LA FACE TRAITEE
C NPIF   : NOMBRE DE POINTS D INTEGRATION SUR LA FACE
C COFADE : COORDONNEES DES POINTS D INTEGRATION DE LA FACE
C N1     : DEGRE DES POLYNOMES + 1
C N3     : N1 ** 3
C POLY   : POLY(I,J,K) = COEFFICIENT DE X**(I-1) YY**(J-1) DU K-EME
C          POLYNOMES LE 1-ER INDICE CONDENSE LES 2 CI-DESSUS
C DP     : TABLEAU AUXILIAIRE(NDEGRE+1,NDEGRE+1) REEL DOUBLE PRECISION
C
C PARAMETRE RESULTAT :
C --------------------
C PFACE  : VALEURS DES POLYNOMES AUX POINTS D INTEGRATION DE LA
C          K-EME FACE
C DPFACE : DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION DE LA
C          K-EME FACE
C NVADPF : NOMBRE DE VARIABLES DE PFACE + DPFACE POUR CETTE FACE K
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C ......................................................................
      DOUBLE PRECISION  PFACE(NBPOLY,NPIF),
     %                  DPFACE(3,NBPOLY,NPIF),COFADE(2,NPIF),S,
     %                  X(3),XX(3),ALFA(3,3),POLY(N3,NBPOLY),DP(1)
C
      DO 100 L=1,NPIF
C
C        RECHERCHE DES COEFFICIENTS DE PARAMETRAGE DE G
C        PAR RAPPORT A F (G = F RESTREINTE A LA FACE)
C        ==============================================
         CALL GFACE(3,NCOGEL,K,ALFA)
C
C        CALCUL DES COORDONNEES DES POINTS SUR LA FACE K
C        ===============================================
         XX(1) = 1D0
         XX(2) = COFADE(1,L)
         XX(3) = COFADE(2,L)
         DO 20 I=1,3
            S = 0D0
            DO 10 J=1,3
               S = S + ALFA(I,J) * XX(J)
   10       ENDDO
            X(I) = S
   20    ENDDO
C
C        BOUCLE SUR LES POLYNOMES
C        ========================
         DO 40 J=1,NBPOLY
C
C           VALEURS DES POLYNOMES
            CALL PN3DVA(N1,POLY(1,J),X(1),X(2),X(3),PFACE(J,L))
C
C           BOUCLE SUR LA DERIVATION
            DO 30 I=1,3
               CALL PN3DDE(I,N1,POLY(1,J),DP)
               CALL PN3DVA(N1,DP,X(1),X(2),X(3),DPFACE(I,J,L))
   30       ENDDO
   40    ENDDO
  100 ENDDO

      NVADPF = 4 * NBPOLY * NPIF
      RETURN
      END


      SUBROUTINE DERFAC( NOMTAB,NCOGEL,NBPOLY,NDEGRE,POLY,DP,
     %                   NPIT,COTRDE,NPIQ,COQUDE,NFPOBA,MOPAGE,NOPAGE)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL ET STOCKAGE DU TABLEAU DES DERIVEES DES POLYNOMES AUX
C ----- POINTS D INTEGRATION DE CHAQUE FACE DE L ELEMENT DE REFERENCE
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMTAB : 4 CARACTERES DU NOM DU TABLEAU A STOCKER
C NCOGEL : CODE DE GEOMETRIE DE L ELEMENT (3 TRIANGLE,...)
C NBPOLY : NOMBRE DE POLYNOMES
C NDEGRE : DEGRE DES POLYNOMES
C POLY   : POLY(I,J,K) = COEFFICIENT DE X**(I-1) YY**(J-1) DU K-EME POLY
C DP     : TABLEAU AUXILIAIRE(NDEGRE+1,NDEGRE+1) REEL DOUBLE PRECISION
C NPIT   : NOMBRE DE POINTS D INTEGRATION SUR LE   TRIANGLE UNITE
C COTRDE : COORDONNEES DES POINTS D INTEGRATION DU TRIANGLE UNITE
C NPIQ   : NOMBRE DE POINTS D INTEGRATION SUR LE   QUADRANGLE UNITE
C COQUDE : COORDONNEES DES POINTS D INTEGRATION DU QUADRANGLE UNITE
C NFPOBA : FICHIER DE STOCKAGE DES TABLEAUX
C MOPAGE  : NOMBRE DE MOTS DE CHAQUE PAGE D ACCES DIRECT DU FICHIER POBA
C
C PARAMETRE MODIFIE :
C -------------------
C NOPAGE : NO DE LA PAGE QUI PRECEDE CELLE DE STOCKAGE DU TABLEAU
C          EN SORTIE NO DE LA DERNIERE PAGE SUPPORT DU TABLEAU SUR POBA
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C ......................................................................
      include"./incl/pppoba.inc"
      COMMON             MCN( MOTMCN )
      DOUBLE PRECISION  DMCN( MOTMCN/2 )
      EQUIVALENCE      (DMCN(1),MCN(1))
      CHARACTER*(*)     NOMTAB
C
      DOUBLE PRECISION  POLY(1),DP(1),COTRDE(1),COQUDE(1)
      COMMON / ELFINI / NBTAST,NODEPA,ITA(4,127),NZ(2)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
10000 FORMAT(/' TABLEAU ',A4,' ADRESSE DANS M ',I8,3X,' NBRE DE MOTS ',
     %I7,3X,' NO DE SA 1-ERE PAGE ',I7/1X,130(1H=)/
     %' ADRESSE DE LA 1-ERE VARIABLE DE CHAQUE FACE =',6I5/
     %(10G10.3))
C
C     INITIALISATIONS
C     ===============
      N1 = NDEGRE + 1
      N3 = N1 ** 3
C
C     ADRESSAGE DU TABLEAU NOMTAB
C     ===========================
      IAT   = 0
      NPITQ = MAX(NPIT,NPIQ)
      LTT   = 6 + 24 * NBPOLY * NPITQ
      CALL TNMCDC( 'REEL2' , LTT , IAT )
C
C     RECHERCHE DU NOMBRE DE FACES DE L ELEMENT
C     =========================================
      NBFAC = NBFACE( NCOGEL )
C
C     BOUCLE SUR LES FACES DE L ELEMENT
C     =================================
      IA0 = (IAT + 6) / 2
      IA  = 1
      DO 100 K=1,NBFAC
C
C          RECHERCHE DU TYPE DE LA FACE (TRIANGLE OU QUADRANGLE)
C          -----------------------------------------------------
           NSOFA = NBSOFA( NCOGEL , K )
           IF( NSOFA .EQ. 3 ) THEN
C
C             TRIANGLE
              IAD = IA + NBPOLY * NPIT
              CALL DERFA1( NCOGEL,NBPOLY,K,NPIT,COTRDE,N1,N3,
     %                     POLY,DP,DMCN(IA0+IA),DMCN(IA0+IAD),NVADPF )
           ELSE
C
C             QUADRANGLE
              IAD = IA + NBPOLY * NPIQ
              CALL DERFA1( NCOGEL,NBPOLY,K,NPIQ,COQUDE,N1,N3,
     %                     POLY,DP,DMCN(IA0+IA),DMCN(IA0+IAD),NVADPF )
           ENDIF
C
C          MISE A JOUR DE IA
           MCN(IAT-1+K) = IA
           IA           = IA + NVADPF
  100 ENDDO
C
C     MISE A JOUR DU TABLEAU POINTEUR DE NBFAC+1 A 6
C     ==============================================
      IF( NBFAC .EQ. 6 ) GOTO 120
      DO 110 K=NBFAC , 5
           MCN( IAT + K ) = 0
  110 ENDDO
C
C     TRANSFERT SUR LE FICHIER
C     ========================
  120 IA = IA - 1
      K  = 6 + IA * 2
      CALL TRPOBA( NOMTAB,IAT,K,NFPOBA,MOPAGE,NOPAGE )
C
C     IMPRESSION DES VALEURS STOCKEES
C     ===============================
      WRITE(IMPRIM,10000) (ITA(K,NBTAST),K=1,4),(MCN(IAT-1+K),K=1,6),
     %                    (DMCN(IA0+K),K=1,IA)
C
C     DESTRUCTION DU TABLEAU
C     ======================
      CALL TNMCDS( 'REEL2' , LTT , IAT )
      RETURN
      END


      SUBROUTINE SUPOBA( NOMTAB, NFPOBA, MOPAGE, NOPAGE, NDIM,
     %                   IPOIDS, IPOLY,  IDPOLY, IDDPOL,
     %                   NDEGRE, NBPOLY, NPI,    POIDS,  COORPI, POLY,
     %                   DPXYZ,  DDPXYZ,
     %                   POLYPI, DPOLYP, DDPOLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER PUIS PORTER SUR LE FICHIER POBA LE TABLEAU DES POIDS,
C ----- DES VALEURS DES POLYNOMES ET DE LEURS DERIVEES AUX
C       POINTS D INTEGRATION NUMERIQUE
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMTAB : NOM DU TABLEAU A PORTER SUR POBA (4 CARACTERES)
C NFPOBA : NO DU FICHIER D ACCES DIRECT POBA
C MOPAGE : NOMBRE DE MOTS D UNE PAGE DU FICHIER NFPOBA
C NOPAGE : NO DE LA PAGE DU FICHIER POBA DU TABLEAU AVANT CELUI A AJOUTE
C NDIM   : DIMENSION DE L ESPACE DE L ELEMENT ( 1 OU 2 OU 3 )
C IPOIDS : INDICATEUR DE STOCKAGE(1) OU NON(0) DU TABLEAU POIDS SUR POBA
C IPOLY  : INDICATEUR DE STOCKAGE(1) OU NON(0) DU TABLEAU POLYPI
C IDPOLY : INDICATEUR DE STOCKAGE(1) OU NON(0) DU TABLEAU DPOLYP
C IDDPOL : INDICATEUR DE STOCKAGE(1) OU NON(0) DU TABLEAU DDPOLY
C NDEGRE : DEGRE DES POLYNOMES
C NBPOLY : NOMBRE DE POLYNOMES DE L ELEMENT FINI
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE
C POIDS  : POIDS DES POINTS DE LA FORMULE D INTEGRATION NUMERIQUE
C COORPI : COORDONNEES DES POINTS D INTEGRATION
C POLY   : COEFFICIENTS DES POLYNOMES
C          POLY(I,J,K) COEFFICIENT DE X**(I-1) , Y**(J-1) , Z**(K-1)
C          K EXISTE SEULEMENT SI NDIM = 3 , J EXISTE ...
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C DPXYZ  : TABLEAU DES COEFFICIENTS DES POLYNOMES DERIVES
C          DE NDIM  * NBPOLY VARIABLES REELLES DOUBLE PRECISION
C DDPXYZ : TABLEAU DES COEFFICIENTS DES POLYNOMES DERIVES 2 FOIS
C          DE NDIM2 * NBPOLY VARIABLES REELLES DOUBLE PRECISION
C
C PARAMETRES RESULTATS :
C ----------------------
C POLYPI : VALEURS DES POLYNOMES AUX NPI POINTS D INTEGRATION
C DPOLYP : VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION
C DDPOLY : VALEURS DES DERIVEES SECONDES DES POLYNOMES AUX POINTS
C          D INTEGRATION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS       AVRIL 1995
C2345X7..............................................................012
      include"./incl/pppoba.inc"
      COMMON             MCN( MOTMCN )
      DOUBLE PRECISION  DMCN( MOTMCN/2 )
      EQUIVALENCE      (DMCN(1),MCN(1))
      DOUBLE PRECISION  POIDS(NPI),COORPI(NDIM,NPI),POLY(*),
     %                  DPXYZ(*), DDPXYZ(*),
     %                  POLYPI(*),DPOLYP(*),DDPOLY(*)
      CHARACTER*(*)     NOMTAB
      INTEGER           ITA(4,127)
      COMMON / ELFINI / NELFI(512)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      EQUIVALENCE       (NELFI(1),NBTAST) , (NELFI(2),NODEPA) ,
     %                  (NELFI(3),ITA(1,1))
10100 FORMAT(/' TABLEAU ',A4,' ADRESSE DANS M ',I9,3X,' NBRE DE MOTS ',
     %I7,3X,' NO DE SA 1-ERE PAGE ',I6/1X,130(1H=)/
     %' DIMENSION=',I3,3X,' NB POLYNOMES=',I3,3X,' NB POINTS INTEGRATION
     %=',I3,3X/' INDICATEUR DE PRESENCE DES POIDS=',I3,3X,' DES VALEURS
     %DES POLYNOMES=',I3,3X,' DES DERIVEES=',I3,3X,
     %' DES DERIVEES SECONDES=',I3,3X/)
10110 FORMAT(' LES',I5,' POIDS'/1X,130(1H-)/(10G13.5))
10114 FORMAT(/' LES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION'/
     %1X,130(1H-))
10115 FORMAT(' POINT D INTEGRATION',I6/(10G13.5))
10120 FORMAT(/' LES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D INTE
     %GRATION'/1X,130(1H-))
10130 FORMAT(/' LES VALEURS DES DERIVEES SECONDES DES POLYNOMES AUX POIN
     %TS D INTEGRATION'/1X,130(1H-))
C
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     LE TABLEAU CONTIENT : NDIM  = DIMENSION DE L ESPACE DE L EF
C     SUR LE FICHIER PXYZ   NBPOLY= NBRE DE POLYNOMES DE BASE
C                           NPI   = NBRE DE POINTS D'INTEGRATION
C                           IPOIDS= INDICATEUR DE PRESENCE DES POIDS
C                           IPOLY = INDICATEUR DE PRESENCE DES POLYNOMES
C                           IDPOLY= INDICATEUR DE PRESENCE DES DERIVEES 1
C                           IDDPOL= INDICATEUR DE PRESENCE DES DERIVEES 2
C     SI IPOIDS = 1         POIDS(NPI) = POIDS DES POINTS
C     SI IPOLY  = 1         POLY(NBPOLY,NPI) =   PI EN (XL,YL)
C                           POLY(  I   , L ) =
C     SI IDPOLY = 1         DPOLYP(NDIM,NBPOLY,NPI)=DPJ/DXI (XL,YL)
C                           DPOLYP( I  , J    , L )=
C     SI IDDPOL = 1         DDPOLY(NDIM2,NBPOLY,NPI)=DDPJ/DXXI (XL,YL)
C                           DDPOLY( 1, J, L )=DDPJ/DDXX (XL,YL)
C                           DDPOLY( 2, J, L )=DDPJ/DDYX (XL,YL)
C                           DDPOLY( 3, J, L )=DDPJ/DDYY (XL,YL)
C                           DDPOLY( 4, J, L )=DDPJ/DDZX (XL,YL,ZL)
C                           DDPOLY( 5, J, L )=DDPJ/DDZY (XL,YL,ZL)
C                           DDPOLY( 6, J, L )=DDPJ/DDZZ (XL,YL,ZL)
C     AUCUNE DERIVEE SECONDE EN DIMENSION 6 (NON PROGRAMME)
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      L     = NDEGRE + 1
      NDIM2 = NDIM * (NDIM+1) / 2
C
      GOTO ( 10, 20, 30, 5, 5, 36 ) , NDIM
 5    WRITE(IMPRIM,*) 'ERREUR DIMENSION ',NDIM,' NON PROGRAMMEE'
      RETURN
C
 10   CALL VAPOR1( NBPOLY, NPI, L, POLY, COORPI,
     %             DPXYZ(1), DDPXYZ(1),
     %             POLYPI,   DPOLYP,   DDPOLY )
      GOTO 40
C
 20   L1 = L * L
      CALL VAPOR2( NBPOLY, NPI, L, POLY, COORPI,
     %             DPXYZ(1),  DPXYZ(1+L1),
     %             DDPXYZ(1), DDPXYZ(1+L1), DDPXYZ(1+2*L1),
     %             POLYPI,    DPOLYP,       DDPOLY )
      GOTO 40
C
 30   L1 = L * L * L
      CALL VAPOR3( NBPOLY, NPI, L, L1, POLY, COORPI,
     %             DPXYZ(1),  DPXYZ(1+L1),  DPXYZ(1+2*L1),
     %             DDPXYZ(1), DDPXYZ(1+L1), DDPXYZ(1+2*L1),
     %             DDPXYZ(1+3*L1), DDPXYZ(1+4*L1), DDPXYZ(1+5*L1),
     %             POLYPI,         DPOLYP,         DDPOLY )
      GOTO 40
C
 36   L1 = L ** 6
      CALL VAPOR6( NBPOLY, NPI, L, L1, POLY, COORPI,
     %             DPXYZ(1),       DPXYZ(1+  L1),  DPXYZ(1+2*L1),
     %             DPXYZ(1+3*L1),  DPXYZ(1+4*L1),  DPXYZ(1+5*L1),
     %             POLYPI,         DPOLYP )
C
 40   NBTAST        = NBTAST + 1
C     CODAGE DES 4 CARACTERES DU NOM DANS UN ENTIER
      ITA(1,NBTAST) = ICHARX( NOMTAB )
      L1            = 2 * ( IPOIDS *                   NPI )
      L2            = 2 * ( IPOLY  *          NBPOLY * NPI )
      L3            = 2 * ( IDPOLY * NDIM   * NBPOLY * NPI )
      L4            = 2 * ( IDDPOL * NDIM2  * NBPOLY * NPI )
C     LE 8-EME MOT EST LIBRE (NECESSAIRE POUR L'ALIGNEMENT DOUBLE PRECISION)
      L             = 8 + L1 + L2 + L3 + L4
      ITA(3,NBTAST) = L
      IA            = 0
      CALL TNMCDC( 'REEL2' , L , IA )
      ITA(2,NBTAST) = IA
C
C     INITIALISATION DU TABLEAU
C     -------------------------
      MCN(IA  ) = NDIM
      MCN(IA+1) = NBPOLY
      MCN(IA+2) = NPI
      MCN(IA+3) = IPOIDS
      MCN(IA+4) = IPOLY
      MCN(IA+5) = IDPOLY
      MCN(IA+6) = IDDPOL
C     LE 8-EME MOT EST LIBRE (NECESSAIRE POUR L'ALIGNEMENT DOUBLE PRECISION)
      MCN(IA+7) = 0
C
C     COPIE DES POIDS
C     ---------------
      IA1 = IA + 8
      IF( IPOIDS .EQ. 0 ) GOTO 50
      CALL TRTATA ( POIDS, MCN(IA1), L1 )
      IA1 = IA1 + L1
C
C     COPIE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
C     --------------------------------------------------------
   50 IF( IPOLY  .EQ. 0 ) GOTO 60
      CALL TRTATA ( POLYPI, MCN(IA1), L2 )
      IA1 = IA1 + L2
C
C     COPIE DES DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION
C     ---------------------------------------------------------
   60 IF( IDPOLY .EQ. 0 ) GOTO 70
      CALL TRTATA ( DPOLYP, MCN(IA1), L3 )
      IA1 = IA1 + L3
C
C     COPIE DES DERIVEES SECONDES DES POLYNOMES AUX POINTS D INTEGRATION
C     ------------------------------------------------------------------
   70 IF( IDDPOL .EQ. 0 ) GOTO 100
      CALL TRTATA ( DDPOLY, MCN(IA1), L4 )
C
C     MISE SUR LE FICHIER
C     -------------------
  100 ITA(4,NBTAST) = NOPAGE + 1
      CALL ESTASF(ITA(2,NBTAST),ITA(3,NBTAST),-1,NFPOBA,MOPAGE,NOPAGE)
C
C     IMPRESSION
C     ----------
      IA1 = ( IA + 7 ) / 2
      WRITE (IMPRIM,10100) (ITA(I,NBTAST),I=1,4),(MCN(IA-1+I),I=1,7)
C
      IF( IPOIDS .EQ. 0 ) GOTO 110
      WRITE(IMPRIM,10110) NPI,(DMCN(IA1+I),I=1,NPI)
      IA1 = IA1 + NPI
C
  110 IF( IPOLY .EQ. 0 ) GOTO 120
      WRITE(IMPRIM,10114)
      DO 115 L=1,NPI
         WRITE(IMPRIM,10115) L,(DMCN(IA1+I),I=1,NBPOLY)
         IA1 = IA1 + NBPOLY
  115 ENDDO
C
  120 IF( IDPOLY .EQ. 0 ) GOTO 130
      WRITE(IMPRIM,10120)
      L1 = NDIM * NBPOLY
      DO 125 L=1,NPI
         WRITE(IMPRIM,10115) L,(DMCN(IA1+I),I=1,L1)
         IA1 = IA1 + L1
  125 ENDDO
C
  130 IF( IDDPOL .EQ. 0 ) GOTO 200
      WRITE(IMPRIM,10130)
      L1 = NDIM2 * NBPOLY
      DO 135 L=1,NPI
         WRITE(IMPRIM,10115) L,(DMCN(IA1+I),I=1,L1)
         IA1 = IA1 + L1
  135 ENDDO
C
  200 CALL TNMCDS( 'REEL2' , L , IA )
      RETURN
      END


      SUBROUTINE TAEQ06( B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERE LES MATRICES CONSTANTES UTILISEES PAR
C       LES MODULES DE LA METHODE EQUILIBRE
C       D'ORDRE 2 POUR LES TRIANGLES (TREQ06,TSEQ06,TCEQ06)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : M.LEGENDRE ANALYSE NUMERIQUE PARIS    FEVRIER 1981
C ................................................................
      DOUBLE PRECISION B(126),BR(30),RG(60),DG(36)
      DATA BR/ 0.D0,-1.D0,0.D0,0.7320508D0,2.7320508D0,
     %         0.D0,-1.D0,0.D0,-2.7320508D0,-0.7320508D0,
     %         1.D0,1.D0,-0.7320508D0,2.7320508D0,0.D0,
     %         1.D0,1.D0,2.7320508D0,-0.7320508D0,0.D0,
     %        -1.D0,0.D0,-2.7320508D0,0.D0,0.7320508D0,
     %        -1.D0,0.D0,0.7320508D0,0.D0,-2.7320508D0/
      DATA RG/ 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,
     %         1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,
     %         0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0,
     %         0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0,
     %        -2.0D0, 0.0D0,-2.0D0, 0.0D0, 4.0D0, 0.0D0,
     %        -2.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,
     %         0.0D0,-2.0D0, 0.0D0, 4.0D0, 0.0D0,-2.0D0,
     %         0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,-2.0D0,
     %         4.0D0,-4.0D0,-2.0D0, 2.0D0,-2.0D0, 2.0D0,
     %         1.0D0,-1.0D0,-2.0D0, 2.0D0, 1.0D0,-1.0D0/
      DATA DG/ 0.0D0, 5.0D0, 7.0D0,-7.0D0, 0.0D0, 0.0D0,
     %         2.5D0,-1.0D0, 3.5D0,-3.5D0, 0.0D0, 4.5D0,
     %         0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0,
     %        -0.5D0, 0.0D0, 0.5D0, 0.5D0, 0.0D0,-0.5D0,
     %         5.0D0, 0.0D0, 0.0D0, 0.0D0,-7.0D0, 7.0D0,
     %         4.5D0, 0.0D0,-3.5D0, 3.5D0,-1.0D0, 2.5D0/
      DO I=1,30
         B(I)=BR(I)/2.D0
      ENDDO
      DO I=1,60
         B(I+30)=RG(I)
      ENDDO
      DO I=1,36
         B(I+90)=-DG(I)/6.D0
      ENDDO
      RETURN
      END


      SUBROUTINE TAHD06( B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERE LES MATRICES CONSTANTES UTILISEES PAR
C       LES MODULES DE LA METHODE HYBRIDE DUALE
C       D'ORDRE 2 POUR LES TRIANGLES (TRHD06,TSHD06,TCHD06)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : M.LEGENDRE ANALYSE NUMERIQUE PARIS    FEVRIER 1981
C ................................................................
      DOUBLE PRECISION B(144),BR(36),RG(72),DG(36)
      DATA BR/-0.5D0,-0.5D0, 1.0D0, 1.0D0, 0.0D0, 0.8D0,
     %         0.5D0, 0.0D0,-1.0D0, 0.0D0,-1.0D0, 0.8D0,
     %         0.0D0, 0.5D0, 0.0D0,-1.0D0, 1.0D0, 0.8D0,
     %         0.0D0,-2.0D0, 0.0D0,-2.0D0, 2.0D0,-0.8D0,
     %         2.0D0, 2.0D0, 2.0D0, 2.0D0, 0.0D0,-0.8D0,
     %        -2.0D0, 0.0D0,-2.0D0, 0.0D0,-2.0D0,-0.8D0/
      DATA RG/ 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,
     %         1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,
     %         0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0,
     %         0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0,
     %        -2.0D0, 0.0D0,-2.0D0, 0.0D0, 4.0D0, 0.0D0,
     %        -2.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,
     %         0.0D0,-2.0D0, 0.0D0, 4.0D0, 0.0D0,-2.0D0,
     %         0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0,-2.0D0,
     %         4.0D0,-4.0D0,-2.0D0, 2.0D0,-2.0D0, 2.0D0,
     %         1.0D0,-1.0D0,-2.0D0, 2.0D0, 1.0D0,-1.0D0,
     %        -2.0D0,-2.0D0, 4.0D0,-2.0D0,-2.0D0, 4.0D0,
     %        -0.5D0, 1.0D0,-0.5D0,-0.5D0, 1.0D0,-0.5D0/
      DATA DG/ 0.0D0, 5.0D0, 7.0D0,-7.0D0, 0.0D0, 0.0D0,
     %         2.5D0,-1.0D0, 3.5D0,-3.5D0, 0.0D0, 4.5D0,
     %         0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0,
     %        -0.5D0, 0.0D0, 0.5D0, 0.5D0, 0.0D0,-0.5D0,
     %         5.0D0, 0.0D0, 0.0D0, 0.0D0,-7.0D0, 7.0D0,
     %         4.5D0, 0.0D0,-3.5D0, 3.5D0,-1.0D0, 2.5D0/
      DO I=1,36
         B(I)=BR(I)/3.D0
      ENDDO
      DO I=1,72
         B(I+36)=RG(I)
      ENDDO
      DO I=1,36
         B(I+108)=-DG(I)/6.D0
      ENDDO
      RETURN
      END


      SUBROUTINE TMX3H1(PFB,PB1,PB2,PB3,PB4,PB5,PB6)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DE TABLEAUX NECESSAIRES AUX CALCULS DES TABLEAUX
C -----  ELEMENTAIRES DE L ELEMENT MIXTE HEXA M3H1
C
C PARAMETRES RESULTATS :
C ---------------------
C PFB,PB1,PB2,PB3,PB4,PB5,PB6 : TABLEAUX NECESSAIRES AUX CALCULS DES
C                               TABLEAUX ELEMENTAIRES DE L ELEMENT MIXTE
C                               HEXA M3H1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : DUPUY J M ; LAN 189 , PARIS , JANVIER 1981
C.......................................................................
      DOUBLE PRECISION PFB(8,21),PB1(8,10),PB2(8,10),PB3(8,10),PB4(8,16)
     %                ,PB5(8,16),PB6(8,16),PM(8),PN(8),PK(8)
      DATA PM/0.211324865405187D0,0.788675134594813D0
     %,0.788675134594813D0,0.211324865405187D0
     %,0.211324865405187D0,0.788675134594813D0
     %,0.788675134594813D0,0.211324865405187D0/
      DATA PN/0.211324865405187D0,0.211324865405187D0
     %,0.788675134594813D0,0.788675134594813D0
     %,0.211324865405187D0,0.211324865405187D0
     %,0.788675134594813D0,0.788675134594813D0/
      DATA PK/0.211324865405187D0,0.211324865405187D0
     %,0.211324865405187D0,0.211324865405187D0
     %,0.788675134594813D0,0.788675134594813D0
     %,0.788675134594813D0,0.788675134594813D0/
C
      DO I = 1,8
      PFB(I,1) = PM(I) * PM(I)
      PFB(I,2) = PM(I) * PN(I)
      PFB(I,3) = PM(I) * PK(I)
      PFB(I,4) = PM(I) * (PM(I) - 1.D0)
      PFB(I,5) = PM(I) * (PN(I) - 1.D0)
      PFB(I,6) = PM(I) * (PK(I) - 1.D0)
      PFB(I,7) = PN(I) * PN(I)
      PFB(I,8) = PN(I) * PK(I)
      PFB(I,9) = PN(I) * (PM(I) - 1.D0)
      PFB(I,10) = PN(I) * (PN(I) - 1.D0)
      PFB(I,11) = PN(I) * (PK(I) - 1.D0)
      PFB(I,12) = PK(I) * PK(I)
      PFB(I,13) = PK(I) * (PM(I) - 1.D0)
      PFB(I,14) = PK(I) * (PN(I) - 1.D0)
      PFB(I,15) = PK(I) * (PK(I) - 1.D0)
      PFB(I,16) = (PM(I) - 1.D0) * (PM(I) - 1.D0)
      PFB(I,17) = (PM(I) - 1.D0) * (PN(I) - 1.D0)
      PFB(I,18) = (PM(I) - 1.D0) * (PK(I) - 1.D0)
      PFB(I,19) = (PN(I) - 1.D0) * (PN(I) - 1.D0)
      PFB(I,20) = (PN(I) - 1.D0) * (PK(I) - 1.D0)
      PFB(I,21) = (PK(I) - 1.D0) * (PK(I) - 1.D0)
C
      PB1(I,1) = PFB(I,20) * PFB(I,20) * 0.125D0
      PB1(I,2) = - PFB(I,20) * PFB(I,11) * 0.25D0
      PB1(I,3) = PFB(I,20) * PFB(I,8) * 0.25D0
      PB1(I,4) = - PFB(I,20) * PFB(I,14) * 0.25D0
      PB1(I,5) = PFB(I,11) * PFB(I,11) * 0.125D0
      PB1(I,6) = - PFB(I,11) * PFB(I,8) * 0.25D0
      PB1(I,7) = PFB(I,11) * PFB(I,14) * 0.25D0
      PB1(I,8) = PFB(I,8) * PFB(I,8) * 0.125D0
      PB1(I,9) = - PFB(I,8) * PFB(I,14) * 0.25D0
      PB1(I,10) = PFB(I,14) * PFB(I,14) * 0.125D0
C
      PB2(I,1) = PFB(I,18) * PFB(I,18) * 0.125D0
      PB2(I,2) = - PFB(I,6) * PFB(I,18) * 0.25D0
      PB2(I,3) = PFB(I,18) * PFB(I,3) * 0.25D0
      PB2(I,4) = - PFB(I,18) * PFB(I,13) * 0.25D0
      PB2(I,5) = PFB(I,6) * PFB(I,6) * 0.125D0
      PB2(I,6) = - PFB(I,6) * PFB(I,3) * 0.25D0
      PB2(I,7) = PFB(I,6) * PFB(I,13) * 0.25D0
      PB2(I,8) = PFB(I,3) * PFB(I,3) * 0.125D0
      PB2(I,9) = - PFB(I,3) * PFB(I,13) * 0.25D0
      PB2(I,10) = PFB(I,13) * PFB(I,13) * 0.125D0
C
      PB3(I,1) = PFB(I,17) * PFB(I,17) * 0.125D0
      PB3(I,2) = - PFB(I,17) * PFB(I,5) * 0.25D0
      PB3(I,3) = PFB(I,17) * PFB(I,2) * 0.25D0
      PB3(I,4) = - PFB(I,17) * PFB(I,9) * 0.25D0
      PB3(I,5) = PFB(I,5) * PFB(I,5) * 0.125D0
      PB3(I,6) = - PFB(I,5) * PFB(I,2) * 0.25D0
      PB3(I,7) = PFB(I,5) * PFB(I,9) * 0.25D0
      PB3(I,8) = PFB(I,2) * PFB(I,2) * 0.125D0
      PB3(I,9) = - PFB(I,2) * PFB(I,9) * 0.25D0
      PB3(I,10) = PFB(I,9) * PFB(I,9) * 0.125D0
C
      PB4(I,1) = PFB(I,20) * PFB(I,18) * 0.125D0
      PB4(I,2) = - PFB(I,20) * PFB(I,6) * 0.125D0
      PB4(I,3) = PFB(I,20) * PFB(I,3) * 0.125D0
      PB4(I,4) = - PFB(I,20) * PFB(I,13) * 0.125D0
      PB4(I,5) = - PFB(I,11) * PFB(I,18) * 0.125D0
      PB4(I,6) = PFB(I,11) * PFB(I,6) * 0.125D0
      PB4(I,7) = - PFB(I,11) * PFB(I,3) * 0.125D0
      PB4(I,8) = PFB(I,11) * PFB(I,13) * 0.125D0
      PB4(I,9) = PFB(I,8) * PFB(I,18) * 0.125D0
      PB4(I,10) = - PFB(I,8) * PFB(I,6) * 0.125D0
      PB4(I,11) = PFB(I,8) * PFB(I,3) * 0.125D0
      PB4(I,12) = - PFB(I,8) * PFB(I,13) * 0.125D0
      PB4(I,13) = - PFB(I,14) * PFB(I,18) * 0.125D0
      PB4(I,14) = PFB(I,14) * PFB(I,6) * 0.125D0
      PB4(I,15) = - PFB(I,14) * PFB(I,3) * 0.125D0
      PB4(I,16) = PFB(I,14) * PFB(I,13) * 0.125D0
C
      PB5(I,1) = PFB(I,20) * PFB(I,17) * 0.125D0
      PB5(I,2) = - PFB(I,20) * PFB(I,5) * 0.125D0
      PB5(I,3) = PFB(I,20) * PFB(I,2) * 0.125D0
      PB5(I,4) = - PFB(I,20) * PFB(I,9) * 0.125D0
      PB5(I,5) = - PFB(I,11) * PFB(I,17) * 0.125D0
      PB5(I,6) = PFB(I,11) * PFB(I,5) * 0.125D0
      PB5(I,7) = - PFB(I,11) * PFB(I,2) * 0.125D0
      PB5(I,8) = PFB(I,11) * PFB(I,9) * 0.125D0
      PB5(I,9) = PFB(I,8) * PFB(I,17) * 0.125D0
      PB5(I,10) = - PFB(I,8) * PFB(I,5) * 0.125D0
      PB5(I,11) = PFB(I,8) * PFB(I,2) * 0.125D0
      PB5(I,12) = - PFB(I,8) * PFB(I,9) * 0.125D0
      PB5(I,13) = - PFB(I,14) * PFB(I,17) * 0.125D0
      PB5(I,14) = PFB(I,14) * PFB(I,5) * 0.125D0
      PB5(I,15) = - PFB(I,14) * PFB(I,2) * 0.125D0
      PB5(I,16) = PFB(I,14) * PFB(I,9) * 0.125D0
C
      PB6(I,1) = PFB(I,18) * PFB(I,17) * 0.125D0
      PB6(I,2) = - PFB(I,18) * PFB(I,5) * 0.125D0
      PB6(I,3) = PFB(I,18) * PFB(I,2) * 0.125D0
      PB6(I,4) = - PFB(I,18) * PFB(I,9) * 0.125D0
      PB6(I,5) = - PFB(I,6) * PFB(I,17) * 0.125D0
      PB6(I,6) = PFB(I,6) * PFB(I,5) * 0.125D0
      PB6(I,7) = - PFB(I,6) * PFB(I,2) * 0.125D0
      PB6(I,8) = PFB(I,6) * PFB(I,9) * 0.125D0
      PB6(I,9) = PFB(I,3) * PFB(I,17) * 0.125D0
      PB6(I,10) = - PFB(I,3) * PFB(I,5) * 0.125D0
      PB6(I,11) = PFB(I,3) * PFB(I,2) * 0.125D0
      PB6(I,12) = - PFB(I,3) * PFB(I,9) * 0.125D0
      PB6(I,13) = - PFB(I,13) * PFB(I,17) * 0.125D0
      PB6(I,14) = PFB(I,13) * PFB(I,5) * 0.125D0
      PB6(I,15) = - PFB(I,13) * PFB(I,2) * 0.125D0
      PB6(I,16) = PFB(I,13) * PFB(I,9) * 0.125D0
C
      ENDDO
      RETURN
      END


      SUBROUTINE TMX3H2(PFB,PB1,PB2,PB3,PB4,PB5,PB6,AA1,AA2,AA3,B1,B2,
     %                  PFU)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES TABLEAUX ELEMENTAIRES INTERMEDIAIRES POUR
C ----- L ELEMENT MIXTE HEXA M3H2
C
C PARAMETRES RESULTATS :
C ---------------------
C AA1    : TABLEAU CONTENANT LES VALEURS AUX POINTS D INTEGRATION DE
C          GAUSS LEGENDRE ( FORMULE A 27 POINTS ) DES POLYNOMES ,
C          ETANT LES PREMIERES COMPOSANTES DES FONCTIONS DE BASE
C          ASSOCIEES A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C          TEMPERATURE
C AA2    : TABLEAU CONTENANT LES VALEURS AUX POINTS D INTEGRATION DE
C          GAUSS LEGENDRE ( FORMULE A 27 POINTS ) DES POLYNOMES ,
C          ETANT LES DEUXIEMES COMPOSANTES DES FONCTIONS DE BASE
C          ASSOCIEES A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C          TEMPERATURE
C AA3    : TABLEAU CONTENANT LES VALEURS AUX POINTS D INTEGRATION DE
C          GAUSS LEGENDRE ( FORMULE A 27 POINTS ) DES POLYNOMES ,
C          ETANT LES TROISIEMES COMPOSANTES DES FONCTIONS DE BASE
C          ASSOCIEES A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C          TEMPERATURE
C B1     : B1(I,J) , VALEUR DE L INTEGRALE SUR LE CUBE UNITE DE
C          U(J) * DIV(P(I) ; P(I): I-IEME FONCTION DE BASE ASSOCIEE
C          AU GRADIENT DE LA TEMPERATURE ; U(J) J-IEME FONCTION DE BASE
C          ASSOCIEE A LA TEMPERATURE , POUR I = 1 , 24 ET J = 1 , 8
C B2     : B2(I,J) , VALEUR DE L INTEGRALE SUR LE CUBE UNITE DE
C          U(J) * DIV(P(I) ; P(I): I-IEME FONCTION DE BASE ASSOCIEE
C          AU GRADIENT DE LA TEMPERATURE ; U(J) J-IEME FONCTION DE BASE
C          ASSOCIEE A LA TEMPERATURE , POUR I = 25 , 36 ET J = 1 , 8
C PB1 , PB2 , PB3 , PB4 , PB5 , PB6 , PFB , PFU :
C TABLEAUX NECESSAIRES AUX CALCULS DES DIFFERENTES MATRICES ELEMENTAIRES
C DE L ELEMENT HEXA M3H2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : DUPUY J M ; LAN189 , PARIS , JANVIER 1981
C.......................................................................
      DOUBLE PRECISION PFB(27,21),PB1(27,10),PB2(27,10),PB3(27,10),
     %                 PB4(27,16),PB5(27,16),PB6(27,16),P1(12,36),
     %                 P2(12,36),P3(12,36),AA1(27,36),AA2(27,36),
     %                 AA3(27,36),PM(27),PN(27),PK(27),PL(27),A(8,36),
     %                 U(8,8),B1(24,8),B2(12,8),PFU(8,21),RM(8),RN(8),
     %                 RK(8),U1(32),U2(32)
      COMMON / UNITES / LECTEU,NUNITE(31)
      EQUIVALENCE (U1(1) , U(1,1)) , (U2(1) , U(1,5))
C
      DATA PM/
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.5D0,0.8872983346207415D0/
      DATA PN/
     %0.1127016653792585D0,0.1127016653792585D0,0.1127016653792585D0,
     %0.5D0,0.5D0,0.5D0,
     %0.8872983346207415D0,0.8872983346207415D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.1127016653792585D0,0.1127016653792585D0,
     %0.5D0,0.5D0,0.5D0,
     %0.8872983346207415D0,0.8872983346207415D0,0.8872983346207415D0,
     %0.1127016653792585D0,0.1127016653792585D0,0.1127016653792585D0,
     %0.5D0,0.5D0,0.5D0,
     %0.8872983346207415D0,0.8872983346207415D0,0.8872983346207415D0/
      DATA PK/
     %0.1127016653792585D0,0.1127016653792585D0,0.1127016653792585D0,
     %0.1127016653792585D0,0.1127016653792585D0,0.1127016653792585D0,
     %0.1127016653792585D0,0.1127016653792585D0,0.1127016653792585D0,
     %0.5D0,0.5D0,0.5D0,
     %0.5D0,0.5D0,0.5D0,
     %0.5D0,0.5D0,0.5D0,
     %0.8872983346207415D0,0.8872983346207415D0,0.8872983346207415D0,
     %0.8872983346207415D0,0.8872983346207415D0,0.8872983346207415D0,
     %0.8872983346207415D0,0.8872983346207415D0,0.8872983346207415D0/
      DATA PL/
     %0.02143347050754457D0,0.03429355281207132D0,0.02143347050754457D0,
     %0.03429355281207132D0,0.05486968449931411D0,0.03429355281207132D0,
     %0.02143347050754457D0,0.03429355281207132D0,0.02143347050754457D0,
     %0.03429355281207132D0,0.05486968449931411D0,0.03429355281207132D0,
     %0.05486968449931411D0,0.08779149519890257D0,0.05486968449931411D0,
     %0.03429355281207132D0,0.05486968449931411D0,0.03429355281207132D0,
     %0.02143347050754457D0,0.03429355281207132D0,0.02143347050754457D0,
     %0.03429355281207132D0,0.05486968449931411D0,0.03429355281207132D0,
     %0.02143347050754457D0,0.03429355281207132D0,0.02143347050754457D0/
      DATA U1/
     %2.549038105676657D0,-3.232050807568877D0,-3.232050807568877D0,
     %-3.232050807568877D0,4.098076211353314D0,4.098076211353314D0,
     %4.098076211353314D0,-5.196152422706625D0,
     %-0.6830127018922185D0,3.232050807568877D0,0.8660254037844386D0,
     %0.8660254037844375D0,-4.098076211353315D0,-4.098076211353314D0,
     %-1.098076211353314D0,5.196152422706626D0,
     %0.1830127018922191D0,-0.8660254037844382D0,-0.8660254037844388D0,
     %-0.2320508075688765D0,4.098076211353315D0,1.098076211353314D0,
     %1.098076211353314D0,-5.196152422706627D0,
     %-0.6830127018922182D0,0.8660254037844379D0,3.232050807568877D0,
     %0.8660254037844382D0,-4.098076211353315D0,-1.098076211353314D0,
     %-4.098076211353315D0,5.196152422706627D0/
      DATA U2/
     %-0.6830127018922190D0,0.8660254037844379D0,0.8660254037844386D0,
     %3.232050807568877D0,-1.098076211353314D0,-4.098076211353315D0,
     %-4.098076211353315D0,5.196152422706626D0,
     %0.1830127018922187D0,-0.8660254037844379D0,-0.2320508075688767D0,
     %-0.8660254037844379D0,1.098076211353314D0,4.098076211353314D0,
     %1.098076211353314D0,-5.196152422706626D0,
     %-0.04903810567665776D0,0.2320508075688770D0,0.2320508075688770D0,
     %0.2320508075688767D0,-1.098076211353314D0,-1.098076211353314D0,
     %-1.098076211353314D0,5.196152422706626D0,
     %0.1830127018922189D0,-0.2320508075688767D0,-0.8660254037844384D0,
     %-0.8660254037844383D0,1.098076211353314D0,1.098076211353314D0,
     %4.098076211353315D0,-5.196152422706627D0/
      DATA RM/0.211324865405187D0,0.788675134594813D0
     %,0.788675134594813D0,0.211324865405187D0
     %,0.211324865405187D0,0.788675134594813D0
     %,0.788675134594813D0,0.211324865405187D0/
      DATA RN/0.211324865405187D0,0.211324865405187D0
     %,0.788675134594813D0,0.788675134594813D0
     %,0.211324865405187D0,0.211324865405187D0
     %,0.788675134594813D0,0.788675134594813D0/
      DATA RK/0.211324865405187D0,0.211324865405187D0
     %,0.211324865405187D0,0.211324865405187D0
     %,0.788675134594813D0,0.788675134594813D0
     %,0.788675134594813D0,0.788675134594813D0/
C
C     LECTURE DES TABLEAUX P1 P2 P3
      NCVALS = 0
      J      = 12 * 36
      CALL LIRTRD( NCVALS , J , P1 )
C
      NCVALS = 0
      CALL LIRTRD( NCVALS , J , P2 )
C
      NCVALS = 0
      CALL LIRTRD( NCVALS , J , P3 )
C
      DO I = 1 , 36
         DO K = 1 , 27
            AA1(K,I) = P1(1,I) +P1(2,I) * PM(K) + P1(3,I) *
     %                 PN(K) + P1(4,I) * PK(K) + P1(5,I) * PM(K) * PN(K)
     %               + P1(6,I) * PM(K) * PK(K) + P1(7,I) * PN(K) * PK(K)
     %               + P1(8,I) * PM(K) * PN(K) * PK(K) + PM(K) * PM(K) *
     %               (P1(9,I) + P1(10,I) * PN(K) + P1(11,I) * PK(K) +
     %                P1(12,I) * PK(K) * PN(K))
C
            AA2(K,I) = P2(1,I) +P2(2,I) * PM(K) + P2(3,I) *
     %               PN(K) + P2(4,I) * PK(K) + P2(5,I) * PN(K) * PM(K)
     %               + P2(6,I) * PN(K) * PK(K) + P2(7,I) * PM(K) * PK(K)
     %               + P2(8,I) * PN(K) * PM(K) * PK(K) + PN(K) * PN(K) *
     %               (P2(9,I) + P2(10,I) * PM(K) + P2(11,I) * PK(K) +
     %               P2(12,I) * PK(K) * PM(K))
C
            AA3(K,I) = P3(1,I) +P3(2,I) * PM(K) + P3(3,I) *
     %               PN(K) + P3(4,I) * PK(K) + P3(5,I) * PK(K) * PM(K)
     %               + P3(6,I) * PK(K) * PN(K) + P3(7,I) * PM(K) * PN(K)
     %               + P3(8,I) * PK(K) * PM(K) * PN(K) + PK(K) * PK(K) *
     %               (P3(9,I) + P3(10,I) * PM(K) + P3(11,I) * PN(K) +
     %               P3(12,I) * PN(K) * PM(K))
         ENDDO
      ENDDO
C
      DO 4 I = 1 , 24
         A(1,I) = P1(2,I) + P2(3,I) + P3(4,I)
         A(2,I) = 2.D0 * P1(9,I) + P2(5,I) + P3(5,I)
         A(3,I) = P1(5,I) + 2.D0 * P2(9,I) + P3(6,I)
         A(4,I) = P1(6,I) + P2(6,I) + 2.D0 * P3(9,I)
         A(5,I) = 2.D0 * (P1(10,I) + P2(10,I)) + P3(8,I)
         A(6,I) = 2.D0 * (P1(11,I) + P3(10,I)) + P2(8,I)
         A(7,I) = 2.D0 * (P2(11,I) + P3(11,I)) + P1(8,I)
         A(8,I) = 2.D0 * (P1(12,I) +P2(12,I) + P3(12,I))
         DO 5 J = 1 , 8
               B1(I,J) =A(1,I) * U(1,J) + (A(1,I) * U(2,J) + A(1,I) *
     %         U(3,J) + A(1,I) * U(4,J) + A(2,I) * U(1,J) + A(3,I) *
     %         U(1,J) + A(4,I) * U(1,J)) * 0.5D0
     %         + (A(2,I) * U(3,J) + A(2,I) * U(4,J) + A(3,I) * U(2,J)
     %         + A(3,I) * U(4,J) + A(4,I) * U(2,J) + A(4,I) * U(3,J)
     %         + A(1,I) * U(5,J) + A(1,I) * U(6,J) + A(1,I) * U(7,J)
     %         + A(5,I) * U(1,J) + A(6,I) * U(1,J) + A(7,I) * U(1,J)) *
     %         0.25D0
     %         + (A(2,I) * U(2,J) + A(3,I) * U(3,J) + A(4,I) * U(4,J))
     %         * 0.333333333333333D0
     %         + (A(2,I) * U(5,J) + A(2,I) * U(6,J) + A(3,I) * U(5,J)
     %         + A(3,I) * U(7,J) + A(4,I) * U(6,J) + A(4,I) * U(7,J)
     %         + A(5,I) * U(2,J) + A(5,I) * U(3,J) + A(6,I) * U(2,J)
     %         + A(6,I) * U(4,J) + A(7,I) * U(3,J) + A(7,I) * U(4,J))
     %         * 0.166666666666666D0
     %         + (A(5,I) * U(5,J) + A(6,I) * U(6,J) + A(7,I) * U(7,J))
     %         * 0.111111111111111D0
               B1(I,J) = B1(I,J) +
     %         (A(2,I) * U(7,J) + A(3,I) * U(6,J) + A(4,I) * U(5,J)
     %         + A(5,I) * U(4,J) + A(6,I) * U(3,J) + A(7,I) * U(2,J)
     %         + A(1,I) * U(8,J) + A(8,I) * U(1,J)) * 0.125D0
     %         + (A(8,I) * U(2,J) + A(8,I) * U(3,J) + A(8,I) * U(4,J)
     %         + A(2,I) * U(8,J) + A(3,I) * U(8,J) + A(4,I) * U(8,J)
     %         + A(5,I) * U(6,J) + A(5,I) * U(7,J) + A(6,I) * U(5,J)
     %         + A(6,I) * U(7,J) + A(7,I) * U(5,J) + A(7,I) * U(6,J))
     %         * 0.083333333333333D0
     %         + (A(8,I) * U(5,J) + A(8,I) * U(6,J) + A(8,I) * U(7,J)
     %         + A(5,I) * U(8,J) + A(6,I) * U(8,J) + A(7,I) * U(8,J))
     %         * 0.055555555555555D0
     %         + A(8,I) * U(8,J) * 0.037037037037037D0
    5    ENDDO
    4 ENDDO
C
      DO 6 I = 25 , 36
         A(1,I) = P1(2,I) + P2(3,I) + P3(4,I)
         A(2,I) = 2.D0 * P1(9,I) + P2(5,I) + P3(5,I)
         A(3,I) = P1(5,I) + 2.D0 * P2(9,I) + P3(6,I)
         A(4,I) = P1(6,I) + P2(6,I) + 2.D0 * P3(9,I)
         A(5,I) = 2.D0 * (P1(10,I) + P2(10,I)) + P3(8,I)
         A(6,I) = 2.D0 * (P1(11,I) + P3(10,I)) + P2(8,I)
         A(7,I) = 2.D0 * (P2(11,I) + P3(11,I)) + P1(8,I)
         A(8,I) = 2.D0 * (P1(12,I) +P2(12,I) + P3(12,I))
         II = I - 24
         DO 7 J = 1 , 8
            B2(II,J) =A(1,I) * U(1,J) + (A(1,I) * U(2,J) + A(1,I) *
     %         U(3,J) + A(1,I) * U(4,J) + A(2,I) * U(1,J) + A(3,I) *
     %         U(1,J) + A(4,I) * U(1,J)) * 0.5D0
     %         + (A(2,I) * U(3,J) + A(2,I) * U(4,J) + A(3,I) * U(2,J)
     %         + A(3,I) * U(4,J) + A(4,I) * U(2,J) + A(4,I) * U(3,J)
     %         + A(1,I) * U(5,J) + A(1,I) * U(6,J) + A(1,I) * U(7,J)
     %         + A(5,I) * U(1,J) + A(6,I) * U(1,J) + A(7,I) * U(1,J)) *
     %         0.25D0
     %         + (A(2,I) * U(2,J) + A(3,I) * U(3,J) + A(4,I) * U(4,J))
     %         * 0.333333333333333D0
     %         + (A(2,I) * U(5,J) + A(2,I) * U(6,J) + A(3,I) * U(5,J)
     %         + A(3,I) * U(7,J) + A(4,I) * U(6,J) + A(4,I) * U(7,J)
     %         + A(5,I) * U(2,J) + A(5,I) * U(3,J) + A(6,I) * U(2,J)
     %         + A(6,I) * U(4,J) + A(7,I) * U(3,J) + A(7,I) * U(4,J))
     %         * 0.166666666666666D0
     %         + (A(5,I) * U(5,J) + A(6,I) * U(6,J) + A(7,I) * U(7,J))
     %         * 0.111111111111111D0
            B2(II,J) = B2(II,J) +
     %         (A(2,I) * U(7,J) + A(3,I) * U(6,J) + A(4,I) * U(5,J)
     %         + A(5,I) * U(4,J) + A(6,I) * U(3,J) + A(7,I) * U(2,J)
     %         + A(1,I) * U(8,J) + A(8,I) * U(1,J)) * 0.125D0
     %         + (A(8,I) * U(2,J) + A(8,I) * U(3,J) + A(8,I) * U(4,J)
     %         + A(2,I) * U(8,J) + A(3,I) * U(8,J) + A(4,I) * U(8,J)
     %         + A(5,I) * U(6,J) + A(5,I) * U(7,J) + A(6,I) * U(5,J)
     %         + A(6,I) * U(7,J) + A(7,I) * U(5,J) + A(7,I) * U(6,J))
     %         * 0.083333333333333D0
     %         + (A(8,I) * U(5,J) + A(8,I) * U(6,J) + A(8,I) * U(7,J)
     %         + A(5,I) * U(8,J) + A(6,I) * U(8,J) + A(7,I) * U(8,J))
     %         * 0.055555555555555D0
     %         + A(8,I) * U(8,J) * 0.037037037037037D0
    7    ENDDO
    6 ENDDO
C
      DO I = 1 , 27
      PFB(I,1) = PM(I) * PM(I)
      PFB(I,2) = PM(I) * PN(I)
      PFB(I,3) = PM(I) * PK(I)
      PFB(I,4) = PM(I) * (PM(I) - 1.D0)
      PFB(I,5) = PM(I) * (PN(I) - 1.D0)
      PFB(I,6) = PM(I) * (PK(I) - 1.D0)
      PFB(I,7) = PN(I) * PN(I)
      PFB(I,8) = PN(I) * PK(I)
      PFB(I,9) = PN(I) * (PM(I) - 1.D0)
      PFB(I,10) = PN(I) * (PN(I) - 1.D0)
      PFB(I,11) = PN(I) * (PK(I) - 1.D0)
      PFB(I,12) = PK(I) * PK(I)
      PFB(I,13) = PK(I) * (PM(I) - 1.D0)
      PFB(I,14) = PK(I) * (PN(I) - 1.D0)
      PFB(I,15) = PK(I) * (PK(I) - 1.D0)
      PFB(I,16) = (PM(I) - 1.D0) * (PM(I) - 1.D0)
      PFB(I,17) = (PM(I) - 1.D0) * (PN(I) - 1.D0)
      PFB(I,18) = (PM(I) - 1.D0) * (PK(I) - 1.D0)
      PFB(I,19) = (PN(I) - 1.D0) * (PN(I) - 1.D0)
      PFB(I,20) = (PN(I) - 1.D0) * (PK(I) - 1.D0)
      PFB(I,21) = (PK(I) - 1.D0) * (PK(I) - 1.D0)
C
      PB1(I,1) = PFB(I,20) * PFB(I,20) * PL(I)
      PB1(I,2) = - PFB(I,20) * PFB(I,11) * 2.D0 * PL(I)
      PB1(I,3) = PFB(I,20) * PFB(I,8) * 2.D0 * PL(I)
      PB1(I,4) = - PFB(I,20) * PFB(I,14) * 2.D0 * PL(I)
      PB1(I,5) = PFB(I,11) * PFB(I,11) * PL(I)
      PB1(I,6) = - PFB(I,11) * PFB(I,8) * 2.D0 * PL(I)
      PB1(I,7) = PFB(I,11) * PFB(I,14) * 2.D0 * PL(I)
      PB1(I,8) = PFB(I,8) * PFB(I,8) * PL(I)
      PB1(I,9) = - PFB(I,8) * PFB(I,14) * 2.D0 * PL(I)
      PB1(I,10) = PFB(I,14) * PFB(I,14) * PL(I)
C
      PB2(I,1) = PFB(I,18) * PFB(I,18) * PL(I)
      PB2(I,2) = - PFB(I,6) * PFB(I,18) * 2.D0 * PL(I)
      PB2(I,3) = PFB(I,18) * PFB(I,3) * 2.D0 * PL(I)
      PB2(I,4) = - PFB(I,18) * PFB(I,13) * 2.D0 * PL(I)
      PB2(I,5) = PFB(I,6) * PFB(I,6) * PL(I)
      PB2(I,6) = - PFB(I,6) * PFB(I,3) * 2.D0 * PL(I)
      PB2(I,7) = PFB(I,6) * PFB(I,13) * 2.D0 * PL(I)
      PB2(I,8) = PFB(I,3) * PFB(I,3) * PL(I)
      PB2(I,9) = - PFB(I,3) * PFB(I,13) * 2.D0 * PL(I)
      PB2(I,10) = PFB(I,13) * PFB(I,13) * PL(I)
C
      PB3(I,1) = PFB(I,17) * PFB(I,17) * PL(I)
      PB3(I,2) = - PFB(I,17) * PFB(I,5) * 2.D0 * PL(I)
      PB3(I,3) = PFB(I,17) * PFB(I,2) * 2.D0 * PL(I)
      PB3(I,4) = - PFB(I,17) * PFB(I,9) * 2.D0 * PL(I)
      PB3(I,5) = PFB(I,5) * PFB(I,5) * PL(I)
      PB3(I,6) = - PFB(I,5) * PFB(I,2) * 2.D0 * PL(I)
      PB3(I,7) = PFB(I,5) * PFB(I,9) * 2.D0 * PL(I)
      PB3(I,8) = PFB(I,2) * PFB(I,2) * PL(I)
      PB3(I,9) = - PFB(I,2) * PFB(I,9) * 2.D0 * PL(I)
      PB3(I,10) = PFB(I,9) * PFB(I,9) * PL(I)
C
      PB4(I,1) = PFB(I,20) * PFB(I,18) * PL(I)
      PB4(I,2) = - PFB(I,20) * PFB(I,6) * PL(I)
      PB4(I,3) = PFB(I,20) * PFB(I,3) * PL(I)
      PB4(I,4) = - PFB(I,20) * PFB(I,13) * PL(I)
      PB4(I,5) = - PFB(I,11) * PFB(I,18) * PL(I)
      PB4(I,6) = PFB(I,11) * PFB(I,6) * PL(I)
      PB4(I,7) = - PFB(I,11) * PFB(I,3) * PL(I)
      PB4(I,8) = PFB(I,11) * PFB(I,13) * PL(I)
      PB4(I,9) = PFB(I,8) * PFB(I,18) * PL(I)
      PB4(I,10) = - PFB(I,8) * PFB(I,6) * PL(I)
      PB4(I,11) = PFB(I,8) * PFB(I,3) * PL(I)
      PB4(I,12) = - PFB(I,8) * PFB(I,13) * PL(I)
      PB4(I,13) = - PFB(I,14) * PFB(I,18) * PL(I)
      PB4(I,14) = PFB(I,14) * PFB(I,6) * PL(I)
      PB4(I,15) = - PFB(I,14) * PFB(I,3) * PL(I)
      PB4(I,16) = PFB(I,14) * PFB(I,13) * PL(I)
C
      PB5(I,1) = PFB(I,20) * PFB(I,17) * PL(I)
      PB5(I,2) = - PFB(I,20) * PFB(I,5) * PL(I)
      PB5(I,3) = PFB(I,20) * PFB(I,2) * PL(I)
      PB5(I,4) = - PFB(I,20) * PFB(I,9) * PL(I)
      PB5(I,5) = - PFB(I,11) * PFB(I,17) * PL(I)
      PB5(I,6) = PFB(I,11) * PFB(I,5) * PL(I)
      PB5(I,7) = - PFB(I,11) * PFB(I,2) * PL(I)
      PB5(I,8) = PFB(I,11) * PFB(I,9) * PL(I)
      PB5(I,9) = PFB(I,8) * PFB(I,17) * PL(I)
      PB5(I,10) = - PFB(I,8) * PFB(I,5) * PL(I)
      PB5(I,11) = PFB(I,8) * PFB(I,2) * PL(I)
      PB5(I,12) = - PFB(I,8) * PFB(I,9) * PL(I)
      PB5(I,13) = - PFB(I,14) * PFB(I,17) * PL(I)
      PB5(I,14) = PFB(I,14) * PFB(I,5) * PL(I)
      PB5(I,15) = - PFB(I,14) * PFB(I,2) * PL(I)
      PB5(I,16) = PFB(I,14) * PFB(I,9) * PL(I)
C
      PB6(I,1) = PFB(I,18) * PFB(I,17) * PL(I)
      PB6(I,2) = - PFB(I,18) * PFB(I,5) * PL(I)
      PB6(I,3) = PFB(I,18) * PFB(I,2) * PL(I)
      PB6(I,4) = - PFB(I,18) * PFB(I,9) * PL(I)
      PB6(I,5) = - PFB(I,6) * PFB(I,17) * PL(I)
      PB6(I,6) = PFB(I,6) * PFB(I,5) * PL(I)
      PB6(I,7) = - PFB(I,6) * PFB(I,2) * PL(I)
      PB6(I,8) = PFB(I,6) * PFB(I,9) * PL(I)
      PB6(I,9) = PFB(I,3) * PFB(I,17) * PL(I)
      PB6(I,10) = - PFB(I,3) * PFB(I,5) * PL(I)
      PB6(I,11) = PFB(I,3) * PFB(I,2) * PL(I)
      PB6(I,12) = - PFB(I,3) * PFB(I,9) * PL(I)
      PB6(I,13) = - PFB(I,13) * PFB(I,17) * PL(I)
      PB6(I,14) = PFB(I,13) * PFB(I,5) * PL(I)
      PB6(I,15) = - PFB(I,13) * PFB(I,2) * PL(I)
      PB6(I,16) = PFB(I,13) * PFB(I,9) * PL(I)
      ENDDO

      DO I = 1 , 8
      PFU(I,1) = RM(I) * RM(I)
      PFU(I,2) = RM(I) * RN(I)
      PFU(I,3) = RM(I) * RK(I)
      PFU(I,4) = RM(I) * (RM(I) - 1.D0)
      PFU(I,5) = RM(I) * (RN(I) - 1.D0)
      PFU(I,6) = RM(I) * (RK(I) - 1.D0)
      PFU(I,7) = RN(I) * RN(I)
      PFU(I,8) = RN(I) * RK(I)
      PFU(I,9) = RN(I) * (RM(I) - 1.D0)
      PFU(I,10) = RN(I) * (RN(I) - 1.D0)
      PFU(I,11) = RN(I) * (RK(I) - 1.D0)
      PFU(I,12) = RK(I) * RK(I)
      PFU(I,13) = RK(I) * (RM(I) - 1.D0)
      PFU(I,14) = RK(I) * (RN(I) - 1.D0)
      PFU(I,15) = RK(I) * (RK(I) - 1.D0)
      PFU(I,16) = (RM(I) - 1.D0) * (RM(I) - 1.D0)
      PFU(I,17) = (RM(I) - 1.D0) * (RN(I) - 1.D0)
      PFU(I,18) = (RM(I) - 1.D0) * (RK(I) - 1.D0)
      PFU(I,19) = (RN(I) - 1.D0) * (RN(I) - 1.D0)
      PFU(I,20) = (RN(I) - 1.D0) * (RK(I) - 1.D0)
      PFU(I,21) = (RK(I) - 1.D0) * (RK(I) - 1.D0)
      ENDDO

      RETURN
      END


      SUBROUTINE TMX3T2(A1,A2,A3,B1,B2,B3,BB1,BB2)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES TABLEAUX INTERMEDIAIRES NECESSAIRES AU CALCUL DES
C ----- TABLEAUX ELEMENTAIRES DE L ELEMENT MIXTE TETR M3T2
C
C PARAMETRES RESULTATS :
C ---------------------
C A1    : VALEURS DES INTEGRALES SUR LE TETRAEDRE UNITE DE P1(I) * P1(J)
C         P1(I) : PREMIERE COMPOSANTE DE LA I-IEME FONCTION DE BASE
C         ASSOCIEE A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C         TEMPERATURE ; TABLEAU STOCKE DE FACON SYMETRIQUE : 1. 2. 4.
C                                                               3. 5.
C                                                                  6.
C A2    : VALEURS DES INTEGRALES SUR LE TETRAEDRE UNITE DE P2(I) * P2(J)
C         P2(I) : DEUXIEME COMPOSANTE DE LA I-IEME FONCTION DE BASE
C         ASSOCIEE A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C         TEMPERATURE ; TABLEAU STOCKE DE FACON SYMETRIQUE : 1. 2. 4.
C                                                               3. 5.
C                                                                  6.
C A3    : VALEURS DES INTEGRALES SUR LE TETRAEDRE UNITE DE P3(I) * P3(J)
C         P3(I) : TROISIEME COMPOSANTE DE LA I-IEME FONCTION DE BASE
C         ASSOCIEE A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C         TEMPERATURE ; TABLEAU STOCKE DE FACON SYMETRIQUE : 1. 2. 4.
C                                                               3. 5.
C                                                                  6.
C B1    : VALEURS DES INTEGRALES SUR LE TETRAEDRE UNITE DE
C        ( P1(I) * P2(J) + P1(J) * P2(I) ) ; P1(I) 1-ERE COMPOSANTE ,
C         P2(I) 2-IEME COMPOSANTE DE LA I-IEME FONCTION DE BASE
C         ASSOCIEE A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C         TEMPERATURE ; TABLEAU STOCKE DE FACON SYMETRIQUE : 1. 2. 4.
C                                                               3. 5.
C                                                                  6.
C B2    : VALEURS DES INTEGRALES SUR LE TETRAEDRE UNITE DE
C        ( P1(I) * P3(J) + P1(J) * P3(I) ) ; P1(I) 1-ERE COMPOSANTE ,
C         P3(I) 3-IEME COMPOSANTE DE LA I-IEME FONCTION DE BASE
C         ASSOCIEE A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C         TEMPERATURE ; TABLEAU STOCKE DE FACON SYMETRIQUE : 1. 2. 4.
C                                                               3. 5.
C                                                                  6.
C B3    : VALEURS DES INTEGRALES SUR LE TETRAEDRE UNITE DE
C        ( P2(I) * P3(J) + P2(J) * P3(I) ) ; P2(I) 2-IEME COMPOSANTE ,
C         P3(I) 3-IEME COMPOSANTE DE LA I-IEME FONCTION DE BASE
C         ASSOCIEE A L INCONNUE VARIATIONNELLE : GRADIENT DE LA
C         TEMPERATURE ; TABLEAU STOCKE DE FACON SYMETRIQUE : 1. 2. 4.
C                                                               3. 5.
C                                                                  6.
C BB1    : VALEUR SUR LE TETRAEDRE UNITE DE L INTEGRALE :
C        : DIV(P(I) * U(J) ; P(I) I-FCT DE BASE ASSOCIEE AU GRADIENT
C          , U(J) J-IEME FCT DE BASE ASSOCIEE A LA TEMPERATURE
C          POUR I = 1 , 12 ; ET J = 1 , 4
C BB2    : VALEUR SUR LE TETRAEDRE UNITE DE L INTEGRALE :
C        : DIV(P(I) * U(J) ; P(I) I-FCT DE BASE ASSOCIEE AU GRADIENT
C          , U(J) J-IEME FCT DE BASE ASSOCIEE A LA TEMPERATURE
C          POUR I = 13 , 15 ; ET J = 1 , 4
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : DUPUY J M ; LAN 189 , PARIS , JANVIER 1981
C.......................................................................
      DOUBLE PRECISION A1(120),A2(120),A3(120),B1(120),B2(120),B3(120),
     %                 BB1(12,4),BB2(3,4),A,B,C,D,P1(5,15),P2(5,15),
     %                 P3(5,15),U(4,4),XP
      DATA P1/0.D0,8.D0,0.D0,0.D0,-10.D0,0.D0,-8.D0,0.D0,0.D0,10.D0,
     %0.D0,-8.D0,0.D0,0.D0,10.D0,-2.D0,12.D0,4.D0,0.D0,-10.D0,
     %2.D0,8.D0,-4.D0,-4.D0,-10.D0,-2.D0,12.D0,0.D0,4.D0,-10.D0,
     %0.D0,-8.D0,0.D0,0.D0,10.D0,0.D0,-8.D0,0.D0,0.D0,10.D0,
     %0.D0,8.D0,0.D0,0.D0,-10.D0,0.D0,-8.D0,0.D0,0.D0,10.D0,
     %0.D0,8.D0,0.D0,0.D0,-10.D0,0.D0,-8.D0,0.D0,0.D0,10.D0,
     %0.D0,120.D0,0.D0,0.D0,-120.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,0.D0,0.D0/
      DATA P2/0.D0,0.D0,-8.D0,0.D0,10.D0,0.D0,0.D0,-8.D0,0.D0,10.D0,
     %0.D0,0.D0,8.D0,0.D0,-10.D0,0.D0,0.D0,8.D0,0.D0,-10.D0,
     %0.D0,0.D0,-8.D0,0.D0,10.D0,0.D0,0.D0,-8.D0,0.D0,10.D0,
     %-2.D0,0.D0,12.D0,4.D0,-10.D0,2.D0,-4.D0,8.D0,-4.D0,-10.D0,
     %-2.D0,4.D0,12.D0,0.D0,-10.D0,0.D0,0.D0,-8.D0,0.D0,10.D0,
     %0.D0,0.D0,-8.D0,0.D0,10.D0,0.D0,0.D0,8.D0,0.D0,-10.D0,
     %0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,120.D0,0.D0,-120.D0,
     %0.D0,0.D0,0.D0,0.D0,0.D0/
      DATA P3/-2.D0,4.D0,0.D0,12.D0,-10.D0,2.D0,-4.D0,-4.D0,8.D0,-10.D0,
     %-2.D0,0.D0,4.D0,12.D0,-10.D0,0.D0,0.D0,0.D0,-8.D0,10.D0,
     %0.D0,0.D0,0.D0,-8.D0,10.D0,0.D0,0.D0,0.D0,8.D0,-10.D0,
     %0.D0,0.D0,0.D0,8.D0,-10.D0,0.D0,0.D0,0.D0,-8.D0,10.D0,
     %0.D0,0.D0,0.D0,-8.D0,10.D0,0.D0,0.D0,0.D0,8.D0,-10.D0,
     %0.D0,0.D0,0.D0,-8.D0,10.D0,0.D0,0.D0,0.D0,-8.D0,10.D0,
     %0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,120.D0,-120.D0/
      DATA U/1.927050983124842D0,-2.236067977499789D0,-2.236067977499789
     %D0,-2.236067977499789D0,-0.3090169943749473D0,2.236067977499789D0,
     %0.D0,0.D0,-0.3090169943749473D0,0.D0,2.236067977499789D0,
     %0.D0,-0.3090169943749473D0,0.D0,0.D0,2.236067977499789D0/
C
      DO 31 I = 1 , 5
         XP = P1(I,2)
         P1(I,2) = P1(I,3)
         P1(I,3) = XP
         XP = P2(I,2)
         P2(I,2) = P2(I,3)
         P2(I,3) = XP
         XP = P3(I,2)
         P3(I,2) = P3(I,3)
         P3(I,3) = XP
   31 ENDDO
C
      DO 1 J=1,15
           DO 2 I=1,J
               II=(J-1)*J/2+I
               A1(II)=(4.D0*P1(1,I)*P1(1,J)+P1(1,I)*P1(2,J)+
     %         P1(2,I)*P1(1,J)+P1(1,I)*P1(3,J)+P1(3,I)*P1(1,J)
     %         +P1(1,I)*P1(4,J)+P1(4,I)*P1(1,J))
     %          *0.041666666666666D0
     %         +(2.D0*(P1(1,I)*P1(5,J)+P1(5,I)*P1(1,J))+
     %         P1(2,I)*P1(2,J)+P1(3,I)*P1(3,J)+
     %         P1(4,I)*P1(4,J))*0.0166666666666666D0
     %         +(P1(2,I)*P1(3,J)+P1(3,I)*P1(2,J)+P1(2,I)*P1(4,J)+
     %         P1(4,I)*P1(2,J)+P1(3,I)*P1(4,J)+P1(4,I)*P1(3,J))
     %          *0.008333333333333D0
     %         +(P1(2,I)*P1(5,J)+P1(5,I)*P1(2,J))
     %          *0.013888888888888D0
     %         +(P1(3,I)*P1(5,J)+P1(5,I)*P1(3,J)+P1(4,I)*P1(5,J)
     %         +P1(5,I)*P1(4,J))*0.006944444444444D0
     %         +(P1(5,I)*P1(5,J))*0.0119047619047619D0
C
               A2(II)=(4.D0*P2(1,I)*P2(1,J)+P2(1,I)*P2(3,J)+
     %         P2(3,I)*P2(1,J)+P2(1,I)*P2(2,J)+P2(2,I)*P2(1,J)
     %         +P2(1,I)*P2(4,J)+P2(4,I)*P2(1,J))
     %          *0.041666666666666D0
     %         +(2.D0*(P2(1,I)*P2(5,J)+P2(5,I)*P2(1,J))+
     %         P2(3,I)*P2(3,J)+P2(2,I)*P2(2,J)+
     %         P2(4,I)*P2(4,J))*0.0166666666666666D0
     %         +(P2(2,I)*P2(3,J)+P2(3,I)*P2(2,J)+P2(2,I)*P2(4,J)+
     %         P2(4,I)*P2(2,J)+P2(3,I)*P2(4,J)+P2(4,I)*P2(3,J))
     %          *0.008333333333333D0
     %         +(P2(3,I)*P2(5,J)+P2(5,I)*P2(3,J))
     %          *0.013888888888888D0
     %         +(P2(2,I)*P2(5,J)+P2(5,I)*P2(2,J)+P2(4,I)*P2(5,J)
     %         +P2(5,I)*P2(4,J))*0.006944444444444D0
     %         +(P2(5,I)*P2(5,J))*0.0119047619047619D0
C
               A3(II)=(4.D0*P3(1,I)*P3(1,J)+P3(1,I)*P3(4,J)+
     %         P3(4,I)*P3(1,J)+P3(1,I)*P3(3,J)+P3(3,I)*P3(1,J)
     %         +P3(1,I)*P3(2,J)+P3(2,I)*P3(1,J))
     %          *0.041666666666666D0
     %         +(2.D0*(P3(1,I)*P3(5,J)+P3(5,I)*P3(1,J))+
     %         P3(4,I)*P3(4,J)+P3(3,I)*P3(3,J)+
     %         P3(2,I)*P3(2,J))*0.0166666666666666D0
     %         +(P3(2,I)*P3(3,J)+P3(3,I)*P3(2,J)+P3(2,I)*P3(4,J)+
     %         P3(4,I)*P3(2,J)+P3(3,I)*P3(4,J)+P3(4,I)*P3(3,J))
     %          *0.008333333333333D0
     %         +(P3(4,I)*P3(5,J)+P3(5,I)*P3(4,J))
     %          *0.013888888888888D0
     %         +(P3(3,I)*P3(5,J)+P3(5,I)*P3(3,J)+P3(2,I)*P3(5,J)
     %         +P3(5,I)*P3(2,J))*0.006944444444444D0
     %         +(P3(5,I)*P3(5,J))*0.0119047619047619D0
C
               B1(II)=(4.D0*(P1(1,I)*P2(1,J)+P1(1,J)*P2(1,I))+
     %         P1(1,I)*P2(2,J)+P1(1,J)*P2(2,I)+P1(2,I)*P2(1,J)
     %         +P1(2,J)*P2(1,I)+P1(1,I)*P2(3,J)+P1(1,J)*P2(3,I)+
     %         P1(3,I)*P2(1,J)+P1(3,J)*P2(1,I)+P1(1,I)*P2(4,J)+
     %         P1(1,J)*P2(4,I)+P1(4,I)*P2(1,J)+P1(4,J)*P2(1,I))
     %          *0.041666666666666D0
     %         +(P1(2,I)*P2(2,J)+P1(2,J)*P2(2,I)+P1(3,I)*P2(3,J)+
     %         P1(3,J)*P2(3,I)+P1(4,I)*P2(4,J)+P2(4,I)*P1(4,J)+
     %         2.D0*(P1(1,I)*P2(5,J)+P1(1,J)*P2(5,I)+P1(5,I)*P2(1,J)
     %         +P1(5,J)*P2(1,I)))*0.016666666666666D0
     %         +(P1(2,I)*P2(3,J)+P1(2,J)*P2(3,I)+P1(3,I)*P2(2,J)+
     %         P1(3,J)*P2(2,I)+P1(2,I)*P2(4,J)+P1(2,J)*P2(4,I)+
     %         P1(4,I)*P2(2,J)+P1(4,J)*P2(2,I)+P1(3,I)*P2(4,J)+
     %         P1(3,J)*P2(4,I)+P1(4,I)*P2(3,J)+P1(4,J)*P2(3,I))
     %          *0.008333333333333D0+
     %         (P1(2,I)*P2(5,J)+P1(2,J)*P2(5,I)+P1(4,I)*P2(5,J)+
     %         P1(4,J)*P2(5,I)+P2(3,J)*P1(5,I)+P2(3,I)*P1(5,J)+
     %         P2(4,J)*P1(5,I)+P2(4,I)*P1(5,J))
     %          *0.006944444444444D0
               B1(II)=B1(II)+(P1(5,I)*P2(2,J)+P1(5,J)*P2(2,I)+P1(3,I)
     %         *P2(5,J)+
     %         P1(3,J)*P2(5,I))*0.013888888888888D0
     %         +(P1(5,I)*P2(5,J)+P1(5,J)*P2(5,I))*0.005952380952380952D0
C
               B2(II)=(4.D0*(P1(1,I)*P3(1,J)+P1(1,J)*P3(1,I))+
     %         P1(1,I)*P3(2,J)+P1(1,J)*P3(2,I)+P1(2,I)*P3(1,J)+
     %         P1(2,J)*P3(1,I)+P1(1,I)*P3(3,J)+P1(1,J)*P3(3,I)
     %         +P1(4,I)*P3(1,J)+P1(4,J)*P3(1,I)+P1(1,I)*P3(4,J)+
     %         P1(1,J)*P3(4,I)+P1(3,I)*P3(1,J)+P1(3,J)*P3(1,I))
     %          *0.041666666666666D0
     %         +(P1(2,I)*P3(2,J)+P1(2,J)*P3(2,I)+P1(3,I)*P3(3,J)+
     %         P1(3,J)*P3(3,I)+P1(4,I)*P3(4,J)+P3(4,I)*P1(4,J)+
     %         2.D0*(P1(1,I)*P3(5,J)+P1(1,J)*P3(5,I)+P1(5,I)*P3(1,J)
     %         +P1(5,J)*P3(1,I)))*0.016666666666666D0
     %         +(P1(2,I)*P3(3,J)+P1(2,J)*P3(3,I)+P1(4,I)*P3(2,J)+
     %         P1(4,J)*P3(2,I)+P1(2,I)*P3(4,J)+P1(2,J)*P3(4,I)+
     %         P1(3,I)*P3(2,J)+P1(3,J)*P3(2,I)+P1(3,I)*P3(4,J)+
     %         P1(3,J)*P3(4,I)+P1(4,I)*P3(3,J)+P1(4,J)*P3(3,I))
     %          *0.008333333333333D0+
     %         (P1(2,I)*P3(5,J)+P1(2,J)*P3(5,I)+P1(3,I)*P3(5,J)+
     %         P1(3,J)*P3(5,I)+P3(3,J)*P1(5,I)+P3(3,I)*P1(5,J)+
     %         P3(4,J)*P1(5,I)+P3(4,I)*P1(5,J))
     %          *0.006944444444444D0
               B2(II)=B2(II)+(P1(5,I)*P3(2,J)+P1(5,J)*P3(2,I)+P1(4,I)
     %         *P3(5,J)+
     %         P1(4,J)*P3(5,I))*0.013888888888888D0
     %         +(P1(5,I)*P3(5,J)+P1(5,J)*P3(5,I))*0.005952380952380952D0
C
               B3(II)=(4.D0*(P3(1,I)*P2(1,J)+P3(1,J)*P2(1,I))+
     %         P3(1,I)*P2(4,J)+P3(1,J)*P2(4,I)+P3(2,I)*P2(1,J)
     %         +P3(2,J)*P2(1,I)+P3(1,I)*P2(3,J)+P3(1,J)*P2(3,I)+
     %         P3(3,I)*P2(1,J)+P3(3,J)*P2(1,I)+P3(1,I)*P2(2,J)+
     %         P3(1,J)*P2(2,I)+P3(4,I)*P2(1,J)+P3(4,J)*P2(1,I))
     %          *0.041666666666666D0
     %         +(P3(2,I)*P2(2,J)+P3(2,J)*P2(2,I)+P3(3,I)*P2(3,J)+
     %         P3(3,J)*P2(3,I)+P3(4,I)*P2(4,J)+P2(4,I)*P3(4,J)+
     %         2.D0*(P3(1,I)*P2(5,J)+P3(1,J)*P2(5,I)+P3(5,I)*P2(1,J)
     %         +P3(5,J)*P2(1,I)))*0.016666666666666D0
     %         +(P3(2,I)*P2(3,J)+P3(2,J)*P2(3,I)+P3(3,I)*P2(4,J)+
     %         P3(3,J)*P2(4,I)+P3(4,I)*P2(3,J)+P3(4,J)*P2(3,I)+
     %         P3(2,I)*P2(4,J)+P3(2,J)*P2(4,I)+P3(3,I)*P2(2,J)+
     %         P3(3,J)*P2(2,I)+P3(4,I)*P2(2,J)+P3(4,J)*P2(2,I))
     %          *0.008333333333333D0+
     %         (P3(2,I)*P2(5,J)+P3(2,J)*P2(5,I)+P3(4,I)*P2(5,J)+
     %         P3(4,J)*P2(5,I)+P2(3,J)*P3(5,I)+P2(3,I)*P3(5,J)+
     %         P2(2,J)*P3(5,I)+P2(2,I)*P3(5,J))
     %          *0.006944444444444D0
               B3(II)=B3(II)+(P3(5,I)*P2(4,J)+P3(5,J)*P2(4,I)+P3(3,I)
     %        *P2(5,J)+
     %         P3(3,J)*P2(5,I))*0.013888888888888D0
     %         +(P3(5,I)*P2(5,J)+P3(5,J)*P2(5,I))*0.005952380952380952D0
    2       ENDDO
    1 ENDDO
C
      DO3 I=1,12
         A=P1(2,I)+P2(3,I)+P3(4,I)
         B=2.D0*P1(5,I)+P2(5,I)+P3(5,I)
         C=P1(5,I)+2.D0*P2(5,I)+P3(5,I)
         D=P1(5,I)+P2(5,I)+2.D0*P3(5,I)
         DO 4 J = 1 , 4
            BB1(I,J) = U(1,J) * A / 6.D0 + (A * U(2,J) + B * U(1,J) +
     %      A * U(3,J) + C * U(1,J) + A * U(4,J) + D * U(1,J)) / 24.D0
     %      + (B * U(2,J) + C * U(3,J) + D * U(4,J)) /60.D0 +
     %      (B * U(3,J) + B * U(4,J) + C * U(2,J) + C * U(4,J) +
     %      D * U(2,J) + D * U(3,J)) / 120.D0
    4    ENDDO
    3 ENDDO
C
      DO5 I=13,15
         A=P1(2,I)+P2(3,I)+P3(4,I)
         B=2.D0*P1(5,I)+P2(5,I)+P3(5,I)
         C=P1(5,I)+2.D0*P2(5,I)+P3(5,I)
         D=P1(5,I)+P2(5,I)+2.D0*P3(5,I)
         II=I-12
         DO 6 J = 1 , 4
            BB2(II,J) = U(1,J) * A / 6.D0 + (A * U(2,J) + B * U(1,J) +
     %      A * U(3,J) + C * U(1,J) + A * U(4,J) + D * U(1,J)) / 24.D0
     %      + (B * U(2,J) + C * U(3,J) + D * U(4,J)) /60.D0 +
     %      (B * U(3,J) + B * U(4,J) + C * U(2,J) + C * U(4,J) +
     %      D * U(2,J) + D * U(3,J)) / 120.D0
    6    ENDDO
    5 ENDDO
      END


      SUBROUTINE TMXQ10( PFB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULE DES TABLEAUX INTERMEDIAIRES POUR L ELEMENT
C ----- QUAD MQ10
C
C TABLEAU PFB:
C VALEURS AUX POINTS D INTEGRATION DES FONCTIONS (DANS L ORDRE):
C (1-Y)*(1-Y),Y*(1-Y),Y*Y,(1-X)*(1-Y),X*(1-Y),Y*(1-X),X*Y
C (1-X)*(1-X),X*(1-X),X*X
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR: DUPUY J M ,LAB:LAN189,PARIS,JANVIER 1981
C.......................................................................
      DOUBLE PRECISION  PM,PN,PFB
      DIMENSION         PM(4),PN(4),PFB(4,10)
      DATA PM/0.211324865405187D0,0.788675134594813D0,
     %0.788675134594813D0,0.211324865405187D0/
      DATA PN/0.211324865405187D0,0.211324865405187D0,
     %0.788675134594813D0,0.788675134594813D0/
C
      DO1 I= 1,4
         PFB(I,1)=(1.D0-PN(I))*(1.D0-PN(I))
         PFB(I,2)=(1.D0-PN(I))*PN(I)
         PFB(I,3)=PN(I)*PN(I)
         PFB(I,4)=(1.D0-PM(I))*(1.D0-PN(I))
         PFB(I,5)=PM(I)*(1.D0-PN(I))
         PFB(I,6)=PN(I)*(1.D0-PM(I))
         PFB(I,7)=PN(I)*PM(I)
         PFB(I,8)=(1.D0-PM(I))*(1.D0-PM(I))
         PFB(I,9)=(1.D0-PM(I))*PM(I)
         PFB(I,10)=PM(I)*PM(I)
    1 ENDDO
      END


      SUBROUTINE TMXQ21(AA1,AA2,E,PFB)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULE DES TABLEAUX INTERMEDIAIRES POUR L ELEMENT
C ----- QUAD MQ21
C
C PL: TABLEAU DES POIDS DE LA FORMULE D INTEGRATION DE GAUSS LEGENDRE
C A NEUF POINTS
C
C PM:TABLEAU DES ABSCISSES DES POINTS D INTEGRATION
C
C PN:TABLEAU DES ORDONNEES DES POINTS D INTEGRATION
C
C TABLEAU AA1:TABLEAU CONTENANT LES VALEURS AUX POINTS D INTEGRATION
C(FORMULE D INTEGRATION A NEUF POINTS,INTEGRANT EXACTEMENT LES POLYNOMES
C DE DEGRE < OU = A 5),DU POLYNOME,ETANT LA PREMIERE COMPOSANTE DE LA
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE:GRADIENT DE LA
C TEMPERATURE
C
C TABLEAU AA2:TABLEAU CONTENANT LES VALEURS AUX POINTS D INTEGRATION
C(FORMULE D INTEGRATION A NEUF POINTS,INTEGRANT EXACTEMENT LES POLYNOMES
C DE DEGRE < OU = A 5),DU POLYNOME,ETANT LA DEUXIEME COMPOSANTE DE LA
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE:GRADIENT DE LA
C TEMPERATURE
C
C TABLEAU E:TABLEAU E(I,J);VALEUR DE L INTEGRALE SUR LE TRIANGLE DE
C REFERENCE DE U(J)*DIVP(I);P(I):I IEME FONCTION VECTORIELLE DE BASE
C ASSOCIEE A L INCONNUE VARIATIONNELLE:GRADIENT DE LA TEMPERATURE
C U(J):J IEME FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE
C LA TEMPERATURE
C
C TABLEAU PFB:
C VALEURS AUX POINTS D INTEGRATION DES FONCTIONS (DANS L ORDRE):
C (1-Y)*(1-Y),Y*(1-Y),Y*Y,(1-X)*(1-Y),X*(1-Y),Y*(1-X),X*Y
C (1-X)*(1-X),X*(1-X),X*X
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: DUPUY Jean Marie ,LAB:LAN189, PARIS,              JANVIER 1981
C ......................................................................
      DOUBLE PRECISION PP1,PP2,PP3,P1,P2,PM,PN,PL,U,A,B,C,D,E
     %                ,AA1,AA2,PFB
      DIMENSION PP1(12,4),PP2(12,4),PP3(12,4),P1(6,12),P2(6,12)
     %,PM(9),PN(9),PL(9),U(4,4),AA1(9,12),AA2(9,12),E(12,4),PFB(9,10)
      DATA U/
     %1.86602540378444D0,-2.36602540378444D0,-2.36602540378444D0,
     %0.300000000000000D+01,
     %-0.499999999999999D0,0.236602540378444D+01,0.633974596215560D0,
     %-0.300000000000000D+01,
     %0.133974596215561D0,-0.633974596215560D0,-0.633974596215560D+0,
     %0.300000000000000D+01,
     %-0.499999999999999D0,0.633974596215560D0,0.236602540378444D+01,
     %-0.300000000000000D+01/
       DATA PP1/
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %-1.36602540378444D0,1.73205080756888D0,5.46410161513776D0,
     %-6.92820323027551D0,-4.09807621135332D0,5.19615242270663D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.366025403784438D0,-1.73205080756888D0,-1.46410161513776D0,
     %6.92820323027551D0,1.09807621135332D0,-5.19615242270663D0,
     %0.D0,-2.73205080756888D0,0.D0,
     %3.46410161513776D0,4.09807621135332D0,-5.19615242270663D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.732050807568877D0,0.D0,
     %-3.46410161513775D0,-1.09807621135332D0,5.19615242270663D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0/
      DATA PP2/
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.732050807568878D0,
     %-3.46410161513776D0,-1.09807621135332D0,5.19615242270663D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,-2.73205080756888D0,
     %3.46410161513776D0,4.09807621135332D0,-5.19615242270663D0,
     %0.366025403784438D0,-1.46410161513776D0,-1.73205080756888D0,
     %6.92820323027551D0,1.09807621135332D0,-5.19615242270663D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %-1.36602540378444D0,5.46410161513776D0,1.73205080756888D0,
     %-6.92820323027551D0,-4.09807621135332D0,5.19615242270663D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0/
       DATA PP3/
     %0.D0,24.D0,0.D0,
     %-36.D0,-24.D0,36.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,24.D0,
     %-36.D0,-24.D0,36.D0,
     %0.D0,-36.D0,0.D0,
     %72.D0,36.D0,-72.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,0.D0,
     %0.D0,0.D0,-36.D0,
     %72.D0,36.D0,-72.D0/
      DATA PM/0.5D0,0.8872983346207415D0,0.1127016653792585D0
     %,0.5D0,0.5D0,0.8872983346207415D0
     %,0.1127016653792585D0,0.8872983346207415D0,0.1127016653792585D0/
      DATA PN/0.5D0,0.5D0,0.5D0,
     %0.8872983346207415D0,0.1127016653792585D0,0.8872983346207415D0
     %,0.1127016653792585D0,0.1127016653792585D0,0.8872983346207415D0/
      DATA PL/0.1975308641975309D0,0.1234567901234568D0,
     %0.1234567901234568D0,0.1234567901234568D0,0.1234567901234568D0
     %,0.07716049382716049D0,0.07716049382716049D0,0.07716049382716049D0
     %,0.07716049382716049D0/
C
      DO 1 J=1,4
         J1=J+4
         J2=J+8
            DO 2 I=1,6
               II=I+6
               P1(I,J)=PP1(I,J)
               P2(I,J)=PP1(II,J)
               P1(I,J1)=PP2(I,J)
               P2(I,J1)=PP2(II,J)
               P1(I,J2)=PP3(I,J)
               P2(I,J2)=PP3(II,J)
    2       ENDDO
    1  ENDDO
      DO 3 I=1,12
           DO 4 K=1,9
               AA1(K,I)=P1(1,I)+P1(2,I)*PM(K)
     %        +P1(3,I)*PN(K)+P1(4,I)*PM(K)*PN(K)
     %        +P1(5,I)*PM(K)*PM(K)
     %        +P1(6,I)*PM(K)*PM(K)*PN(K)
               AA2(K,I)=P2(1,I)+P2(2,I)*PM(K)
     %        +P2(3,I)*PN(K)+P2(4,I)*PM(K)*PN(K)
     %        +P2(5,I)*PN(K)*PN(K)
     %        +P2(6,I)*PN(K)*PM(K)*PN(K)
    4      ENDDO
    3 ENDDO
      DO 9 I=1,12
      A=P1(2,I)+P2(3,I)
      B=2.D0*P1(5,I)+P2(4,I)
      C=P1(4,I)+2.D0*P2(5,I)
      D=2.D0*(P1(6,I)+P2(6,I))
C
      DO10 J=1,4
      E(I,J)=A*U(1,J)+(U(1,J)*(B+C)+A*(U(2,J)+U(3,J)))*0.5D0
     %+(U(1,J)*D+U(2,J)*C+U(3,J)*B+U(4,J)*A)*0.25D0
     %+(U(2,J)*B+U(3,J)*C)*0.3333333333333333D0
     %+(U(2,J)*D+U(3,J)*D+(B+C)*U(4,J))*0.1666666666666666D0
     %+D*U(4,J)*0.1111111111111111D0
   10 ENDDO
    9 ENDDO
      DO20 I=1,9
      PFB(I,1)=PL(I)*(1.D0-PN(I))*(1.D0-PN(I))
      PFB(I,2)=PL(I)*PN(I)*(1.D0-PN(I))
      PFB(I,3)=PL(I)*PN(I)*PN(I)
      PFB(I,4)=PL(I)*(1.D0-PM(I))*(1.D0-PN(I))
      PFB(I,5)=PL(I)*PM(I)*(1.D0-PN(I))
      PFB(I,6)=PL(I)*PN(I)*(1.D0-PM(I))
      PFB(I,7)=PL(I)*PM(I)*PN(I)
      PFB(I,8)=PL(I)*(1.D0-PM(I))*(1.D0-PM(I))
      PFB(I,9)=PL(I)*PM(I)*(1.D0-PM(I))
      PFB(I,10)=PL(I)*PM(I)*PM(I)
   20 ENDDO
      RETURN
      END


      SUBROUTINE TMXT21(AA,BB,CC,E)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULE DES TABLEAUX INTERMEDIAIRES POUR L ELEMENT
C ----- TRIA MT21
C
C LE TABLEAU AA CONTIENT LES VALEURS DES INTEGRALES,SUR LE TRIANGLE
C UNITE,DE P1(I)*P1(J);P1(I),ETANT LA PREMIERE COMPOSANTE DE LA I-IEME
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE: GRADIENT
C DE LA TEMPERATURE;P1(J),ETANT LA PREMIERE COMPOSANTE DE LA J-IEME
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE: GRADIENT
C DE LA TEMPERATURE
C
C LE TABLEAU BB CONTIENT LES VALEURS DES INTEGRALES,SUR LE TRIANGLE
C UNITE,DE P1(I)*P2(J)+P2(I)*P1(J);P1(I),ETANT LA PREMIERE COMPOSANTE DE
C LA I-IEME FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE :
C GRADIENT DE LA TEMPERATURE;P2(J),ETANT LA DEUXIEME COMPOSANTE DE LA
C J-IEME FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE: GRADIENT
C DE LA TEMPERATURE;P1(J),ETANT LA PREMIERE COMPOSANTE DE LA J-IEME
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE: GRADIENT
C DE LA TEMPERATURE;P2(I),ETANT LA DEUXIEME COMPOSANTE DE LA I-IEME
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE :GRADIENT
C DE LA TEMPERATURE
C
C LE TABLEAU CC CONTIENT LES VALEURS DES INTEGRALES,SUR LE TRIANGLE
C UNITE,DE P2(I)*P2(J);P2(I),ETANT LA DEUXIEME COMPOSANTE DE LA I-IEME
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE: GRADIENT
C DE LA TEMPERATURE;P2(J),ETANT LA DEUXIEME COMPOSANTE DE LA J-IEME
C FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE: GRADIENT
C DE LA TEMPERATURE
C
C TABLEAU E:TABLEAU E(I,J);VALEUR DE L INTEGRALE SUR LE TRIANGLE DE
C REFERENCE DE U(J)*DIVP(I);P(I):I IEME FONCTION VECTORIELLE DE BASE
C ASSOCIEE A L INCONNUE VARIATIONNELLE:GRADIENT DE LA TEMPERATURE
C U(J):J JEME FONCTION DE BASE ASSOCIEE A L INCONNUE VARIATIONNELLE
C LA TEMPERATURE
C POUR J VARIANT DE 1 A 3 ET I DE 1 A 6
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR:DUPUY J M ,LAB:LAN189,PARIS,JANVIER 1981
C ......................................................................
      DOUBLE PRECISION AA,BB,CC,A,E,P1,P2,P,D,B,B0,B1,B2
      DIMENSION AA(36),BB(36),CC(36),A(36),B(36),P1(4,8),P2(4,8)
     %,P(8,4),D(8,4),E(6,3),B0(6),B1(6),B2(6)
      DATA P/
     % 0.D0       , 1.098076211353315D0, 0.D0       ,
     %-1.464101615137753D0,-1.366025403784438D0, 1.732050807568877D0,
     % 5.366025403784440D0,-4.000000000000002D0,
     % 0.D0       ,-4.098076211353314D0, 0.D0       ,
     % 5.464101615137752D0, 0.3660254037844385D0,-1.732050807568877D0,
     % 3.633974596215561D0,-3.999999999999999D0,
     % 0.D0       ,-4.098076211353316D0, 0.D0       ,
     % 5.464101615137755D0, 0.D0       , 0.D0       ,
     % 1.098076211353315D0,-1.464101615137753D0,
     % 0.D0       , 1.098076211353312D0, 0.D0       ,
     %-1.464101615137751D0, 0.D0       , 0.D0       ,
     %-4.098076211353315D0, 5.464101615137753D0/
      DATA D/
     % 0.3660254037844392D0, 3.633974596215557D0,-1.732050807568877D0,
     %-3.999999999999996D0, 0.D0       , 0.D0       ,
     %-4.098076211353315D0, 5.464101615137753D0,
     %-1.366025403784438D0, 5.366025403784439D0, 1.732050807568877D0,
     %-4.000000000000002D0, 0.D0       , 0.D0       ,
     % 1.098076211353315D0,-1.464101615137753D0,
     % 0.D0       , 23.99999999999999D0, 0.D0       ,
     %-24.00000000000000D0, 0.D0       , 0.D0       ,
     % 0.D0       , 0.D0       ,
     % 0.D0       , 0.D0       , 0.D0       ,
     % 0.D0       , 0.D0       , 0.D0       ,
     % 23.99999999999999D0,-24.00000000000000D0/
      DO1 J=1,4
      JJ=J+4
              DO2 I=1,4
      P1(I,J)=P(I,J)
      P1(I,JJ)=D(I,J)
    2         ENDDO
              DO3 I=5,8
      II=I-4
      P2(II,J)=P(I,J)
      P2(II,JJ)=D(I,J)
    3         ENDDO
    1 ENDDO
      DO4 I=1,8
               DO5 J=I,8
      II=J*(J-1)/2+I
      AA(II)=P1(1,J)*P1(1,I)*0.5D0+
     %(P1(1,I)*P1(2,J)+P1(2,I)*P1(1,J)+P1(1,I)*P1(3,J)+P1(3,I)*
     %P1(1,J))*0.1666666666666666D0
     %+(P1(1,I)*P1(4,J)+P1(4,I)*P1(1,J))*0.125D0
     %+(P1(2,I)*P1(2,J)+P1(3,I)*P1(3,J))*0.0833333333333333D0
     %+(P1(2,I)*P1(3,J)+P1(3,I)*P1(2,J))*0.0416666666666666D0
     %+(P1(2,I)*P1(4,J)+P1(4,I)*P1(2,J))*0.0666666666666666D0
     %+(P1(3,I)*P1(4,J)+P1(4,I)*P1(3,J))*0.0333333333333333D0
     %+P1(4,I)*P1(4,J)*0.0555555555555555D0
C
C
      CC(II)=P2(1,J)*P2(1,I)*0.5D0+
     %(P2(1,I)*P2(3,J)+P2(3,I)*P2(1,J)+P2(1,I)*P2(2,J)+P2(2,I)*
     %P2(1,J))*0.1666666666666666D0
     %+(P2(1,I)*P2(4,J)+P2(4,I)*P2(1,J))*0.125D0
     %+(P2(3,I)*P2(3,J)+P2(2,I)*P2(2,J))*0.0833333333333333D0
     %+(P2(3,I)*P2(2,J)+P2(2,I)*P2(3,J))*0.0416666666666666D0
     %+(P2(3,I)*P2(4,J)+P2(4,I)*P2(3,J))*0.0666666666666666D0
     %+(P2(2,I)*P2(4,J)+P2(4,I)*P2(2,J))*0.0333333333333333D0
     %+P2(4,I)*P2(4,J)*0.0555555555555555D0
C
C
      A(II)=P1(1,I)*P2(1,J)*0.5D0+
     %(P1(1,I)*P2(2,J)+P1(2,I)*P2(1,J)+P1(1,I)*P2(3,J)+P1(3,I)
     %*P2(1,J))*0.1666666666666666D0
     %+(P1(1,I)*P2(4,J)+P1(4,I)*P2(1,J))*0.125D0
     %+(P1(2,I)*P2(3,J)+P1(3,I)*P2(2,J))*0.0416666666666666D0
     %+(P1(2,I)*P2(2,J)+P1(3,I)*P2(3,J))*0.0833333333333333D0
     %+(P1(2,I)*P2(4,J)+P1(4,I)*P2(3,J))*0.0333333333333333D0
     %+(P1(3,I)*P2(4,J)+P1(4,I)*P2(2,J))*0.0666666666666666D0
     %+P1(4,I)*P2(4,J)*0.0277777777777777D0
C
C
      B(II)=P1(1,J)*P2(1,I)*0.5D0+
     %(P1(1,J)*P2(2,I)+P1(2,J)*P2(1,I)+P1(1,J)*P2(3,I)+P1(3,J)
     %*P2(1,I))*0.1666666666666666D0
     %+(P1(1,J)*P2(4,I)+P1(4,J)*P2(1,I))*0.125D0
     %+(P1(2,J)*P2(3,I)+P1(3,J)*P2(2,I))*0.0416666666666666D0
     %+(P1(2,J)*P2(2,I)+P1(3,J)*P2(3,I))*0.0833333333333333D0
     %+(P1(2,J)*P2(4,I)+P1(4,J)*P2(3,I))*0.0333333333333333D0
     %+(P1(3,J)*P2(4,I)+P1(4,J)*P2(2,I))*0.0666666666666666D0
     %+P1(4,J)*P2(4,I)*0.0277777777777777D0
C
      BB(II)=A(II)+B(II)
    5 ENDDO
    4 ENDDO
      DO20 I=1,6
      B0(I)=P1(2,I)+P2(3,I)
      B1(I)=2.D0*P1(4,I)+P2(4,I)
      B2(I)=P1(4,I)+2.D0*P2(4,I)
      E(I,1)=B0(I)*0.1666666666666666D0+B1(I)*0.0833333333333333D0
      E(I,2)=B0(I)*0.1666666666666666D0+(B1(I)+B2(I))
     %*0.0833333333333333D0
      E(I,3)=B0(I)*0.1666666666666666D0+B2(I)*0.0833333333333333D0
   20 ENDDO
      RETURN
      END


      SUBROUTINE TRPOBA( NTABL, IATABL, LTABL, NFPOBA, MOPAGE, NOPAGE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFERT SUR LE FICHIER D ACCES DIRECT NFPOBA DU CONTENU DU
C ----- TABLEAU DE NOM NTABL, D ADRESSE IATABL DANS LE SUPER-TABLEAU MCN
C       DE LTABL MOTS A PARTIR DE LA NOPAGE+1 -EME PAGE , ELLE-MEME
C       DE MOPAGE MOTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1989
C.......................................................................
      include"./incl/pppoba.inc"
      COMMON MCN(MOTMCN)
      CHARACTER*(*)   NTABL
      INTEGER         ITA(4,127)
      COMMON /ELFINI/ NELFI(512)
      EQUIVALENCE    (NELFI(1),NBTAST),(NELFI(2),NODEPA),
     %               (NELFI(3),ITA(1,1))
C
C     NO DU TABLEAU SUR POBA ET DANS LE COMMON /ELFINI/
      NBTAST = NBTAST + 1
C
C     SON NOM DE 4 CARACTERES EST CODE DANS UN ENTIER
      ITA(1,NBTAST) = ICHARX( NTABL )
C
C     SON ADRESSE DANS MCN
      ITA(2,NBTAST) = IATABL
C
C     SON NOMBRE DE MOTS
      ITA(3,NBTAST) = LTABL
C
C     SON NO DE 1-ERE PAGE SUR LE FICHIER POBA
      ITA(4,NBTAST) = NOPAGE + 1
C
C     LE TRANSFERT EFFECTIF SUR LE FICHIER POBA
      CALL ESTASF(IATABL,LTABL,-1,NFPOBA,MOPAGE,NOPAGE)
      END


      SUBROUTINE VAPOR1( NBPOLY, NPI,    N1, POLY, COORPO,
     %                   DPX,    DDPX,
     %                   VAPOLY, DAPOLY, DDPOLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DES NBPOLY POLYNOMES A UNE VARIABLE AUX
C ----- NPI POINTS D INTEGRATION DE COORDONNEE COORPO
C
C PARAMETRES D ENTREE :
C ---------------------
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS DE L ESPACE REEL
C N1     : DEGRE + 1 DES POLYNOMES POLY
C POLY   : POLY(I+1,J) COEFFICIENT DE X ** I    POUR I=0,...,N1-1
C          DU J-EME POLYNOME
C COORPO : COORPO(J) COORDONNEE DU J-EME POINT
C DPX    : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1 VARIABLES)
C DDPX   : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1 VARIABLES)
C
C PARAMETRE RESULTAT :
C --------------------
C VAPOLY : VAPOLY(I,J) VALEUR DU I-EME POLYNOME AU J-EME POINT
C DAPOLY : DAPOLY(I,J) VALEUR DE LA DERIVEE DU I-EME POLYNOME AU
C                        J-EME POINT
C DDPOLY : DAPOLY(I,J) VALEUR DE LA DERIVEE DU I-EME POLYNOME AU
C                        J-EME POINT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1989
C ......................................................................
      DOUBLE PRECISION POLY(N1,NBPOLY),COORPO(NPI),VAPOLY(NBPOLY,NPI),
     %                 DAPOLY(NBPOLY,NPI),DDPOLY(NBPOLY,NPI),
     %                 DPX(N1),DDPX(N1)
C
      DO 10 I=1,NBPOLY
         CALL PN1DDE( N1, POLY(1,I), DPX  )
         CALL PN1DDE( N1, DPX      , DDPX )
         DO 1 J=1,NPI
            CALL PN1DVA( N1, POLY(1,I), COORPO(J), VAPOLY(I,J) )
            CALL PN1DVA( N1, DPX(1)   , COORPO(J), DAPOLY(I,J) )
            CALL PN1DVA( N1, DDPX(1)  , COORPO(J), DDPOLY(I,J) )
    1    ENDDO
   10 ENDDO
      END


      SUBROUTINE VAPOR2( NBPOLY, NPI,    N1,    POLY,  COORPO,
     %                   DPX,    DPY,    DDPXX, DDPYX, DDPYY,
     %                   VAPOLY, DAPOLY, DDPOLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DES NBPOLY POLYNOMES A 2 VARIABLES AUX
C ----- NPI POINTS D INTEGRATION DE COORDONNEES COORPO
C
C PARAMETRES D ENTREE :
C ---------------------
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS DE L ESPACE REEL
C N1     : DEGRE + 1 DES POLYNOMES POLY
C POLY   : POLY(I+1,J+1,K)  COEFFICIENT DE  X**I  Y**J  I,J=0,...,N1-1
C          DU K-EME POLYNOME
C COORPO : COORPO(I,J) I-EME COORDONNEE DU J-EME POINT
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C DPX    : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DPY    : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPXX  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPYX  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPYY  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C
C PARAMETRE RESULTAT :
C --------------------
C VAPOLY : VAPOLY(I,J) VALEUR DU I-EME POLYNOME AU J-EME POINT
C DAPOLY : DAPOLY(K,I,J) VALEUR DE LA DERIVEE K-EME DU I-EME POLYNOME AU
C                        J-EME POINT
C DDPOLY : DDPOLY(K,I,J) VALEUR DE LA DERIVEE SECONDE K-EME DU I-EME POLYNOME AU
C                        J-EME POINT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C ......................................................................
      DOUBLE PRECISION COORPO(2,NPI),VAPOLY(NBPOLY,NPI),
     %                 POLY(N1,N1,NBPOLY),
     %                 DAPOLY(2,NBPOLY,NPI),
     %                 DDPOLY(3,NBPOLY,NPI),
     %                 DPX(1),DPY(1), DDPXX(1),DDPYX(1),DDPYY(1),
     %                 X,Y
C
       DO 12 I=1,NBPOLY
         CALL PN2DDE( 1, N1, POLY(1,1,I), DPX )
         CALL PN2DDE( 2, N1, POLY(1,1,I), DPY )
C        D/X  DE D/X P
         CALL PN2DDE( 1, N1,DPX, DDPXX)
C        D/Y  DE D/X P
         CALL PN2DDE( 2, N1,DPX, DDPYX)
C        D/Y  DE D/Y P
         CALL PN2DDE( 2, N1,DPY, DDPYY)
C
         DO 2 J=1,NPI
            X = COORPO(1,J)
            Y = COORPO(2,J)
C
            CALL PN2DVA( N1,POLY(1,1,I),X,Y, VAPOLY(I,J) )
C
            CALL PN2DVA( N1,DPX,X,Y, DAPOLY(1,I,J) )
            CALL PN2DVA( N1,DPY,X,Y, DAPOLY(2,I,J) )
C
            CALL PN2DVA( N1,DDPXX,X,Y, DDPOLY(1,I,J) )
            CALL PN2DVA( N1,DDPYX,X,Y, DDPOLY(2,I,J) )
            CALL PN2DVA( N1,DDPYY,X,Y, DDPOLY(3,I,J) )
    2    ENDDO
C
   12 ENDDO
      END


      SUBROUTINE VAPOR3(NBPOLY, NPI, N1, N3, POLY, COORPO,
     %                  DPX,DPY,DPZ,DDPXX,DDPYX,DDPYY,DDPZX,DDPZY,DDPZZ,
     %                  VAPOLY, DAPOLY, DDPOLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DES NBPOLY POLYNOMES A 3 VARIABLES AUX
C ----- NPI POINTS D INTEGRATION DE COORDONNEES COORPO
C
C PARAMETRES D ENTREE :
C ---------------------
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS DE L ESPACE REEL
C N1     : DEGRE + 1 DES POLYNOMES POLY
C N3     : N1 ** 3
C POLY   : POLY(I+1,J+1,K)  COEFFICIENT DE  X**I  Y**J  I,J=0,...,N1-1
C          DU K-EME POLYNOME
C COORPO : COORPO(I,J) I-EME COORDONNEE DU J-EME POINT
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C DPX    : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DPY    : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DPZ    : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPXX  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPYX  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPYY  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPZX  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPZY  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C DDPZZ  : TABLEAU AUXILIAIRE DOUBLE PRECISION (N1,N1)
C
C PARAMETRE RESULTAT :
C --------------------
C VAPOLY : VAPOLY(I,J) VALEUR DU I-EME POLYNOME AU J-EME POINT
C DAPOLY : DAPOLY(K,I,J) VALEUR DE LA DERIVEE K-EME DU I-EME POLYNOME AU
C                        J-EME POINT
C DDPOLY : DDPOLY(K,I,J) VALEUR DE LA DERIVEE SECONDE K-EME DU I-EME
C                        POLYNOME AU J-EME POINT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1995
C ......................................................................
      DOUBLE PRECISION COORPO(3,NPI),VAPOLY(NBPOLY,NPI),
     %                 POLY(N3,NBPOLY),
     %                 DAPOLY(3,NBPOLY,NPI),
     %                 DDPOLY(6,NBPOLY,NPI),
     %                 DPX(1),DPY(1),DPZ(1),
     %                 DDPXX(1),DDPYX(1),DDPYY(1),DDPZX(1),DDPZY(1),
     %                 DDPZZ(1),
     %                 X,Y,Z
C
      DO 12 I=1,NBPOLY
C
         CALL PN3DDE( 1,N1,POLY(1,I), DPX )
         CALL PN3DDE( 2,N1,POLY(1,I), DPY )
         CALL PN3DDE( 3,N1,POLY(1,I), DPZ )
C
C        D/X  DE D/X P
         CALL PN3DDE( 1, N1,DPX, DDPXX )
C        D/Y  DE D/X P
         CALL PN3DDE( 2, N1,DPX, DDPYX )
C        D/Y  DE D/Y P
         CALL PN3DDE( 2, N1,DPY, DDPYY )
C        D/Z  DE D/X P
         CALL PN3DDE( 3, N1,DPX, DDPZX )
C        D/Z  DE D/X P
         CALL PN3DDE( 3, N1,DPY, DDPZY )
C        D/Z  DE D/Y P
         CALL PN3DDE( 3, N1,DPZ, DDPZZ )
C
         DO 2 J=1,NPI
            X = COORPO(1,J)
            Y = COORPO(2,J)
            Z = COORPO(3,J)
C           VALEUR
            CALL PN3DVA( N1,POLY(1,I),X,Y,Z, VAPOLY(I,J) )
C           VALEUR DE LA DERIVEE PREMIERE
            CALL PN3DVA( N1,DPX,X,Y,Z, DAPOLY(1,I,J) )
            CALL PN3DVA( N1,DPY,X,Y,Z, DAPOLY(2,I,J) )
            CALL PN3DVA( N1,DPZ,X,Y,Z, DAPOLY(3,I,J) )
C           VALEUR DE LA DERIVEE SECONDE
            CALL PN3DVA( N1,DDPXX,X,Y,Z, DDPOLY(1,I,J) )
            CALL PN3DVA( N1,DDPYX,X,Y,Z, DDPOLY(2,I,J) )
            CALL PN3DVA( N1,DDPYY,X,Y,Z, DDPOLY(3,I,J) )
            CALL PN3DVA( N1,DDPZX,X,Y,Z, DDPOLY(4,I,J) )
            CALL PN3DVA( N1,DDPZY,X,Y,Z, DDPOLY(5,I,J) )
            CALL PN3DVA( N1,DDPZZ,X,Y,Z, DDPOLY(6,I,J) )
    2    ENDDO
C
   12 ENDDO

      RETURN
      END


      SUBROUTINE VAPOR6( NBPOLY, NPI, N1, N6, POLY, COORPO,
     %                   DPX, DPY, DPZ, DPU, DPV, DPW,
     %                   VAPOLY, DAPOLY )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DES NBPOLY POLYNOMES A 6 VARIABLES AUX
C ----- NPI POINTS D INTEGRATION DE COORDONNEES COORPO
C
C ENTREES:
C --------
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D'INTEGRATION
C N1     : DEGRE + 1 DES POLYNOMES POLY
C N6     : N1 ** 6
C POLY   : POLY(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY(I,J,K,L,M,N,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                             U**(L-1) V**(M-1) W**(N-1)
C COORPO : COORPO(I,J) I-EME COORDONNEE DU J-EME POINT
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C DPX    : TABLEAU AUXILIAIRE DOUBLE PRECISION (NBPOLY)
C DPY    : TABLEAU AUXILIAIRE DOUBLE PRECISION (NBPOLY)
C DPZ    : TABLEAU AUXILIAIRE DOUBLE PRECISION (NBPOLY)
C DPU    : TABLEAU AUXILIAIRE DOUBLE PRECISION (NBPOLY)
C DPV    : TABLEAU AUXILIAIRE DOUBLE PRECISION (NBPOLY)
C DPW    : TABLEAU AUXILIAIRE DOUBLE PRECISION (NBPOLY)
C
C SORTIES:
C --------
C VAPOLY : VAPOLY(I,J)   VALEUR DU I-EME POLYNOME AU J-EME POINT
C DAPOLY : DAPOLY(K,I,J) VALEUR DE LA DERIVEE /XK-EME DU I-EME POLYNOME
C                        AU J-EME POINT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C ......................................................................
      DOUBLE PRECISION COORPO(6,NPI),
     %                 POLY( N6,NBPOLY),
     %                 VAPOLY(  NBPOLY,NPI),
     %                 DAPOLY(6,NBPOLY,NPI),
     %                 DPX(NBPOLY), DPY(NBPOLY), DPZ(NBPOLY),
     %                 DPU(NBPOLY), DPV(NBPOLY), DPW(NBPOLY),
     %                 X, Y, Z, U, V, W
C
      DO 20 I=1,NBPOLY
C
C        LES COEFFICIENTS DES POLYNOMES DERIVES UNE FOIS
         CALL PN6DDE( 1, N1, POLY(1,I), DPX )
         CALL PN6DDE( 2, N1, POLY(1,I), DPY )
         CALL PN6DDE( 3, N1, POLY(1,I), DPZ )
         CALL PN6DDE( 4, N1, POLY(1,I), DPU )
         CALL PN6DDE( 5, N1, POLY(1,I), DPV )
         CALL PN6DDE( 6, N1, POLY(1,I), DPW )
C
         DO 10 J=1,NPI
C
C           LES 6 COORDONNEES DU POINT J D'INTEGRATION
            X = COORPO(1,J)
            Y = COORPO(2,J)
            Z = COORPO(3,J)
            U = COORPO(4,J)
            V = COORPO(5,J)
            W = COORPO(6,J)
C
C           VALEUR DU POLYNOME I AU POINT D'INTEGRATION J
            CALL PN6DVA( N1, POLY(1,I), X,Y,Z,U,V,W, VAPOLY(I,J) )
C
C           VALEUR DES 6 DERIVEES PREMIERES
            CALL PN6DVA( N1, DPX, X,Y,Z,U,V,W, DAPOLY(1,I,J) )
            CALL PN6DVA( N1, DPY, X,Y,Z,U,V,W, DAPOLY(2,I,J) )
            CALL PN6DVA( N1, DPZ, X,Y,Z,U,V,W, DAPOLY(3,I,J) )
            CALL PN6DVA( N1, DPU, X,Y,Z,U,V,W, DAPOLY(4,I,J) )
            CALL PN6DVA( N1, DPV, X,Y,Z,U,V,W, DAPOLY(5,I,J) )
            CALL PN6DVA( N1, DPW, X,Y,Z,U,V,W, DAPOLY(6,I,J) )

 10      ENDDO

 20   ENDDO

      RETURN
      END

      SUBROUTINE VALXYZ( NBPI, POIDPI, NBCOOR, COORPI, VALX2Y2Z2 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA VALEUR DE L'INTEGRALE DE X2 Y2 Z2 SUR UN EF 3D
C ----- A PARTIR DES NBPI POINTS D INTEGRATION DE COORDONNEES COORPI

C ENTREES:
C --------
C NBPI   : NOMBRE DE POINTS D'INTEGRATION
C POIDPI : POIDS                  DE LA FORMULE D'INTEGRATION NUMERIQUE
C COORPI : COORDONNEES DES POINTS DE LA FORMULE D'INTEGRATION NUMERIQUE

C SORTIE :
C --------
C VOLUME : APPROCHE DE L'EF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L LIONS UPMC Paris Octobre 2007
C ......................................................................
      DOUBLE PRECISION  POIDPI(NBPI), COORPI(NBCOOR,NBPI), VALX2Y2Z2

      VALX2Y2Z2 = 0D0
      DO L=1,NBPI
         VALX2Y2Z2 = VALX2Y2Z2 + POIDPI(L)
     %                         * COORPI(1,L)**2
     %                         * COORPI(2,L)**2
     %                         * COORPI(3,L)**2
      ENDDO

      RETURN
      END
