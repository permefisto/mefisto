      SUBROUTINE ELINTE( NOMPB , NOTYEF, NDIMF , NUINTF,
     &                   NBINVA, NUINVA, NUINTI, NBNDIN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER SELON LE PROBLEME TRAITE , LE NOM DE L ELEMENT FINI
C ----- LE NOMBRE DE COMPOSANTES DE L APPLICATION
C       ELEMENT REFERENCE --> ELEMENT COURANT
C       LE NO D INTERPOLATION COMMUN DE SES COMPOSANTES
C       LE NOMBRE D INCONNUES VARIATIONNELLES DE L ELEMENT
C       LEUR NO ET LE NO DE LEUR INTERPOLATION
C       LE NOMBRE DE NOEUDS NECESSAIRE A CHAQUE INTERPOLATION
C       EVENTUELLEMENT LA LISTE DES NOEUDS ELEMENTAIRES DE CHACUNE
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMPB  : NOM DE 4 CARACTERES DE LA CLASSE DU PROBLEME TRAITE
C NOTYEF : NUMERO DU TYPE DE L'EF
C
C PARAMETRES RESULTATS :
C ----------------------
C NDIMF  : NOMBRE DE COMPOSANTES DE L APPLICATION
C          ELEMENT REFERENCE --> ELEMENT COURANT
C NUINTF : NO DE L INTERPOLATION COMMUNE AUX NDIMF COMPOSANTES
C          ( CF SP INTERP )
C NBINVA : NOMBRE D INCONNUES VARIATIONNELLES SUR CET ELEMENT
C NUINVA : NO DE CES NBINVA INCONNUES VARIATIONNELLES
C          ( CF SP INCVAR )
C NUINTI : NO D INTERPOLATION DE CES NBINVA INCONNUES VARIATIONNELLES
C          ( CF SP INTERP )
C NBNDIN : NOMBRE DE NOEUDS NECESSAIRES A CHAQUE INTERPOLATION
C          DE CHAQUE INCONNUE VARIATIONNELLE
C          -1 SIGNIFIE QUE TOUS LES NOEUDS DE L ELEMENT SERVENT
C             ET QUE NUNOIN DE CETTE INCONNUE N EST PAS DONNE
C NUNOIN : NO DES NBNDIN NOEUDS DANS LA NUMEROTATION ELEMENTAIRE
C          NECESSAIRE A CHAQUE INCONNUE VARIATIONNELLE
C          NUNOIN(I,J) J=1,...,NBINVA . I=1,...,NBNDIN(J) SI NON = -1
C                      DANS TOUS LES CAS NBINVA<6 ET NBNDIN(I)<30
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1991
C2345X7..............................................................012
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*10      LISTE(3)
      CHARACTER*(*)     NOMPB
      CHARACTER*4       NOMELE(2)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DATA LLISTE,LISTE /3,'THERMIQUE ','ELASTICITE','FLUIDE    '/
C
C     CALCUL DE K NO DE LA CLASSE DE PROBLEME TRAITE
C     ----------------------------------------------
      DO 5 I=1,LLISTE
         K = I
         IF( NOMPB .EQ. LISTE(I) ) GOTO 8
    5 CONTINUE
C
C     CLASSE DE TRAVAIL NON RETROUVEE
C     -------------------------------
      NBLGRC(NRERR) = 1
      KERR(1) = 'ELINTE:LA CLASSE '//NOMPB//' EST INCONNUE'
      CALL LEREUR
      RETURN
C
C     ******************************************************************
  8   GOTO ( 1000 , 10 , 10 ) , K
C     ******************************************************************
 10   NBLGRC(NRERR) = 2
      CALL ELNUNM( NOTYEF, NOMELE )
      KERR(1) = 'ELINTE: ELEMENT INCONNU '
      KERR(2) = NOMELE(1) // ' '// NOMELE(2)
      CALL LEREUR
      RETURN
C
C     ==================================================================
C     LE PROBLEME THERMIQUE
C     ==================================================================
C
C     ******************************************************************
 1000 GOTO ( 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1100,
     &       1110, 1120, 1010, 1140, 1020, 1030, 1170, 1040, 1190, 1200,
     &       1210, 1220, 1230, 1240, 1140, 1250,   10, 1320, 1010, 1260,
     &       1270, 1280, 1330, 1340   ) ,NOTYEF
C     ******************************************************************
C
C     ==================================================================
C     TRIA AP1D ET TRIA 2P1D ET TRIA 2P1C
C     ==================================================================
 1010 NUINTF    = 2001
      NUINTI(1) = 2001
C
 1011 NDIMF     =  2
C
 1012 NBINVA    =  1
      NUINVA(1) =  1
      NBNDIN(1) = -1
      RETURN
C
C     ==================================================================
C     TRIA AP2C ET TRIA 2P2C
C     ==================================================================
 1020 NUINTF    = 2002
      NUINTI(1) = 2002
      GOTO 1011
C
C     ==================================================================
C     QUAD AQ1C ET QUAD 2Q1C
C     ==================================================================
 1030 NUINTF    = 2031
      NUINTI(1) = 2031
      GOTO 1011
C
C     ==================================================================
C     QUAD AQ2C ET QUAD 2Q2C
C     ==================================================================
 1040 NUINTF    = 2032
      NUINTI(1) = 2032
      GOTO 1011
C
C     ==================================================================
C     TRIA MT10
C     ==================================================================
 1050 NUINTF    = 2001
      NUINTI(1) = 2004
C
 1051 NDIMF     = 2
      NBINVA    = 1
      NUINVA(1) = 2
      NBNDIN(1) = -1
      RETURN
C
C     ==================================================================
C     TRIA MT21
C     ==================================================================
 1060 NUINTF    = 2001
      NUINTI(1) = 2005
      GOTO 1051
C
C     ==================================================================
C     QUAD MQ10
C     ==================================================================
 1070 NUINTF    = 2031
      NUINTI(1) = 2034
      GOTO 1051
C
C     ==================================================================
C     QUAD MQ21
C     ==================================================================
 1080 NUINTF    = 2031
      NUINTI(1) = 2035
      GOTO 1051
C
C     ==================================================================
C     TETR M3T1
C     ==================================================================
 1090 NUINTF    = 3001
      NUINTI(1) = 3004
C
 1091 NDIMF     = 3
      NBINVA    = 1
      NUINVA(1) = 2
      NBNDIN(1) = -1
      RETURN
C
C     ==================================================================
C     TETR M3T2
C     ==================================================================
 1100 NUINTF    = 3001
      NUINTI(1) = 3005
      GOTO 1091
C
C     ==================================================================
C     HEXA M3H1
C     ==================================================================
 1110 NUINTF    = 3061
      NUINTI(1) = 3064
      GOTO 1091
C
C     ==================================================================
C     HEXA M3H2
C     ==================================================================
 1120 NUINTF    = 3061
      NUINTI(1) = 3065
      GOTO 1091
C
C     ==================================================================
C     TRIA 2P2D    TRIA HD06
C     ==================================================================
 1140 NUINTF    = 2001
      NUINTI(1) = 2002
      GOTO 1011
C
C     ==================================================================
C     QUAD 2Q2D
C     ==================================================================
 1170 NUINTF    = 2031
      NUINTI(1) = 2032
      GOTO 1011
C
C     ==================================================================
C     TETR 3P1D
C     ==================================================================
 1190 NUINTF    = 3001
C
 1191 NDIMF     = 3
      NUINTI(1) = NUINTF
      GOTO 1012
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
 1200 NUINTF    = 3002
      GOTO  1191
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
 1210 NUINTF    = 3031
      GOTO  1191
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
 1220 NUINTF    = 3032
      GOTO  1191
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
 1230 NUINTF    = 3061
      GOTO  1191
C
C     ==================================================================
C     HEXA 3Q2C
C     ==================================================================
 1240 NUINTF    = 3062
      GOTO  1191
C
C     ==================================================================
C     TRIA EQ06
C     ==================================================================
 1250 NUINTF=2001
      NUINTI(1)=2006
      GOTO 1011
C
C     ==================================================================
C     6CUB 6Q1C
C     ==================================================================
 1260 NUINTF   = 3061
      NDIMF    = 6
      NUINTI(1)= NUINTF
      GOTO 1012
C
C     ==================================================================
C     PYRA 3PY1
C     ==================================================================
 1270 NUINTF = 3091
      GOTO  1191
C
C     ==================================================================
C     PYRA 3PY2
C     ==================================================================
 1280 NUINTF = 3092
      GOTO  1191
C
C     ==================================================================
C     SEGM 1P1D
C     ==================================================================
 1320 NUINTF    = 1001
C
 1325 NDIMF     =  1
      NUINTI(1) = NUINTF
      NBINVA    =  1
      NUINVA(1) =  1
      NBNDIN(1) = -1
      RETURN
C
C     ==================================================================
C     SEGM 1P2D
C     ==================================================================
 1330 NUINTF = 1002
      GOTO 1325
C
C     ==================================================================
C     SEGM 1P2D
C     ==================================================================
 1340 NUINTF = 1003
      GOTO 1325
      END
