      SUBROUTINE ELINT1( NOMPB,  NO,
     &                   NDIM,   NOINTF, NBPOF,
     &                   NBINVA, NOINVA, NOINTI, NBNOIN, NONOIN)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFIN DE TRAITER LES INCONNUES SECONDAIRES DES ELEMENTS
C ----- APRES RESOLUTION DU SYSTEME LINEAIRE CE SP
C       RETOURNE SELON LE PROBLEME TRAITE , LE NOM DE L ELEMENT FINI
C       LE NOMBRE DE COMPOSANTES DE L APPLICATION
C       ELEMENT REFERENCE --> ELEMENT COURANT
C       LE NO D INTERPOLATION COMMUN DE SES COMPOSANTES
C       LE NOMBRE D INCONNUES VARIATIONNELLES DE L ELEMENT FINI
C       LEUR NO ET LE NO DE LEUR INTERPOLATION
C       LE NOMBRE DE NOEUDS NECESSAIRE A CHAQUE INTERPOLATION
C       EVENTUELLEMENT LA LISTE DES NOEUDS ELEMENTAIRES DE CHACUNE
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMPB  : NOM DE LA CLASSE DU PROBLEME TRAITE
C NO     : NUMERO DU TYPE DE L'EF
C
C PARAMETRES RESULTATS :
C ----------------------
C NDIM   : NOMBRE DE COMPOSANTES DE L APPLICATION
C          ELEMENT REFERENCE --> ELEMENT COURANT
C NOINTF : NO DE L INTERPOLATION COMMUNE AUX NDIM COMPOSANTES
C          ( CF SP INTERP )
C NBPOF  : NOMBRE DE POINTS DE L'APPLICATION F
C          <0 ENTRAINE LA NON UTILISATION DE NONOIN(I,.)=I
C NBINVA : NOMBRE D INCONNUES VARIATIONNELLES SUR CET ELEMENT
C NOINVA : NO DE CES NBINVA INCONNUES VARIATIONNELLES
C          ( CF SP INCVAR )
C NOINTI : NO D INTERPOLATION DE CES NBINVA INCONNUES VARIATIONNELLES
C          ( CF SP INTERP )
C NBNOIN : NOMBRE DE NOEUDS NECESSAIRES A CHAQUE INTERPOLATION
C          DE CHAQUE INCONNUE VARIATIONNELLE
C          <0 ENTRAINE LA NON FOURNITURE DE NONOIN(I,.)=I
C          >0 NONOIN(I,J)=NO DU I-EME NOEUD DE L INTERPOLATION J
C                         0 SI LE NOEUD EST INUTILE A CETTE INTERPOLATIO
C NONOIN : NO DES NBNOIN NOEUDS DANS LA NUMEROTATION ELEMENTAIRE
C          NECESSAIRE A CHAQUE INCONNUE VARIATIONNELLE
C          NONOIN(I,J) J=1,...,NBINVA . I=1,...,NBNOIN(J) SI NON = -1
C                      DANS TOUS LES CAS NBINVA<6 ET NBNOIN(I)<30
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JANVIER 1991
C2345X7..............................................................012
      INTEGER           NOINVA(5),NOINTI(5),NBNOIN(5),NONOIN(30,5)
      CHARACTER*4       NOMPB,LISTE(3)
      CHARACTER*4       NOMELE(2)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DATA LLISTE,LISTE /3,'THER','ELAS','FLUI'/
10005 FORMAT('ERREUR ELINT1 : LA CLASSE ',A4,' NON RETROUVEE PARMI'/
     &20(1X,A4))
10010 FORMAT('0ERREUR ELINT1:PROBLEME:',A4,' ELEMENT:',2(1X,A4),
     &' NON DEFINI'/)
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
      WRITE(IMPRIM,10005) NOMPB,LISTE
      CALL ARRET(100)
C
C     ******************************************************************
    8 GOTO ( 1000 , 10 , 10 ) , K
C     ******************************************************************
C
   10 CALL ELNUNM( NO , NOMELE )
      WRITE(IMPRIM,10010) NOMPB,NOMELE
      CALL ARRET(100)
C
C     ==================================================================
C     LE PROBLEME THERMIQUE
C     ==================================================================
C     ******************************************************************
 1000 GOTO ( 1020, 1030, 1040, 1050 ,   10,   10,   10,  10,   10 ,  10,
     &        10 ,   10, 1010, 1030 , 1030, 1040, 1050, 1050, 1310,1320,
     &      1330 , 1340, 1350, 1360 , 1060, 1060,   10, 10, 1010) ,  NO
C     ******************************************************************
C
C     ==================================================================
C     TRIA 2P1D ET TRIA 2P1C
C     ==================================================================
 1010 NOINTF    = 2001
      NBNOIN(1) = -1
C
 1015 NDIM      = 2
C
 1018 NBINVA    = 1
      NOINVA(1) = 1
      NOINTI(1) = NOINTF
      NBPOF=NBNOIN(1)
      RETURN
C
C     ==================================================================
C     TRIA AP1D
C     ==================================================================
 1020 NOINTF    = 2001
      NBNOIN(1) = -3
      GOTO 1015
C
C     ==================================================================
C     TRIA AP2C ET 2P2C
C     ==================================================================
 1030 NOINTF    = 2002
      NBNOIN(1) = -6
      GOTO 1015
C
C     ==================================================================
C     QUAD AQ1C ET 2Q1C
C     ==================================================================
C
 1040 NOINTF    = 2031
      NBNOIN(1) =  4
      NONOIN(1,1) = 1
      NONOIN(2,1) = 2
      NONOIN(3,1) = 4
      NONOIN(4,1) = 3
      GOTO 1015
C
C     ==================================================================
C     QUAD AQ2C ET 2Q2C
C     ==================================================================
 1050 NOINTF    = 2032
      NBNOIN(1) = 8
      NONOIN(1,1) = 1
      NONOIN(2,1) = 3
      NONOIN(3,1) = 9
      NONOIN(4,1) = 7
      NONOIN(5,1) = 2
      NONOIN(6,1) = 6
      NONOIN(7,1) = 8
      NONOIN(8,1) = 4
      GOTO 1015
C
C     ==================================================================
C     TRIA 2P2D ET HD06 ET EQ06
C     ==================================================================
 1060 NDIM       = 2
      NOINTF     = 2001
      NBPOF      = -3
      NBINVA     = 1
      NOINVA(1)  = 1
      NOINTI(1)  = 2002
      NBNOIN(1)  = -6
      RETURN
C
C     =================================================================
C     TETR 3P1D
C     ==================================================================
 1310 NOINTF    = 3001
      NBNOIN(1) = -1
C
 1315 NDIM      = 3
      GOTO 1018
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
 1320 NOINTF    = 3002
      NBNOIN(1) = 10
      NONOIN( 1,1) =  2
      NONOIN( 2,1) =  3
      NONOIN( 3,1) =  4
      NONOIN( 4,1) =  5
      NONOIN( 5,1) = 12
      NONOIN( 6,1) = 13
      NONOIN( 7,1) = 15
      NONOIN( 8,1) = 11
      NONOIN( 9,1) = 14
      NONOIN(10,1) = 10
      GOTO 1315
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
 1330 NOINTF    = 3031
      NBNOIN(1) = -6
      GOTO 1315
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
 1340 NOINTF    = 3032
      NBNOIN(1) = 15
      NONOIN( 1,1) = 1
      NONOIN( 2,1) = 2
      NONOIN( 3,1) = 3
      NONOIN( 4,1) = 15
      NONOIN( 5,1) = 16
      NONOIN( 6,1) = 17
      NONOIN( 7,1) = 4
      NONOIN( 8,1) = 5
      NONOIN( 9,1) = 6
      NONOIN(10,1) = 8
      NONOIN(11,1) = 9
      NONOIN(12,1) = 10
      NONOIN(13,1) = 18
      NONOIN(14,1) = 19
      NONOIN(15,1) = 20
      GOTO 1315
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
 1350 NOINTF    = 3061
      NBNOIN(1) =  8
      NONOIN(1,1) = 1
      NONOIN(2,1) = 2
      NONOIN(3,1) = 4
      NONOIN(4,1) = 3
      NONOIN(5,1) = 5
      NONOIN(6,1) = 6
      NONOIN(7,1) = 8
      NONOIN(8,1) = 7
      GOTO 1315
C
C     ==================================================================
C     HEXA 3Q2C
C     ==================================================================
 1360 NOINTF    = 3062
      NBNOIN(1) = 20
      NONOIN( 1,1) =  1
      NONOIN( 2,1) =  3
      NONOIN( 3,1) =  9
      NONOIN( 4,1) =  7
      NONOIN( 5,1) = 19
      NONOIN( 6,1) = 21
      NONOIN( 7,1) = 27
      NONOIN( 8,1) = 25
      NONOIN( 9,1) =  2
      NONOIN(10,1) =  6
      NONOIN(11,1) =  8
      NONOIN(12,1) =  4
      NONOIN(13,1) = 10
      NONOIN(14,1) = 12
      NONOIN(15,1) = 18
      NONOIN(16,1) = 16
      NONOIN(17,1) = 20
      NONOIN(18,1) = 24
      NONOIN(19,1) = 26
      NONOIN(20,1) = 22
      GOTO 1315
      END
