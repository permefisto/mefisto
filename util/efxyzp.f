      SUBROUTINE EFXYZP( NDIM, MNXYZP, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                   COOPEF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES NDIM COORDONNEES DES NBPOE POINTS DE L'EF
C -----    NUELEM PARMI LES NBELEM EF DE CE TYPE
C
C ENTREES:
C --------
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 1 ou 2 ou 3 ou 6
C          (SI AXISYMETRIE NDIM=2 et X=>R>=0 Y=>Z Z=0)
C MNXYZP : ADRESSE MCN DU TMS XYZPOINT DE L'OBJET
C NBELEM : NOMBRE D'EF DE CE TYPE D'EF
C NUELEM : NUMERO DE L'EF A TRAITER
C MNPGEL : ADRESSE MCN DES NUMEROS DES POINTS GEOMETRIQUES DES EF
C NBPOE  : NOMBRE DE POINTS DE L'EF
C
C SORTIE :
C --------
C COOPEF : LES 3 ou 6 COORDONNEES DES POINTS DE l'EF NUELEM
C          COOPEF(I,J)=J-EME COORDONNEE DU I-EME POINT DE L'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS     Fevrier 1999
C MODIFS: ALAIN PERRONNET  Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      include"./incl/a___xyzpoint.inc"
C
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))
C
      REAL           COOPEF( NBPOE, NDIM )
C
C     NBC NOMBRE DE COORDONNEES D'UN POINT DU TABLEAU XYZPOI
      IF( NDIM .LE. 3 ) THEN
         NBC = 3
      ELSE
         NBC = NDIM
      ENDIF
C
C     ADRESSES MCN DE DEBUT DE TABLEAUX
      MNP = MNPGEL - 1 + NUELEM - NBELEM
      MNC = MNXYZP + WYZPOI - NBC
C
C     BOUCLE SUR LES POINTS DE CET EF NUELEM
      DO 10 I=1,NBPOE
C
C        NO NUMERO DU I-EME POINT DE L'ELEMENT FINI NUELEM
         MNP = MNP + NBELEM
         NOP = MCN( MNP )
C
C        ADRESSE DE LA PREMIERE COORDONNEE DU POINT I DANS XYZPOINT
         MNCE = MNC + NBC * NOP
C
C        LA PREMIERE COORDONNEE
         COOPEF( I, 1 ) = RMCN( MNCE )
C
C        LA SECONDE COORDONNEE
         IF( NDIM .GE. 2 ) COOPEF( I, 2 ) = RMCN( MNCE + 1 )
C
C        LA TROISIEME COORDONNEE
         IF( NDIM .GE. 3 ) COOPEF( I, 3 ) = RMCN( MNCE + 2 )
C
C        CAS DE LA DIMENSION 6
         IF( NDIM .EQ. 6 ) THEN
C           LES 3 DERNIERES COORDONNEES
            COOPEF( I, 4 ) = RMCN( MNCE + 3 )
            COOPEF( I, 5 ) = RMCN( MNCE + 4 )
            COOPEF( I, 6 ) = RMCN( MNCE + 5 )
         ENDIF
C
 10   CONTINUE
C
      RETURN
      END
