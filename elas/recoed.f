      SUBROUTINE RECOED( NYOBJT, NUOBJT, NDIM, XPI,YPI,ZPI, MNCOED,
     %                   COEFDE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU COEDEP DES COEFFICIENTS DU DEPLACEMENT
C -----
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NDIM   : 2 OU 3 DIMENSION DE L'ESPACE  ( 2 SI PB AXISYMETRIQUE )
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNCOED : ADRESSE MCN DU TABLEAU 'COEFFICIENT DEPLACEMENT'
C
C SORTIES:
C --------
C COEFDE : TABLEAU (1:NDIM) DU COEFFICIENT DU DEPLACEMENT
C
C ATTENTION: SEULES NBCOED<=NDIM COMPOSANTES DU COEFFICIENT A DONNER
C ---------- MAIS RESULTAT SUR NDIM COMPOSANTES!
C            LES AUTRES VALEURS SONT MISES A ZERO
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1999
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___coefdeplacement.inc"
      include"./incl/ctemps.inc"
      DOUBLE PRECISION  XPI,YPI,ZPI,XYZ(7)
      DOUBLE PRECISION  COEFDE(NDIM)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     MISE A ZERO DES NDIM COMPOSANTES DU COEFFICIENT
      DO 5 J=1,NDIM
         COEFDE( J ) = 0D0
 5    CONTINUE
C
C     LE NOMBRE DE COMPOSANTES DE COEFFICIENTS A DEFINIR
      NBCOED = MCN( MNCOED + WBCOED )
C
C     REMPLISSAGE SELON LE TYPE DES DONNEES DES COEFFICIENTS
      MN = MNCOED + WUCOED + NBCOED
      IF( MCN( MNCOED + WTCOED ) .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CET OBJET
C        =================================
         DO 10 J=1,NBCOED
C           LA VALEUR DE LA COMPOSANTE DU COEFFICIENT EST EMPILEE
            NOCOMP = MCN( MNCOED + WUCOED + J - 1 )
            COEFDE( NOCOMP ) = RMCN( MN - 1 + J )
 10      CONTINUE
C
      ELSE
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XPI
         XYZ(3) = YPI
         XYZ(4) = ZPI
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         DO 20 J=1,NBCOED
C           LE NUMERO DE LA COMPOSANTE
            NOCOMP = MCN( MNCOED + WUCOED + J - 1 )
            XYZ(7) = NOCOMP
            CALL FONVAL( MCN(MN), 7, XYZ,
     %                   NCODEV, COEFDE(NOCOMP) )
 20      CONTINUE
      ENDIF
      END
