      SUBROUTINE COBABF( NBSOM, MNNPEF, XYZBEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU BARYCENTRE DE TOUS LES EF BREZZI-FORTIN
C -----
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C MNNPEF : ADRESSE MCN DU TABLEAU NPEF"2P1D ou NPEF"3P1D DU FLUIDE
C
C MODIFIES:
C --------
C XYZBEF : AJOUT DES 3 COORDONNEES DU BARYCENTRE DES NBELEM TRIANGLES
C          OU TETRAEDRES
C          DECLARE XYZBEF( 3, NBSOM+NBELEM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray     Mars 2010
C23456---------------------------------------------------------------012
      include"./incl/a___npef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/ponoel.inc"
      include"./incl/pp.inc"
      COMMON    MCN(MOTMCN)
C
      REAL              XYZBEF(3,*)
      INTEGER           NONOEF(4)
      DOUBLE PRECISION  S
      INTRINSIC         REAL
C
C     MNELE : ADRESSE DU TABLEAU NPEF"2P1D ou NPEF"3P1D
      MNELE = MCN( MNNPEF )
C
C     NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELEM = MCN( MNELE + WBELEM )
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )
C     NUTYEL : 13 TRIANGLE 2P1D, 19 TETRAEDRE 3P1D de BREZZI-FORTIN
C     NDIM   : 2 ou 3 DIMENSION DE L'ESPACE DES EF
      IF( NUTYEL .EQ. 13 ) THEN
         NDIM = 2
      ELSE
         NDIM = 3
      ENDIF
C
C     NOMBRE DE SOMMETS DE L'EF
      NBSTEF = NDIM + 1
C
C     LA BOUCLE SUR LES ELEMENTS FINIS
C     ================================
      DO 50 NUELEM = 1, NBELEM
C
C        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C
C        XYZ DU BARYCENTRE DE L'EF NUELEM
         NUBARY = NBSOM + NUELEM
         DO 20 K=1,3
            S = 0D0
            DO 10 I = 1, NBSTEF
               S = S + XYZBEF( K, NONOEF(I) )
 10         CONTINUE
C           LA COORDONNEE K DU BARYCENTRE
            XYZBEF( K, NUBARY ) = REAL( S / NBSTEF )
 20      CONTINUE
C
 50   CONTINUE
C
      RETURN
      END
