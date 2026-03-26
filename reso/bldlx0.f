      SUBROUTINE BLDLX0( NTDL, NDSM, NBDLFX, NODLFX, VADLFX,  XG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :   PRISE EN COMPTE DES C.L. SUR LE VECTEUR XG
C  -----   LA VALEUR DE VADLFX EST IMPOSEE AU VECTEUR XG
C
C         | XG 1 |          | XG 1     |
C         |      | DEVIENT  |          |
C         | XG 2 |          | VADLFX 2 |
C
C ENTREES:
C --------
C NTDL   : NOMBRE D'INCONNUES DU SYSTEME
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NBDLFX : NOMBRE DE CONDITIONS AUX LIMITES
C VADLFX : LES VALEURS IMPOSEES
C XG     : LE VECTEUR INITIAL
C
C SORTIE :
C --------
C XG     : EST MODIFIE EN SORTIE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS AOUT 1998
C23456---------------------------------------------------------------012
      INTEGER           NODLFX(NBDLFX)
      DOUBLE PRECISION  VADLFX(NDSM,NBDLFX), XG(NTDL,NDSM)
C
      DO 10 N=1,NBDLFX
C
C        LE NUMERO DU DL FIXE
         NODL = NODLFX(N)
C
C        BOUCLE SUR SECONDS MEMBRES
         DO 5 L=1,NDSM
            XG( NODL, L ) = VADLFX( L, N )
  5      CONTINUE
C
 10   CONTINUE
      RETURN
      END
