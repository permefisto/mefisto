      SUBROUTINE BLDLFX( NTDL, NDSM, NBDLFX, NODLFX, VADLFX, STGV, XG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :   PRISE EN COMPTE DES C.L. SUR LE VECTEUR XG
C  -----   LA VALEUR DE VADLFX EST IMPOSEE AU VECTEUR XG
C
C         | XG 1 |          | XG 1            |
C         |      | DEVIENT  |                 |
C         | XG 2 |          | VADLFX * STGV 2 |     DL FIXES

C ENTREES:
C --------
C NTDL   : NOMBRE D'INCONNUES DU SYSTEME
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NBDLFX : NOMBRE DE CONDITIONS AUX LIMITES
C NODLFX : LES NUMEROS DES DEGRES DE LIBERTE FIXES
C VADLFX : LES VALEURS DES DEGRES DE LIBERTE FIXES
C STGV   : LA VALEUR TRES GRANDE (1D30) SI MATRICE PROFIL UTILISEE
C          ou 1D0 SI MATRICE MORSE CONDENSEE UTILISEE ou ...
C XG     : LE VECTEUR INITIAL

C SORTIE :
C --------
C XG     : LE VECTEUR FINAL MODIFIE EN SORTIE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS AOUT 1998
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY             Avril 2022
C23456---------------------------------------------------------------012
      INTEGER           NODLFX(NBDLFX)
      DOUBLE PRECISION  VADLFX(NDSM,NBDLFX), STGV, XG(NTDL,NDSM)

      IF( NBDLFX .GT. 0 ) THEN
         DO N=1,NBDLFX

C           LE NUMERO DU DL FIXE
            NODL = NODLFX( N )

C           BOUCLE SUR LES NDSM SECONDS MEMBRES
            DO L=1,NDSM
               XG( NODL, L ) = VADLFX( L, N ) * STGV
            ENDDO

         ENDDO
      ENDIF

      RETURN
      END
