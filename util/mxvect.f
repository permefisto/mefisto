      SUBROUTINE MXVECT( NTDL, NDSM,   VECTEUR,
     %                   VMIN, NDLMIN, NCAMIN, VMAX, NDLMAX, NCAMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE MIN ET MAX DES NDSM VECTEURS DE NTDL COMPOSANTES
C -----    REELLES DOUBLE PRECISION
C
C ENTREES :
C ---------
C NDSM    : NOMBRE TOTAL DE VECTEURS
C NTDL    : NOMBRE TOTAL DE COMPOSANTES
C VECTEUR : TABLEAU (NTDL,NDSM) DES VECTEURS
C
C SORTIES :
C ---------
C VMIN   : VALEUR MINIMALE REELLE DOUBLE PRECISION DES NDSM VECTEURS
C NDLMIN : NUMERO DE LA COMPOSANTE DE VALEUR MINIMALE
C NCAMIN : NUMERO DU CAS           DE VALEUR MINIMALE
C VMAX   : VALEUR MAXIMALE REELLE DOUBLE PRECISION DES NDSM VECTEURS
C NDLMAX : NUMERO DE LA COMPOSANTE DE VALEUR MAXIMALE
C NCAMAX : NUMERO DU CAS           DE VALEUR MAXIMALE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MARS 1998
C MODIFS: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VECTEUR(NTDL,NDSM), VMIN, VMAX, V
      INTRINSIC         ABS
C
      VMIN   = VECTEUR( 1, 1 )
      VMAX   = VMIN
      NDLMIN = 1
      NCAMIN = 1
      NDLMAX = 1
      NCAMAX = 1
C
      DO 40 NCAS = 1, NDSM
         DO 20 I = 1, NTDL
            V = VECTEUR( I, NCAS )
            IF( V .LT. VMIN ) THEN
               VMIN   = V
               NDLMIN = I
               NCAMIN = NCAS
            ELSE IF( V .GT. VMAX ) THEN
               VMAX   = V
               NDLMAX = I
               NCAMAX = NCAS
            ENDIF
 20      CONTINUE
 40   CONTINUE
C
C     PROTECTION CONTRE LE CAS CONSTANT => LES DIVISIONS PAR (VMAX-VMIN)
      IF( ABS(VMAX-VMIN) .LE. 1D-7*ABS(VMAX) ) THEN
         IF( VMAX .GT. 0D0 ) THEN
            VMAX = VMAX * 1.0000001D0
         ELSE IF( VMAX .EQ. 0 ) THEN
            VMAX = 1D-7
         ELSE
            VMAX = VMAX * 0.9999999D0
         ENDIF
      ENDIF
C
      RETURN
      END
