      SUBROUTINE DSBOUL( XYZ, NB, NBBOUL, CENTRE, RAYON,
     %                   NUBOIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NUMERO DE LA PREMIERE BOULE SAUF LA BOULE NB
C -----    SOIT CONTENENT OU SUR ET 0 SI LE POINT XYZ EST EXTERNE
C
C ENTREES:
C --------
C XYZ    : 3 COORDONNEES DU POINT
C NB     : NUMERO DE LA BOULE A EVITER
C NBBOUL : NOMBRE DE BOULES
C CENTRE : 3 COORDONNEES DU CENTRE DE CHACUNES DES BOULES
C RAYON  : RAYON DE CHACUNES DES BOULES
C
C SORTIES:
C --------
C NUBOIN :  NUMERO DE LA PREMIERE BOULE CONTENANT LE POINT XYZ
C          -NUMERO DE LA PREMIERE BOULE SUR LAQUELLE EST LE POINT XYZ
C           0 SI XYZ EST EXTERNE A TOUTES LES BOULES SAUF NB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET  ANALYSE NUMERIQUE UPMC  PARIS         AOUT 1998
C2345X7--------------------------------------------------------------012
      PARAMETER  (PRECBO=1E-5)
C     PRECBO: PRECISION RELATIVE POUR DECIDER SI UN POINT EST SUR LA BOULE
C
      REAL  XYZ(3), CENTRE(3,NBBOUL), RAYON(NBBOUL)
C
      DO 10 I=1,NBBOUL
C
         IF( I .NE. NB ) THEN
C
C           LE CARRE DU RAYON
            RAY = RAYON(I) ** 2
C
C           DISTANCE DU POINT AU CENTRE
            D = ( XYZ(1) - CENTRE(1,I) ) ** 2
     %        + ( XYZ(2) - CENTRE(2,I) ) ** 2
     %        + ( XYZ(3) - CENTRE(3,I) ) ** 2
C
            IF( ABS( D-RAY ) .LE. PRECBO*RAY ) THEN
C
C              POINT SUR LA BOULE I
               NUBOIN = -I
               RETURN
C
            ELSE IF( D .LE. RAY*(1.0-PRECBO) ) THEN
C
C              POINT INTERNE A LA BOULE I
               NUBOIN = I
               RETURN
C
            ENDIF
         ENDIF
C
 10   CONTINUE
C
C     LE POINT EST EXTERNE A TOUTES LES BOULES SAUF NB
      NUBOIN = 0
      RETURN
      END
