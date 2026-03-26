      SUBROUTINE MATAXO6
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CONSTRUCTION DE LA MATRICE DE LA PROJECTION AXONOMETRIQUE
C -----  A PARTIR DES 6 COORDONNEES DU POINT VISE PTV
C                 DES 6 COORDONNEES POSITION DE L'OEIL
C
C SORTIE :
C --------
C AXOMAT(6,6) DANS LE COMMON TRVAR4 DE ./incl/trvari.inc"
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C2345X789............................................................012
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     IDENTITE
      CALL MATAXI6
C
C     PTV EST SUPPOSE DIFFERENT DE OEIL
C     ---------------------------------
C     LA DIRECTION e3 d'AXONOMETRIE
      AXODIS = 0
      DO 10 N=1,6
         AXOMAT(N,3) = AXOEIL(N) - AXOPTV(N)
         AXODIS = AXODIS +  AXOMAT(N,3) ** 2
 10   CONTINUE
      AXODIS = SQRT( AXODIS )
      DO 20 N=1,6
         AXOMAT(N,3) = AXOMAT(N,3) / AXODIS
 20   CONTINUE
C
C     LA DIRECTION e1 d'AXONOMETRIE
      AXOMAT(1,1) =-AXOMAT(2,3)
      AXOMAT(2,1) = AXOMAT(1,3)
      R = SQRT( AXOMAT(1,1)**2 + AXOMAT(2,1)**2 )
      IF( R .GE. 1E-3 ) THEN
C
C        z N'EST PAS COLINEAIRE A PTV - OEIL
C        e2 = e3 * e1
C        -----------------------------------
         AXOMAT(1,1) = AXOMAT(1,1) / R
         AXOMAT(2,1) = AXOMAT(2,1) / R
C
         AXOMAT(1,2) = -AXOMAT(2,1) * AXOMAT(3,3)
         AXOMAT(2,2) =  AXOMAT(1,1) * AXOMAT(3,3)
         AXOMAT(3,2) =  AXOMAT(2,1) * AXOMAT(1,3)
     %                - AXOMAT(1,1) * AXOMAT(2,3)
c
         AXOMAT(4,2) = 1.0
ccc         AXOMAT(5,2) = 1.0
ccc         AXOMAT(6,2) = 1.0
C
C        ROTATION DE PI/4 de e4 e5 e6
         CALL MATAXI6
         AXOMAT(1,1) = SQRT(2.0) /2
         AXOMAT(2,1) = SQRT(2.0) /2
c
         AXOMAT(2,2) = 1.0
         AXOMAT(3,2) = 1.0
c
         AXOMAT(3,3) = 1.0
         AXOMAT(4,3) = 1.0
c
         AXOMAT(4,4) = 1.0
         AXOMAT(5,4) = 1.0
C
         AXOMAT(5,5) = 1.0
         AXOMAT(6,5) = 1.0
C
         AXOMAT(6,6) = 1.0
         AXOMAT(3,6) = 1.0
         AXOMAT(1,6) = 1.0
C
C        ORTHONORMALISATION DES 6 DIRECTIONS AXONOMETRIQUES
         DO 70 K=2,6
            DO 50 I=1,K-1
C              (ei,ek)
               S=0
               DO 30 L=1,6
                  S = S + AXOMAT(L,I) * AXOMAT(L,K)
 30            CONTINUE
C              ek = ek - (ei,ek) ei
               DO 40 L=1,6
                  AXOMAT(L,K) = AXOMAT(L,K) - S * AXOMAT(L,I)
 40            CONTINUE
 50         CONTINUE
C           NORME A 1
            R = 0
            DO 55 L=1,6
               R = R + AXOMAT(L,K)**2
 55         CONTINUE
            R = SQRT(R)
ccc            print*,'||e',k,'||=',R
            DO 60 L=1,6
               AXOMAT(L,K) = AXOMAT(L,K) / R
 60         CONTINUE
 70      CONTINUE
      ELSE IF( AXOEIL(3) .GE. AXOPTV(3) ) THEN
C
C        z EST COLINEAIRE A PTV - OEIL   OEIL AU DESSUS DE PTV
C        -----------------------------------------------------
C        e1=x  e2=y  e3=z  e4=u  e5=v  e6=w
         CALL MATAXI6
C
      ELSE
C
C        z EST COLINEAIRE A PTV - OEIL   OEIL AU DESSOUS DE PTV
C        ------------------------------------------------------
C        e1=-x  e2=y  e3=-z  e4=u  e5=v  e6=w
         CALL MATAXI6
         AXOMAT(1,1) = -1.0
         AXOMAT(3,3) = -1.0
      ENDIF
C
CCC      WRITE(IMPRIM,*) 'MATRICE AXOMAT EN SORTIE DU SP MATAXO6'
CCC      DO 10 I=1,6
CCC         WRITE(IMPRIM,10010) (AXOMAT(I,J),J=1,6)
CCC 10   CONTINUE
CCC10010 FORMAT(6E15.6)
      END
