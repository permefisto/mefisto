      SUBROUTINE ANRAL1( A , B , R , PHI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER L'ANGLE EN RADIANS DE -PI A PI
C ----- DU SEGMENT AB AVEC OX ET SA LONGUEUR R
C
C ENTREES:
C --------
C A,B    : LES EXTREMITES DU SEGMENT DANS R ** 2
C
C SORTIES:
C --------
C R      : LONGUEUR DU SEGMENT AB
C PHI    : ANGLE EN RADIANS DE -PI A PI DE AB AVEC L'AXE OX
c          0 SI A ET B SONT CONFONDUS AVEC IMPRESSION D'UN DIAGNOSTIC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS JANVIER 1986
C..............................................................................
      PARAMETER ( PI2 = 3.14159265358979312 * 2. )
      include"./incl/gsmenu.inc"
      COMMON    / UNITES/ LECTEU,IMPRIM,NUNITE(30)
      REAL      A(2),B(2)
C
C     LE SEGMENT AB EST-IL PARALLELE  A L'UN DES AXES ?
      IF( B(2) - A(2) .NE. 0. ) THEN
C        ORDONNEES A ET B NON CONFONDUES
         IF( B(1) - A(1) .NE. 0. ) THEN
C           ABCISSES A ET B NON CONFONDUES : AB NON PARALLELE AUX AXES
            R = ( B(1) - A(1) ) ** 2 + ( B(2) - A(2) ) ** 2
            R = SQRT( R )
         ELSE
C           ABCISSES CONFONDUES
            R = ABS( B(2) - A(2) )
         ENDIF
      ELSE
C        ORDONNEES CONFONDUES
         R = ABS( B(1) - A(1) )
      ENDIF
C
C     L'ANGLE PHI EN RADIANS
      IF( R .GT. 0 ) THEN
         PHI = ACOS(  MIN( (B(1)-A(1))/R , 1.0 )  )
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(2)(1:15),'(G15.7)') A(1)
         WRITE(KERR(2)(16:30),'(G15.7)') A(2)
         KERR(1) = 'ANRAL1:POINTS A ET B CONFONDUS X='
     %           //KERR(MXLGER)(1:15)//' Y=' //KERR(MXLGER)(16:30)
         CALL LEREUR
         PHI = 0.
      ENDIF
C
      IF( B(2) - A(2) .LT. 0. ) PHI = - PHI
      END
