      SUBROUTINE ZEROS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  MODIFIER LA PRECISION AUTOUR DE POINTS OU DE L'ORIGINE
C -----  OU REINITIALISER LEURS VALEURS INITIALES
C        EPSXYZ EPZERO  PRECISION POUR L'IDENTIFICATION DE 2 POINTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1910
C2345X7..............................................................012
      COMMON / EPSSSS /  EPZERO, EPSXYZ

C     LECTURE DU MOT CLE A TRAITER
C     ----------------------------
 10   CALL LIMTCL( 'zeros', NMTCL )
      IF( NMTCL .LE. 0 ) GOTO 9000
      GOTO( 100, 200, 300 ), NMTCL

 100  CALL INVITE( 126 )
      NCVALS = 5
      R      = EPSXYZ
      CALL LIRRSP( NCVALS, R )
      IF( R .GE. 0 ) THEN
          EPSXYZ = R
      ENDIF
      GOTO 10

 200  CALL INVITE( 125 )
      NCVALS = 5
      R      = EPZERO
      CALL LIRRSP( NCVALS, R )
      IF( R .GE. 0 ) THEN
          EPZERO = R
      ENDIF
      GOTO 10

C     ICI CE DOIT ETRE LES MEMES VALEURS QUE DANS util/initia.f
 300  EPZERO = 1E-6
      EPSXYZ = 1E-4
      GOTO 10

 9000 RETURN
      END
