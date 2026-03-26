      SUBROUTINE EFLEGPOS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EFFACER LA LEGENDE SUR LE FICHIER POSTSCRIPT DU TRACE
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & Saint PIERRE du PERRAY   AOUT 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
C
C     EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
        IF ( LASOPS .EQ. 1 ) THEN
          LASOPS = -11
        ELSE
          IF ( LASOPS .EQ. 2 ) THEN
            LASOPS = -12
          ELSE
            LASOPS = 0
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'eflegpos: MAUVAISE VALEUR de LASOPS'
               KERR(2) = '          ARRET du TRACE POSTSCRIPT'
            ELSE
               KERR(1) = 'eflegpos: BAD VALUE of LASOPS'
               KERR(2) = '          STOP of the POSTSCRIPT DRAWING'
            ENDIF
            CALL LEREUR
            GOTO 9000
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT(LASOPS)
        LASOPS = - LASOPS
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
 9000 RETURN
      END
