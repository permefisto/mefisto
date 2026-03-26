      SUBROUTINE VRSION( NMVERS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  NOM de la DATE de la VERSION ACTUELLE de MEFISTO
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint Pierre du Perray  &  Veulettes sur mer
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      CHARACTER*32   NMVERS

      IF( LANGAG .EQ. 0 ) THEN
         NMVERS = 'Version du 29 janvier 2026'
      ELSE
         NMVERS = 'January 29, 2026 version'
      ENDIF

      RETURN
      END
