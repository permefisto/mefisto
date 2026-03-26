      SUBROUTINE PROJ6C( NOPROJ, KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUT DU TYPE DE LA PROJECTION SI OBJET 6D A KNOM
C -----
C
C ENTREES:
C --------
C NOPROJ : SI OBJET EN 6D
C          TYPE DE PROJECTION 0 CI-DESSOUS FIXE LES COORDONNEES A ZERO
C          <0 PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C           0 : 'KDIRPRO' CALCULE DANS LE SP c6i123000.f cf dirpro.inc
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C
C MODIFIES:
C ---------
C KNOM : NOM A COMPLETER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L LIONS UPMC PARIS Novembre 2006
C23456---------------------------------------------------------------012
      include"./incl/dirpro.inc"
      CHARACTER*(*) KNOM
      CHARACTER*100 NOMAUX
C
      IF( NOPROJ .LE. 0 .OR. NOPROJ .GT. 4 ) RETURN
C
C     NO DANS KNOM DU DERNIER CARACTERE NON BLANC
      N = NUDCNB( KNOM )
C
      GOTO( 1, 10, 20, 30, 40 ),NOPROJ+1
C
 1    NOMAUX = KNOM(1:N) // '  Projection ' // KDIRPRO
      GOTO 100
C
 10   NOMAUX = KNOM(1:N) // '  Projection XYZ000'
      GOTO 100
C
 20   NOMAUX = KNOM(1:N) // '  Projection XY0U00'
      GOTO 100
C
 30   NOMAUX = KNOM(1:N) // '  Projection X00UV0'
      GOTO 100
C
 40   NOMAUX = KNOM(1:N) // '  Projection 000UVW'
C
 100  KNOM = NOMAUX
      RETURN
      END
