      SUBROUTINE COEXOB
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SAISIR LE CADRE EXTREME DES COORDONNEES DES OBJETS
C -----
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1992
C ......................................................................
      include"./incl/gsmenu.inc"
      include"./incl/xyzext.inc"

C     CADRE ENGLOBANT DE TOUS LES POINTS DE TOUS LES OBJETS
C     =====================================================
      PRINT*
      KERR(1) = 'REMPLACER LES XYZ min MAX du CADRE ACTUEL:'
11200 FORMAT(A1,' MIN:',G13.6,T22,A1,' MAX:',G13.6)
      WRITE(KERR(2),11200) 'X',COOEXT(1,1),'X',COOEXT(1,2)
      WRITE(KERR(3),11200) 'Y',COOEXT(2,1),'Y',COOEXT(2,2)
      WRITE(KERR(4),11200) 'Z',COOEXT(3,1),'Z',COOEXT(3,2)
      NBLGRC(NRERR) = 4
      CALL LERESU

      CALL INVITE( 107 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(1,1) )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL INVITE( 105 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(1,2) )
      IF( NCVALS .LE. 0 ) GOTO 10

      CALL INVITE( 116 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(2,1) )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL INVITE( 114 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(2,2) )
      IF( NCVALS .LE. 0 ) GOTO 10

      CALL INVITE( 124 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(3,1) )
      IF( NCVALS .LE. 0 ) GOTO 10
      CALL INVITE( 123 )
      NCVALS = 5
      CALL LIRRSP( NCVALS , COOEXT(3,2) )

 10   KERR(1) = 'PAR les XYZ min MAX du NOUVEAU CADRE:'
      WRITE(KERR(2),11200) 'X',COOEXT(1,1),'X',COOEXT(1,2)
      WRITE(KERR(3),11200) 'Y',COOEXT(2,1),'Y',COOEXT(2,2)
      WRITE(KERR(4),11200) 'Z',COOEXT(3,1),'Z',COOEXT(3,2)
      NBLGRC(NRERR) = 4
      CALL LERESU
      PRINT*

      RETURN
      END
