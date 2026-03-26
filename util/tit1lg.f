      SUBROUTINE TIT1LG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE DE LA 1-ERE LIGNE EN HAUT A GAUCHE DE LA FENETRE DE TRACE
C -----  AVEC 'MEFISTO', NomUtilisateur, DATE, ...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      INTEGER         LHPXCA
      PARAMETER      (LHPXCA=20)
C
      include"./incl/langue.inc"
      include"./incl/nmproj.inc"
      include"./incl/typnoobj.inc"
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"
      CHARACTER*80   KINFO, CARAC
      CHARACTER*192  LIGNE
C
C     LA LIGNE EST TRACEE EN NOIR
      CALL XVCOULEUR( NCNOIR )
C
C     CHANGEMENT DE POLICE DE CARACTERES POUR UNE DE LHPXCA PIXELS DE HAUT
ccc      print*,'tit1lg 1: cooext X=',cooext(1,1),cooext(1,2),
ccc     %       '  Y=',cooext(2,1),cooext(2,2),
ccc     %       '  Z=',cooext(3,1),cooext(3,2)
      CALL CHOIXFONTE( LHPXCA )
ccc      print*,'tit1lg 2: cooext X=',cooext(1,1),cooext(1,2),
ccc     %       '  Y=',cooext(2,1),cooext(2,2),
ccc     %       '  Z=',cooext(3,1),cooext(3,2)
C
C     Mefisto et Projet
      IF( LANGAG .EQ. 0 ) THEN
         LIGNE = 'Mefisto Projet: ' // NMPROJ
      ELSE
         LIGNE = 'Mefisto Project: ' // NMPROJ
      ENDIF
C
C     AJOUT DE LA DATE
      CARAC = KINFO( 'DATE' )
      IF( CARAC(4:4) .EQ. ' ' ) CARAC(4:4) = '0'
      IF( CARAC(7:7) .EQ. ' ' ) CARAC(7:7) = '0'
      L = NUDCNB( LIGNE )
      IF( LANGAG .EQ. 0 ) THEN
         LIGNE = LIGNE(1:L) // ' ' // CARAC(1:6) // '20' // CARAC(7:8)
      ELSE
         LIGNE = LIGNE(1:L) // ' 20' // CARAC(1:8)
      ENDIF
C
C     AJOUT DE L'HEURE
      CARAC = KINFO( 'HEURE' )
      L     = NUDCNB( LIGNE )
      LIGNE = LIGNE(1:L) // '  ' // CARAC(1:2) // 'h '
     %                           // CARAC(3:4) // 'm '
     %                           // CARAC(5:6) // 's '
C
C     AJOUT DU NOM DE L'UTILISATEUR
      CARAC = KINFO( 'UTILISATEUR' )
      LL    = NUDCNB( CARAC )
      L     = NUDCNB( LIGNE )
      IF( LANGAG .EQ. 0 ) THEN
         LIGNE = LIGNE(1:L) // '  Auteur: ' // CARAC(1:LL)
      ELSE
         LIGNE = LIGNE(1:L) // '  Author: ' // CARAC(1:LL)
      ENDIF

      IF( NUMTYPOBJ .GT. 0 ) THEN

C        AJOUT DU NOM DU TYPE DE L'OBJET et NOM DE L'OBJET
         L = NUDCNB( LIGNE )
         LIGNE = LIGNE(1:L) // '  ' // KNMTYPOBJ

         L = NUDCNB( LIGNE )
         LIGNE = LIGNE(1:L) // ': ' // KNMOBJLX

      ENDIF

C     TRACE DE LA 1-ERE LIGNE EN HAUT A GAUCHE DE LA FENETRE DE TRACE
      CALL SANSDBL( LIGNE, L )
ccc      L = NUDCNB( LIGNE )
ccc      print *, ligne(1:L)
      CALL XVTEXTE( LIGNE(1:L), L, 20, LHPXCA+3 )
c
ccc      print*,'tit1lg 3: cooext X=',cooext(1,1),cooext(1,2),
ccc     %       '  Y=',cooext(2,1),cooext(2,2),
ccc     %       '  Z=',cooext(3,1),cooext(3,2)
C
      RETURN
      END
