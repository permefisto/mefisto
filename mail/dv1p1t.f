      SUBROUTINE DV1P1T( NP,     PXYD,   NT,     NOCOTE,
     %                   N1TRVI, NOTRIA, CETRIA, NOTRSO,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE TRIANGLE DE SOMMET NP ET D'ARETE OPPOSEE NOCOTE
C -----    DU TRIANGLE NT
C
C ENTREES:
C --------
C NP     : NUMERO PXYD DU POINT
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C NT     : NUMERO NOTRIA DU TRIANGLE
C NOCOTE : NUMERO DU COTE DU TRIANGLE NT
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : NUMERO DU 1 PREMIER TRIANGLE VIDE DANS LE TABLEAU NOTRIA
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOTRIA(4,.)
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE i
C
C CETRIA : COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT ET
C          CARRE DU RAYON
C                         ------- ------- --------
C          PAR TRIANGLE : XCENTRE YCENTRE RAYON**2
C                         ------- ------- --------
C          TABLEAU REEL(3,MXTRIA)
C
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1>0 SINON
C          10  SI TRIANGLE SOLUTION ENCORE DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    OCTOBRE 1994
C....................................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,*),
     %                  NOTRSO(*)
      DOUBLE PRECISION  PXYD(3,*),
     %                  CETRIA(3,*)
C
C     VERIFICATION ARETE FRONTALIERE
      IF( NOTRIA(3+NOCOTE,NT) .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'DV1P1T: ARETE NON FRONTALIERE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     MISE A JOUR DES TABLEAUX
      IF( N1TRVI .LE. 0 ) THEN
C        SATURATION DES TRIANGLES
         NBLGRC(NRERR) = 1
         KERR(1) = 'SATURATION DES TRIANGLES'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     MISE A JOUR DU 1-ER TRIANGLE VIDE
      NT1    = N1TRVI
      N1TRVI = NOTRIA(4,N1TRVI)
C
      IF( NOCOTE .LT. 3 ) THEN
         N = NOCOTE + 1
      ELSE
         N = 1
      ENDIF
C
C     LES SOMMETS
      NS1 = NOTRIA( N, NT )
      NS2 = NOTRIA( NOCOTE, NT )
      NOTRIA(1,NT1) = NS1
      NOTRIA(2,NT1) = NS2
      NOTRIA(3,NT1) = NP
C
      NOTRIA(4,NT1) = NT
      NOTRIA(5,NT1) = 0
      NOTRIA(6,NT1) = 0
C
      NOTRIA(3+NOCOTE,NT) = NT1
C
      NOTRSO( NP ) = NT1
C
C     LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C     AU TRIANGLE NT ET CARRE DE SON RAYON
      N = 1
      CALL CENCED( PXYD(1,NS1), PXYD(1,NS2) ,
     %             PXYD(1,NP),  CETRIA(1,NT1), N )
      IF( N .NE. 0 ) IERR=10
      END
