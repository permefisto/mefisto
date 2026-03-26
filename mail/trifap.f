      SUBROUTINE TRIFAP( NBFAPE, NOFAPE, LEFACO, PTXYZD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRIER SELON LES SURFACES DECROISSANTES LES FACES PERDUES
C -----
C
C
C ENTREES:
C --------
C NBFAPE : NOMBRE DE FACES PERDUES
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C
C MODIFIE:
C --------
C NOFAPE : NUMERO DES FACES PERDUES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    OCTOBRE 1991
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      DOUBLE PRECISION  PTXYZD(1:4,1:*), SURTRD
      INTEGER           LEFACO(11,0:*),NOFAPE(1:NBFAPE)
C
      IF( NBFAPE .LE. 1 ) GOTO 9900
C
C     DECLARATION DE 2 TABLEAUX
      MNSUFP = 0
      CALL TNMCDC( 'REEL',   NBFAPE,   MNSUFP )
      MNNOAN = 0
      CALL TNMCDC( 'ENTIER', NBFAPE*2, MNNOAN )
C
C     CALCUL DE LA SURFACE DE CHAQUE FACE PERDUE
      MN = MNNOAN - 1
      DO 10 I=1,NBFAPE
         NOFAPE(I) = ABS( NOFAPE(I) )
         NF = NOFAPE( I )
C        LA SURFACE DE NF
         RMCN( MNSUFP-1+I ) = REAL( SURTRD( PTXYZD(1,LEFACO(1,NF)),
     %                                      PTXYZD(1,LEFACO(2,NF)),
     %                                      PTXYZD(1,LEFACO(3,NF)) ) )
C        IDENTITE DANS NO ANCIEN
         MCN( MN+I ) = I
 10   CONTINUE
C
C     TRI CROISSANT
      CALL TRITRP( NBFAPE, RMCN(MNSUFP),  MCN(MNNOAN) )
C
C     TRI DECROISSANT DE NOFAPE
      MN = MNNOAN + NBFAPE
      DO 20 I=NBFAPE,1,-1
C        LE NUMERO DANS NOFAPE DU PLUS GRAND ACTUEL
         NF = MCN( MNNOAN-1+I )
C        MISE DANS LA 2-EME PARTIE DE NOAN
         MCN( MN + NBFAPE - I ) = NOFAPE( NF )
 20   CONTINUE
C
C     RECOPIE DANS NOFAPE
      MN = MN - 1
      DO 30 I=1,NBFAPE
         NOFAPE(I) = MCN(MN+I)
 30   CONTINUE
C
C     DESTRUCTION DES 2 TABLEAUX
      CALL TNMCDS( 'REEL',   NBFAPE,   MNSUFP )
      CALL TNMCDS( 'ENTIER', NBFAPE*2, MNNOAN )
C
 9900 RETURN
      END
