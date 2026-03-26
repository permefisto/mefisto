      SUBROUTINE REELREEL2( NBV, REEL, DREEL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRANSFORMER LES NBV REELS EN REELS DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC& St Pierre du Perray SEPTEMBRE 2011
C2345X7..............................................................012
      REAL               REEL(NBV)
      DOUBLE PRECISION  DREEL(NBV)
      INTRINSIC         DBLE
C
      DO K = 1, NBV
         DREEL( K ) = DBLE( REEL(K) )
      ENDDO
C
      RETURN
      END
