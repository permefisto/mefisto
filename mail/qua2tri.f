      SUBROUTINE QUA2TRI( NBEFSU, MCNUSO,  NBTRI, MCTRI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   DECOUPAGE DES QUADRANGLES EN 2 TRIANGLES
C -----   DANS UN NOUVEAU TABLEAU D'ADRESSE MCN MCTRI
C
C
C ENTREES:
C --------
C NBEFSU: NOMBRE DE TRIANGLES OU QUADRANGLES
C MCNUSO: ADRESSE MCN DU TABLEAU DES NO DES 4 SOMMETS DES NBEFSU T ou Q
C
C SORTIES:
C --------
C NBTRI : NOMBRE DE TRIANGLES APRES DECOUPAGE
C         =0 SI PAS DE QUADRANGLE DECOUPE
C MCTRI : ADRESSE MCN DU TABLEAU DES NO DES 4 SOMMETS DES TRIANGLES
C         =0 SI PAS DE QUADRANGLE DECOUPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET Alain LJLL UPMC & St Pierre du Perray  NOVEMBRE 2011
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      NBTRI = 0
      MCTRI = 0
C
C     EXISTE-T-IL DES QUADRANGLES?
      MN     = MCNUSO + 3
      NBQUAD = 0
      DO NEF=1,NBEFSU
         IF( MCN(MN) .GT. 0 ) THEN
            NBQUAD = NBQUAD + 1
         ENDIF
         MN = MN + 4
      ENDDO
      WRITE(IMPRIM,*) 'QUA2TRI:',NBQUAD,' QUADRANGLES PARMI',
     %                 NBEFSU,' EF'
C
      IF( NBQUAD .GT. 0 ) THEN
C
C        OUI: DECOUPAGE DES QUADRANGLES EN 2 TRIANGLES
         CALL TNMCDC( 'ENTIER', 4*(NBEFSU+NBQUAD), MCTRI )
         MN  = MCNUSO - 1
         MNT = MCTRI  - 1
         DO NEF=1,NBEFSU
C
            IF( MCN(MN+4) .GT. 0 ) THEN
C
C              LE QUADRANGLE NEF EST DECOUPE EN 2 TRIANGLES
C              LE PREMIER TRIANGLE
               NBTRI = NBTRI + 1
               DO K=1,3
                  MCN(MNT+K) = MCN(MN+K)
               ENDDO
               MCN(MNT+4) = 0
               MNT = MNT + 4
C
C              LE SECOND TRIANGLE
               NBTRI = NBTRI + 1
               MCN(MNT+1) = MCN(MN+1)
               MCN(MNT+2) = MCN(MN+3)
               MCN(MNT+3) = MCN(MN+4)
               MCN(MNT+4) = 0
               MNT = MNT + 4
C
            ELSE
C
C              COPIE SIMPLE DU TRIANGLE NEF
               NBTRI = NBTRI + 1
               DO K=1,4
                  MCN(MNT+K) = MCN(MN+K)
               ENDDO
               MNT = MNT + 4
C
            ENDIF
            MN = MN + 4
         ENDDO
C
      ENDIF
C
      RETURN
      END
