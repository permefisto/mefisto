         SUBROUTINE VOEXH7(MNTSMA,MNSOFA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES NUMEROS DES SOMMETS DES NSEF
C -----    POUR UN HEXAEDRE ALGEBRIQUE STRUCTURE
C
C ENTREES :
C ---------
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOFA : ADRESSE MCN DU TABLEAU DES SOMMETS DES CUBES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         FEVRIER 1989
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     NOMBRE DE SOMMETS DES NSEF
      NBARXH = MCN( MNTSMA + WBARXH )
      NBARYH = MCN( MNTSMA + WBARYH )
      NBARZH = MCN( MNTSMA + WBARZH )
      NBX    = NBARXH + 1
      NBY    = NBARYH + 1
      IA     = MNSOFA - 1
      NBXNBY = NBX * NBY
      DO NZ=1,NBARZH
         DO NY=1,NBARYH
            NUMSO =  (NY-1) * NBX + (NZ-1) * NBXNBY
            DO NX=1,NBARXH
               NUMSO     = NUMSO + 1
               MCN(IA+1) = NUMSO
               MCN(IA+2) = NUMSO + 1
               MCN(IA+3) = NUMSO + NBX + 1
               MCN(IA+4) = NUMSO + NBX
               MCN(IA+5) = NUMSO + NBXNBY
               MCN(IA+6) = NUMSO + NBXNBY + 1
               MCN(IA+7) = NUMSO + NBX + NBXNBY + 1
               MCN(IA+8) = NUMSO + NBX + NBXNBY
               IA = IA + 8
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
