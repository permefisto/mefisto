         SUBROUTINE  VOEXP5(MNTSMA,MNSOFA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES NUMEROS DES SOMMETS DES NSEF
C -----    POUR UN PENTAEDRE STRUCTURE
C ENTREES :
C ---------
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOFA : ADRESSE MCN DU TABLEAU DES SOMMETS DES CUBES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS FEVRIER 1989
C23456---------------------------------------------------------------012
C
      IMPLICIT INTEGER (W)
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     NOMBRE D'ARETES SUR CHACUN DES COTES DES TRIANGLES
      NBAT = MCN( MNTSMA + WBARTP )
      NBSA = NBAT + 1
C     NOMBRE D'ARETES SUR CHACUN DES COTES DES QUADRILATERES
      NBAQ = MCN( MNTSMA + WBARZP )
      NSP  = NBSA * (NBSA+1) / 2
      NTP  = NBAT * NBAT
      NSPH = 0
      NSPB = NSPH + NSP
C
      DO 1 NBH=1,NBAQ
         DO 2 N1=1,NBAT
            NBVO = (N1-1)*(N1-1) + NTP * (NBH-1)
            NEL1 = NBVO - 1
            I1   = (N1*N1-N1)/2 + 1
            I2   = (N1*N1+N1)/2 + 1
            DO 3 N2=1,N1
               NEL1 = NEL1 + 2
               IA   = MNSOFA - 1 + (NEL1-1) * 8
               MCN(IA+1) = I1 + NSPH
               MCN(IA+2) = I2 + NSPH
               MCN(IA+3) = I2 + 1 + NSPH
               MCN(IA+4) = I1 + NSPB
               MCN(IA+5) = I2 + NSPB
               MCN(IA+6) = I2 + 1 + NSPB
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               I1 = I1 + 1
               I2 = I2 + 1
 3          CONTINUE
            NEL2 = NBVO
            I1   = (N1*N1+N1)/2 + 1
            I2   = (N1*N1-N1)/2 + 1
            DO 4 N2=1,N1-1
               NEL2 = NEL2 + 2
               IA   = MNSOFA - 1 + (NEL2-1) * 8
               MCN(IA+1) = I1 + 1 + NSPH
               MCN(IA+2) = I2 + 1 + NSPH
               MCN(IA+3) = I2 + NSPH
               MCN(IA+4) = I1 + 1 + NSPB
               MCN(IA+5) = I2 + 1 + NSPB
               MCN(IA+6) = I2 + NSPB
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               I1 = I1 + 1
               I2 = I2 + 1
 4          CONTINUE
 2       CONTINUE
         NSPH = NSPB
         NSPB = NSPB + NSP
 1    CONTINUE
C
      RETURN
      END
