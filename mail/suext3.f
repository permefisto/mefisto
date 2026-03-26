         SUBROUTINE  SUEXT3(MNTSMA,MNSOFA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES NUMEROS DES SOMMETS DES NSEF
C -----    POUR UN TRIANGLE ALGEBRIQUE STRUCTURE
C ENTREES :
C ---------
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOFA : ADRESSE MCN DU TABLEAU DES SOMMETS DES FACES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS FEVRIER 1989
C23456---------------------------------------------------------------012
C
      IMPLICIT INTEGER (W)
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     NOMBRE DE SEGMENTS PAR COTE
      NBAR = MCN( MNTSMA + WBARTR )
C
      DO 1 N1=1,NBAR
         NEL  = (N1-1)**2
         NEL1 =  NEL - 1
         I1   = (N1*N1-N1)/2 + 1
         I2   = (N1*N1+N1)/2 + 1
         DO 2 N2=1,N1
            NEL1 = NEL1 + 2
            IA   = MNSOFA - 1 + (NEL1-1) * 4
            MCN(IA+1) = I1
            MCN(IA+2) = I2
            MCN(IA+3) = I2+1
            MCN(IA+4) = 0
            I1 = I1 + 1
            I2 = I2 + 1
 2       CONTINUE
         NEL2 = NEL
         I1   = (N1*N1+N1)/2 + 1
         I2   = (N1*N1-N1)/2 + 1
         DO 3 N2=1,N1-1
            NEL2 = NEL2 + 2
            IA   = MNSOFA - 1 + (NEL2-1) * 4
            MCN(IA+1) = I1+1
            MCN(IA+2) = I2+1
            MCN(IA+3) = I2
            MCN(IA+4) = 0
            I1 = I1 + 1
            I2 = I2 + 1
 3       CONTINUE
 1    CONTINUE
C
      RETURN
      END
