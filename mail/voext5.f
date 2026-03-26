         SUBROUTINE  VOEXT5(MNTSMA,MNSOFA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES NUMEROS DES SOMMETS DES NSEF
C -----    POUR UN TETRAEDRE STRUCTURE
C ENTREES :
C ---------
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C MNSOFA : ADRESSE MCN DU TABLEAU DES SOMMETS DES CUBES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS FEVRIER 1989
C23456---------------------------------------------------------------012
C
      IMPLICIT INTEGER (W)
      INTEGER NSOM(6)
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     NOMBRE DE SEGMENTS PAR COTE
      NBAR = MCN( MNTSMA + WBARTE)
      NBVO = 0
      NSPH = 0
      NSPB = 1
      IA = MNSOFA - 1
      DO 1 NB=1,NBAR

C        ** 1) LES ELEMENTS AYANT LA BASE EN BAS
         DO NB1=1,NB
            DO NB2=1,NB1
C              LE NUMERO DU VOLUME
               NBVO = NBVO + 1
C              LES SOMMETS
               MCN(IA+1) = NB1*(NB1-1)/2+NB2+NSPH
               MCN(IA+2) = NB1*(NB1-1)/2+NB2+NSPB
               MCN(IA+3) = NB1*(NB1+1)/2+NB2+NSPB
               MCN(IA+4) = NB1*(NB1+1)/2+NB2+1+NSPB
               MCN(IA+5) = 0
               MCN(IA+6) = 0
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               IA = IA + 8
            ENDDO
         ENDDO

C        ** 2) LES ELEMENTS INTERMEDIAIRES (BASE EN HAUT OU EN BAS)
         DO NB1=1,NB-1
            DO NB2=1,NB1
C              LES SOMMETS DU HAUT
               NSOM(1) = NB1*(NB1-1)/2+NB2+NSPH
               NSOM(2) = NB1*(NB1+1)/2+NB2+NSPH
               NSOM(3) = NB1*(NB1+1)/2+NB2+1+NSPH
C              LES SOMMETS DU BAS
               NSOM(4) = NB1*(NB1+1)/2+NB2+NSPB
               NSOM(5) = (NB1+1)*(NB1+2)/2+NB2+1+NSPB
               NSOM(6) = NB1*(NB1+1)/2+NB2+1+NSPB
C              LES 4 TETRAEDRES
               NBVO = NBVO + 1
               MCN(IA+1) = NSOM(1)
               MCN(IA+2) = NSOM(2)
               MCN(IA+3) = NSOM(3)
               MCN(IA+4) = NSOM(4)
               MCN(IA+5) = 0
               MCN(IA+6) = 0
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               IA = IA + 8
               NBVO = NBVO + 1
               MCN(IA+1) = NSOM(2)
               MCN(IA+2) = NSOM(3)
               MCN(IA+3) = NSOM(4)
               MCN(IA+4) = NSOM(5)
               MCN(IA+5) = 0
               MCN(IA+6) = 0
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               IA = IA + 8
               NBVO = NBVO + 1
               MCN(IA+1) = NSOM(3)
               MCN(IA+2) = NSOM(4)
               MCN(IA+3) = NSOM(5)
               MCN(IA+4) = NSOM(6)
               MCN(IA+5) = 0
               MCN(IA+6) = 0
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               IA = IA + 8
               NBVO = NBVO + 1
               MCN(IA+1) = NSOM(4)
               MCN(IA+2) = NSOM(6)
               MCN(IA+3) = NSOM(3)
               MCN(IA+4) = NSOM(1)
               MCN(IA+5) = 0
               MCN(IA+6) = 0
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               IA = IA + 8
            ENDDO
         ENDDO

C        ** 3) LES ELEMENTS AYANT LA BASE EN HAUT
         DO NB1=1,NB-2
            DO NB2=1,NB1
C              LE NUMERO DU VOLUME
               NBVO = NBVO + 1
C              LES SOMMETS
               MCN(IA+1) = NB1*(NB1+1)/2+NB2+NSPH
               MCN(IA+2) = (NB1+1)*(NB1+2)/2+NB2+1+NSPH
               MCN(IA+3) = NB1*(NB1+1)/2+NB2+1+NSPH
               MCN(IA+4) = (NB1+1)*(NB1+2)/2+NB2+1+NSPB
               MCN(IA+5) = 0
               MCN(IA+6) = 0
               MCN(IA+7) = 0
               MCN(IA+8) = 0
               IA = IA + 8
            ENDDO
         ENDDO
         NSPH = NSPB
         NSPB = NSPB + (NB+1)*(NB+2)/2

 1    ENDDO
C
      RETURN
      END
