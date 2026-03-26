      SUBROUTINE NSEFCUSP( NANOYAU, NC3CUB, NBSTCH, NUSTDF, NBEF1C,
     %                     NSEFDF, NSEFCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 8 NUMEROS DE SOMMETS DES HEXAEDRES
C -----    FORMES PAR NANOYAU DIFFERENCES FINIES DU NOYAU ET
C          NC3CUB COUCHES HOMOTHETIQUES EN PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C NANOYAU: NOMBRE D'ARETES DANS UNE DIRECTION DES DIFFERENCES FINIES
C NC3CUB : NOMBRE DE COUCHES HOMOTHETIQUES
C NBSTCH : NOMBRE DE SOMMETS D'UNE COUCHE
C NUSTDF : NO GLOBAL DES SOMMETS FRONTALIERS POUR LES DF
C NBEF1C : NOMBRE D'EF D'UNE COUCHE
C
C SORTIES:
C --------
C NSEFDF : NO DES 8 SOMMETS DES EF DIFFERENCES FINIES
C NSEFCH : NO DES 8 SOMMETS DES EF DES NC3CUB COUCHES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L LIONS UMPC PARIS Decembre 2006
C.......................................................................
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      INTEGER  NUSTDF(*)
CCC                   0:NANOYAU,0:NANOYAU,0:NANOYAU)
      INTEGER  NSEFDF(8,*)
CCC                     1:NANOYAU,1:NANOYAU,1:NANOYAU)
      INTEGER  NSEFCH(8, 1:NBEF1C,1:NC3CUB)
C
      INTEGER  NIVO(3)
      EQUIVALENCE  (NIVO(1),IA),(NIVO(2),JA),(NIVO(3),KA)
C
      include"./incl/nusc2c3.inc"
C     LES NUMEROS DES 4 SOMMETS DES 6 FACES=2-CUBES DANS UN 3-CUBE
C     SONT DANS include"./incl/nusc2c3.inc"
C
C     NOST1 = NO DU PREMIER SOMMET DE L'EF A TRAITER
C     LES AUTRES NO DE SOMMETS SE DEDUISENT PAR TRANSLATION DE NX1**j
      NX1  = NANOYAU + 1
      NX12 = NX1  * NX1
C
      NEF = 0
      DO 230 KA=0, NANOYAU-1
         DO 220 JA=0, NANOYAU-1
            DO 210 IA=0, NANOYAU-1
C              NUMERO DE L'EF
               NEF   = NEF + 1
C              NUMERO DU PREMIER SOMMET DE L'EF
               NOST1 = 1 + IA + JA * NX1  + KA * NX12
C              NUMERO DU SOMMET ELEMENTAIRE
               NSE = 0
               DO 203 K=0,1
                  DO 202 J=0,1
                     DO 201 I=0,1
                        NSE = NSE + 1
                        NSEFDF(NSE,NEF) = NOST1 + I + J * NX1 + K * NX12
 201                 CONTINUE
 202              CONTINUE
 203           CONTINUE
 210        CONTINUE
 220     CONTINUE
 230  CONTINUE
C
C     LES EF DE LA 1-ERE COUCHE AUTOUR DES HEXAEDRES DIFFERENCES FINIES
      NEF  = 0
      NEF0 = 0
      DO 430 KA=0, NANOYAU-1
         DO 420 JA=0, NANOYAU-1
            DO 410 IA=0, NANOYAU-1
C
               NEF0 = NEF0 + 1
C              LES 3 INDICES DIFFERENCES FINIES DU 3-CUBE NEF0
C              SONT IA, JA, KA  cf EQUIVALENCE sur NIVO
C              NUMERO GLOBAL DES 8 SOMMETS DE NEF0
C              SONT DANS NSEFDF(.,NEF0)
C
C              PARCOURS DES 6 FACES DE NEF0
               DO 405 NV=1,3
C
                  IF( NIVO(NV) .EQ. 0 ) THEN
C                    CONSTRUCTION DE NEF A PARTIR DE LA FACE 1 DE NEF0
C                    QUI EST LA FACE 2 DE NEF
                     NEF = NEF + 1
                     NF  = 2 * NV - 1
                     DO 402 NSF=1,4
C                       NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                        NSE0 = NUS2C3C(NSF,NF)
C                       NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF+1 DE NEF
                        NSE1 = NUS2C3C(NSF,NF+1)
                        NSEFCH(NSE1,NEF,1)=NSEFDF(NSE0,NEF0)
C                       NUMERO GLOBAL DU SOMMET I DE LA FACE NF DE NEF
                        NSEFCH(NSE0,NEF,1)=NUSTDF(NSEFDF(NSE0,NEF0))
 402                 CONTINUE
                  ENDIF
C
                  IF( NIVO(NV) .EQ. NANOYAU-1 ) THEN
C                    CONSTRUCTION DE NEF A PARTIR DE LA FACE 2 DE NEF0
C                    QUI EST LA FACE 1 DE NEF
                     NEF = NEF + 1
                     NF  = 2 * NV
                     DO 404 NSF=1,4
C                       NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                        NSE0 = NUS2C3C(NSF,NF)
C                       NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF-1 DE NEF
                        NSE1 = NUS2C3C(NSF,NF-1)
                        NSEFCH(NSE1,NEF,1)=NSEFDF(NSE0,NEF0)
C                       NUMERO GLOBAL DU SOMMET I DE LA FACE NF DE NEF
                        NSEFCH(NSE0,NEF,1)=NUSTDF(NSEFDF(NSE0,NEF0))
 404                 CONTINUE
                  ENDIF
C
 405           CONTINUE
 410        CONTINUE
 420     CONTINUE
 430  CONTINUE
C
C     LES EF DE LA 2-EME COUCHE AUTOUR DES 3-CUBES DIFFERENCES FINIES
      NEF = 0
      DO 530 KA=0, NANOYAU-1
         DO 520 JA=0, NANOYAU-1
            DO 510 IA=0, NANOYAU-1
C
               NEF0 = NEF0 + 1
C              LES 3 INDICES DIFFERENCES FINIES DU 3-CUBE NEF0
C              SONT IA, JA, KA  cf EQUIVALENCE sur NIVO
C              NUMERO GLOBAL DES 8 SOMMETS DE NEF0 SONT DANS NSEFDF(.,NEF0)
C              PARCOURS DES 6 FACES DE NEF0
               DO 505 NV=1,3
C
                  IF( NIVO(NV) .EQ. 0 ) THEN
C                    CONSTRUCTION DE NEF A PARTIR DE LA FACE NF DE NEF0
C                    QUI EST LA FACE 2 DE NEF
                     NEF = NEF + 1
                     NF  = 2 * NV - 1
                     DO 502 NSF=1,4
C                       NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                        NSE0 = NUS2C3C(NSF,NF)
C                       NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF
                        NSG0 = NSEFCH(NSE0,NEF,1)
C                       NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF+1 DE NEF
                        NSE1 = NUS2C3C(NSF,NF+1)
                        NSEFCH(NSE1,NEF,2)=NSG0
C                       NUMERO GLOBAL DU SOMMET DE LA FACE NF DE NEF
                        NSEFCH(NSE0,NEF,2)=NSG0 + NBSTCH
 502                 CONTINUE
                  ENDIF
C
                  IF( NIVO(NV) .EQ. NANOYAU-1 ) THEN
C                    CONSTRUCTION DE NEF A PARTIR DE LA FACE NF DE NEF0
C                    QUI EST LA FACE 1 DE NEF
                     NEF = NEF + 1
                     NF  = 2 * NV
                     DO 504 NSF=1,4
C                       NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                        NSE0 = NUS2C3C(NSF,NF)
C                       NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF
                        NSG0 = NSEFCH(NSE0,NEF,1)
C                       NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF-1 DE NEF
                        NSE1 = NUS2C3C(NSF,NF-1)
                        NSEFCH(NSE1,NEF,2)=NSG0
C                       NUMERO GLOBAL DU SOMMET DE LA FACE NF DE NEF
                        NSEFCH(NSE0,NEF,2)=NSG0 + NBSTCH
 504                 CONTINUE
                  ENDIF
C
 505           CONTINUE
 510        CONTINUE
 520     CONTINUE
 530  CONTINUE
C
C     LES EF AU DELA DE LA 2-EME COUCHE AUTOUR DES 3-CUBES DIFFERENCES FINIES
C     TOUS LES NUMEROS DE SOMMET SE DEDUISENT PAR UNE TRANSLATION DE NBSTCH
      DO 690 NCH = 3, NC3CUB
         DO 680 NEF0 = 1, NBEF1C
            DO 670 NSE0 = 1, 8
               NSEFCH(NSE0,NEF0,NCH)= NSEFCH(NSE0,NEF0,NCH-1) + NBSTCH
 670        CONTINUE
 680     CONTINUE
 690  CONTINUE
      RETURN
      END
