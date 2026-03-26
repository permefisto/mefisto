      SUBROUTINE NSEF6C( NA6CUB, NC6CUB, NBSTCH, NUSTDF, NBEF1C,
     %                   NSEFDF, NSEFCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 64 NUMEROS DE SOMMETS DES 6-CUBES
C -----    FORMES PAR NA6CUB DIFFERENCES FINIES ET NC6CUB COUCHES
C          HOMOTHETIQUES EN PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C NA6CUB : NOMBRE D'ARETES DANS UNE DIRECTION DES DIFFERENCES FINIES
C NC6CUB : NOMBRE DE COUCHES HOMOTHETIQUES
C NBSTCH : NOMBRE DE SOMMETS D'UNE COUCHE
C NUSTDF : NO GLOBAL DES SOMMETS FRONTALIERS POUR LES DF
C NBEF1C : NOMBRE D'EF D'UNE COUCHE
C
C SORTIES:
C --------
C NSEFDF : NO DES 64 SOMMETS DES EF DIFFERENCES FINIES
C NSEFCH : NO DES 64 SOMMETS DES EF DES NC6CUB COUCHES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JLL UMPC PARIS       OCTOBRE 2005
C.......................................................................
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      INTEGER  NUSTDF(*)
CCC     %             0:NA6CUB,0:NA6CUB,0:NA6CUB,
CCC     %             0:NA6CUB,0:NA6CUB,0:NA6CUB)
      INTEGER  NSEFDF(64,*)
CCC     %                1:NA6CUB,1:NA6CUB,1:NA6CUB,
CCC     %                1:NA6CUB,1:NA6CUB,1:NA6CUB)
      INTEGER  NSEFCH(64, 1:NBEF1C,1:NC6CUB)
C
      INTEGER  NIVO(6)
      EQUIVALENCE  (NIVO(1),IA),(NIVO(2),JA),(NIVO(3),KA),
     %             (NIVO(4),LA),(NIVO(5),MA),(NIVO(6),NA)
C
      include"./incl/nusc5c6.inc"
C     LES NUMEROS DES 32 SOMMETS DES 12 FACES=5-CUBES DANS UN 6-CUBE
C     SONT DANS include"./incl/nusc5c6.inc"
C
C     NOST1 = NO DU PREMIER SOMMET DE L'EF A TRAITER
C     LES AUTRES NO DE SOMMETS SE DEDUISENT PAR TRANSLATION DE NX1**j
      NX1  = NA6CUB + 1
      NX12 = NX1  * NX1
      NX13 = NX12 * NX1
      NX14 = NX13 * NX1
      NX15 = NX14 * NX1
C
      NEF = 0
      DO 260 NA=0, NA6CUB-1
         DO 250 MA=0, NA6CUB-1
            DO 240 LA=0, NA6CUB-1
               DO 230 KA=0, NA6CUB-1
                  DO 220 JA=0, NA6CUB-1
                     DO 210 IA=0, NA6CUB-1
C                       NUMERO DE L'EF
                        NEF   = NEF + 1
C                       NUMERO DU PREMIER SOMMET DE L'EF
                        NOST1 = 1 + IA + JA * NX1  + KA * NX12
     %                            + LA * NX13 + MA * NX14 + NA * NX15
C                       NUMERO DU SOMMET ELEMENTAIRE
                        NSE = 0
                        DO 206 N=0,1
                           DO 205 M=0,1
                              DO 204 L=0,1
                                 DO 203 K=0,1
                                    DO 202 J=0,1
                                       DO 201 I=0,1
                                          NSE = NSE + 1
                                          NSEFDF(NSE,NEF) = NOST1
     %                                  + I        + J * NX1  + K * NX12
     %                                  + L * NX13 + M * NX14 + N * NX15
 201                                   CONTINUE
 202                                CONTINUE
 203                             CONTINUE
 204                          CONTINUE
 205                       CONTINUE
 206                    CONTINUE
 210                 CONTINUE
 220              CONTINUE
 230           CONTINUE
 240        CONTINUE
 250     CONTINUE
 260  CONTINUE
      IF( NC6CUB .LE. 0 ) GOTO 9999
C
C     LES EF DE LA 1-ERE COUCHE AUTOUR DES 6-CUBES DIFFERENCES FINIES
      NEF  = 0
      NEF0 = 0
      DO 460 NA=0, NA6CUB-1
         DO 450 MA=0, NA6CUB-1
            DO 440 LA=0, NA6CUB-1
               DO 430 KA=0, NA6CUB-1
                  DO 420 JA=0, NA6CUB-1
                     DO 410 IA=0, NA6CUB-1
C
C                       L'EF DE LA DERNIERE COUCHE DU NOYAU DIFFERENCES FINIES
                        NEF0 = NEF0 + 1
C
C                       LES 6 INDICES DIFFERENCES FINIES DU 6-CUBE NEF0
C                       SONT IA, JA, KA, LA, MA, NA en EQUIVALENCE avec NIVO()
C                       LE NUMERO GLOBAL DES 64 SOMMETS DE NEF0 SONT DANS NSEFDF
C
C                       PARCOURS DES 12 FACES DE NEF0
                        DO 405 NV=1,6
C
                           IF( NIVO(NV) .EQ. 0 ) THEN
C                             CONSTRUCTION DE NEF A PARTIR DE LA FACE 1 DE NEF0
C                             QUI EST LA FACE 2 DE NEF DE LA COUCHE 1
                              NEF = NEF + 1
                              NF  = 2 * NV - 1
                              DO 402 NSF=1,32
C                                NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                                 NSE0 = NUS5C6C(NSF,NF)
C                                NUMERO GLOBAL DU SOMMET NSE0 DE NEF0
                                 NSG0 = NSEFDF(NSE0,NEF0)
C                                NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF+
                                 NSE1 = NUS5C6C(NSF,NF+1)
C                                NUMERO GLOBAL DU SOMMET NSE1 DE NEF de la COUCH
                                 NSEFCH(NSE1,NEF,1) = NSG0
C                                NUMERO GLOBAL DU SOMMET NSE0 DE NEF de la COUCH
                                 NSEFCH(NSE0,NEF,1) = NUSTDF(NSG0)
 402                          CONTINUE
                           ENDIF
C
                           IF( NIVO(NV) .EQ. NA6CUB-1 ) THEN
C                             CONSTRUCTION DE NEF A PARTIR DE LA FACE 2 DE NEF0
C                             QUI EST LA FACE 1 DE NEF
                              NEF = NEF + 1
                              NF  = 2 * NV
                              DO 404 NSF=1,32
C                                NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                                 NSE0 = NUS5C6C(NSF,NF)
C                                NUMERO GLOBAL DU SOMMET NSE0 DE NEF0
                                 NSG0  = NSEFDF(NSE0,NEF0)
C                                NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF-
                                 NSE1 = NUS5C6C(NSF,NF-1)
C                                NUMERO GLOBAL DU SOMMET NSE1 DE NEF de la COUCH
                                 NSEFCH(NSE1,NEF,1) = NSG0
C                                NUMERO GLOBAL DU SOMMET NSE0 DE NEF de la COUCH
                                 NSEFCH(NSE0,NEF,1) = NUSTDF(NSG0)
 404                          CONTINUE
                           ENDIF
C
 405                    CONTINUE
 410                 CONTINUE
 420              CONTINUE
 430           CONTINUE
 440        CONTINUE
 450     CONTINUE
 460  CONTINUE
      IF( NC6CUB .EQ. 1 ) GOTO 9999
C
C     LES EF DE LA 2-EME COUCHE AUTOUR DES 6-CUBES DIFFERENCES FINIES
      NEF = 0
      DO 560 NA=0, NA6CUB-1
         DO 550 MA=0, NA6CUB-1
            DO 540 LA=0, NA6CUB-1
               DO 530 KA=0, NA6CUB-1
                  DO 520 JA=0, NA6CUB-1
                     DO 510 IA=0, NA6CUB-1
C
                        NEF0 = NEF0 + 1
C                       LES 6 INDICES DIFFERENCES FINIES DU 6-CUBE NEF0
C                       SONT IA, JA, KA, LA, MA, NA cf EQUIVALENCE sur NIVO
C                       NUMERO GLOBAL DES 64 SOMMETS DE NEF0
C                       SONT DANS NSEFDF(.,NEF0)
C                       PARCOURS DES 12 FACES DE NEF0
                        DO 505 NV=1,6
C
                           IF( NIVO(NV) .EQ. 0 ) THEN
C                             CONSTRUCTION DE NEF A PARTIR DE LA FACE NF DE NEF0
C                             QUI EST LA FACE 2 DE NEF
                              NEF = NEF + 1
                              NF  = 2 * NV - 1
                              DO 502 NSF=1,32
C                                NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                                 NSE0 = NUS5C6C(NSF,NF)
C                                NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF sur C
                                 NSG0 = NSEFCH(NSE0,NEF,1)
C                                NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF+1 DE
                                 NSE1 = NUS5C6C(NSF,NF+1)
C                                NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF sur C
                                 NSEFCH(NSE1,NEF,2) = NSG0
C                                NUMERO GLOBAL DU SOMMET DE LA FACE NF DE NEF
                                 NSEFCH(NSE0,NEF,2) = NSG0 + NBSTCH
 502                          CONTINUE
                           ENDIF
C
                           IF( NIVO(NV) .EQ. NA6CUB-1 ) THEN
C                             CONSTRUCTION DE NEF A PARTIR DE LA FACE NF DE NEF0
C                             QUI EST LA FACE 1 DE NEF
                              NEF = NEF + 1
                              NF  = 2 * NV
                              DO 504 NSF=1,32
C                                NUMERO ELEMENTAIRE DU SOMMET NSF DE LA FACE NF
                                 NSE0 = NUS5C6C(NSF,NF)
C                                NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF sur C
                                 NSG0 = NSEFCH(NSE0,NEF,1)
C                                NUMERO GLOBAL DU SOMMET NSF DE LA FACE NF-1 DE
                                 NSE1 = NUS5C6C(NSF,NF-1)
                                 NSEFCH(NSE1,NEF,2) = NSG0
C                                NUMERO GLOBAL DU SOMMET DE LA FACE NF DE NEF
                                 NSEFCH(NSE0,NEF,2) = NSG0 + NBSTCH
 504                          CONTINUE
                           ENDIF
C
 505                    CONTINUE
 510                 CONTINUE
 520              CONTINUE
 530           CONTINUE
 540        CONTINUE
 550     CONTINUE
 560  CONTINUE
C
C     LES EF AU DELA DE LA 2 EME COUCHE AUTOUR DES 6-CUBES DIFFERENCES FINIES
C     TOUS LES NUMEROS DES SOMMETS SE DEDUISENT PAR UNE TRANSLATION DE NBSTCH
      DO 690 NCH = 3, NC6CUB
         DO 680 NEF0 = 1, NBEF1C
            DO 670 NSE0 = 1, 64
               NSEFCH(NSE0,NEF0,NCH) = NSEFCH(NSE0,NEF0,NCH-1) + NBSTCH
 670        CONTINUE
 680     CONTINUE
 690  CONTINUE
C
 9999 RETURN
      END
