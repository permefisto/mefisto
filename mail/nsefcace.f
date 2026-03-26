      SUBROUTINE NSEFCACE( NANOYAU, NC2QUA, NBSTCH, NUSTDF, NBEF1C,
     %                     NSEFDF,  NSEFCH )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 4 NUMEROS DE SOMMETS DES QUADRANGLES
C -----    FORMES PAR NANOYAU DIFFERENCES FINIES DU NOYAU ET
C          NC2QUA COUCHES HOMOTHETIQUES EN PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C NANOYAU: NOMBRE D'ARETES DANS UNE DIRECTION DES DIFFERENCES FINIES
C NC2QUA : NOMBRE DE COUCHES HOMOTHETIQUES DE QUADRANGLES
C NBSTCH : NOMBRE DE SOMMETS D'UNE COUCHE
C NUSTDF : NUMERO GLOBAL DES SOMMETS FRONTALIERS POUR LES DF
C NBEF1C : NOMBRE DE QUADRANGLES D'UNE COUCHE
C
C SORTIES:
C --------
C NSEFDF : NO DES 4 SOMMETS DES CARRES DIFFERENCES FINIES
C NSEFCH : NO DES 4 SOMMETS DES QUADRANGLES DES NC2QUA COUCHES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & LJUBLJANA SLOVENIE   Novembre 2011
C.......................................................................
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      INTEGER  NUSTDF(*)
CCC                   0:NANOYAU,0:NANOYAU)
      INTEGER  NSEFDF(4,*)
CCC                   1:NANOYAU,1:NANOYAU)
      INTEGER  NSEFCH(4,1:NBEF1C,1:NC2QUA)
C
      INTEGER  NIVO(2)
      EQUIVALENCE  (NIVO(1),IA),(NIVO(2),JA)
C
      INTEGER  NUS2C3C(2,4)
      DATA     NUS2C3C/ 1,3, 2,4, 1,2, 3,4 /
C     LES NUMEROS DES 2 SOMMETS DES ARETES DU CARRE
C
C     NOST1 = NO DU PREMIER SOMMET DE L'EF A TRAITER
C     LES AUTRES NO DE SOMMETS SE DEDUISENT PAR TRANSLATION DE NX1**j
      NX1  = NANOYAU + 1
      NX12 = NX1 * NX1
C
C     LES CARRES DU NOYAU
      NEF = 0
      DO JA=0, NANOYAU-1
         DO IA=0, NANOYAU-1
C           NUMERO DE L'EF
            NEF   = NEF + 1
C           NUMERO DU PREMIER SOMMET DE L'EF
            NOST1 = 1 + IA + JA * NX1
C           NUMERO DU SOMMET ELEMENTAIRE
            NSE = 0
            DO J=0,1
               DO I=0,1
                  NSE = NSE + 1
                  NSEFDF(NSE,NEF) = NOST1 + I + J * NX1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     LES EF DE LA 1-ERE COUCHE AUTOUR DES CARRES DIFFERENCES FINIES
      NEF  = 0
      NEF0 = 0
      DO JA=0, NANOYAU-1
         DO IA=0, NANOYAU-1
C
            NEF0 = NEF0 + 1
C           LES 2 INDICES DIFFERENCES FINIES DU CARRE NEF0
C           SONT IA, JA  cf EQUIVALENCE sur NIVO
C           NUMERO GLOBAL DES 4 SOMMETS DE NEF0 SONT DANS NSEFDF(.,NEF0)
C
C           PARCOURS DES 4 ARETES DE NEF0
            DO NV=1,2
C
               IF( NIVO(NV) .EQ. 0 ) THEN
C                 CONSTRUCTION DE NEF A PARTIR DE L'ARETE 1 DE NEF0
C                 QUI EST L'ARETE 2 DE NEF
                  NEF = NEF + 1
                  NA  = 2 * NV - 1
                  DO NSF=1,2
C                    NUMERO ELEMENTAIRE DU SOMMET NSF DE L'ARETE NA
                     NSE0 = NUS2C3C(NSF,NA)
C                    NUMERO GLOBAL DU SOMMET NSF DE L'ARETE NA+1 DE NEF
                     NSE1 = NUS2C3C(NSF,NA+1)
                     NSEFCH(NSE1,NEF,1)=NSEFDF(NSE0,NEF0)
C                    NUMERO GLOBAL DU SOMMET I DE L'ARETE NA DE NEF
                     NSEFCH(NSE0,NEF,1)=NUSTDF(NSEFDF(NSE0,NEF0))
                  ENDDO
               ENDIF
C
               IF( NIVO(NV) .EQ. NANOYAU-1 ) THEN
C                 CONSTRUCTION DE NEF A PARTIR DE L'ARETE 2 DE NEF0
C                 QUI EST L'ARETE 1 DE NEF
                  NEF = NEF + 1
                  NA  = 2 * NV
                  DO NSF=1,2
C                    NUMERO ELEMENTAIRE DU SOMMET NSF DE L'ARETE NA
                     NSE0 = NUS2C3C(NSF,NA)
C                    NUMERO GLOBAL DU SOMMET NSF DE L'ARETE NA-1 DE NEF
                     NSE1 = NUS2C3C(NSF,NA-1)
                     NSEFCH(NSE1,NEF,1)=NSEFDF(NSE0,NEF0)
C                    NUMERO GLOBAL DU SOMMET I DE L'ARETE NA DE NEF
                     NSEFCH(NSE0,NEF,1)=NUSTDF(NSEFDF(NSE0,NEF0))
                  ENDDO
               ENDIF
C
            ENDDO
         ENDDO
      ENDDO
C
C     LES EF DE LA 2-EME COUCHE AUTOUR DES 3-CUBES DIFFERENCES FINIES
      NEF = 0
      DO JA=0, NANOYAU-1
         DO IA=0, NANOYAU-1
C
            NEF0 = NEF0 + 1
C           LES 2 INDICES DIFFERENCES FINIES DU CARRE NEF0
C           SONT IA, JA  cf EQUIVALENCE sur NIVO
C           NUMERO GLOBAL DES 4 SOMMETS DE NEF0 SONT DANS NSEFDF(.,NEF0)
C           PARCOURS DES 4 ARETES DE NEF0
            DO NV=1,2
C
               IF( NIVO(NV) .EQ. 0 ) THEN
C                 CONSTRUCTION DE NEF A PARTIR DE L'ARETE NA DE NEF0
C                 QUI EST L'ARETE 2 DE NEF
                  NEF = NEF + 1
                  NA  = 2 * NV - 1
                  DO NSF=1,2
C                    NUMERO ELEMENTAIRE DU SOMMET NSF DE L'ARETE NA
                     NSE0 = NUS2C3C(NSF,NA)
C                    NUMERO GLOBAL DU SOMMET NSF DE L'ARETE NA
                     NSG0 = NSEFCH(NSE0,NEF,1)
C                    NUMERO GLOBAL DU SOMMET NSF DE L'ARETE NA+1 DE NEF
                     NSE1 = NUS2C3C(NSF,NA+1)
                     NSEFCH(NSE1,NEF,2)=NSG0
C                    NUMERO GLOBAL DU SOMMET DE L'ARETE NA DE NEF
                     NSEFCH(NSE0,NEF,2)=NSG0 + NBSTCH
                  ENDDO
               ENDIF
C
               IF( NIVO(NV) .EQ. NANOYAU-1 ) THEN
C                 CONSTRUCTION DE NEF A PARTIR DE L'ARETE NA DE NEF0
C                 QUI EST L'ARETE 1 DE NEF
                  NEF = NEF + 1
                  NA  = 2 * NV
                  DO NSF=1,2
C                    NUMERO ELEMENTAIRE DU SOMMET NSF DE L'ARETE NA
                     NSE0 = NUS2C3C(NSF,NA)
C                    NUMERO GLOBAL DU SOMMET NSF DE L'ARETE NA
                     NSG0 = NSEFCH(NSE0,NEF,1)
C                    NUMERO GLOBAL DU SOMMET NSF DE L'ARETE NA-1 DE NEF
                     NSE1 = NUS2C3C(NSF,NA-1)
                     NSEFCH(NSE1,NEF,2)=NSG0
C                    NUMERO GLOBAL DU SOMMET DE L'ARETE NA DE NEF
                     NSEFCH(NSE0,NEF,2)=NSG0 + NBSTCH
                  ENDDO
               ENDIF
C
            ENDDO
         ENDDO
      ENDDO
C
C     LES EF AU DELA DE LA 2-EME COUCHE AUTOUR DES CARRES DIFFERENCES FINIES
C     TOUS LES NUMEROS DE SOMMET SE DEDUISENT PAR UNE TRANSLATION DE NBSTCH
      DO NCH = 3, NC2QUA
         DO NEF0 = 1, NBEF1C
            DO NSE0 = 1, 4
               NSEFCH(NSE0,NEF0,NCH)= NSEFCH(NSE0,NEF0,NCH-1) + NBSTCH
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
