       SUBROUTINE IDENT2SM( NUMST1, NUMST2, MNXYZ1, MNNSEF1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER 2 SOMMETS EN UN SEUL POUR FAIRE UNE COUTURE DU
C -----    MAILLAGE, SA POSITION DEVIENT LE MILIEU DE CES 2 SOMMETS

C ENTREES:
C --------
C NUMST1 : NUMERO DU SOMMET A IDENTIFIER A NUMST2 # NUMST1
C NUMST2 : NUMERO DU SOMMET A IDENTIFIER A NUMST1
C MNXYZ1 : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF1: ADRESSE DU TABLEAU NSEF      DE LA SURFACE

C SORTIES:
C --------
C NUMST1 : NUMERO DU SOMMET FINAL (PLUS PETIT DES 2 NUMEROS DE SOMMETS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & St PIERRE du PERRAY NOVEMBRE 2015
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NUMST1, NUMST2

C     NBTQ NOMBRE DES EF
      NBTQ   = MCN(MNNSEF1+WBEFOB)
      MNSOEL = MNNSEF1 + WUSOEF

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZ1+WNBSOM)

C     LE SOMMET NUMST1 DEVIENT LE MILIEU DE SOMMET1-SOMMET2 (NUMST1-NUMST2)
      IF( NUMST1 .GT. NUMST2 ) THEN
         K      = NUMST1
         NUMST1 = NUMST2
         NUMST2 = K
      ENDIF

      NTEMP1 = MNXYZ1+WYZSOM+3*(NUMST1-1)
      NTEMP2 = MNXYZ1+WYZSOM+3*(NUMST2-1)
      DO K = 0,2
         RMCN(NTEMP1+K) = ( RMCN(NTEMP1+K) + RMCN(NTEMP2+K) ) / 2
      ENDDO

C     DANS CHAQUE QUADRANGLE-TRIANGLE LE NO NUMST2 EST REMPLACE PAR NUMST1
      DO J = 1,NBTQ
         NTEMP = MNSOEL+4*(J-1)-1
         DO K = 1,4
            IF( MCN(NTEMP+K) .EQ. NUMST2 ) THEN
                MCN(NTEMP+K) = NUMST1
            ENDIF
         ENDDO

C        SI UN QUADRANGLE A 2 SOMMETS NUMST1, IL DEVIENT UN TRIANGLE
         IF( MCN(NTEMP+4) .GT. 0 ) THEN
            K1 = 0
            K2 = 0
            NB = 0
            DO K = 1,4
               IF( MCN(NTEMP+K) .EQ. NUMST1 ) THEN
                  NB = NB + 1
                  IF( NB .EQ. 1 ) THEN
                     K1 = K
                  ELSE
                     K2 = K
                  ENDIF
               ENDIF
            ENDDO
            IF( NB .GE. 2 ) THEN
C              LE SECOND NUMERO NUMST1 EST MIS EN POSITION 4
               DO K = K2+1,4
                  MCN(NTEMP+K-1) = MCN(NTEMP+K)
               ENDDO
               MCN(NTEMP+4) = 0
            ENDIF
         ENDIF
      ENDDO

C     SUPPRESSION DES TRIANGLES AYANT 2 SOMMETS EGAUX A NUMST1
      NBTQAS = 0
      NTEMP2 = MNSOEL-1
      NTEMP  = MNSOEL-1
      DO J = 1,NBTQ

         NB = 0
         DO K = 1,4
            IF( MCN(NTEMP+K) .EQ. NUMST1 ) THEN
               NB = NB + 1
            ENDIF
         ENDDO

         IF( NB .GE. 2 ) THEN
C           LES 4 SOMMETS DU QUADRANGLE J SONT SUPPRIMES DU TABLEAU NSEF
            WRITE(IMPRIM,*) 'ident2sm: SUPPRESSION DU TRIANGLE',
     %                      (MCN(NTEMP+K),K=1,4)
            DO K = 1,4
               MCN(NTEMP+K) = 0
            ENDDO
            NBTQAS = NBTQAS + 1
         ELSE
C           LE NO DES 4 SOMMETS SANS DOUBLE NUMST1
            DO K = 1,4
               MCN(NTEMP2+K) = MCN(NTEMP+K)
            ENDDO
            NTEMP2 = NTEMP2 + 4
         ENDIF

         NTEMP = NTEMP + 4

      ENDDO

      IF( NBTQAS .LE. 0 ) GOTO 9900

C     LE TABLEAU NSEF EST MIS A JOUR
      MCN(MNNSEF1+WBEFOB) = NBTQ-NBTQAS
      MCN(MNNSEF1+WUTFMA) = 0
      MCN(MNNSEF1+WBEFTG) = 0
      MCN(MNNSEF1+WBEFAP) = 0
      CALL ECDATE( MCN(MNNSEF1) )
      MCN(MNNSEF1+MOTVAR(6)) = NONMTD('~>>>NSEF')
C
C     LE TABLEAU XYZSOMMET EST MIS A JOUR
C     ON A SUPPRIME LE SOMMET NUMST2
C     ON VA LE PERMUTER AVEC LE DERNIER SOMMET PUIS LE SUPPRIMER
C     EST-CE LE DERNIER SOMMET???
      IF( NBSOM .NE. NUMST2 ) THEN

C        PAS LE DERNIER
         NTEMP1 = MNXYZ1+WYZSOM+3*(NUMST2-1)
         NTEMP2 = MNXYZ1+WYZSOM+3*(NBSOM-1)
         RMCN(NTEMP1  ) = RMCN(NTEMP2)
         RMCN(NTEMP1+1) = RMCN(NTEMP2+1)
         RMCN(NTEMP1+2) = RMCN(NTEMP2+2)

C        ON BOUCLE DANS NSEF POUR REMPLACER NBSOM PAR NUMST2
         DO J = 1,NBTQ
            NTEMP = MNNSEF1+WUSOEF+4*J-5
            DO L = 1,4
               IF (MCN(NTEMP+L) .EQ. NBSOM ) MCN(NTEMP+L) = NUMST2
            ENDDO
         ENDDO

      ENDIF
C
C     MAINTENANT ON SUPPRIME LE DERNIER SOMMET
      NTEMP = MNXYZ1+WYZSOM+3*(NBSOM-1)
      RMCN(NTEMP  ) = 0
      RMCN(NTEMP+1) = 0
      RMCN(NTEMP+2) = 0
      MCN(MNXYZ1+WNBSOM) = NBSOM-1
      CALL ECDATE( MCN(MNXYZ1) )
      MCN(MNXYZ1+MOTVAR(6)) = NONMTD('~>>>XYZSOMMET')

 9900 RETURN
      END
