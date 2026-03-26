      SUBROUTINE VOLPLUS( MNSOCU, MNCUVO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRANSFORMER LES EF 3D DE VOLUME NEGATIF EN VOLUME POSITIF
C -----    PAR PERMUTATION DE CERTAINS SOMMETS
C
C ENTREES:
C --------
C MNSOCU : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF '~td/d/a___xyzsommet'
C MNCUVO : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ELEMENTS
C          CF '~td/d/a___nsef'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J.L. LIONS UPMC Paris   Mars 2007
C23456+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))
      REAL          V1(3),V2(3),V3(3),V4(3)
C
C     VERIFICATION DU VOLUME POSITIF DES TETRAEDRES PENTAEDRES
C                                        PYRAMIDES ET HEXAEDRES
C     ET SINON PERMUTATION DE SOMMETS POUR L'OBTENIR
      MNC = MNCUVO + WUSOEF - 1
      MNX = MNSOCU + WYZSOM - 4
C
C     NOMBRE D'ELEMENTS FINIS 3D
      NBEF = MCN(MNCUVO+WBEFOB)
      DO 400 NBT=1,NBEF
C
         IF( MCN(MNC+5) .EQ. 0 ) THEN
C           TETRAEDRE
            NS = 4
            NCOGEL = 5
         ELSE IF( MCN(MNC+6) .EQ. 0 ) THEN
C           PYRAMIDE DE BASE CARREE
            NS = 5
            NCOGEL = 9
         ELSE IF( MCN(MNC+7) .EQ. 0 ) THEN
C           PENTAEDRE
            NS = 4
            NCOGEL = 6
         ELSE
C           HEXAEDRE
            NS = 5
            NCOGEL = 7
         ENDIF
C
C        12 PRODUIT VECTORIEL 13
         DO 350 I=1,3
            V     = RMCN(MNX+3*MCN(MNC+ 1)+I)
            V1(I) = RMCN(MNX+3*MCN(MNC+ 2)+I) - V
            V2(I) = RMCN(MNX+3*MCN(MNC+ 3)+I) - V
            V4(I) = RMCN(MNX+3*MCN(MNC+NS)+I) - V
 350     CONTINUE
         CALL PROVER( V1, V2, V3 )
C
C        ( 12 PRODUIT VECTORIEL 13 ) PRODUIT SCALAIRE 14 ou 15
         V = PROSCR( V3, V4, 3 )
C
         IF( V .LT. 0 ) THEN
C           VOLUME NEGATIF
C           PERMUTATION DES SOMMETS
            IF( NCOGEL .EQ. 5 ) THEN
C              TETRAEDRE  PERMUTATION des SOMMETS 2 et 3
               NS         = MCN(MNC+2)
               MCN(MNC+2) = MCN(MNC+3)
               MCN(MNC+3) = NS
            ELSE IF( NCOGEL .EQ. 9 ) THEN
C              PYRAMIDE  PERMUTATION des SOMMETS 2 et 4
               NS         = MCN(MNC+2)
               MCN(MNC+2) = MCN(MNC+4)
               MCN(MNC+4) = NS
            ELSE IF( NCOGEL .EQ. 6 ) THEN
C              PENTAEDRE  PERMUTATION des SOMMETS 2-3 et 5-6
               NS         = MCN(MNC+2)
               MCN(MNC+2) = MCN(MNC+3)
               MCN(MNC+3) = NS
               NS         = MCN(MNC+5)
               MCN(MNC+5) = MCN(MNC+6)
               MCN(MNC+6) = NS
            ELSE
C              HEXAEDRE  PERMUTATION des SOMMETS 2-4 et 6-8
               NS         = MCN(MNC+2)
               MCN(MNC+2) = MCN(MNC+4)
               MCN(MNC+4) = NS
               NS         = MCN(MNC+6)
               MCN(MNC+6) = MCN(MNC+8)
               MCN(MNC+8) = NS
            ENDIF
         ENDIF
         MNC = MNC + 8
 400  CONTINUE
C
C     LA DATE DE MODIFICATION
      CALL ECDATE( MCN(MNSOCU) )
      CALL ECDATE( MCN(MNCUVO) )
C
      RETURN
      END
