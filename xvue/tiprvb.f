      SUBROUTINE TIPRVB( TEINTE, INTENS, PURETE,
     &                   ROUGE , VERT  , BLEU  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : A PARTIR DE LA COULEUR EXPRIMEE EN : TEINTE INTENSITE PURETE ,
C ----- DETERMINER LES VALEUR CORRESPONDANTES EN : ROUGE VERT BLEU
C         ( TEINTE , INTENSITE , PURETE  -->  ROUGE , VERT , BLEU )
C
C ENTREES :
C ---------
C TEINTE  : LA TEINTE DE LA COULEUR EN DEGRES ( DE 0 A 360 )
C                                             ( 0 DEGRES --> ROUGE )
C N.B : TEINTE = -1. SI LA TEINTE N'EST PAS DEFINIE ( CAS DES GRIS )
C -----
C INTENS  : L'INTENSITE DE LA COULEUR ENTRE 0. ET 1. ,
C           = 0. COULEUR FONCEE.  = 1. COULEUR CLAIRE.
C PURETE  : LA PURETE DE LA COULEUR ENTRE 0. ET 1. ,
C           = 0. COULEUR GRISE (LAVEE DE BLANC)
C           = 1. COULEUR PURE  (SATUREE)
C
C SORTIES :
C ---------
C ROUGE  : LA VALEUR DE ROUGE ENTRE 0. ET 1.
C VERT   : LA VALEUR DE VERT  ENTRE 0. ET 1.
C BLEU   : LA VALEUR DE BLEU  ENTRE 0. ET 1.
C
C N.B: LA PRECISION EST DE 1 DEGRES POUR LA TEINTE,
C ----                  DE 0.01 POUR L'INTENSITE ET LA PURETE.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : SVIGA J.M.  ANALYSE NUMERIQUE UPMC  PARIS        OCTOBRE 1985
C2345X7..............................................................012
      REAL   TEINTE, INTENS, PURETE, ROUGE, VERT, BLEU
C
C     LA COULEUR EST-ELLE ACHROMATIQUE ( C.A.D. EST-CE UN GRIS ) ?
C     ============================================================
      IF( TEINTE .LT. 0. ) THEN
C
C         C'EST UN GRIS.
C         --------------
          ROUGE = INTENS
          VERT  = INTENS
          BLEU  = INTENS
C
      ELSE
C
C        C'EST UNE COULEUR CHROMATIQUE ( DIFFERENTE D'UN GRIS )
C        ------------------------------------------------------
         T = TEINTE
         IF ( T .GE. 359. )  T = 0.
C
C        LA TEINTE DANS (0. , 6.(        6. EXCLU
C        ----------------------------------------
         T = T / 60.
C
C        DETERMINATION DES POUCENTAGES.
C        ------------------------------
         I  = INT(T)
         F  = T - FLOAT(I)
         P1 = INTENS * ( 1. - PURETE )
         P2 = INTENS * ( 1. - ( PURETE * F ) )
         P3 = INTENS * ( 1. - ( PURETE * ( 1. - F ) ) )
C
C        AFFECTATION .
C        -------------
         I = I + 1
         GOTO ( 11, 12, 13, 14, 15, 16 ), I
C
C                  0 <= TEINTE <= 60     ( DE ROUGE A JAUNE )
C                  ------------------------------------------
 11                ROUGE = INTENS
                   VERT  = P3
                   BLEU  = P1
         GOTO 10
C                  60 <= TEINTE <= 120   ( DE JAUNE A VERT )
C                  -----------------------------------------
 12                ROUGE = P2
                   VERT  = INTENS
                   BLEU  = P1
         GOTO 10
C                  120 <= TEINTE <= 180  ( DE VERT  A CYAN )
C                  -----------------------------------------
 13                ROUGE = P1
                   VERT  = INTENS
                   BLEU  = P3
         GOTO 10
C                  180 <= TEINTE <= 240  ( DE CYAN  A MAGENTA )
C                  --------------------------------------------
 14                ROUGE = P1
                   VERT  = P2
                   BLEU  = INTENS
         GOTO 10
C                  240 <= TEINTE <= 300  ( DE MAGENTA A BLEU )
C                  -------------------------------------------
 15                ROUGE = P3
                   VERT  = P1
                   BLEU  = INTENS
         GOTO 10
C                  300 <= TEINTE < 360   ( DE BLEU  A ROUGE )
C                  ------------------------------------------
 16                ROUGE = INTENS
                   VERT  = P1
                   BLEU  = P2
10      CONTINUE
      END IF
      END
