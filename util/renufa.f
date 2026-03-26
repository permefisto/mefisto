      SUBROUTINE RENUFA( LESENS, NOTYFA, NBNOFA, NONOFA,  NONO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NO DES POINTS ET NOEUDS D UNE FACE SELON LA
C ----- PERMUTATION DE STOCKAGE INDIQUEE PAR LESENS
C
C PARAMETRES D ENTREE :
C --------------------
C LESENS : NO DU 1-ER SOMMET DE LA FACE APRES LA PERMUTATION AMENANT
C          LE PLUS PETIT NO EN TETE
C          LE SIGNE EST POSITIF SI LE 2-EME SOMMET EST LE MINIMUM
C          DU 2-EME SOMMET ET DU DERNIER DE LA FACE
C          LE SIGNE EST NEGATIF SINON
C NOTYFA : NO DU TYPE DE CETTE FACE
C NBNOFA : NOMBRE DE NOEUDS DE LA FACE
C NONOFA : NO DES NOEUDS DE LA FACE DANS LA NUMEROTATION DES NOEUDS
C          DE L ELEMENT
C
C PARAMETRES RESULTATS :
C ---------------------
C NONO   : NO DES NOEUDS APRES PERMUTATION ET PRISE EN COMPTE DU SIGNE
C          DE LESENS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET   ANALYSE NUMERIQUE PARIS        MARS 1981
C.......................................................................
      include"./incl/gsmenu.inc"
      INTEGER  NONOFA(NBNOFA), NONO(NBNOFA)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
      GOTO ( 10 , 20 , 1 , 10 , 20 , 1 ) , NOTYFA
 1    NBLGRC(NRERR) =1
      KERR(1) = 'RENUFA:TYPE DE FACE '//KERR(MXLGER)(1:4)
     %        //' NON PROGRAMME'
      CALL LEREUR
      RETURN
C
C     ==================================================================
C     FACE TRIANGULAIRE : 0 POINT INTERNE  1 NOEUD INTERNE
C     FACE QUADRANGULAIRE 0 POINT INTERNE  1 NOEUD INTERNE
C     ==================================================================
C
   10 NONO(1) = NONOFA(1)
      RETURN
C
C     ==================================================================
C     FACE TRIANGULAIRE : 0 POINT INTERNE  3 NOEUDS INTERNES
C     FACE QUADRANGULAIRE 0 POINT INTERNE  4 NOEUDS INTERNES
C     ==================================================================
C
   20 IF( LESENS .LT. 0 ) GOTO 25
C
C     FACE DIRECTE
C     ------------
      K = LESENS
      DO 22 I=1,NBNOFA
           NONO(I) = NONOFA(K)
           K       = K + 1
           IF( K .GT. NBNOFA ) K = K - NBNOFA
   22 CONTINUE
      RETURN
C
C     FACE INDIRECTE
C     --------------
C25   K = NBNOFA + 1 + LESENS
   25 K = - LESENS - 1
      DO 27 I=1,NBNOFA
           IF( K .LE. 0 ) K = K + NBNOFA
           NONO(I) = NONOFA(K)
           K       = K - 1
   27 CONTINUE
      RETURN
      END
