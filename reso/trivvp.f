      SUBROUTINE TRIVVP( NBVAPR, VALPR, NBCOVP, VECPR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CLASSEMENT DES VALEURS ET VECTEURS PROPRES DANS L ORDRE CROISSANT
C -----
C
C ENTREES:
C --------
C NBVAPR : NOMBRE DE VALEURS PROPRES
C NBCOVP : NOMBRE DE COMPOSANTES D'UN VECTEUR PROPRE
C          NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C
C MODIFIES:
C ---------
C VALPR  : NBVAPR VALEURS PROPRES
C VECPR  : NBVAPR VECTEURS PROPRES DE NBCOVP COMPOSANTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Paris  St Pierre du Perray Mai 2014
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)

      REAL              VALPR(NBVAPR), VP
      DOUBLE PRECISION  VECPR(NBCOVP,NBVAPR), V

      IF( NBVAPR .GT. 1 ) THEN

C        TRI CROISSANT DES VALEURS PROPRES
 10      IS = 0
         DO I=1,NBVAPR-1

            IF( VALPR(I+1) .LT. VALPR(I) ) THEN
               IS         = IS + 1
               VP         = VALPR(I+1)
               VALPR(I+1) = VALPR(I)
               VALPR(I)   = VP

               DO N=1,NBCOVP
                  V            = VECPR(N,I+1)
                  VECPR(N,I+1) = VECPR(N,I)
                  VECPR(N,I)   = V
               ENDDO

            ENDIF

         ENDDO
         IF( IS .GT. 0 ) GOTO 10

      ENDIF
C
C     AFFICHAGE DES VALEURS PROPRES ACTUELLES
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,10020) (K,VALPR(K),K=1,NBVAPR)
      ELSE
         WRITE (IMPRIM,20020) (K,VALPR(K),K=1,NBVAPR)
      ENDIF
10020 FORMAT('VALEURS PROPRES TRIEES:'/5(I5,':',G14.6))
20020 FORMAT('SORTED EIGENVALUES:'/5(I5,':',G14.6))

      RETURN
      END
