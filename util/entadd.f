      SUBROUTINE ENTADD( BASE , A , B , R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  R EST LA SOMME DE 2 ENTIERS RELATIFS MULTI MOTS A ET B
C -----
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C           ATTENTION ELLE DOIT ETRE CALCULEE DE TELLE SORTE QUE
C                     BASE + (BASE-1)  NE DEBORDE PAS UN ENTIER MACHINE
C
C A       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NA
C           A = SOMME ( A(I) * BASE ** I )   ET  A(-1) = SIGNE( A )
C                I=0                             A(-2) = NA
C
C B       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NB
C           B = SOMME ( B(I) * BASE ** I )   ET  B(-1) = SIGNE( B )
C                I=0                             B(-2) = NB
C
C SORTIES :
C ---------
C R       : ENTIER MULTI-MOTS  SOMME DE A ET B  CODE SOUS LA FORME
C                I=NR
C           R = SOMME ( R(I) * BASE ** I )   ET  R(-1) = SIGNE( R )
C                I=0                             R(-2) = NR
C
C ATTENTION : A B ET R PEUVENT ETRE CONFONDUS A L'APPEL
C
C             A - B SE CALCULE AVEC ENTADD EN INVERSANT LE SIGNE DE B
C                   AUPARAVANT     B(-1)=-B(-1)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER A(-2:*),B(-2:*),R(-2:*)
      INTEGER BASE
C
C     LE TRAITEMENT SELON LE SIGNE DES ENTIERS
      IF( A(-1) .EQ. B(-1) ) THEN
C
C        LES SIGNES SONT IDENTIQUES
C        ==========================
C
C        R = SIGNE(A) * ( ABS(A) + ABS(B) )
C        ----------------------------------
         R(-1) = A(-1)
         CALL ENTSOM( BASE , A , B , R )
C
      ELSE
C
C        LES SIGNES SONT OPPOSES
C        =======================
C
C        LEQUEL DE A ET B EST LE PLUS GRAND ?
         CALL ENTMAX( A , B , MAXIMA )
C
         IF( MAXIMA .GT. 0 ) THEN
C
C           ABS(A) > ABS(B) => R =  SIGNE(A) * (ABS(A) - ABS(B))
C           ----------------------------------------------------
            R(-1) = A(-1)
            CALL ENTDIF( BASE , A , B , R )
C
         ELSE IF( MAXIMA .LT. 0 ) THEN
C
C           ABS(A) < ABS(B) => R =  SIGNE(B) * (ABS(B) - ABS(A))
C           ----------------------------------------------------
            R(-1) = B(-1)
            CALL ENTDIF( BASE , B , A , R )
C
         ELSE
C
C           A = B => R = A-B = 0
C           --------------------
            R(-2) = 0
            R(-1) = 1
            R( 0) = 0
         ENDIF
      ENDIF
      END
