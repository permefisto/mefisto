      SUBROUTINE ENTSOM( BASE , A , B , R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  R EST LA SOMME DE 2 ENTIERS MULTI MOTS A ET B TELS QUE
C -----  A>=0 B>=0  . EN FAIT LE SIGNE N'INTERVIENT PAS ICI.
C        POUR DES ENTIERS RELATIFS VOIR LE SOUS-PROGRAMME ENTADD
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
C R       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NR
C           R = SOMME ( R(I) * BASE ** I )   ET  R(-1) = SIGNE( R )
C                I=0                             R(-2) = NR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER A(-2:*),B(-2:*),R(-2:*)
      INTEGER BASE,SOM,RETENU
C
C     LE NOMBRE DE MOTS DE R = A + B EVENTUELLEMENT RETOUCHE ENSUITE
      NA = A(-2)
      NB = B(-2)
      NR = NA
C
C     LA RETENUE EST NULLE AU DEBUT
      RETENU = 0
C
C         I=MIN(NA,NB)
C     R = SOMME ( A(I) + B(I) )
C         I=0
C
      N = MIN( NA , NB )
      DO 10 I=0,N
         SOM = A(I) + B(I) + RETENU
C
C        0 + 0 + 0 =< SOM =< (BASE-1) + (BASE-1) + 1
C                            (BASE-1) + BASE
C        LA RETENUE EST AU PLUS EGALE A 1
C
         IF( SOM .GE. BASE ) THEN
            RETENU = 1
            SOM    = SOM - BASE
         ELSE
            RETENU = 0
         ENDIF
C
         R(I) = SOM
 10   CONTINUE
C
      IF( NA .GT. NB ) THEN
C                 I=NA
C        R = R + SOMME ( A(I) )
C                 I=N+1
C
         DO 20 I=N+1,NA
            SOM = A(I) + RETENU
C
C           0 + 0 =< SOM =< (BASE-1) + 1 = BASE
C           LA RETENUE EST AU PLUS EGALE A 1
C
            IF( SOM .GE. BASE ) THEN
               RETENU = 1
               SOM    = SOM - BASE
            ELSE
               RETENU = 0
            ENDIF
C
            R(I) = SOM
 20      CONTINUE
C
      ELSE IF( NA .LT. NB ) THEN
C                 I=NB
C        R = R + SOMME ( B(I) )
C                 I=N+1
C
         NR = NB
         DO 30 I=N+1,NB
            SOM = B(I) + RETENU
C
C           0 + 0 =< SOM =< (BASE-1) + 1
C                            BASE
C           LA RETENUE EST AU PLUS EGALE A 1
C
            IF( SOM .GE. BASE ) THEN
               RETENU = 1
               SOM    = SOM - BASE
            ELSE
               RETENU = 0
            ENDIF
C
            R(I) = SOM
 30      CONTINUE
      ENDIF
C
C     LE TRAITEMENT DE LA RETENUE FINALE
      IF( RETENU .GT. 0 ) THEN
         NR = NR + 1
         R( NR ) = RETENU
      ENDIF
      R(-2) = NR
      IF( NR.EQ.0 .AND. R(0).EQ.0 ) THEN
         R(-1 ) = 1
      ENDIF
      END
