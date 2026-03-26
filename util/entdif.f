      SUBROUTINE ENTDIF( BASE , A , B , R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  R EST LA DIFFERENCE DE 2 ENTIERS MULTI MOTS A + B TELS QUE
C -----  A>=0 B>=0 A>=B . EN FAIT LE SIGNE N'INTERVIENT PAS ICI.
C        LA SOMME DE 2 ENTIERS RELATIFS EST FAITE DANS LE SP ENTADD
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
      INTEGER BASE,RETENU,DIF
C
C     LE NOMBRE DE MOTS DE R = A + B EVENTUELLEMENT RETOUCHE ENSUITE
      NA = A(-2)
      NB = B(-2)
C
C     LA RETENUE EST NULLE AU DEBUT
      RETENU = 0
C
C         I=MIN(NA,NB)
C     R = SOMME ( A(I) - B(I) )
C         I=0
C
      DO 10 I=0,NB
C
         DIF = A(I) - B(I) - RETENU
C
C        0 - (BASE-1) - 1 =< DIF =< (BASE-1) - 0 - 0
C                - BASE   =< DIF =< (BASE-1)
C        LA RETENUE EST AU PLUS EGALE A 1
C
         IF( DIF .LT. 0 ) THEN
            RETENU = 1
            DIF    = DIF + BASE
         ELSE
            RETENU = 0
         ENDIF
C
         R(I) = DIF
 10   CONTINUE
C
C             I=NA
C     R = R + SOMME( A(I) )
C             I=NB+1
C
      DO 20 I=NB+1,NA
C
         DIF = A(I) - RETENU
C
C        0 - 1 =< DIF =< (BASE-1) - 0
C
C        LA RETENUE EST AU PLUS EGALE A 1
C
         IF( DIF .LT. 0 ) THEN
            RETENU = 1
            DIF    = DIF + BASE
         ELSE
            RETENU = 0
         ENDIF
C
         R(I) = DIF
 20   CONTINUE
C
C     LA DERNIERE RETENUE EST NULLE CAR A(NA) >0 ET A>= B
C     LE NOMBRE DE MOTS - 1 DU RESULTAT
      NR = NA
C
C     LES ZEROS EN TETE DE R SONT RESORBES
 30   IF( R(NR) .EQ. 0 ) THEN
         IF( NR .GT. 0 ) THEN
            NR = NR - 1
            GOTO 30
         ENDIF
      ENDIF
      R(-2) = NR
      IF( NR.EQ.0 .AND. R(0).EQ.0 ) THEN
         R(-1 ) = 1
      ENDIF
      END
