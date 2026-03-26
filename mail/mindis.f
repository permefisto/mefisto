         SUBROUTINE MINDIS( NBD, XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETERMINER LE POINT M QUI MINIMISE LA SOMME DES
C -----    CARRES DES DISTANCES A NBD DROITES DE L'ESPACE
C
C          SI NBD>6 METTRE A JOUR LE PARAMETER MXD
C
C ENTREES :
C ---------
C XYZ(.,1),  XYZ(.,1+NBD)  : COORDONNEES DES SOMMETS DE LA DROITE (D1)
C XYZ(.,2),  XYZ(.,2+NBD)  : COORDONNEES DES SOMMETS DE LA DROITE (D2)
C XYZ(.,NBD),XYZ(.,NBD+NBD): COORDONNEES DES SOMMETS DE LA DROITE (DNBD)
C
C SORTIES :
C ---------
C XYZ(.,0) : COORDONNEES DU POINT CHERCHE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         JANVIER 1989
C23456---------------------------------------------------------------012
      PARAMETER  (MXD=6)
      REAL        XYZ(3,0:NBD+NBD)
      REAL        V(MXD,3),VN(MXD),A(3,3),B(3),
     %            ORI(MXD,3),AP(3,3),BP(3)
C
C     LES VECTEURS DIRECTEURS DES NBD DROITES
      DO 1 J=1,3
         DO 14 I=1,NBD
            V(I,J)=XYZ(J,I+NBD)-XYZ(J,I)
            ORI(I,J)=XYZ(J,I)
 14      CONTINUE
 1    CONTINUE
C
      DO 12 I=1,NBD
         VN(I) = 0
         DO 11 J=1,3
            VN(I)=VN(I)+V(I,J)**2
 11      CONTINUE
 12   CONTINUE
C
C     LE SYSTEME LINEAIRE A RESOUDRE
C
      DO 2 I=1,3
         B(I)=0
         DO 22 J=1,3
            A(I,J)=0
22       CONTINUE
2     CONTINUE
C
      DO 3 ND=1,NBD
         DO 4 I=1,3
            BP(I)=0
            DO 24 J=1,3
               AP(I,J)=0
24          CONTINUE
4        CONTINUE
         DO 5 I=1,3
            DO 25 J=1,3
               AP(I,I)=AP(I,I)+V(ND,J)**2
25          CONTINUE
5        CONTINUE
         DO 6 I=1,3
            DO 26 J=1,3
               AP(I,J)=AP(I,J)-V(ND,I)*V(ND,J)
26          CONTINUE
6        CONTINUE
         DO 7 I=1,3
            DO 27 J=1,3
               BP(I)=BP(I)+AP(I,J)*ORI(ND,J)/VN(ND)
27          CONTINUE
7        CONTINUE
         DO 8 I=1,3
            B(I)=B(I)+BP(I)
            DO 28 J=1,3
               A(I,J)=A(I,J)+AP(I,J)/VN(ND)
28          CONTINUE
8        CONTINUE
3     CONTINUE
C
C     RESOLUTION DU SYSTEME LINEAIRE
      DO 9 I=1,3
         DO 39 J=I+1,3
            B(J)=B(J)-A(J,I)*B(I)/A(I,I)
               DO 29 K=I+1,3
                  A(J,K)=A(J,K)-A(J,I)*A(I,K)/A(I,I)
29             CONTINUE
39          CONTINUE
9     CONTINUE
C
      DO 10 I=3,1,-1
         XYZ(I,0)=B(I)/A(I,I)
         DO 30 J=I-1,1,-1
            B(J)=B(J)-A(J,I)*XYZ(I,0)
30       CONTINUE
10    CONTINUE
C
C     LA SOLUTION
C     WRITE(IMPRIM,100) (XYZ(I,0),I=1,3)
C 100 FORMAT(' XYZ',3E15.6)
      RETURN
      END
