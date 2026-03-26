         SUBROUTINE VOEXH6( XYZ, CUXYZ, HEXYZ, NBX, NBY, NBZ, NBXYZ, F )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :         CALCULER LES COORDONNEES DE L'IMAGE D'UN SOMMET
C -----         INTERNE DU CUBE UNITE SUR CHAQUE FACE DE L'HEXAEDRE
C ENTREES :
C ---------
C XYZ         : LES 3 COORDONNEES DU SOMMET COURANT
C CUXYZ       : LES 3 COORDONNEES DES SOMMETS DES 6 FACES DU CUBE UNITE
C HEXYZ       : LES 3 COORDONNEES DES SOMMETS DES 6 FACES DE L'HEXAEDRE
C NBX,NBY,NBZ : NOMBRE DE POINTS DANS CHAQUE DIRECTION
C NBXYZ       : NOMBRE DE POINTS DANS CHAQUE FACE (MAJORATION)
C
C SORTIES :
C ---------
C F       : COORDONNEES DES SOMMETS SUR LES FACES DE L'HEXAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         JANVIER 1989
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      REAL     XYZ(3),F(3,6),T(3,2),PXYZ(2),SXYZ(2,3)
      REAL     Q(4),X(4),Y(4),Z(4)
      REAL     HEXYZ(3,6,NBXYZ),CUXYZ(3,6,NBXYZ)
      INTEGER  NUMERO(4),SUIV(4),PREC(4),NBS(2,6),JJ(2,6)
      DATA     SUIV/2,3,4,1/, PREC/4,1,2,3/
C
C-----------------------------------------------------------------------
C     LES PARAMETRES EN FONCTION DE LA FACE
C-----------------------------------------------------------------------
      NBS(1,1) = NBX
      NBS(2,1) = NBY
      JJ(1,1)  = 1
      JJ(2,1)  = 2
C
      NBS(1,2) = NBX
      NBS(2,2) = NBZ
      JJ(1,2)  = 1
      JJ(2,2)  = 3
C
      NBS(1,3) = NBY
      NBS(2,3) = NBZ
      JJ(1,3)  = 2
      JJ(2,3)  = 3
C
      NBS(1,4) = NBX
      NBS(2,4) = NBY
      JJ(1,4)  = 1
      JJ(2,4)  = 2
C
      NBS(1,5) = NBX
      NBS(2,5) = NBZ
      JJ(1,5)  = 1
      JJ(2,5)  = 3
C
      NBS(1,6) = NBY
      NBS(2,6) = NBZ
      JJ(1,6)  = 2
      JJ(2,6)  = 3
C
C-----------------------------------------------------------------------
C     RECHERCHE DE LA FACETTE QUI CONTIENT LE POINT
C-----------------------------------------------------------------------
      DO 1 NF=1,6

         NBS1 = NBS(1,NF)
         NBS2 = NBS(2,NF)
         J1   = JJ(1,NF)
         J2   = JJ(2,NF)
         DO 2 N1=1,NBS1-1
            DO 3 N2=1,NBS2-1
C              LA FACETTE EST DIVISEE EN QUATRE TRIANGLES
               NUMERO(1)=(N2-1)*NBS1+N1
               NUMERO(2)=NUMERO(1)+1
               NUMERO(3)=NUMERO(2)+NBS1
               NUMERO(4)=NUMERO(1)+NBS1
               DO 40 K=1,4
                  DO 50 J=1,3
                     T(J,1)=CUXYZ(J,NF,NUMERO(K))-XYZ(J)
                     T(J,2)=CUXYZ(J,NF,NUMERO(SUIV(K)))-XYZ(J)
 50               ENDDO
                  SURF = T(J1,1)*T(J2,2)-T(J2,1)*T(J1,2)
                  IF ( SURF .LT. - EPSXYZ ) GO TO 3
 40            ENDDO
               GO TO 4
 3          ENDDO
 2       ENDDO
         NBLGRC(NRERR) = 1
         KERR(1) ='UN POINT NE PEUT ETRE PROJETE SUR UNE FACE'
         CALL LEREUR
         WRITE(IMPRIM,2000) (XYZ(J),J=1,3),NF
 2000    FORMAT(' ERREUR VOEXH6 : LE POINT ',3E15.6,
     %   ' NE PEUT ETRE PROJETE DANS LA FACE',I4)
         RETURN
C
C-----------------------------------------------------------------------
C     CALCUL DANS L'HEXAEDRE PAR INTERPOLATION
C-----------------------------------------------------------------------
C        COORDONNEES BARYCENTRIQUES DANS CHACUN DES QUATRE TRIANGLES
 4       PXYZ(1) = XYZ(J1)
         PXYZ(2) = XYZ(J2)
         DO K=1,4
            SXYZ(1,1) = CUXYZ(J1,NF,NUMERO(K))
            SXYZ(2,1) = CUXYZ(J2,NF,NUMERO(K))
            SXYZ(1,2) = CUXYZ(J1,NF,NUMERO(SUIV(K)))
            SXYZ(2,2) = CUXYZ(J2,NF,NUMERO(SUIV(K)))
            SXYZ(1,3) = CUXYZ(J1,NF,NUMERO(PREC(K)))
            SXYZ(2,3) = CUXYZ(J2,NF,NUMERO(PREC(K)))
            CALL  COBARY(PXYZ,SXYZ,X(K),Y(K),Z(K))
         ENDDO
C        LA VALEUR MOYENNE
         XYZ1 = (Y(1)+X(2)+Y(2)+X(3)+Z(3)+Z(4))*.25
         XYZ2 = (Z(1)+Y(2)+X(3)+Y(3)+X(4)+Z(4))*.25
C        INTERPOLATION DE TYPE Q1
         Q(1) = (1.-XYZ1) * (1.-XYZ2)
         Q(2) = XYZ1 * (1.-XYZ2)
         Q(3) = XYZ1 * XYZ2
         Q(4) = (1.-XYZ1) * XYZ2
         DO J=1,3
            F(J,NF) = 0.
            DO K=1,4
               F(J,NF) = F(J,NF) + Q(K) * HEXYZ(J,NF,NUMERO(K))
            ENDDO
         ENDDO

 1    ENDDO
C
      RETURN
      END
