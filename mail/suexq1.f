       SUBROUTINE SUEXQ1( NUCOTE, XYZI, XYZF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LE CONTOUR D'UN QUADRANGLE DEFINI PAR LES
C -----    LIGNES STRUCTUREES DE CHACUN DES 4 COTES
C
C ENTREES:
C --------
C XYZI : LES COORDONNEES DES POINTS INITIAUX DES COTES
C XYZF : LES COORDONNEES DES POINTS FINAUX   DES COTES
C
C SORTIES:
C --------
C NUCOTE : RECLASSEMENT DES COTES DANS L'ORDRE DIRECT
C          NUCOTE(I)=+-NO INITIAL DE LA LIGNE DE COTE I
C          ORDRE ET SENS DES 4 COTES:  S1->S2 S2->S3 S3->S4 S4->S1
C          + SI LA LIGNE A LE SENS DU COTE
C          - SI LA LIGNE A LE SENS INVERSE A CELUI DU COTE
C
C                S4----<-------S3
C                |             |
C                |             |
C               \/             /\
C                |             |
C                |             |
C                S1---->-------S2
C
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS        NOVEMBRE 1988
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS        OCTOBRE  1996
C....................................................................012
      include"./incl/gsmenu.inc"
      INTEGER   NUCOTE(4)
      REAL      XYZI(3,4), XYZF(3,4)
C
C     INITIALISATION
      IERR=0
      DO 1 N=1,4
         NUCOTE(N)=N
 1    CONTINUE
C
C     RECHERCHE DU COTE 2 S2->S3
      DO 2 N=2,4
         CALL XYZIDE(XYZF(1,1),XYZI(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NUCOTE(2)=N
            GO TO 4
         ENDIF
         CALL XYZIDE(XYZF(1,1),XYZF(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NUCOTE(2)=-N
            GO TO 4
         ENDIF
 2    CONTINUE
      GOTO 9900
C
C     RECHERCHE DU COTE 4 S4->S1
 4    DO 6 N=2,4
         CALL XYZIDE(XYZI(1,1),XYZF(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NUCOTE(4)=N
            GO TO 8
         ENDIF
         CALL XYZIDE(XYZI(1,1),XYZI(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NUCOTE(4)=-N
            GO TO 8
         ENDIF
 6    CONTINUE
      GOTO 9900
C
C     RECHERCHE DU COTE 3  S3->S4
 8    IF ( NUCOTE(2) .GE. 0 ) THEN
         DO 10 N=2,4
            CALL XYZIDE(XYZF(1,NUCOTE(2)),XYZI(1,N),IDENT)
            IF (IDENT.EQ.1) THEN
               NUCOTE(3)=N
               GO TO 12
            ENDIF
            IF ( N .EQ. NUCOTE(2) ) GO TO 10
            CALL XYZIDE(XYZF(1,NUCOTE(2)),XYZF(1,N),IDENT)
            IF (IDENT.EQ.1) THEN
               NUCOTE(3)=-N
               GO TO 12
            ENDIF
 10      CONTINUE
         GOTO 9900
      ELSE
         DO 14 N=2,4
            IF (N .EQ. -NUCOTE(2) ) GO TO 14
            CALL XYZIDE(XYZI(1,-NUCOTE(2)),XYZI(1,N),IDENT)
            IF (IDENT.EQ.1) THEN
               NUCOTE(3)=N
               GO TO 12
            ENDIF
            CALL XYZIDE(XYZI(1,-NUCOTE(2)),XYZF(1,N),IDENT)
            IF (IDENT.EQ.1) THEN
               NUCOTE(3)=-N
               GO TO 12
            ENDIF
 14      CONTINUE
         GOTO 9900
      ENDIF
C     LE CONTOUR EST COMPLET
 12   CONTINUE
      RETURN
C
C     ERREUR QUADRANGLE NON FERME
 9900 NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR: QUADRANGLE NON FERME'
      CALL LEREUR
      IERR = 6
      END
