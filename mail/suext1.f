       SUBROUTINE SUEXT1(NCOTE,NSENS,XYZI,XYZF,IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LE CONTOUR D'UN TRIANGLE DEFINI PAR LES
C -----    LIGNES STRUCTUREES DE CHACUN DES 3 COTES
C
C ENTREES:
C --------
C XYZI : LES COORDONNEES DES POINTS INITIAUX
C XYZF : LES COORDONNEES DES POINTS FINAUX
C
C SORTIES:
C --------
C NCOTE : RECLASSEMENT DES COTES DANS L'ORDRE USUEL
C NSENS : ORIENTATION DES COTES DANS LE SENS DIRECT
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS   FEVRIER 1989
C.......................................................................
      include"./incl/gsmenu.inc"
      INTEGER  NCOTE(3),NSENS(3)
      REAL     XYZI(3,3),XYZF(3,3)
C
C     INITIALISATION
C
      IERR=0
      DO 1 N=1,3
         NCOTE(N)=N
         NSENS(N)=1
 1    CONTINUE
C
C     RECHERCHE DU COTE 2
C
      DO 2 N=2,3
         CALL XYZIDE(XYZF(1,1),XYZI(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NCOTE(2)=N
            GO TO 4
         ENDIF
         CALL XYZIDE(XYZF(1,1),XYZF(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NCOTE(2)=N
            NSENS(N)=-1
            GO TO 4
         ENDIF
 2    CONTINUE
      NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR: LE CONTOUR N''EST PAS FERME'
      CALL LEREUR
      IERR = 6
      RETURN
C
C     RECHERCHE DU COTE 3
C
 4    DO 6 N=2,3
         CALL XYZIDE(XYZI(1,1),XYZF(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NCOTE(3)=N
            GO TO 8
         ENDIF
         CALL XYZIDE(XYZI(1,1),XYZI(1,N),IDENT)
         IF (IDENT.EQ.1) THEN
            NCOTE(3)=N
            NSENS(N)=-1
            GO TO 8
         ENDIF
 6    CONTINUE
      NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR: LE CONTOUR N''EST PAS FERME'
      CALL LEREUR
      IERR = 6
      RETURN
C
C     LE CONTOUR EST COMPLET
 8    CONTINUE
C
      RETURN
      END
