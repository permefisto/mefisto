      SUBROUTINE BISSEC( NUFONC , P1 , F1 , P2 , F2 , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TROUVER LE POINT ANNULANT UNE FONCTION UTILISATEUR NUFONC
C -----
C ENTREES:
C --------
C NUFONC : NUMERO DE LA FONCTION UTILISATEUR
C
C ENTREES ET SORTIES :
C --------------------
C P1     : POINT A DROITE OU A GAUCHE DE LA SURFACE
C F1     : VALEUR DE LA FONCTION NUFONC AU POINT P1
C P2     : POINT A DROITE OU A GAUCHE DE LA SURFACE
C F2     : VALEUR DE LA FONCTION NUFONC AU POINT P2
C          P1 P2 SONT TELS QUE F1 * F2 < 0
C
C ATTENTION : P1 EN SORTIE CONTIENT LA PROJECTION SUR LA SURFACE
C
C IERR   : 0 SI PAS D'ERREUR ; 1 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET ANALYSE NUMERIQUE PARIS   DECEMBRE 1980
C ......................................................................
      PARAMETER       (ITEMAX=128)
      include"./incl/gsmenu.inc"
      COMMON /UNITES/  LECTEU , IMPRIM , INTERA , NUNITE(29)
      REAL             P1(3),P2(3),P(3), F
      DOUBLE PRECISION XYZ(3),DF
      INTRINSIC        REAL
C
C     LA FONCTION EST ELLE DESIGNEE ?
      IF( NUFONC .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FONCTION NON DESIGNEE'
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     INITIATIONS
C     -----------
      IERR   = 0
      ITER   = 0
C
C     LE POINT MILIEU DU SEGMENT ET LA VALEUR DE LA FONCTION EN CE POINT
C     ------------------------------------------------------------------
 5    DO 10 I=1,3
         P(I) = ( P1(I) + P2(I) ) * 0.5
 10   CONTINUE
C
C     LE TEST DE CONVERGENCE
C     ----------------------
      IF( (P(1).EQ.P1(1) .AND. P(2).EQ.P1(2) .AND. P(3).EQ.P1(3))
     %.OR.(P(1).EQ.P2(1) .AND. P(2).EQ.P2(2) .AND. P(3).EQ.P2(3))) THEN
C        CONVERGENCE ASSUREE
         RETURN
      ENDIF
C
C     CHOIX DU POINT REMPLACANT LE MILIEU SELON LE SIGNE DE F1*F
C     ----------------------------------------------------------
C     LES PARAMETRES ET LE RESULTAT SONT REELS DOUBLE PRECISION
      XYZ(1) = P(1)
      XYZ(2) = P(2)
      XYZ(3) = P(3)
      CALL FONVAL( NUFONC , 3 , XYZ , NCODEV , DF )
      IF( NCODEV .LE. 0 ) GOTO 9900
C
      F = REAL( DF )
      IF( F .EQ. 0. ) THEN
C        P EST SOLUTION
         DO 20 I=1,3
            P1(I) = P(I)
 20      CONTINUE
         RETURN
      ELSE IF( F1*F .LT. 0. ) THEN
C        P2 DEVIENT P MILIEU DE P1P2
         DO 30 I=1,3
            P2(I) = P(I)
 30      CONTINUE
         F2 = F
      ELSE
C        P1 DEVIENT P MILIEU DE P1P2
         DO 40 I=1,3
            P1(I) = P(I)
 40      CONTINUE
         F1 = F
      ENDIF
C
      ITER = ITER + 1
      IF( ITER .LE. ITEMAX ) GOTO 5
C
C     NON CONVERGENCE
C     ---------------
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') ITEMAX
      KERR(1) =' LE CALCUL DU ZERO DE LA FONCTION NECESSITE PLUS DE'//
     %          KERR(MXLGER)(1:4) //' ITERATIONS'
      CALL LEREUR
C
 9900 IERR = 1
      END
