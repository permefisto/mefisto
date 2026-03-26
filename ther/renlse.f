      SUBROUTINE RENLSE( NYOBJT, NUOBJT, NBCOOR, XYZPI,
     %                   TEMPS,  ONDEPR, ONDEPI, MNCNLSE,
     %                   CONLSE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT : CALCULER LE COEFFICIENT N(PR(t)**2+PI(t)**2) NON LINEAIRE DE NLSE
C ----- DEVANT U=PR(t) +i PI(t) EN UN POINT XYZPI A L'INSTANT TEMPS

C exemple: N(PR(t)**2+PI(t)**2) = BETA ( PR(t)**2 + PI(t)**2 )
C                               + GAMA ( PR(t)**2 + PI(t)**2 )**2 + ...
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET 1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3
C          EN 2D FOURNIR NBCOOR=3 ET LA VALEUR XYZPI(3)=0D0
C TEMPS  : LE TEMPS POUR LE CALCUL DE N
C ONDEPR : PARTIE REELLE     DE L'ONDE
C ONDEPI : PARTIE IMAGINAIRE DE L'ONDE
C XYZPI  : NBCOOR COORDONNEES DU POINT DE CALCUL DU COEFFICIENT
C MNCNLSE: ADRESSE MCN DU TABLEAU a___coefnlse DE L'OBJET

C SORTIE :
C --------
C CONLSE : VALEUR DU COEFFICIENT N(PR(t)**2+PI(t)**2) DU TERME
C          NON LINEAIRE DE NLSE AU POINT XYZPI AU TEMPS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR  Fevrier 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2014
C....6...............................................................012
      include"./incl/a___coefnlse.inc"
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))

      REAL             TEMPS
      DOUBLE PRECISION XYZPI(NBCOOR), ONDEPR, ONDEPI, CONLSE, PARAMF(8)
C
C     LE TYPE DES DONNEES DU COEFFICIENT N DE NLSE
      LTNLSE = MCN( MNCNLSE + WTNLSE )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTNLSE .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE CONSTANT AU COURS DU TEMPS
C        =====================================================
         CONLSE = RMCN( MNCNLSE + WONLSE )
C
      ELSE IF( LTNLSE .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR(t,x,y,z,ntyobj,nuobj,OndePR,OndePI)
C        ====================
         PARAMF(1) = TEMPS
         N = 1
         DO I=1,NBCOOR
            N = N + 1
            PARAMF( N ) = XYZPI(I)
         ENDDO
         PARAMF(N+1) = NYOBJT
         PARAMF(N+2) = NUOBJT
         PARAMF(N+3) = ONDEPR
         PARAMF(N+4) = ONDEPI
C        CALCUL DE LA VALEUR DU COEFFICIENT NLSE
         CALL FONVAL( MCN(MNCNLSE+WFNLSE), N+4, PARAMF,
     %                NCODEV, CONLSE )
C
      ELSE
C
C        ERREUR DE DONNEES
         CONLSE = 0.D0
C
      ENDIF
C
      RETURN
      END
