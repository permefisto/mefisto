      SUBROUTINE SUEY40( NBSOEF, NBQUAD, NOSOEL, XYZINI,
     %                   NBSOMM, INDXEF, XYZFIN,
     %                   NBARX,  NBARY,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RESTRUCTURER UNE QUADRANGULATION NON-STRUCTUREE
C -----
C
C ENTREES:
C --------
C NBSOEF : NOMBRE DE SOMMETS PAR FACE (4)
C NBQUAD : NOMBRE DE QUADRANGLES DE LA SURFACE A RESTRUCTURER
C NOSOEL : LISTE DES SOMMETS DES FACES DE LA SURFACE
C XYZINI : 3 COORDONNEES DES SOMMETS   DE LA SURFACE INITIALE
C NBSOMM : NOMBRE DE SOMMETS DE LA SURFACE NON-STRUCTUREE (XYZINI)
C          ET DE L'EVENTUELLE QUADRANGULATION  STRUCTUREE (XYZFIN)
C
C TRAVAIL:
C --------
C INDXEF : TABLEAU DE TRAVAIL
C
C SORTIES:
C --------
C NBARX  : NOMBRE D'ARETES DANS LA PREMIERE DIRECTION
C NBARY  : NOMBRE D'ARETES DANS LA SECONDE  DIRECTION
C XYZFIN : 3 COORDONNEES DES SOMMETS DE LA QUADRANGULATION FINALE
C IERR   : 0 SI PAS D'ERREUR
C          1 QUADRANGULATION NON STRUCTURABLE
C          2 ERREUR TRIANGLES DANS LA QUADRANGULATION ...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : DOURSAT CHRISTOPHE ANALYSE NUMERIQUE UPMC PARIS  JANVIER 1994
C MODIFS : PERRONNET ALAIN    ANALYSE NUMERIQUE UPMC PARIS DECEMBRE 1999
C2345+7...............................................................72
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN( 1 )
      EQUIVALENCE      (RMCN( 1 ),MCN( 1 ))
C
      REAL      XYZINI(3,NBSOMM),
     %          XYZFIN(3,NBSOMM)
      INTEGER   NOSOEL(NBSOEF,NBQUAD),
     %          INDXEF(NBQUAD)
      INTEGER   NOVOIS(-1:1,4),
     %          NGS(2)
      DATA      NOVOIS/ 4,3,2, 1,4,3, 2,1,4, 3,2,1 /
C
      IERR   = 0
      NBSOM1 = 0
      NBSOM2 = 0
C
C     VERIFICATION QUE TOUS LES EF SONT DES QUADRANGLES
C     =================================================
      DO 10 I=1,NBQUAD
        IF( NOSOEL(4,I) .EQ. 0 ) THEN
C          UN TRIANGLE INTERDIT ICI
           NBLGRC(NRERR) = 1
           WRITE(KERR(5)(1:7),'(I7)') I
           KERR(1)='EF '//KERR(5)(1:7)//' EST UN TRIANGLE INTERDIT ICI'
           CALL LEREUR
           IERR = 2
           RETURN
        ENDIF
   10 CONTINUE
C
C     VERIFICATION DE LA COMPATIBILITE NOMBRE DE SOMMETS-QUADRANGLES
C     ==============================================================
C     POUR UNE QUADRANGULATION STRUCTUREE DE SOMMETS NBS1 NBAS2 ON A
C     NBQUAD = (NBS1-1) * (NBS2-1) = NBS1*NBS2 - (NBS1+NBS2) + 1
C     NBSOMM =  NBS1    *  NBS2    => NBS2 = NBSOMM / NBS1
C     NBQUAD = NBSOMM -(NBS1+NBSOMM/NBS1) + 1 soit encore
C     NBS1**2 + NBS1 * (NBQUAD-NBSOMM-1) + NBSOMM = 0
C     NBQUAD-NBSOMM-1<0 => 2 SOLUTIONS simples OU une double
C     NBS1 = ( -NBQUAD+NBSOMM+1 + SQRT( (NBQUAD-NBSOMM-1)**2 - 4 * NBSOMM ) / 2
C     NBS2 = ( -NBQUAD+NBSOMM+1 - SQRT( (NBQUAD-NBSOMM-1)**2 - 4 * NBSOMM ) / 2
      NSOMME = NBSOMM + 1 - NBQUAD
      DELTA  = NSOMME**2 - 4 * NBSOMM
      IF( DELTA .LT. 0 ) THEN
         IERR = 3
         RETURN
      ENDIF
      DELTA = SQRT( DELTA )
      NBS1  = NINT( ( NSOMME - DELTA ) / 2 )
      NBS2  = NINT( ( NSOMME + DELTA ) / 2 )
      IF( NBS1*NBS2 .NE. NBSOMM ) THEN
        IERR = 3
        RETURN
      ENDIF
      WRITE(IMPRIM,*) 'NOMBRE D''ARETES EN "X" =',NBS1-1
      WRITE(IMPRIM,*) 'NOMBRE D''ARETES EN "Y" =',NBS2-1
C
C     CONSTRUCTION PAR HACHAGE DE LA LISTE DES ARETES DE LA SURFACE
C     =============================================================
      CALL GEARSU( NBSOEF, NBQUAD, NOSOEL,
     %             2,      L1ARET, L2ARET, MNARET, IERR )
      IF( IERR .NE. 0 ) THEN
         WRITE(IMPRIM,*) 'UNE ARETE APPARTIENT A AU MOINS 3 EF'
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE D'UN COIN DU MAILLAGE
C     ===============================
      DO 40 IARET=1,L2ARET
C
         NADR = MNARET + 5*(IARET-1)-1
         IF( (MCN(NADR+4).NE.0) .AND. (MCN(NADR+5).EQ.0) ) THEN
C
C           ARETE FRONTALIERE REPEREE NSOM1-NSOM2
            NUMELT = MCN(NADR+4)
            NSOM1  = MCN(NADR+1)
            NSOM2  = MCN(NADR+2)
C
C           RECHERCHE DANS LE QUADRANGLE NUMELT DES SOMMETS
C           NBSOM1 ET NBSOM2 RESPECTIVEMENT NOVOISINS DE NSOM1 ET NSOM2
            IF( NUMELT .GT. 0 ) THEN
               DO 20 I=1,4
                  IF( NSOM1.EQ.NOSOEL(I,NUMELT)) THEN
                     NBSOM1 = NOSOEL(NOVOIS(-1,I),NUMELT)
                     NBSOM2 = NOSOEL(NOVOIS( 0,I),NUMELT)
                  ENDIF
   20          CONTINUE
            ELSE
               DO 30 I=1,4
                  IF( NSOM1.EQ.NOSOEL(I,-NUMELT)) THEN
                     NBSOM1 = NOSOEL(NOVOIS(1,I),-NUMELT)
                     NBSOM2 = NOSOEL(NOVOIS(0,I),-NUMELT)
                  ENDIF
   30          CONTINUE
            ENDIF
C
C           TEST SUR ARETE NSOM1-NBSOM1 POUR SAVOIR SI ELLE EST FRONTIERE
            NGS(1) = MIN( NSOM1, NBSOM1 )
            NGS(2) = MAX( NSOM1, NBSOM1 )
            CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                   3, L2ARET, NOAR )
            IF( NOAR .LE. 0 ) THEN
               GOTO 990
            ENDIF
            NADR = MNARET + 5*(NOAR-1)-1
            NT1 = ABS( MCN(NADR+4))
            NT2 = MCN(NADR+5)
            IF( NT2.EQ.0) THEN
               IF( NT1 .NE. ABS(NUMELT) ) THEN
                  GOTO 991
               ENDIF
C              COIN TROUVE = NSOM1
               NSOM1 = NSOM1
               NSOM2 = NSOM2
               ISENS = SIGN(1,NUMELT)
               GOTO 50
            ENDIF
C
C           TEST SUR ARETE NSOM2-NBSOM2 POUR SAVOIR SI ELLE EST FRONTIERE
            NGS(1) = MIN( NSOM2, NBSOM2 )
            NGS(2) = MAX( NSOM2, NBSOM2 )
            CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                   3, L2ARET, NOAR )
            IF( NOAR .LE. 0 ) THEN
                GOTO 990
            ENDIF
            NADR = MNARET + 5*(NOAR-1)-1
            NT1 = ABS( MCN(NADR+4))
            NT2 = MCN(NADR+5)
            IF( NT2 .EQ. 0 ) THEN
               IF( NT1 .NE. ABS(NUMELT) ) THEN
                  GOTO 991
               ENDIF
C              COIN TROUVE = NSOM2
               NSOM1 = NSOM2
               NSOM2 = NSOM1
               ISENS = -SIGN(1,NUMELT)
               GOTO 50
            ENDIF
         ENDIF
   40 CONTINUE
      GOTO 995
   50 CONTINUE
C
C     RECHERCHE DE LA PREMIERE LIGNE D ELEMENTS
C     =========================================
      NELT   = 1
      NUMELT = ABS(NUMELT)
      INDXEF(1) = NUMELT
C
C     REORIENTATION DE NUMELT
      DO 100 I=1,4
        IF( NSOM1.EQ.NOSOEL(I,NUMELT)) GOTO 110
  100 CONTINUE
  110 CONTINUE
      N1 = NOSOEL(I,NUMELT)
      N2 = NOSOEL(NOVOIS(ISENS,I),NUMELT)
      N3 = NOSOEL(NOVOIS(0,I),NUMELT)
      N4 = NOSOEL(NOVOIS(-ISENS,I),NUMELT)
      NOSOEL(1,NUMELT) = N1
      NOSOEL(2,NUMELT) = N2
      NOSOEL(3,NUMELT) = N3
      NOSOEL(4,NUMELT) = N4
      NSOM1 = N2
      NSOM2 = N3
C
C     ARETE NSOM1-NSOM2
  120 CONTINUE
      IF(  NSOM1 .LT. NSOM2 ) THEN
         NGS(1) = NSOM1
         NGS(2) = NSOM2
         ISENS  = -1
      ELSE
         NGS(1) = NSOM2
         NGS(2) = NSOM1
         ISENS  = 1
      ENDIF
      CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S             3, L2ARET, NOAR )
      IF( NOAR .LE. 0 ) THEN
          GOTO 990
      ENDIF
      NADR = MNARET + 5*(NOAR-1)-1
      IF( MCN(NADR+5) .EQ. 0 ) THEN
C
C        BOUT DE LA LIGNE D'ELEMENTS FINIS TROUVE
         IF( (NELT+1.NE.NBS1) .AND. (NELT+1.NE.NBS2) ) THEN
            GOTO 995
         ENDIF
         NBS1  = NELT+1
         NBARX = NELT
         NBS2  = NBSOMM / NBS1
         NBARY = NBS2 - 1
         GOTO 150
      ELSE
         NADR = MNARET + 5*(NOAR-1) - 1
         NT1  = ABS( MCN(NADR+4) )
         IF( NT1 .EQ. NUMELT ) THEN
            NUMELT = ABS(MCN(NADR+5))
            ISENS  = ISENS*SIGN(1,MCN(NADR+5))
         ELSE
            NUMELT = NT1
            ISENS  = ISENS*SIGN(1,MCN(NADR+4))
         ENDIF
C
C        REORIENTATION DE NUMELT
         DO 130 I=1,4
            IF( NSOM1 .EQ. NOSOEL(I,NUMELT) ) GOTO 140
  130    CONTINUE
  140    CONTINUE
         N1 = NOSOEL(I,NUMELT)
         N2 = NOSOEL( NOVOIS(ISENS,I),  NUMELT )
         N3 = NOSOEL( NOVOIS(0,I),      NUMELT )
         N4 = NOSOEL( NOVOIS(-ISENS,I), NUMELT )
         NOSOEL(1,NUMELT) = N1
         NOSOEL(2,NUMELT) = N2
         NOSOEL(3,NUMELT) = N3
         NOSOEL(4,NUMELT) = N4
C
C        VERIFICATION QUE L'ARETE N1-N2 EST FRONTIERE
         NGS(1) = MIN( N1, N2 )
         NGS(2) = MAX( N1, N2 )
         CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                3, L2ARET, NOAR )
         IF( NOAR .LE. 0 ) THEN
             GOTO 990
         ENDIF
         NADR = MNARET + 5*(NOAR-1)-1
         IF( MCN(NADR+5) .NE. 0 ) THEN
            GOTO 995
         ENDIF
C
C        AJOUT DANS LA LISTE ET ITERATION
         NELT  = NELT + 1
         INDXEF(NELT) = NUMELT
         NSOM1 = NOSOEL(2,NUMELT)
         NSOM2 = NOSOEL(3,NUMELT)
         GOTO 120
      ENDIF
  150 CONTINUE
C
C     CONSTRUCTION DE INDXEF
C     ======================
      DO 250 J=2,NBARY
C
C        RECHERCHE DU PREMIER ELEMENT DE LA LIGNE
         NELT   = (J-2)*NBARX+1
         NUMELT =INDXEF(NELT)
         NSOM1 = NOSOEL(4,NUMELT)
         NSOM2 = NOSOEL(3,NUMELT)
         IF(  NSOM1 .LT. NSOM2 ) THEN
            NGS(1) = NSOM1
            NGS(2) = NSOM2
            ISENS  = 1
         ELSE
            NGS(1) = NSOM2
            NGS(2) = NSOM1
            ISENS  = -1
         ENDIF
         CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                3, L2ARET, NOAR )
         IF( NOAR .LE. 0 ) THEN
             GOTO 990
         ENDIF
         NADR = MNARET + 5*(NOAR-1)-1
         IF( MCN(NADR+5) .EQ. 0 ) THEN
            GOTO 995
         ENDIF
         IF( ABS(MCN(NADR+4)) .EQ. NUMELT ) THEN
            NUMELT = MCN(NADR+5)
         ELSE
            NUMELT = MCN(NADR+4)
         ENDIF
         ISENS  = ISENS*SIGN(1,NUMELT)
         NUMELT = ABS(NUMELT)
C
C        REORIENTATION DE NUMELT
         DO 200 I=1,4
            IF( NSOM1 .EQ. NOSOEL(I,NUMELT) ) GOTO 210
  200    CONTINUE
  210    CONTINUE
         N1 = NOSOEL( I, NUMELT )
         N2 = NOSOEL( NOVOIS(ISENS,I),  NUMELT )
         N3 = NOSOEL( NOVOIS(0,I),      NUMELT )
         N4 = NOSOEL( NOVOIS(-ISENS,I), NUMELT )
         NOSOEL(1,NUMELT) = N1
         NOSOEL(2,NUMELT) = N2
         NOSOEL(3,NUMELT) = N3
         NOSOEL(4,NUMELT) = N4
         NSOM1 = N2
         NSOM2 = N3
         NELT  = (J-1)*NBARX+1
         INDXEF(NELT) = NUMELT
C
C        VERIFICATION QUE CET ELEMENT EST NOVOISIN DU BORD
         NGS(1) = MIN( N1, N4 )
         NGS(2) = MAX( N1, N4 )
         CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                3, L2ARET, NOAR )
         IF( NOAR .LE. 0 ) THEN
             GOTO 990
         ENDIF
         NADR = MNARET + 5*(NOAR-1)-1
         IF( MCN(NADR+5) .NE. 0 ) THEN
            GOTO 995
         ENDIF
C
C        RECHERCHE DES AUTRES ELEMENTS DE LA LIGNE
         DO 240 I=2,NBARX
            NELT = (J-1)*NBARX+I
            IF(  NSOM1 .LT. NSOM2 ) THEN
               NGS(1) = NSOM1
               NGS(2) = NSOM2
               ISENS  = -1
            ELSE
               NGS(1) = NSOM2
               NGS(2) = NSOM1
               ISENS  = 1
            ENDIF
            CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                   3, L2ARET, NOAR )
            IF( NOAR .LE. 0 ) THEN
               GOTO 990
            ENDIF
            NADR = MNARET + 5*(NOAR-1)-1
            IF( MCN(NADR+5).EQ.0) THEN
               GOTO 995
            ENDIF
            NT1 = ABS( MCN(NADR+4) )
            IF( NT1 .EQ. NUMELT ) THEN
               NUMELT = ABS( MCN(NADR+5) )
               ISENS  = ISENS*SIGN(1,MCN(NADR+5))
            ELSE
               NUMELT = NT1
               ISENS  = ISENS*SIGN(1,MCN(NADR+4))
            ENDIF
C
C           REORIENTATION DE NUMELT
            DO 220 II=1,4
               IF( NSOM1 .EQ. NOSOEL(II,NUMELT) ) GOTO 230
  220       CONTINUE
  230       CONTINUE
            N1 = NOSOEL( II, NUMELT )
            N2 = NOSOEL( NOVOIS(ISENS, II), NUMELT )
            N3 = NOSOEL( NOVOIS(0,II),      NUMELT )
            N4 = NOSOEL( NOVOIS(-ISENS,II), NUMELT )
            NOSOEL(1,NUMELT) = N1
            NOSOEL(2,NUMELT) = N2
            NOSOEL(3,NUMELT) = N3
            NOSOEL(4,NUMELT) = N4
C
C           VERIFICATION QUE L'ARETE N1-N2 EST COMMUNE A L'ELT DESSOUS
            NGS(1) = MIN( N1, N2 )
            NGS(2) = MAX( N1, N2 )
            CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                   3, L2ARET, NOAR )
            IF( NOAR .LE. 0 ) THEN
               GOTO 990
            ENDIF
            NADR = MNARET + 5*(NOAR-1)-1
            NT1  = ABS(MCN(NADR+4))
            NT2  = ABS(MCN(NADR+5))
            IF( (INDXEF((J-2)*NBARX+I) .NE. NT1 ) .AND.
     S          (INDXEF((J-2)*NBARX+I) .NE. NT2 ) ) THEN
               GOTO 995
            ENDIF
C
C           AJOUT DANS LA LISTE ET ITERATION
            INDXEF(NELT) = NUMELT
            NSOM1 = NOSOEL(2,NUMELT)
            NSOM2 = NOSOEL(3,NUMELT)
  240    CONTINUE
C
C        ON EST SUR LE DERNIER ELEMENT DE LA LIGNE
C        VERIFICATION QU'IL EST FRONTIERE
         NGS(1) = MIN( NSOM1, NSOM2 )
         NGS(2) = MAX( NSOM1, NSOM2 )
         CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                3, L2ARET, NOAR )
         IF( NOAR .LE. 0 ) THEN
            GOTO 990
         ENDIF
         NADR = MNARET + 5*(NOAR-1)-1
         IF( MCN(NADR+5) .NE. 0 ) THEN
            GOTO 995
         ENDIF
  250 CONTINUE
C
C     VERIFICATION QUE LE HAUT DE LA DERNIERE LIGNE EST FRONTIERE
C     ===========================================================
      DO 260 I=1,NBARX
         NELT   = (NBS2-2)*NBARX+I
         NUMELT = INDXEF(NELT)
         NSOM1  = NOSOEL(4,NUMELT)
         NSOM2  = NOSOEL(3,NUMELT)
         NGS(1) = MIN( NSOM1, NSOM2 )
         NGS(2) = MAX( NSOM1, NSOM2 )
         CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET),
     S                3, L2ARET, NOAR )
         IF( NOAR .LE. 0 ) THEN
            GOTO 990
         ENDIF
         NADR = MNARET + 5*(NOAR-1)-1
         IF( MCN(NADR+5) .NE. 0 ) THEN
            GOTO 995
         ENDIF
  260 CONTINUE
C
C     LA SURFACE EST BIEN STRUCTUREE = ARRANGEMENT COORDONNEES
C     ========================================================
      DO 330 J=1,NBS2-1
         DO 310 I=1,NBARX
            NELT   = (J-1)*NBARX+I
            NUMELT =INDXEF(NELT)
            NSOM   = (J-1)*NBS1+I
            DO 300 K=1,3
               XYZFIN(K,NSOM) = XYZINI(K,NOSOEL(1,NUMELT))
  300       CONTINUE
  310    CONTINUE
         NSOM = J*NBS1
         DO 320 K=1,3
            XYZFIN(K,NSOM) = XYZINI(K,NOSOEL(2,NUMELT))
  320    CONTINUE
  330 CONTINUE
      DO 350 I=1,NBARX
         NELT   = (NBS2-2)*NBARX+I
         NUMELT = INDXEF(NELT)
         NSOM   = (NBS2-1)*NBS1+I
         DO 340 K=1,3
            XYZFIN(K,NSOM) = XYZINI(K,NOSOEL(4,NUMELT))
  340    CONTINUE
  350 CONTINUE
      NSOM = NBS2*NBS1
      DO 360 K=1,3
         XYZFIN(K,NSOM) = XYZINI(K,NOSOEL(3,NUMELT))
  360 CONTINUE
      GOTO 999
C
C     SORTIES
C     =======
  990 WRITE(IMPRIM,*) 'ARETE NON RETROUVEE DANS LE HACHAGE'
C
  991 NBLGRC(NRERR) = 3
      KERR(1) = 'UNE ARETE N''EST PAS RETROUVEE'
      WRITE(KERR(MXLGER)(1:24),'(2I12)') NGS(1),NGS(2)
      KERR(2) = 'SOMMETS ' // KERR(MXLGER)(1:24)
      KERR(3) = 'QUADRANGULATION NON RESTRUCTURABLE'
      CALL LEREUR
  995 IERR = 1
C
  999 CALL TNMCDS( 'ENTIER', L1ARET*L2ARET, MNARET )
      RETURN
      END
