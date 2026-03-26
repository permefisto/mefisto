      SUBROUTINE dft0t1dt( TPSINI, TPSFIN, DT, NBVESTK, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ENTRER les TEMPS pour une RESOLUTION d'un PROBLEME
C -----    INSTATIONNAIRE D'ORDRE 1 en TEMPS

C SORTIES:
C --------
C TPSINI : TEMPS INITIAL DU CALCUL
C TPSFIN : TEMPS FINAL   DU CALCUL
C DT     : PAS DE TEMPS CONSTANT DU TEMPS ENTRE 2 ITERATIONS
C NBVESTK: NOMBRE DE NOUVEAUX VECTEURS SOLUTION A STOCKER ENTRE
C          TPSINI et TPSFIN
C         (le VECTEUR de TPSINI est AJOUTE SI CE N'EST PAS DEJA FAIT)
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2022
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      REAL           TPSINI, TPSFIN, DT

C     ENTREE DU TEMPS INITIAL : TPSINI
C     ================================
 10   CALL INVITE( 95 )
      NCVALS = 0
      CALL LIRRSP( NCVALS, TPSINI )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        GOTO 9999
      ENDIF

C     ENTREE DU TEMPS FINAL : TPSFIN
C     ==============================
      CALL INVITE( 94 )
      NCVALS = 0
      CALL LIRRSP( NCVALS, TPSFIN )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        GOTO 9999
      ENDIF
      IF( TPSFIN .LE. TPSINI ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: TEMPS FINAL INCORRECT'
            WRITE(KERR(5)(1:15),'(G15.7)') TPSINI
            WRITE(KERR(6)(1:15),'(G15.7)') TPSFIN
            KERR(2) = 'TEMPS INITIAL ' //  KERR(5)(1:15)
            KERR(2) = 'TEMPS FINAL   ' //  KERR(6)(1:15)
         ELSE
            KERR(1) = 'ERROR: INCORRECT FINAL TIME'
            WRITE(KERR(5)(1:15),'(G15.7)') TPSINI
            WRITE(KERR(6)(1:15),'(G15.7)') TPSFIN
            KERR(2) = 'INITIAL TIME ' //  KERR(5)(1:15)
            KERR(2) = 'FINAL   TIME ' //  KERR(6)(1:15)
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     ENTREE DU PAS DE TEMPS CONSTANT INITIAL: DT
C     ===========================================
 15   CALL INVITE( 86 )
      NCVALS = 0
      CALL LIRRSP( NCVALS, DT )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        GOTO 9999
      ENDIF
      IF( DT .LE. 0 ) THEN
        NBLGRC(NRERR) = 2
        WRITE(KERR(5)(1:15),'(G15.7)') DT
        IF( LANGAG .EQ. 0 ) THEN
           KERR(1) = 'ERREUR: PAS DE TEMPS <=0 INCORRECT'
           KERR(2) = 'PAS DE TEMPS ' // KERR(5)(1:15)
        ELSE
           KERR(1) = 'ERROR: INCORRECT TIME STEP <=0'
           KERR(2) = 'STEP TIME ' // KERR(5)(1:15)
        ENDIF
        CALL LEREUR
        GOTO 15
      ENDIF

C     ENTREE DU NOMBRE MAXIMAL DE NOUVEAUX VECTEURS SOLUTION A STOCKER
C     ================================================================
 20   CALL INVITE( 76 )
      NBVESTK = 2
      NCVALS  = 4
      CALL LIRENT( NCVALS, NBVESTK )
      IF( NCVALS .LT. 0 ) THEN
C       ABANDON DE LA LECTURE DES DONNEES
        IERR = 9
        GOTO 9999
      ENDIF
      IF( NBVESTK .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'STOCKAGE MINIMAL VITESSE PRESSION >=2 VECTEURS'
            KERR(2) = 'AU TEMPS MINIMAL ET MAXIMAL'
         ELSE
            KERR(1) = 'MINIMAL STORAGE VELOCITIES-PRESSURES >=2 VECTORS'
            KERR(2) = 'for MINIMAL and MAXIMAL TIME'
         ENDIF
         CALL LEREUR
         GOTO 20
      ENDIF

 9999 RETURN
      END
