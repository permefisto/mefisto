        SUBROUTINE LISTFE( NBLGST, MNLGST,
     %                     LAPERM, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    VERIFIER SI LES NBLGST LIGNES STRUCTUREES CONSTITUENT UNE
C -----    LIGNE FERMEE
C          CREATION DE LA PERMUTATION POUR PARCOURIR LES
C          LIGNES STRUCTUREES DE LA LIGNE
C
C ENTREES:
C --------
C NBLGST : NOMBRE DE LIGNES STRUCTUREES DE LA LIGNE FERMEE
C          CES LIGNES EXISTENT ET SONT SANS ERREUR (CF SP LISTRE)
C MNLGST : ADRESSE MCN DU TMS 'XYZSOMMET' DE CHAQUE LIGNE STRUCTUREE
C
C SORTIES:
C --------
C LAPERM : TABLEAU DES PERMUTATIONS POUR PARCOURIR LA LIGNE FERMEE
C          <0 SIGNIFIE UN SENS CONTRAIRE DE PARCOURS
C IERR   : 0 SI PAS D'ERREUR, 1 SI LIGNE NON FERMEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       SEPTEMBRE 1993
C234567..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      INTEGER           MNLGST(1:NBLGST), LAPERM(1:NBLGST)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     INITIALISATION DES PERMUTATIONS
C     ===============================
      DO 10 I=1,NBLGST
         LAPERM(I) = I
 10   CONTINUE
      IERR = 0
C
C     LE PREMIER POINT DE LA LIGNE 1
      NBS = MCN( MNLGST(1) + WNBSOM )
      MN1 = MNLGST(1) + WYZSOM + NBS * 3 - 3
C
      DO 90 I=1,NBLGST-1
C
C        LE PARCOURS POUR LA LIGNE SUIVANTE DU CONTOUR FERME
         DO 80 J=I+1,NBLGST
            L   = LAPERM( J )
            NBS = MCN( MNLGST(L) + WNBSOM )
            MN2 = MNLGST(L) + WYZSOM
            MN3 = MN2 + NBS * 3 - 3
            CALL XYZIDE( RMCN(MN1), RMCN(MN2), IDENTQ )
            IF( IDENTQ .NE. 0 ) THEN
C              LA LIGNE J EST LA SUITE DE LA LIGNE I
               L           = LAPERM( J )
               LAPERM( J ) = LAPERM(I+1)
               LAPERM(I+1) = L
               MN1 = MN3
               GOTO 90
            ENDIF
C
C           L'AUTRE BOUT DE LA LIGNE
            CALL XYZIDE( RMCN(MN1), RMCN(MN3), IDENTQ )
            IF( IDENTQ .NE. 0 ) THEN
C              LA LIGNE J EST LA SUITE DE LA LIGNE I
               L           = LAPERM( J )
               LAPERM( J ) = LAPERM(I+1)
               LAPERM(I+1) = -L
               MN1 = MN2
               GOTO 90
            ENDIF
 80      CONTINUE
C
C        PAS DE LIGNE SUIVANTE
         GOTO 9999
 90   CONTINUE
C
C     VERIFICATION DE LA FERMETURE DE LA LIGNE
C     ========================================
      MN2 = MNLGST(1) + WYZSOM
      CALL XYZIDE( RMCN(MN1), RMCN(MN2), IDENTQ )
      IF( IDENTQ .NE. 0 ) THEN
C        FERMETURE VERIFIEE
         RETURN
      ENDIF
C
C     ERREUR : LIGNE NON FERMEE
      I = NBLGST
 9999 NBLGRC(NRERR) = 2
      KERR(1) = 'CONTOUR NE FORMANT PAS UNE LIGNE FERMEE'
      WRITE(KERR(MXLGER)(1:6),'(I6)') I
      KERR(2) = 'COTE ' // KERR(MXLGER)(1:6) // ' SANS LIGNE SUIVANTE'
      CALL LEREUR
      IERR = 1
      END
