      SUBROUTINE MT4SQA( NA,  MOARTR, NOARTR, MOSOAR, NOSOAR,
     %                   NS1, NS2, NS3, NS4 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU NUMERO DES 4 SOMMETS DE L'ARETE NA DE NOSOAR
C -----    FORMANT UN QUADRANGLE
C
C ENTREES:
C --------
C NA     : NUMERO DE L'ARETE DANS NOSOAR A TRAITER
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1=0 SI TRIANGLE VIDE => ARETE2=TRIANGLE VIDE SUIVANT
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES, CHAINAGES EN +
C          SOMMET 1 = 0 SI ARETE VIDE => SOMMET 2 = ARETE VIDE SUIVANTE
C
C SORTIES:
C --------
C NS1,NS2,NS3 : LES 3 NUMEROS DES SOMMETS DU TRIANGLE T1 EN SENS DIRECT
C NS1,NS4,NS2 : LES 3 NUMEROS DES SOMMETS DU TRIANGLE T2 EN SENS DIRECT
C
C SI ERREUR RENCONTREE => NS4 = 0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1997
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NOARTR(MOARTR,*), NOSOAR(MOSOAR,*)
C
C     LE NUMERO DE TRIANGLE EST IL CORRECT  ?
C     A SUPPRIMER APRES MISE AU POINT
      IF( NA .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:6),'(I6)') NA
         KERR(1) = KERR(MXLGER)(1:6) //
     %           ' NO INCORRECT ARETE DANS NOSOAR'
         CALL LEREUR
         NS4 = 0
         RETURN
      ENDIF
C
      IF( NOSOAR(1,NA) .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:6),'(I6)') NA
         KERR(1) = KERR(MXLGER)(1:6) //
     %           ' ARETE NON ACTIVE DANS NOSOAR'
         CALL LEREUR
         NS4 = 0
         RETURN
      ENDIF
C
C     RECHERCHE DE L'ARETE NA DANS LE PREMIER TRIANGLE
      NT = NOSOAR(4,NA)
      IF( NT .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:6),'(I6)') NA
         KERR(1) =  'TRIANGLE 1 INCORRECT POUR L''ARETE ' //
     %               KERR(MXLGER)(1:6)
         CALL LEREUR
         NS4 = 0
         RETURN
      ENDIF
C
      DO 5 I=1,3
         IF( ABS( NOARTR(I,NT) ) .EQ. NA ) GOTO 8
 5    CONTINUE
C     SI ARRIVEE ICI => BOGUE AVANT
      WRITE(IMPRIM,*) 'MT4SQA: ARETE',NA,' NON ARETE DU TRIANGLE',NT,
     %' D''ARETES:',(NOARTR(k,NT),k=1,3)
      NS4 = 0
      RETURN
C
C     LES 2 SOMMETS DE L'ARETE NA
 8    IF( NOARTR(I,NT) .GT. 0 ) THEN
         NS1 = 1
         NS2 = 2
      ELSE
         NS1 = 2
         NS2 = 1
      ENDIF
      NS1 = NOSOAR(NS1,NA)
      NS2 = NOSOAR(NS2,NA)
C
C     L'ARETE SUIVANTE
      IF( I .LT. 3 ) THEN
         I = I + 1
      ELSE
         I = 1
      ENDIF
      NAA = ABS( NOARTR(I,NT) )
C
C     LE SOMMET NS3 DU TRIANGLE 123
      NS3 = NOSOAR(1,NAA)
      IF( NS3 .EQ. NS1 .OR. NS3 .EQ. NS2 ) THEN
         NS3 = NOSOAR(2,NAA)
      ENDIF
C
C     LE TRIANGLE DE L'AUTRE COTE DE L'ARETE NA
C     =========================================
      NT = NOSOAR(5,NA)
      IF( NT .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:6),'(I6)') NA
         KERR(1) =  'TRIANGLE 2 INCORRECT POUR L''ARETE ' //
     %               KERR(MXLGER)(1:6)
         CALL LEREUR
         NS4 = 0
         RETURN
      ENDIF
C
C     LE NUMERO DE L'ARETE NAA DU TRIANGLE NT
      NAA = ABS( NOARTR(1,NT) )
      IF( NAA .EQ. NA ) NAA = ABS( NOARTR(2,NT) )
      NS4 = NOSOAR(1,NAA)
      IF( NS4 .EQ. NS1 .OR. NS4 .EQ. NS2 ) THEN
         NS4 = NOSOAR(2,NAA)
      ENDIF
      END
