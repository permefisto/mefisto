        SUBROUTINE LISTRE( NOLI,
     %                     NTLXLI, NTARLI, MNARLI, NTSOLI, MNSOLI,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE LEXIQUE DE LA LIGNE NOLI SUPPOSEE NON FERMEE
C -----    RESTRUCTURER EVENTUELLEMENT CETTE LIGNE
C          RETROUVER L'ADRESSE ET LE TABLEAU NSEF XYZSOMMET
C
C ENTREES:
C --------
C NOLI   : NUMERO DE LA LIGNE DANS LE LEXIQUE DES LIGNES
C
C SORTIES:
C --------
C NTLXLI : NUMERO DU TABLEAU TS DU LEXIQUE DE LA LIGNE
C NTARLI : NUMERO      DU TMS 'NSEF' DES NUMEROS DES EF
C MNARLI : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES EF
C NTSOLI : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       SEPTEMBRE 1993
C234567..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_ligne__definition.inc"
C
      CHARACTER*24      KNOMLG
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      IERR  = 0
      MNADR = 0
C
      IF( NOLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NOLI
         KERR(1) = 'LIGNE DE NUMERO' // KERR(MXLGER)(1:10)
     %          // ' INCORRECT'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE TABLEAU LEXIQUE DE CETTE LIGNE
      CALL LXNLOU( NTLIGN, NOLI, NTLXLI, N )
      IF( NTLXLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) ='LIGNE INCONNUE'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LE NOM DE LA LIGNE
      CALL NMOBNU( 'LIGNE', NOLI, KNOMLG )
C
C     LE TABLEAU 'NSEF' DE CETTE LIGNE
      CALL LXTSOU( NTLXLI, 'NSEF', NTARLI, MNARLI )
      IF( NTARLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS ARETES : ' // KNOMLG
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
      CALL LXTSOU( NTLXLI, 'XYZSOMMET', NTSOLI, MNSOLI )
      IF( NTSOLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS XYZSOMMET : ' // KNOMLG
         CALL LEREUR
         IERR = 8
         RETURN
      ENDIF
C
C     SI LA LIGNE N'EST PAS STRUCTUREE ALORS TENTATIVE DE STRUCTURATION
      CALL LIGSTR( NTLXLI, NTARLI, MNARLI, NTSOLI, MNSOLI, IERR )
      END
