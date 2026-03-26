      SUBROUTINE LIEX40( NTLXLI , LADEFI ,
     %                   NTARLI , MNARLI , NTSOLI , MNSOLI , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    STRUCTURATION D'UNE LIGNE NON STRUCTUREE
C -----
C
C ENTREES:
C --------
C NTLXLI : NUMERO DU TABLEAU TS DU LEXIQUE DE LA LIGNE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA LIGNE
C          CF ~td/d/a_ligne__definition
C
C SORTIES:
C --------
C NTARLI : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C MNARLI : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C          CF ~td/d/a___nsef
C NTSOLI : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C IERR   : 0 SI PAS D'ERREUR
C          1 SI NOMBRE DE SOMMETS INCORRECT <2
C          2 POINT INITIAL OU FINAL NON INITIALISE
C          3 POINT INITIAL ET FINAL CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC  PARIS   DECEMBRE 1988
C...............................................................................
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           LADEFI(0:*)
C
      IERR = 0
C
C     LA LIGNE A STRUCTURER
C     =====================
C     LE NUMERO DE CETTE LIGNE
      NULIIN = LADEFI(WULIIN)
C
C     LE TABLEAU LEXIQUE DE CETTE LIGNE
      CALL LXNLOU( NTLIGN , NULIIN , NTLXLG , MN )
      IF( NTLXLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LIGNE INCONNUE'
         ELSE
            KERR(1) = 'UNKNOWN LINE'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE TABLEAU 'NSEF' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG, 'NSEF', NTARLI, MNARLI )
C
C     LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
      CALL LXTSOU( NTLXLI, 'XYZSOMMET',  NTSOLI, MNSOLI )
C
C     LA LIGNE EST ELLE STRUCTUREE ?
      NUTYLI = MCN( MNARLI + WUTYMA )
      IF( NUTYLI .EQ. 2 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LIGNE DEJA STRUCTUREE'
         ELSE
            KERR(1) = 'ALREADY STRUCTURED LINE'
         ENDIF
         CALL LERESU
         IERR = 0
      ELSE
C        TENTATIVE DE STRUCTURATION DE LA LIGNE
         CALL LIGSTR( NTLXLI, NTARLI, MNARLI, NTSOLI, MNSOLI, IERR )
      ENDIF
      MCN( MNSOLI + WBCOOR ) = 3
C
      RETURN
      END
