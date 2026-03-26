      SUBROUTINE LIEX49( NTLXLI , LADEFI ,
     %                   NTARLI , MNARLI , NTSOLI , MNSOLI , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   DESTRUCTURATION D'UNE LIGNE STRUCTUREE EN LIGNE NON STRUCTUREE
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC  PARIS   DECEMBRE 1988
C.......................................................................
      IMPLICIT INTEGER (W)
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
C     LA LIGNE A DESTRUCTURER
C     =======================
C     LE NOM DE CETTE LIGNE
      NULIIN = LADEFI(WULIIN)
C     LE TABLEAU LEXIQUE DE CETTE LIGNE
      CALL LXNLOU( NTLIGN , NULIIN , NTLXLG , MN )
C     LE TABLEAU 'NSEF' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG , 'NSEF' , NTARLI , MNARLI )
C     LA LIGNE EST ELLE STRUCTUREE ?
      NUTYLI = MCN( MNARLI + WUTYMA )
      IF( NUTYLI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'LIGNE DEJA NON STRUCTUREE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE NOMBRE DE NSEF
      NBEFOB = MCN( MNARLI + WBARSE )
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBSOLI = NBEFOB + 1
C     TEST SUR LE NOMBRE DE SOMMETS
      IF( NBSOLI .LE. 1 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBSOLI
         KERR(1) = 'NOMBRE DE SOMMETS'//KERR(MXLGER)(1:4)
         KERR(2) = 'AU LIEU DE >1'
        IERR = 1
        RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE LIGNE
      CALL LXTSOU( NTLXLG , 'XYZSOMMET' , NTSOLI , MNSOL0 )
C
C     CREATION DE LA LIGNE NON STRUCTUREE
C     ===================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     ---------------------------------
      MOTS = WYZSOM + 3 * NBSOLI
      CALL LXTNDC( NTLXLI , 'XYZSOMMET' , 'MOTS' , MOTS )
      CALL LXTSOU( NTLXLI , 'XYZSOMMET' ,  NTSOLI  , MNSOLI )
C     COPIE DU TABLEAU 'XYZSOMMET' STRUCTURE DANS CELUI NON STRUCTURE
      CALL TRTATA( MCN(MNSOL0) , MCN(MNSOLI) , MOTS )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOLI) )
C
C     CONSTRUCTION DU TABLEAU 'NSEF'
C     ------------------------------
      CALL LXTNDC( NTLXLI , 'NSEF' , 'ENTIER' ,
     S             WUSOEF + 2 * NBEFOB )
      CALL LXTSOU( NTLXLI , 'NSEF' ,  NTARLI  , MNARLI )
C
C     LE TYPE DE L'OBJET : ICI LIGNE
      MCN( MNARLI + WUTYOB ) = 2
C     LE TYPE DU MAILLAGE : ICI SEGMENT NON STRUCTURE
      MCN( MNARLI + WUTYMA ) = 0
C     LE NOMBRE DE SOMMETS PAR SOUS-OBJET
      MCN( MNARLI + WBSOEF ) = 2
C     LE NOMBRE DE NSEF
      MCN( MNARLI + WBEFOB ) = NBEFOB
C
      IA = MNARLI + WUSOEF - 1
      DO 10 I=1,NBSOLI-1
         MCN( IA + 1 ) = I
         MCN( IA + 2 ) = I + 1
         IA = IA + 2
 10   CONTINUE
C
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNARLI + WUTFMA ) = 0
C     PAS DE TANGENTES STOCKEES
      MCN( MNARLI + WBTGEF ) = 0
      MCN( MNARLI + WBEFTG ) = 0
      MCN( MNARLI + WBEFAP ) = 0
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARLI + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      RETURN
      END
