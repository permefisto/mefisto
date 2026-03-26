      SUBROUTINE LIEX07( NTLXLI , LADEFI , RADEFI ,
     %                   NTARLI , MNARLI , NTSOLI , MNSOLI , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES ARETES DE LA LIGNE DEFINIE PAR UNE LISTE DE
C -----    3 COORDONNEES ( 1 ARETE ENTRE 2 SOMMETS SUCCESSIFS )
C
C ENTREES:
C --------
C NTLXLI : NUMERO DU TABLEAU TMS DU LEXIQUE DE LA LIGNE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA LIGNE
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA LIGNE
C          CF ~/td/d/a_ligne__definition
C
C SORTIES:
C --------
C NTARLI : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C MNARLI : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C          CF ~/td/d/a___nsef
C NTSOLI : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C          CF ~/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC  PARIS   DECEMBRE 1989
C...............................................................................
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
C     LE NOMBRE DE POINTS DE LA LIGNE
      NBPOLI = LADEFI( WBPOLI )
C
C     DONNEES CORRECTES?
      IF( NBPOLI .LE. 1 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NBPOLI
         KERR(1) ='NOMBRE INCORRECT DE POINTS'//KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
C     COMPARAISON DES 3 COORDONNEES DU PREMIER ET DERNIER POINT DE LA LIGNE
      MNC = WYZPLI
      MNS = WYZPLI + 3 * NBPOLI - 3
      DO 20 K = 0 , 2
         R1 = ABS( RADEFI( MNC + K ) )
         R2 = ABS( RADEFI( MNS + K ) )
         IF     ( R1 .LE. EPZERO ) THEN
            IF  ( R2 .GT. EPZERO ) GOTO 30
         ELSE IF( R2 .LE. EPZERO ) THEN
            GOTO 30
         ELSE IF( ABS( RADEFI(MNC+K)-RADEFI(MNS+K) ) .GT.
     %              R1 * EPSXYZ    ) THEN
            GOTO 30
         ENDIF
 20   CONTINUE
C     LES PREMIER ET DERNIER POINTS SONT IDENTIQUES : LIGNE FERMEE
      LIGFER = 1
C     LE NOMBRE DE SOMMETS DE LA LIGNE
      NBPOLI = NBPOLI - 1
      NBARLI = NBPOLI
      MOARLI = WUSOEF + 2 * NBARLI
      GOTO 40
C
C     LES PREMIER ET DERNIER POINTS SONT DIFFERENTS : LIGNE OUVERTE
  30  LIGFER = 0
      NBARLI = NBPOLI - 1
      MOARLI = 1 + WBARSE
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     ---------------------------------
  40  CALL LXTNDC( NTLXLI , 'XYZSOMMET' , 'MOTS' , WYZSOM + 3*NBPOLI )
      CALL LXTSOU( NTLXLI , 'XYZSOMMET' ,  NTSOLI  , MNSOLI )
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOLI + WNBSOM ) = NBPOLI
C
C     PAS DE TANGENTES
      MCN( MNSOLI + WNBTGS ) = 0
C
C     GENERATION DES NBPOLI POINTS DANS SOMMETS
      CALL TRTATA( RADEFI(WYZPLI) , RMCN(MNSOLI+WYZSOM) , 3*NBPOLI )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOLI + WBCOOR ) = 3
      MCN( MNSOLI + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' LIGNE STRUCTURE OU NON
C     -----------------------------------------------------------
      CALL LXTNDC( NTLXLI , 'NSEF' , 'MOTS'   , MOARLI )
      CALL LXTSOU( NTLXLI , 'NSEF' ,  NTARLI  , MNARLI )
C
C     LE TYPE DE L'OBJET : ICI LIGNE
      MCN( MNARLI + WUTYOB ) = 2
C
      IF( LIGFER .EQ. 0 ) THEN
C
C        LIGNE OUVERTE => LIGNE STRUCTUREE
         MCN( MNARLI + WUTYMA ) = 2
C        LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
         MCN( MNARLI + WUTFMA ) = 0
C        LE NOMBRE DE SOMMETS ET TANGENTES PAR ARETE
         MCN( MNARLI + WBSOEF ) = 2
         MCN( MNARLI + WBTGEF ) = 0
C        LE NOMBRE D'ARETES DU SEGMENT STRUCTURE
         MCN( MNARLI + WBEFOB ) = NBPOLI - 1
         MCN( MNARLI + WBARSE ) = NBPOLI - 1
C
      ELSE
C
C        LIGNE FERMEE => LIGNE NON STRUCTUREE
         MCN( MNARLI + WUTYMA ) = 0
C        LE TYPE FERME DE FERMETURE DU MAILLAGE
         MCN( MNARLI + WUTFMA ) = 1
C        LE NOMBRE DE SOMMETS DE CHAQUE SOUS OBJET SEGMENT OU ARETE
         MCN( MNARLI + WBSOEF ) = 2
C        PAS DE TANGENTES STOCKEES
         MCN( MNARLI + WBTGEF ) = 0
C        LE NOMBRE D'ARETES DE LA LIGNE NON STRUCTUREE
         MCN( MNARLI + WBEFOB ) = NBARLI
         MNS = MNARLI + WUSOEF
         DO 50 K = 1 , NBARLI
            MCN( MNS     ) = K
            MCN( MNS + 1 ) = K+1
            MNS = MNS + 2
 50      CONTINUE
C        LIGNE FERMEE
         MCN( MNS - 1 ) = 1
      ENDIF
C
C     PAS D'EF A TG
      MCN( MNARLI + WBEFTG ) = 0
      MCN( MNARLI + WBEFAP ) = 0
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARLI + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      IERR = 0
C
C     ERREUR
C     ======
 9999 RETURN
      END
