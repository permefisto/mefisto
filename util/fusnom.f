      SUBROUTINE FUSNOM( KNOMTD, KNOMTS, KNOMFU, NCODER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FUSIONNER LE NOM DE TD AVEC LE NOM TS
C ----- EXEMPLE : TS : ~>POINT>P1>DEFINITION
C                 TD : ~>TRANSFO
C                 FU : ~>TRANSFO
C       EN L'ABSENCE DE >> FU = TD
C
C       OU        TS : ~>POINT>P1>DEFINITION"SUIVANT
C                 TD : ~>POINT>>XYZSOMMET
C                 FU : ~>POINT>P1>XYZSOMMET"SUIVANT
C       ON COMPLETE >> DU TD PAR LE NOM CORRESPONDANT ENTRE >> DU TS
C                 ET AJOUT de "... du TS dans FU
C ENTREE :
C --------
C KNOMTD : NOM DU TABLEAU DESCRIPTEUR ASSOCIE AU TABLEAU MS
C KNOMTS : NOM DU TABLEAU MS
C
C SORTIES :
C ---------
C KNOMFU  : NOM DU TABLEAU APRES FUSION
C NCODER  : CODE D'ERREUR 1 SI KNOMFU = KNOMTD
C                         2 SI KNOMFU = KNOMTD + KNOMTS
C                         0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      CHARACTER*(*)     KNOMTD,KNOMTS,KNOMFU
      CHARACTER*80      KNOMTM
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     PASSAGE AUX MAJUSCULES
      KNOMTM = KNOMTD
      CALL MAJUSC( KNOMTM )
C
C     EXISTE-T-IL >> ?
      LTS = 0
      LTD = INDEX( KNOMTM , '>>' )
      IF( LTD .LE. 0 ) THEN
         KNOMFU = KNOMTM
         NCODER = 1
         GOTO 8000
      ENDIF
C
C     INITIALISATION A BLANC
      KNOMFU = ' '
C
C     LE NOMBRE DE CARACTERES DECLARES DES 3 NOMS
      LTD   = LEN( KNOMTM )
C
C     LA POSITION DU DERNIER CARACTERE NON BLANC DE KNOMTM
      L1BTD = INDEX( KNOMTM , ' ' )
      IF( L1BTD .LE. 0 ) THEN
C        PAS DE BLANC . DERNIER CARACTERE = LONGUEUR
         L1BTD = LTD
      ELSE
C        DERNIER CARACTERE AVANT LE 1-ER BLANC
         L1BTD = L1BTD - 1
      ENDIF
      LTS = LEN( KNOMTS )
      LFU = LEN( KNOMFU )
C
C     LE CARACTERE COURANT DES 3 NOMS
      NCTD = 0
      NCTS = 0
      NCFU = 0
C
 10   NCTD = NCTD + 1
      NCTS = NCTS + 1
C     LE DERNIER CARACTERE DE TD EST IL ATTEINT ?
      IF( NCTD .GT. L1BTD ) THEN
         NCODER = 2
         GOTO 8000
      ENDIF
C
      IF( KNOMTM(NCTD:NCTD) .EQ. KNOMTS(NCTS:NCTS) ) THEN
C        MEME CARACTERE TD ET TS . AJOUT DANS KNOMFU
         NCFU = NCFU + 1
         KNOMFU(NCFU:NCFU) = KNOMTM(NCTD:NCTD)
         GOTO 10
      ENDIF
C
C     LES 2 CARACTERES SONT DIFFERENTS
      IF( (KNOMTM(NCTD-1:NCTD)   .EQ. '>>') .AND.
     %    (KNOMTS(NCTS-1:NCTS-1) .EQ. '>' ) ) THEN
C        DANS TD >> ET DANS TS >NOM> . COPIE DE NOM DANS FU
         LAFIN = INDEX( KNOMTS(NCTS:LTS) , '>' )
         IF( LAFIN .LE. 0 ) GOTO 9000
         KNOMFU(NCFU+1:NCFU+LAFIN) = KNOMTS(NCTS:NCTS-1+LAFIN)
         NCFU = NCFU + LAFIN
         NCTS = NCTS + LAFIN - 1
         GOTO 10
      ELSE
         LAFIN = INDEX( KNOMTS(NCTS:LTS) , '>' )
         IF( LAFIN .LE. 0 ) THEN
C           DERNIER NOM. COPIE DU NOM TD DANS FU
            KNOMFU(NCFU+1:LFU) = KNOMTM(NCTD:L1BTD)
            NCODER = 2
         ELSE
            GOTO 9000
         ENDIF
      ENDIF
C
C     AJOUT DU TERME FINAL EVENTUEL "XXXXX de TS dans FU
 8000 L1BTD = INDEX( KNOMTS , '"' )
      IF( L1BTD .GT. 0 ) THEN
         LFU = INDEX( KNOMFU , ' ' )
         KNOMFU(LFU:LFU+LTS-L1BTD) = KNOMTS(L1BTD:LTS)
      ENDIF
      RETURN
C
C     ERREUR
C     ======
 9000 NBLGRC(NRERR) = 3
      KERR(1) = KNOMTM
      KERR(2) = KNOMTS
      KERR(3) ='NON FUSIONNABLES'
      CALL LEREUR
      NCODER = 0
      RETURN
      END
