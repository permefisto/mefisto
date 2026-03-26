      SUBROUTINE TRCOEF( NMTCL,  NOFORE, KNOMOB, NTLXOB, MODECO, NDIM,
     %                   NBTYEL, MNTOPO, MNNPEF, MNPOGE,
     %                   NCAS,   NOPT,   CONTMX, CMFLEC, CMPCON)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FLECHES REPRESENTANT LES CONTRAINTES ET
C -----    DIRECTIONS PRINCIPALES DANS TOUS LES EF D'UN OBJET
C
C ENTREES :
C ---------
C NMTCL  : NO DU MOT CLE ACTIF
C NOFORE : 0 SI PAS DE FONCTION LU REGION POUR LE TRACE
C         >0 NO DE LA FONCTION REGION
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET A TRAITER
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT DES DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NDIM   : DIMENSION DE L'ESPACE =2 OU 3
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS FINIS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS DEPLACEMENTS
C
C MODIFIES EVENTUELLEMENT :
C -------------------------
C NCAS   : NUMERO DU CAS A TRAITER
C NOPT   : NUMERO DE L'OPTION DE TRAITEMENT
C CONTMX : MAXIMUM EN VALEUR ABSOLUE DES CONTRAINTES
C CMFLEC : NOMBRE DE CM DU DE LA FLECHE MAXIMALE
C CMPCON : NOMBRE DE CM POUR L'UNITE DE CONTRAINTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      REAL              CONTMX
ccc      REAL              HEXSEC(6,2)
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      CHARACTER*4       NOMELE(2)
C
      MOREE2 = MOTVAR(6)
      IF( NMTCL .EQ. 9 ) GOTO 330
C
C     TRACE des CONTRAINTES dans tous les EF
C     --------------------------------------
C     AGRANDIR REDUIRE le TRACE
 300  CALL VISEE(NMTCL)
      IF( NMTCL .LT. 0 ) GOTO 9000
C
      IF( NDIM  .EQ. 2 ) THEN
C        INITIALISATION DU ZOOM DEPLACEMENT
         CALL ZOOM2D0( NOTYEV )
      ELSE
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
      ENDIF
      IF( NOTYEV .EQ. 0 ) GOTO 300
      GOTO 304
C
C     EN 2D: ZOOM OU TRANSLATION ACTIFS
 302  CALL ZOOM2D1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 300
      GOTO 304
C
C     EN 3D: ORBITE OU ZOOM OU TRANSLATION ACTIFS
 303  CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 300
C
C     EXECUTION DU TRACE DES CONTRAINTES
 304  IF( NOPT .NE. 1 .OR. CMFLEC .GT. 0. ) GOTO 307
C
C     ERREUR
 305  WRITE(KERR(MXLGER-3)(1:13),'(I13)')   NOPT
      WRITE(KERR(MXLGER-2)(1:13),'(G13.5)') CMFLEC
      WRITE(KERR(MXLGER-1)(1:13),'(G13.5)') CONTMX
      WRITE(KERR(MXLGER)(1:13),'(G13.5)')   CMPCON
      NBLGRC(NRERR) = 4
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR :  OPTION            '
     &           // KERR(MXLGER-3)(1:13)
         KERR(2) = 'FLECHE MAXIMALE EN CM       '
     &           // KERR(MXLGER-2)(1:13)
         KERR(3) = 'CONTRAINTE D''UN CM DE TRACE'
     &           // KERR(MXLGER-1)(1:13)
         KERR(4) = 'CM / CONTRAINTE             '
     &           // KERR(MXLGER)(1:13)
      ELSE
         KERR(1) = 'ERROR :  OPTION          '
     &           // KERR(MXLGER-3)(1:13)
         KERR(2) = 'MAXIMALE ARROW CM        '
     &           // KERR(MXLGER-2)(1:13)
         KERR(3) = 'STRESS of 1 CM of DRAWING'
     &           // KERR(MXLGER-1)(1:13)
         KERR(4) = 'CM / STRESS              '
     &           // KERR(MXLGER)(1:13)
      ENDIF
      CALL LEREUR
      GOTO 9000
C
 307  IF( NOPT .EQ. 2 .AND. CONTMX .LE. 0. ) GOTO 305
      IF( NOPT .EQ. 3 .AND. CMPCON .LE. 0. ) GOTO 305
ccc 310  IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,10310) NCAS
ccc      ELSE
ccc         WRITE(IMPRIM,20310) NCAS
ccc      ENDIF
ccc10310 FORMAT(' CAS DE CHARGE VISUALISE =',I5/)
ccc20310 FORMAT(' DRAWING of CASE =',I5/)
C
C     MISE A JOUR DE L'ECHELLE
      IF( NOPT .EQ. 1 ) THEN
         CMPCON = CMFLEC / CONTMX
      ELSE IF( NOPT .EQ. 2 ) THEN
         CMPCON = 1. / CONTMX
      ENDIF
C
C     LA PREPARATION DU TRACE
C
C     LE CODE DE TRAITEMENT DE L'OBJET
C     NDPGST : CODE TRAITEMENT DE L ELEMENT
C                 0 : NOEUDS=POINTS=SOMMETS
C                 1 : NOEUDS=POINTS#SOMMETS
C                 2 : NOEUDS#POINTS=SOMMETS
C                 3 : NOEUDS#POINTS#SOMMETS
 330  NDPGST = MCN( MNTOPO + WDPGST )
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS DU MAILLAGE
      DO 400 I = 0, NBTYEL-1
C        L'ADRESSE MCN DU DEBUT DU TABLEAU ELEMENTS
         MNELE = MCN( MNNPEF + I )
C        LE NOM DU TABLEAU CONTRAINTES ASSOCIE
         NUTYEL = MCN( MNELE + WUTYEL )
         IF( NUTYEL .LE. 4 ) THEN
            NOAXIS=1
         ELSE
            NOAXIS=0
         ENDIF
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         KNOM = 'CONTRAINTE"' // NOMELE(2)
C        OUVERTURE DU TABLEAU
         CALL LXTSOU( NTLXOB, KNOM, NTCONT, MNCONT )
         IF( NTCONT .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'OBJET SANS CONTRAINTES CALCULEES'
            ELSE
               WRITE(IMPRIM,*) 'OBJECT WITHOUT COMPUTED STRESSES'
            ENDIF
            GOTO 400
         ENDIF
C        LES CARACTERISTIQUES DE L'ELEMENT
         CALL ELTYCA( NUTYEL )
C
C        TRACE EFFECTIF DES FLECHES
C        LES DIFFERENTS TABLEAUX ET VARIABLES DE 'CONTRAINTE'
C        AFIN DE BENEFICIER DES INDICES FORTRAN EN CLAIR
         NBELFI = MCN(MNCONT+WBELFI)
         NBPIEF = MCN(MNCONT+WBPIEF)
         NBJECA = MCN(MNCONT+WBJECA)
         NDIMES = MCN(MNCONT+WDIMES)
         LDCOPR = MCN(MNCONT+WDCOPR)
         MNDIPR = MNCONT+LDCOPR+NBELFI*NBPIEF*NDIMES*NBJECA*MOREE2
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + MCN(MNELE+WBELEM) * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CONTRAINTES PRINCIPALES SELON LA DIMENSION DE L'OBJET
         IF( NMTCL .EQ. 9 ) THEN
C           AFFICHAGE DES CONTRAINTES POUR CE TYPE D'EF
            CALL AFCONT( NOFORE, NOAXIS, NOMELE,
     %                   NBELFI, NBPIEF, NDIMES, NBJECA ,
     %                   MCN(MNCONT+WOPIEF), MCN(MNCONT+LDCOPR),
     %                   MCN(MNDIPR), NCAS )
         ELSE
            IF( NDIMES .EQ. 2 ) THEN
C              TRACE DES CONTRAINTES D'UN OBJET DE DIMENSION 2 OU AXISYMETRIQUE
               CALL TRAXE2
               CALL TRCON2( NOFORE, NBELFI, NBPIEF, NDIMES, NBJECA ,
     %                      MCN(MNCONT+WOPIEF), MCN(MNCONT+LDCOPR),
     %                      MCN(MNDIPR), MCN(MNPOGE+WYZPOI) ,
     %                      MCN(MNPGEL), NCAS, CMPCON )
               ELSE
C              TRACE DES CONTRAINTES PRINCIPALES D'UN OBJET DE DIMENSION 3
               CALL TRAXE3
                  CALL TRCON3( NOFORE, NBELFI, NBPIEF, NDIMES, NBJECA ,
     %                      MCN(MNCONT+WOPIEF), MCN(MNCONT+LDCOPR),
     %                      MCN(MNDIPR), MCN(MNPOGE+WYZPOI) ,
     %                      NCAS, CMPCON , KNOMOB)
            ENDIF
            ENDIF
 400  CONTINUE
      IF( NMTCL .EQ. 9 ) RETURN
C
C     LE TRACE DU TITRE
      CALL XVCOULEUR( NCNOIR )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'OBJET: ' // KNOMOB
      ELSE
         KNOM = 'OBJECT: ' // KNOMOB
      ENDIF
      I = NUDCNB( KNOM )
      IF( MODECO .EQ. 1 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KNOM(I+1:I+9) = ' au TEMPS'
         ELSE
            KNOM(I+1:I+9) = ' at TIME'
         ENDIF
         WRITE( KNOM(I+10:I+23), '(G14.6)' ) TEMPS
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            KNOM(I+1:I+12) = '  FREQUENCE '
         ELSE
            KNOM(I+1:I+12) = '  FREQUENCY '
         ENDIF
         WRITE( KNOM(I+13:I+26), '(G14.6)' ) TEMPS
         KNOM(I+27:I+30) = ' Hz '
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 50, 60 )
C
ccc      X = (HEXSEC(1,1)+HEXSEC(1,2)) * 0.5
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = '1CM DE FLECHE= '
         WRITE( KNOM(15:27),'(G13.5)') 1.0/CMPCON
         KNOM(29:51)='UNITES DE CONTRAINTES'
      ELSE
         KNOM = '1CM of ARROW= '
         WRITE( KNOM(15:27),'(G13.5)') 1.0/CMPCON
         KNOM(29:51)='UNITIES of STRESS'
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 50, 75 )
C
C     FIN DU TRACE
      IAVTIT = 1
      IF( LANGAG .EQ. 0 ) THEN
         CALL TRFINS( 'Les CONTRAINTES et DIRECTIONS PRINCIPALES' )
      ELSE
         CALL TRFINS( 'The MAIN STRESSES and DIRECTIONS' )
      ENDIF
      IF( LORBITE .NE. 0 ) THEN
         IF( NDIM .EQ. 2 ) THEN
            GOTO 302
         ELSE
            GOTO 303
         ENDIF
      ENDIF
C
      CALL CLICSO
      GOTO 300
C
 9000 RETURN
      END
