      SUBROUTINE CHVATS( NLD , NCD , NLF , NCF ,
     %                   NOVATS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE A PARTIR DU CARACTERE (NLD,NCD) D'UNE VARIABLE
C ----- D'UN TABLEAU TMS
C       CE PEUT NE PAS ETRE UNE TELLE VARIABLE
C       SI CETTE VARIABLE TMS EXISTE DANS LES MS ALORS ELLE EST AJOUTEE
C
C ENTREES :
C ---------
C NLD,NCD : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C
C SORTIES :
C ---------
C NLF,NCF : SI NOVATS>0  POSITION DANS KLG DU DERNIER CARACTERE
C                        DE LA VARIABLE
C           SI NOVATS=-1 <=> FRAPPE DE @ ET NLF=-1
C           SINON NLF=0
C NOVATS : -1 FRAPPE DE @
C          =0 PAS DE VARIABLE RETROUVEE
C          >0 LE NUMERO DE LA VARIABLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/impres.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*80      NOM
C
C     RECHERCHE DU NOM
      CALL LIUTCH( NLD,NCD, NL,NC, NLF,NCF )
      IF( NLF .EQ. -1 ) THEN
C        FRAPPE DE @ => ABANDON
         NOVATS = -1
         RETURN
      ENDIF
      IF( NLF .EQ.  0 ) GOTO 9000
      IF( KLG(NLF)(NCF:NCF) .EQ. ';' ) NCF = NCF - 1
C
C     LE NOM DE L'EVENTUELLE VARIABLE TMS
C     COMPLETE EVENTUELLEMENT PAR ~>
      ND = 0
      IF( KLG(NL)(NC:NC) .EQ. '~' ) THEN
         IF( KLG(NL)(NC+1:NC+1) .EQ. '>' ) THEN
            NOM = KLG(NL)(NC:NCF)
         ELSE
            NOM = '~>' // KLG(NL)(NC+1:NCF)
            ND  = 1
         ENDIF
      ELSE
         IF( KLG(NL)(NC:NC) .EQ. '>' ) THEN
            NOM = '~' // KLG(NL)(NC:NCF)
            ND  = 1
         ELSE
            NOM = '~>' // KLG(NL)(NC:NCF)
            ND  = 2
         ENDIF
      ENDIF
C
C     RECHERCHE DU NOM DANS LES LEXIQUES
C     SUPPRESSION DE L'AFFICHAGE DU NOM DU TMS
      IMNM   = IMNMTS
      IMNMTS = 0
      CALL VATSTD( NOM , NOTYPE , MNTMS , LDTMS , NCFF )
C     SI EN SORTIE NCFF=0 CE N'EST PAS UNE VARIABLE
      IMNMTS = IMNM
      IF( MNTMS .GT. 0 .AND. NCFF .GT. 0 ) THEN
C
C        AFFICHAGE DE LA VALEUR DE LA VARIABLE
         NCFF = MAX(1,NCFF)
         CALL AFVATS( NOM(1:NCFF) , NOTYPE , MNTMS , LDTMS )
C
C        VARIABLE RETROUVEE DANS LES MS
C        CETTE VARIABLE EST ELLE STOCKEE DANS KVATS?
         DO 10 NOVATS=1,NBVATS
            NBC = INDEX( KVATS(NOVATS) , ' ' ) - 1
            IF( NBC .LE. 0 ) NBC = NBCATS
            NBC = NBC - 1
            IF( KVATS(NOVATS)(1:1+NBC) .EQ. NOM(1:NCFF) ) THEN
C              LA VARIABLE NOVATS EST RETROUVEE
               GOTO 15
            ENDIF
 10      CONTINUE
C
C        AJOUT DE LA VARIABLE TMS DANS LA TABLE
         IF( NBVATS .GE. MXVATS ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = 'LU: TABLE DES VARIABLES TMS SATUREE'
            KERR(2) = 'LU: AUGMENTER MXVATS DANS LU.INC'
            CALL LEREUR
            GOTO 9000
         ENDIF
C
C        IL RESTE DE LA PLACE
         NBVATS = NBVATS + 1
         NOVATS = NBVATS
         KVATS(NBVATS) = NOM(1:NCFF)
C
 15      NLF    = NL
         NCF    = NC - 1 - ND + NCFF
         RETURN
      ENDIF
C
C     VARIABLE NON RETROUVEE
C     ======================
 9000 NOVATS = 0
      NLF    = 0
      END
