      SUBROUTINE LXIM( NTLX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IMPRIMER LE LEXIQUE NTLX
C -----
C ENTREE :
C --------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS     OCTOBRE 1984
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/langue.inc"
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / MSKNOM / KNOM
      CHARACTER*160     KNOM
C
C     OUVERTURE DU TABLEAU MS = LEXIQUE NTLX
      CALL TAMSOU( NTLX , MNLX )
      IF( MNLX .LE. 0 ) RETURN
C
C     MNLX : ADRESSE MCN DU LEXIQUE
C     M1LX : NOMBRE D ENTIERS PAR NOM ET ATTRIBUTS DU LEXIQUE
      M1LX   = MCN( MNLX )
C     NBENNM : NOMBRE D'ENTIERS POUR STOCKER LES CARACTERES D'UN NOM
      NBENNM = MCN( MNLX + 2 )
C     NBCANM : NOMBRE DE CARACTERES D'UN NOM DU LEXIQUE
      NBCANM = MCN( MNLX + 3 )
C
C     IMPRESSION DES GENERALITES DU LEXIQUE
CC      WRITE(IMPRIM,10000) NTLX,(MCN(MNLX+N),N=0,6)
CC10000 FORMAT(' LEXIQUE ( TMS=',I6,' )',1X,49(1H-)/
CC     %I5,' ENTIERS PAR NOM',T40,
CC     %I5,' NOMS AU PLUS'/
CC     %I5,' ENTIERS POUR STOCKER UN NOM',T40,
CC     %I5,' CARACTERES PAR NOM'/
CC    %I5,' NUMERO TAMS DU LEXIQUE PERE'/
CC     %I5,' EST LE 1-ER NOM OCCUPE',T40,
CC     %I5,' EST LE 1-ER NOM VIDE'/
CC     %    ' CE LEXIQUE CONTIENT LES NOMS')

ccc      NBLGRC(NRERR) = NBLGRC(NRERR) + 1
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(NBLGRC(NRERR)) = 'CE LEXIQUE CONTIENT LES NOMS'
      ELSE
         KERR(NBLGRC(NRERR)) = 'This LEXICON contains the NAMES'
      ENDIF
C
C     N NUMERO DU 1-ER NOM OCCUPE DU LEXIQUE
      N = MCN( MNLX + 5 )
      IF( N .LE. 0 ) THEN
C        LEXIQUE VIDE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) ' LEXIQUE VIDE'
         ELSE
            WRITE(IMPRIM,*) ' EMPTY LEXICON'
         ENDIF
         RETURN
      ENDIF
C
C     LE NOM EST IL OCCUPE ?
 10   IF( N .GT. 0 ) THEN
C        OUI.TRANSFORMATION DES NBENNM ENTIERS EN CARACTERES
         MN = MNLX + M1LX * N
         CALL ENTNOM( NBENNM , MCN(MN) , KNOM )
C        IMPRESSION DU NOM
         IF( NBLGRC(NRERR) .GE. MXLGER ) THEN
             KERR(MXLGER) = ' ... '
             GOTO 20
         ENDIF
         NBLGRC(NRERR) = NBLGRC(NRERR) + 1
         KERR(NBLGRC(NRERR)) = KNOM(1:NBCANM)
C        PASSAGE AU SUIVANT
         N = MCN( MN + NBENNM )
         GOTO 10
C     ELSE
C        NON.RETOUR
      ENDIF
C
C     AFFICHAGE EFFECTIF
 20   CALL LERESU
      END
