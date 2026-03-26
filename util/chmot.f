      SUBROUTINE CHMOT( MOT , NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE DE MOT
C       A L'EXCEPTION DES MOTS DANS UNE CHAINE
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DU MOT
C NLF,NCF : POSITION DU DERNIER CARACTERE DU MOT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     MOT
C
C     LE 1-ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
C     LE NOMBRE D'APOSTROPHES
      NBAPOS = 0
C
 5    CALL CAR1NB( NLD , NCD )
      NLDD = NLD
      NCDD = NCD
C
C     RECHERCHE DU MOT
 10   NNN = INDEX( KTD(NLD)(NCD:NCKTD) , MOT )
      IF( NNN .LE. 0 ) THEN
         NLD = NLD + 1
         NCD = 1
         IF( NLD .GT. LHKTD ) THEN
CCC            NBLGRC(NRERR) = 1
CCC            KERR(1) = 'ERREUR CHMOT:DEBORDEMENT DU TABLEAU KTD'
CCC            CALL LEREUR
            NLD = 0
            RETURN
         ENDIF
         GOTO 10
      ELSE
         NCD = NCD + NNN - 1
      ENDIF
C
C     LA LONGUEUR DU MOT
      NLF = NLD
      NCF = NCD + LEN( MOT ) - 1
C
C     LE RESULTAT EST CORRECT SI LE MOT N'EST PAS DANS UN COMMENTAIRE
C     LE NOMBRE DE ' EST COMPTE ENTRE LE CARACTERE DE DEPART ET LE
C     PREMIER CARACTERE DU MOT
      NC1    = NCDD
      NC2    = NCKTD
      DO 20 N=NLDD,NLD
         IF( N .EQ. NLD ) NC2 = NCD
 15      NNN = INDEX( KTD(N)(NC1:NC2) , '''' )
         IF( NNN .GT. 0 ) THEN
            NBAPOS = NBAPOS + 1
            NC1    = NC1 + NNN
            GOTO 15
         ENDIF
         NC1 = 1
 20   CONTINUE
C
C     SI LE NOMBRE D'APOSTROPHES EST IMPAIR LE MOT EST DANS UNE CHAINE
      IF( MOD(NBAPOS,2) .EQ. 1 ) THEN
C        LA CHAINE CONTIENT LE MOT
         NLD = NLF
         NCD = NCF
         GOTO 5
      ENDIF
      END
