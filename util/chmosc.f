      SUBROUTINE CHMOSC( MOT , NLMAX , NCMAX , NL , NC ,
     %                   NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE DE MOT SOUS LA CONDITION
C       DE NE PAS DEPASSER LE CARACTERE NLMAX,NCMAX
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C NLMAX   : NUMERO DE LA DERNIERE LIGNE DE RECHERCHE
C NCMAX   : NUMERO DE LA DERNIERE COLONNE DE NLMAX DE RECHERCHE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DU MOT
C           NLD = 0 SI LE MOT N'EST PAS RETROUVE
C NLF,NCF : POSITION DU DERNIER CARACTERE DU MOT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
      CHARACTER*(*)  MOT
C     LE 1-ER CARACTERE APRES
      NLD = NL
      NCD = NC
      CALL CARAPR( NLD , NCD )
      NLL = NLD
      IF( NLD .EQ. NLMAX ) GOTO 20
C
C     RECHERCHE SUR LA PREMIERE LIGNE
      I   = INDEX( KTD(NLL)(NCD:NCKTD) , MOT )
      IF( I .GT. 0 ) GOTO 100
C
C     RECHERCHE SUR LES LIGNES INTERMEDIAIRES
      NCD = 1
      DO 10 NLL = NLD+1 , NLMAX-1
         I = INDEX( KTD(NLL) , MOT )
         IF( I .GT. 0 ) GOTO 100
 10   CONTINUE
C
C     RECHERCHE SUR LA DERNIERE LIGNE
      NLL = NLMAX
 20   I   = INDEX( KTD(NLMAX)(NCD:NCMAX) , MOT )
      IF( I .GT. 0 ) GOTO 100
C
C     MOT NON RETROUVE
      NLD = 0
      RETURN
C
C     MOT RETROUVE LOCALISATION DU DERNIER CARACTERE DE MOT
 100  NLD = NLL
      NLF = NLL
      NCD = NCD + I - 1
      NCF = NCD - 1
      DO 110 I=1,LEN(MOT)
         NCF = NCF + 1
         IF( MOT(I:I) .NE. KTD(NLD)(NCF:NCF) ) THEN
            NCF = NCF - 1
            RETURN
         ENDIF
 110  CONTINUE
      END
