      SUBROUTINE LIJUMO( NBMOT, MOT1, MOT2, NL , NC , NLD , NCD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOURNER (NLD,NCD)  DEBUT DANS KLG
C -----     DU 1-ER CARACTERE DE MOT1 FRACAIS OU MOT2 ANGLAIS
C           LE PLUS PROCHE EN DEMANDANT A L'UTILISATEUR
C           AUTANT DE LIGNES QUE NECESSAIRE POUR TROUVER MOT
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KLG DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DU MOT
C           NLD= 0 SI LE MOT N'EST PAS RETROUVE
C           NLD=-1 SI LE CARACTERE D'ABANDON @ A ETE FRAPPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     MOT1, MOT2
C
C     SAUVEGARDE DE NL NC
      NLD = NL
      NCD = NC
C
C     RECHERCHE DU MOT SUR LA PREMIERE LIGNE
 10   I = INDEX( KLG(NLD)(NCD:NBCALI) , MOT1 )
      IF( NBMOT .GT. 1 ) THEN
         I2 = INDEX( KLG(NLD)(NCD:NBCALI), MOT2 )
         I  = MAX( I, I2 )
      ENDIF
      IF( I .LE. 0 ) THEN
C        AUCUN MOT N'EST RETROUVE   PASSAGE A LA LIGNE SUIVANTE
         IF( NLD .LE. LHKLG ) THEN
C           LECTURE D'UNE LIGNE DEMANDEE A L'UTILISATEUR
            IF( INTERA .GE. 3 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'En ATTENTE de ' // MOT1
               ELSE
                  KERR(1) = 'WAITING for ' // MOT2
               ENDIF
               CALL LERESU
            ENDIF
            CALL LIRLIG( I )
            IF( I .EQ. -1 ) THEN
               NLD = -1
               RETURN
            ELSE IF( I .NE.  0 ) THEN
C              PROBLEME DANS LA LECTURE
               NLD = 0
               RETURN
            ELSE
               NLD = LHKLG - 1
            ENDIF
         ENDIF
         NLD = NLD + 1
         NCD = 1
         GOTO 10
      ENDIF
C
C     ICI LE MOT A ETE RETROUVE
      NCD = NCD - 1 + I
      END
