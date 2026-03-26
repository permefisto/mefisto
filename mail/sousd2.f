      SUBROUTINE SOUSD2( NUTYOB , NUOBJE ,
     S                   MNOBPR , NBOBPR , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CREER LA LISTE DES OBJETS PREMIERS D'UN OBJET DONNE
C -----
C
C ENTREES :
C ---------
C NUTYOB : LE NUMERO DU TYPE DE L'OBJET
C NUOBJE : LE NUMERO DE L'OBJET DANS LE LEXIQUE
C
C SORTIES :
C ---------
C MNOBPR : L'ADRESSE MCN DU TABLEAU DES OBJETS AUX LIMITES
C NBOBPR : LE NOMBRE DES OBJETS AUX LIMITES
C IERR   : 0 SI PAS D'ERREUR , > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1990
C23456---------------------------------------------------------------012
      include"./incl/a_objet__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      INTEGER          NUCOM(6)
      CHARACTER*24     KNOM
      CHARACTER*10     NOTOB,NMTYOB
      EQUIVALENCE     (RMCN(1),MCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
      IERR  = 0
      NTYCO = 0
      NBCOM = 0
C
C     LA PILE DU TYPE ET NUMERO DES OBJETS EST DECLAREE
      MOPILE = 64
      CALL TNMCDC( 'ENTIER' , MOPILE , MNPILE )
C
C     LE TABLEAU DU TYPE ET NUMERO DES OBJETS PREMIERS
      MOOBPR = 64
      CALL TNMCDC( 'ENTIER' , MOOBPR , MNOBPR )
      NBOBPR = 0
C
C     L'OBJET INITIAL EST EMPILE
      LHPILE = 2
      MCN( MNPILE )     = NUTYOB
      MCN( MNPILE + 1 ) = NUOBJE
C
C     TANT QUE LA PILE EST NON VIDE DECOMPOSER EN OBJETS PREMIERS
 10   IF( LHPILE .GT. 0 ) THEN
C        LE TYPE ET NUMERO DE L'OBJET
         MN     = MNPILE - 2 + LHPILE
C        LE TYPE DE L'OBJET
         NUTYOB = MCN( MN )
         ICAS = 1
         IF (NUTYOB.LE.0) ICAS = -1
         NUTYOB = ICAS * NUTYOB
         NOTOB=NMTYOB(NUTYOB)
C        SON NUMERO
         NUOBJE = MCN( MN + 1 )
C        LE NOM DE L'OBJET
         CALL NMOBNU( NOTOB , NUOBJE , KNOM )
C        CET OBJET EST DEPILE
         LHPILE = LHPILE - 2
C
C        L'OBJET EST IL UN OBJET PREMIER ?
         IF (ICAS.EQ.1) THEN
         IF( NUTYOB .GT. 0 .AND. NUTYOB .LT. 5 ) THEN
C           OUI : STOCKAGE DE CET OBJET PREMIER
            IF( NBOBPR + NBOBPR .GE. MOOBPR ) THEN
C              LE TABLEAU TROP PETIT . SA TAILLE EST AUGMENTEE
               CALL TNMCAU( 'ENTIER' , MOOBPR , MOOBPR+MOOBPR ,
     %                       MOOBPR , MNOBPR )
               MOOBPR = MOOBPR + MOOBPR
            ENDIF
            MN = MNOBPR + NBOBPR + NBOBPR
            MCN( MN     ) = NUTYOB
            MCN( MN + 1 ) = NUOBJE
            NBOBPR = NBOBPR + 1
         ELSE IF( NUTYOB .NE. 5 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYOB
            KERR(1)=' OBJET '//KNOM//' TYPE '//KERR(MXLGER)(1:4)
            KERR(2)=' NON ACCEPTE EN SOUS-DOMAINES'
            CALL LEREUR
            IERR = 1
            RETURN
         ENDIF
         ENDIF
C
C        LE LEXIQUE DE CE TYPE D'OBJETS
         NTLX = NTMN( NUTYOB )
C        LE LEXIQUE DE CET OBJET
         CALL LXNLOU( NTLX , NUOBJE , NTLXOB , NT )
         IF( NTLXOB .LE. 0 ) THEN
C           LE NOM DE L'OBJET
            CALL NMOBNU( NOTOB , NUOBJE , KNOM )
            NBLGRC(NRERR) = 1
            KERR(1) =  NOTOB // KNOM // ' SANS LEXIQUE'
            CALL LEREUR
            IERR = 1
            GOTO 10
         ENDIF
C
C        RECHERCHE DU TABLEAU DEFINITION DE L'OBJET
         CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C        S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
         IF( NTDFOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  NOTOB // KNOM // ' SANS DEFINITION'
            CALL LEREUR
            IERR = 1
            GOTO  10
         ENDIF
C        RECHERCHE DU TABLEAU XYZSOMMET DE L'OBJET
         CALL LXTSOU( NTLXOB , 'XYZSOMMET' , NTXYZ , MNXYZ )
C        S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
         IF( NTXYZ .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  NOTOB // KNOM // ' SANS TABLEAU XYZSOMMET'
            CALL LEREUR
            IERR = 1
            GOTO  10
         ENDIF
C
C        TRAITEMENT SUIVANT LE TYPE DE L'OBJET
C        -------------------------------------
C
         IF (NUTYOB.LE.1) THEN
            GO TO 10
         ELSE IF (NUTYOB.EQ.2) THEN
C           CAS D'UNE LIGNE
            NUTYLI=MCN(MNDFOB+WUTYLI)
C           TEST SUR LE TYPE DE LIGNE
            IF (NUTYLI .EQ. 1) THEN
               NBCOM  = 2
               NTYCO  = 1
               NBPOI  = MCN(MNDFOB+WBPOLI)
               NUCOM(1) = MCN(MNDFOB+WUPOLI)
               NUCOM(2) = MCN(MNDFOB+WUPOLI+NBPOI-1)
               GO TO 11
            ELSE IF (NUTYLI .EQ. 2) THEN
               NBCOM  = 2
               NTYCO  = 1
               NUCOM(1) = MCN(MNDFOB+WUPTIN)
               NUCOM(2) = MCN(MNDFOB+WUPTFI)
               GO TO 11
            ELSE IF (NUTYLI .EQ. 3) THEN
               NBCOM  = 2
               NTYCO  = 1
               NUCOM(1) = MCN(MNDFOB+WUPTIN)
               NUCOM(2) = MCN(MNDFOB+WUPTFI)
               GO TO 11
            ELSE IF (NUTYLI .GE. 11 .AND. NUTYLI.LE. 14) THEN
               NBCOM  = 2
               NTYCO  = 1
               NBPCBL = MCN(MNDFOB+WBPCBL)
               NUCOM(1) = MCN(MNDFOB+WUPCBL)
               NUCOM(2) = MCN(MNDFOB+WUPCBL+NBPCBL-1)
               GO TO 11
            ELSE IF (NUTYLI .GE. 40 .AND. NUTYLI.LE. 41) THEN
               NBCOM = 1
               NTYCO = -22
               NUCOM(1) = MCN(MNDFOB+WULIIN)
               GO TO 11
            ELSE
C           CAS NON TRAITE : ON N'EMPILE PAS
               NBLGRC(NRERR) = 2
               WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYLI
               KERR(1)=' LIGNE '//KNOM//' TYPE '//KERR(MXLGER)(1:4)
               KERR(2)=' CAS EN COURS DE PROGRAMMATION'
CCC               CALL LEREUR
CCC               WRITE(IMPRIM,1000) KNOM,NUTYLI
CCC 1000          FORMAT(/,' ATTENTION : LIGNE ',A24,/
CCC     %         ' TYPE ',I4,' CAS EN COURS DE PROGRAMMATION')
               GO TO 10
CCC               RETURN
            ENDIF
         ELSE IF (NUTYOB.EQ.3) THEN
C           CAS D'UNE SURFACE
            NUTYSU=MCN(MNDFOB+WUTYSU)
            IF (NUTYSU.LE.2) THEN
               NBCOM = 4
               NTYCO = 2
               DO 13 L=1,4
                  NUCOM(L)=MCN(MNDFOB+WU4COT+L-1)
 13            CONTINUE
               GO TO 11
            ELSE IF (NUTYSU.EQ.6) THEN
               NBCOM = 3
               NTYCO = 2
               DO 14 L=1,NBCOM
                  NUCOM(L) = MCN(MNDFOB+WU3COT+L-1)
 14            CONTINUE
               GO TO 11
            ELSE IF (NUTYSU.EQ.9) THEN
               NBCOM = MCN(MNDFOB+WBLFTR)
               NTYCO = 2
               DO 15 L=1,NBCOM
                  NUCOM(L) = MCN(MNDFOB+WULFTR+L-1)
 15            CONTINUE
               GO TO 11
C           cf ~/td/f/xxtstd.ftn
            ELSE IF (NUTYSU .GE. 29 .AND. NUTYLI.LE. 31) THEN
               NBCOM = 1
               NTYCO = -3
               NUCOM(1) = MCN(MNDFOB+WUSUQU)
               GO TO 11
            ELSE
C           CAS NON TRAITE : ON N'EMPILE PAS
               NBLGRC(NRERR) = 2
               WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYSU
               KERR(1)=' SURFACE '//KNOM//' TYPE '//KERR(MXLGER)(1:4)
               KERR(2) =  ' CAS EN COURS DE PROGRAMMATION'
CCC               CALL LEREUR
CCC               WRITE(IMPRIM,2000) KNOM,NUTYSU
CCC 2000          FORMAT(/,' ATTENTION : SURFACE ',A24,/
CCC     %         ' TYPE ',I4,' CAS EN COURS DE PROGRAMMATION')
               GO TO 10
CCC               RETURN
            ENDIF
         ELSE IF (NUTYOB.EQ.4) THEN
C           CAS D'UN VOLUME
            NUTYVO=MCN(MNDFOB+WUTYVO)
C           TEST SUR LE TYPE DE VOLUME
            IF (NUTYVO .LE. 2) THEN
               NBCOM = 6
               NTYCO = 3
               MNCOM = WU6FAC
               GO TO 12
            ELSE IF (NUTYVO .EQ. 3) THEN
               NBCOM = 4
               NTYCO = 3
               MNCOM = WU4FAC
               GO TO 12
            ELSE IF (NUTYVO .EQ. 4) THEN
               NBCOM = 5
               NTYCO = 3
               MNCOM = WU5FAC
               GO TO 12
            ELSE IF (NUTYVO .GE. 30 .AND. NUTYVO .LE. 31) THEN
               NBCOM = 1
               NTYCO = -4
               NUCOM(1) = MCN(MNDFOB+WUVOIN)
               GO TO 11
            ELSE IF (NUTYVO .EQ. 41 ) THEN
               NBCOM = 1
               NTYCO = -4
               NUCOM(1) = MCN(MNDFOB+WUVOIN)
               GO TO 11
            ELSE
C           CAS NON TRAITE : ON N'EMPILE PAS
               NBLGRC(NRERR) = 2
               WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYVO
               KERR(1)=' VOLUME '//KNOM//' TYPE '//KERR(MXLGER)(1:4)
               KERR(2) =  ' CAS EN COURS DE PROGRAMMATION'
CCC               CALL LEREUR
CCC               WRITE(IMPRIM,3000) KNOM,NUTYVO
CCC 3000          FORMAT(/,' ATTENTION : VOLUME ',A24,/
CCC     %         ' TYPE ',I4,' CAS EN COURS DE PROGRAMMATION')
               GO TO 10
C               RETURN
            END IF
 12         DO 16 L=1,NBCOM
               NUCOM(L)=MCN(MNDFOB+MNCOM+L-1)
 16         CONTINUE
C           CAS D'UN OBJET
         ELSE IF ( NUTYOB .EQ. 5) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYOB
            KERR(1)=' OBJET '//KNOM//' TYPE '//KERR(MXLGER)(1:4)
            KERR(2)=' NON ACCEPTE EN SOUS-DOMAINES'
            CALL LEREUR
CCC            WRITE(IMPRIM,4000) KNOM,NUTYOB
CCC 4000       FORMAT(/,' ERREUR SOUSD2 : OBJET ',A24,/
CCC     %            ' TYPE ',I4, ' CAS NON PROGRAMME')
            RETURN
         ENDIF
C
C        LES NBCOM COMPOSANTES DE CET OBJET SONT EMPILEES
 11      MN = NBCOM * 2
         IF( LHPILE + MN .GT. MOPILE ) THEN
C           LE TABLEAU TROP PETIT . SA TAILLE EST AUGMENTEE
            CALL TNMCAU( 'ENTIER' , MOPILE , MOPILE+MOPILE ,
     %                     MOPILE , MNPILE )
            MOPILE = MOPILE + MOPILE
         ENDIF
         DO 17 NB = 1 , NBCOM
            MCN(MNPILE+LHPILE)   =  NTYCO
            MCN(MNPILE+LHPILE+1) =  NUCOM(NB)
            LHPILE = LHPILE + 2
 17      CONTINUE
         GOTO 10
      ENDIF
C
C     DESTRUCTION DE LA PILE DEVENUE INUTILE
      CALL TNMCDS( 'ENTIER' , MOPILE , MNPILE )
C
C     REAGENCEMENT DU TABLEAU
C     -----------------------
      MN = MNOBPR
      DO 20 NB = 1 , NBOBPR
         NUTY = MCN( MN )
         NUOB = MCN( MN + 1 )
         MN = MN + 2
         IF (NUTY .EQ. 0) GO TO 20
         MNL = MN
         DO 21 NBL = NB + 1 , NBOBPR
            NUTYL = MCN( MNL )
            NUOBL = MCN( MNL + 1 )
            IF ( NUTYL .EQ. NUTY .AND.
     %          NUOBL .EQ. NUOB ) THEN
                MCN( MNL ) = 0
                MCN( MNL+ 1 ) = 0
            END IF
            MNL = MNL + 2
 21      CONTINUE
 20   CONTINUE
      MN   = MNOBPR
      MNOB = MNOBPR
      NBOB = 0
      DO 22 NB = 1 , NBOBPR
         NUTY = MCN( MN )
         NUOB = MCN( MN + 1 )
         MN = MN + 2
         IF (NUTY .NE. 0) THEN
            MCN( MNOB   ) = NUTY
            MCN( MNOB+1 ) = NUOB
            NBOB = NBOB + 1
            MNOB = MNOB + 2
         END IF
 22   CONTINUE
      NBOBPR = NBOB
      CALL TNMCRA('ENTIER' , MOOBPR , NBOBPR*2 , MNOBPR )
C
      RETURN
      END
