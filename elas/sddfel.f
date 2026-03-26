      SUBROUTINE SDDFEL( KNOMST , LISTE1 , LLIST1 , LISTE2 , LLIST2 ,
     %                   KPB    , KPBMTC , KPBMCL , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LES CARACTERISTIQUES PHYSIQUES DE L'ELASTICITE
C -----    DE L'OBJET KNOMSD ET LES FORCES EXERCEES DANS LE CAS DE
C          RESOLUTION PAR SOUS-DOMAINES
C
C ENTREE :
C --------
C KNOMST : NOM DE L'OBJET A ELASTICITER
C
C SORTIES :
C ---------
C IERR   : 0 SI PAS D'ERREUR , NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  FEVRIER 1990
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*10      NMTYOB,KNOMTY
      CHARACTER*24      KNOMST,KNOMOB,KNOMSD,KNOMLI,KNOM
      CHARACTER*(*)     LISTE1(LLIST1)
      CHARACTER*(*)     LISTE2(LLIST2)
      CHARACTER*(*)     KPB,KPBMTC,KPBMCL
C
      IERR = 0
C
C     LE NUMERO DE L'OBJET KNOMST DANS LE LEXIQUE
      CALL LXNMNO( NTOBJE , KNOMST , NUOBST , MNLXST )
C     NOM DE L'OBJET "SOUS-DOMAINES"
      L = INDEX( KNOMST , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMST )
      ENDIF
C
C     CREATION D'UN OBJET SOUS-DOMAINE
C     ================================
C
      KNOMSD = KNOMST(1:L) // '_SD'
C     RECHERCHE DE L'OBJET KNOMSD DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE , KNOMSD , NTLXSD , MNLXSD )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXSD .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOMSD
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU DEFINITION DE L'OBJET KNOMSD
      CALL LXTSOU( NTLXSD , 'DEFINITION' , NTDFSD , MNDFSD )
C     CALCUL DU NOMBRE DE SOUS-DOMAINES ET DE LIGNES
      NBDOBJ = MCN(MNDFSD+WBDOBJ)
      MNSD   = MNDFSD + WTYOBJ
      NBSURF = 0
      NBLIGN = 0
      NBPOIN = 0
      DO 99 NO = 1 , NBDOBJ
         NUTYOB = MCN(MNSD)
         MNSD   = MNSD + 2
         IF (NUTYOB .EQ. 3 ) THEN
            NBSURF =  NBSURF + 1
         ELSE IF (NUTYOB .EQ. 2 ) THEN
            NBLIGN =  NBLIGN + 1
         ELSE IF (NUTYOB .EQ. 1 ) THEN
            NBPOIN =  NBPOIN + 1
         END IF
 99   CONTINUE
C
CCCC     LECTURE DES DONNEES INTERNES AU NIVEAU DE L'OBJET
CCC      CALL SDDEF1( KNOMST , LISTE1 , KPB , KPBMTC )
C
C     LECTURE DES DONNEES INTERNES
C     AU NIVEAU DE CHAQUE SOUS-DOMAINE
C     --------------------------------
      MNSD = MNDFSD + WTYOBJ
      DO 100 NO = 1 , NBDOBJ
         NUTYSD = MCN(MNSD)
         NUOBSD = MCN(MNSD+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNSD   = MNSD + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           TRAITEMENT DES SOUS-DOMAINES
C           LE NOM DU SOUS-DOMAINE
            CALL NMOBNU( KNOMTY , NUOBSD , KNOM )
C           RECHERCHE DE SON NOM DANS LE LEXIQUE DES OBJETS
            CALL LXLXOU( NTOBJE , KNOM , NTLXOB , MNLXOB )
            IF( NTLXOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOM
               CALL LEREUR
               GOTO 9900
            ENDIF
C           LE TABLEAU TOPOLOGIE EST RECHERCHE
            CALL LXTSOU( NTLXOB , 'TOPOLOGIE' , NTTOPO , MNTOPO )
            IF( NTTOPO .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = ' ERREUR : OBJET ' // KNOM //
     %                   ' SANS TOPOLOGIE'
               CALL LEREUR
               GOTO 9900
            ENDIF
C           BOUCLE SUR LES OBJETS INTERNES
            NBTYEL = MCN( MNTOPO + WBTYEL )
            NBOBIN = MCN( MNTOPO + WBOBIN )
            MNOB   = MNTOPO + WMTYEL + NBTYEL - 2
            DO 101 I=0,NBOBIN-1
C              LE TYPE DE L'OBJET
               MNOB   = MNOB + 2
               NUTYOB = MCN( MNOB )
               NUMOBJ = MCN( MNOB + 1 )
               CALL LXNLOU(NTMN(NUTYOB),NUMOBJ,NTOB,MNMNOB)
               KNOMTY = NMTYOB( NUTYOB )
               CALL NMOBNU( KNOMTY , NUMOBJ , KNOMOB )
               IF( NTOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOMOB
               CALL LEREUR
                   GOTO 9900
               ENDIF
               LLISTE = 4
               CALL SDDEF2( NUOBST , NUTYOB , NUMOBJ ,
     %                      LISTE1 , LLIST1 ,
     %                      KPB    , KPBMTC , 1 )
 101        CONTINUE
         END IF
 100  CONTINUE
C
C     LECTURE DES DONNEES AUX LIMITES
C     AU NIVEAU DE L'OBJET
C     --------------------------------
C     LES LIGNES DE L'OBJET
      MNSD = MNDFSD + WTYOBJ
      MNLIMI = 0
      NBLIMI = NBSURF + NBLIGN + NBPOIN
      MOLIMI = NBLIMI * 2
      CALL TNMCDC( 'ENTIER' , MOLIMI , MNLIMI )
      NULIMI = 0
      DO 200 NO = 1 , NBDOBJ
         NUTYOB = MCN(MNSD)
         NUMEOB = MCN(MNSD+1)
         KNOMTY = NMTYOB(NUTYOB)
         MNSD   = MNSD + 2
C        SEUL LE CAS LIGNE EST TRAITE
         IF (NUTYOB .EQ. 2 ) THEN
C           LE NOM DE LA LIGNE
            CALL NMOBNU( KNOMTY , NUMEOB , KNOMLI )
C           SON LEXIQUE
            CALL LXLXOU( NTLIGN , KNOMLI , NTLXOB , MNLXOB )
            IF( NTLXOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = ' ERREUR : LIGNE INCONNUE ' // KNOMLI
               GOTO 9900
            ENDIF
            MCN(MNLIMI+NULIMI)=NUTYOB
            MCN(MNLIMI+NULIMI+1)=NTLXOB
            NULIMI = NULIMI + 2
            CALL SDDEF3( KNOMLI , KNOMTY , LISTE2 , KPBMCL )
         ENDIF
 200  CONTINUE
      NBLIMI = NULIMI / 2
C
C     LECTURE DES DONNEES AUX LIMITES
C     AU NIVEAU DE CHAQUE SOUS-DOMAINE
C     --------------------------------
      MNSD = MNDFSD + WTYOBJ
      DO 300 NO = 1 , NBDOBJ
         NUTYSD = MCN(MNSD)
         NUOBSD = MCN(MNSD+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNSD   = MNSD + 2
         IF (NUTYSD .EQ. 5 ) THEN
C           TRAITEMENT DES SOUS-DOMAINES
C           LE NOM DU SOUS-DOMAINE
            CALL NMOBNU( KNOMTY , NUOBSD , KNOM )
C           SON LEXIQUE
            CALL LXLXOU( NTOBJE , KNOM , NTLXOB , MNLXOB )
            IF( NTLXOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOM
               GOTO 9900
            ENDIF
C           LE TABLEAU TOPOLOGIE EST RECHERCHE
            CALL LXTSOU( NTLXOB , 'TOPOLOGIE' , NTTOPO , MNTOPO )
            IF( NTTOPO .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = ' ERREUR : OBJET ' // KNOM //
     %                   ' SANS TOPOLOGIE'
               GOTO 9900
            ENDIF
C           BOUCLE SUR LES OBJETS AUX LIMITES
            NBOBCL = MCN( MNTOPO + WBOBCL )
            NBTYEL = MCN( MNTOPO + WBTYEL )
            NBOBIN = MCN( MNTOPO + WBOBIN )
            MN     = MNTOPO + WMTYEL + NBTYEL - 2 + MOTVAR(13)*NBOBIN
            DO 301 I=0,NBOBCL-1
C              LE TYPE DE L'OBJET
               MN   = MN + 2
               NUTYOB = MCN( MN )
C              LE CAS DES POINTS ET DES LIGNES
               IF (NUTYOB .GE. 3) GO TO 301
               NUMOBJ = MCN( MN + 1 )
               CALL LXNLOU( NTMN(NUTYOB) , NUMOBJ , NTOBRE , MNOB )
               KNOMTY = NMTYOB( NUTYOB )
               CALL NMOBNU( KNOMTY , NUMOBJ , KNOM )
               IF( NTOBRE .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOM
                  CALL LEREUR
                  GOTO 9900
               ENDIF
C              RECHERCHE DES OBJETS AUX LIMITES QUI LE CONTIENNENT
C              ** ATTENTION L'IDENTIFICATION EST FAITE SUR LES SOMMETS ***
               CALL SDDEF5( NBLIMI, MCN(MNLIMI), NUTYOB, NTOBRE, NTOBTR)
               IF (NTOBTR .GT. 0) THEN
C              L'OBJET EST RETROUVE : TRANSMISSION DES DONNEES AUX LIMITES
                  IF (NUTYOB .EQ. 1 ) THEN
                     LISTE0 = 2
                  ELSE
                     LISTE0 = 1
                  ENDIF
                  LLISTE = 2
                  CALL SDDEF4(NTOBTR,NTOBRE,LISTE2,LISTE0,LLIST2)
               END IF
 301        CONTINUE
         END IF
 300  CONTINUE
      CALL TNMCDS( 'ENTIER' , MOLIMI , MNLIMI )
C
C     ERREUR
 9900 IF( INTERA .LE. 1 ) CALL ARRET(100)
      RETURN
      END
