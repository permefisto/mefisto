      SUBROUTINE TRTOBJ( NMTOBJ, UNPARUN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER TOUS LES PLSVO ACTUELS DE TYPE NMTOBJ
C -----
C
C ENTREE :
C --------
C NMTOBJ : NOM DE CE TYPE D'OBJET
C          'POINT' 'LIGNE' 'SURFACE' 'VOLUME' 'OBJET'
C UNPARUN: 1 SI TRACE DEMANDE UN PAR UN
C          0 SI TRACE DE TOUS EN UNE SEULE FOIS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1988
C....................................................................012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      include"./incl/xyzext.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*24      NOMOBJ
      CHARACTER*(*)     NMTOBJ
      INTEGER           UNPARUN
C
C     NOM DU PLSV TRAITE EN DERNIER
      CHARACTER*26      KNOMA
      COMMON / KNMDSV / KNOMA
C
C     LE NUMERO DE TYPE DE CE TYPE D'OBJET
      NUTYOB = NTYOBJ( NMTOBJ )
C
C     RECHERCHE DE L'ADRESSE MCN DU LEXIQUE DE CE TYPE D'OBJETS
      CALL TAMSOU( NTMN(NUTYOB) , MNOBJT )
C
C     LE DEBUT DU CHAINAGE DES OBJETS OCCUPES DANS LE LEXIQUE
      NUOBJT = MCN( MNOBJT + 5 )
C
C     NOMBRE D'ENTIERS POUR STOCKER UN NOM DE OBJET
      NBENNM = MCN( MNOBJT + 2 )
C
cccC     TRACE DES AXES INITIAL
ccc      CALL TRAXES
C
C     LA BOUCLE SUR LES OBJETS OCCUPES
C     ================================
 10   IF( NUOBJT .GT. 0 ) THEN
C
C        ADRESSE MCN DU DEBUT DU OBJET DANS LE LEXIQUE
         MNOBJ = MNOBJT + MCN(MNOBJT) * NUOBJT
C
C        LE LEXIQUE DE CE OBJET EXISTE-T-IL ?
         NTOBJ = MCN( MNOBJ + NBENNM + 2 )
C
         IF( NTOBJ .GT. 0 ) THEN
C
C           CET OBJET EXISTE : RECHERCHE DE SON NOM
            CALL ENTNOM( NBENNM , MCN( MNOBJ ) , NOMOBJ )
C
C           TRACE DE CE PLSVO DE NOM NOMOBJ
            NBC = INDEX( NOMOBJ , ' ' ) - 1
            IF( NBC .LE. 0 ) NBC = LEN( NOMOBJ )
C
cccC           TRACE DU NOM DU PLSVO
ccc            IF( LORBITE .EQ. 0 ) THEN
ccc               CALL RECTEF( NRHIST )
ccc               NBLGRC(NRHIST) = 1
ccc               NB = INDEX( NMTOBJ , ' ' ) - 1
ccc               IF( NB .LE. 0 ) NB = LEN( NMTOBJ )
ccc               KHIST(1) = NMTOBJ(1:NB) // ': ' // NOMOBJ(1:NBC)
ccc               CALL LHISTO
ccc            ENDIF
C
C           NOM DU PLSV TRAITE EN DERNIER
            KNOMA = NOMOBJ
C
            IF( NMTOBJ(1:1) .EQ. 'O' ) THEN
C              OBJET
               CALL T1OBJE( NOMOBJ(1:NBC) )
            ELSE
C              POINT LIGNE SURFACE VOLUME
               CALL T1MOBJ( NMTOBJ, NOMOBJ(1:NBC), NUOBJT )
            ENDIF
C
cccC           RE-TRACE DU NOM DU PLSVO
ccc            IF( LORBITE .EQ. 0 ) THEN
ccc               NBLGRC(NRHIST) = 1
ccc               KHIST(1) = NMTOBJ(1:NB) // ': ' // NOMOBJ(1:NBC)
ccc               CALL LHISTO
cccccC              UN CLIC SOURIS OU FRAPPE CLAVIER POUR CONTINUER
ccccc               IF( UNPARUN .GT. 0 ) THEN 1dwall => 4 clics pour 4 materiaux
ccccc                  CALL CLICSO
ccccc               ENDIF
ccc            ENDIF
C
            IF( UNPARUN .GT. 0 ) THEN
               NB = INDEX( NMTOBJ , ' ' ) - 1
               IF( NB .LE. 0 ) NB = LEN( NMTOBJ )
               IF( LANGAG .EQ. 0 ) THEN
                  print *,'Trace ',NMTOBJ(1:NB),' : ',NOMOBJ(1:NBC)
               ELSE
                  print *,'Drawing  ',NMTOBJ(1:NB),' : ',NOMOBJ(1:NBC)
               ENDIF
            ENDIF
C
         ENDIF
C
C        PASSAGE A L' OBJET SUIVANT
         NUOBJT = MCN( MNOBJ + NBENNM )
         GOTO 10
C
      ENDIF
C
      RETURN
      END
