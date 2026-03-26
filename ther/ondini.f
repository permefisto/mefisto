      SUBROUTINE ONDINI( KNOMOB, NTLXOB, MOREE2, NTDL, TPSINI, IEONIN,
     %                   NBTYEL, MNNPEF, NDPGST,
     %                   MNXYZN, NUMIOB, MNDOEL,
     %                   NTVECT, MNVECT, NBVECT, MNUG0, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE VECTEUR DE L'ONDE(NBNOMA,2) COMPLEXE INITIALE
C -----    DE NLSE AUX NBNOMA NOEUDS DU MAILLAGE DE L'OBJET
C          L'ONDE COMPLEXE INITIALE A ETE DONNEE EVENTUELLEMENT SUR
C          L'OBJET soit deja CALCULEE, soit CONSTANTE, soit FONCTION
C          LES VOLUMES  soit Constante, soit Fonction
C          LES SURFACES soit Constante, soit Fonction
C          LES LIGNES   soit Constante, soit Fonction
C          LES POINTS   soit Constante, soit Fonction
C
C          LES POSSIBILITES selon le TMS a___ondeinit
C          3 : 'VECTEUR PROPRE EXISTANT'
C          2 : 'Onde INITIALE EXISTANTE'
C          1 : 'Onde INITIALE CONSTANTE'
C         -1 : 'FONCTION Onde INITIALE'
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C NTLXOB : NUMERO DU TMS DU LEXIQUE DE L'OBJET KNOMOB
C MOREE2 : NOMBRE DE MOTS   D'UNE VARIABLE REELLE DOUBLE PRECISION
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE SUR LE MAILLAGE
C          ICI ONDE COMPLEXE NLSE  NTDL=2 FOIS LE NOMBRE DE NOEUDS
C          POUR STOCKER LA PARTIE REELLE ET IMAGINAIRE
C
C TPSINI : TEMPS INITIAL DU CALCUL (PEUT ETRE MODIFIE!)
C IEONIN :=0 SI AUCUN PLSV DE L'OBJET SUPPORTE UN TMS ONDEINIT
C         >0 SINON
C
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES PLSV
C MNDOEL : ADRESSE MCN DES DONNEES THERMIQUES DE L'OBJET
C
C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"ONDENLSE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"ONDENLSE DE L'OBJET
C NBVECT : NUMERO DU DERNIER VECTEUR TEMPERATURE STOCKE = VECTEUR A TPSINI
C MNUG0  : ADRESSE MCN DE UG0(NBNOMA,2) ONDE NLSE A L'INSTANT INITIAL
C IERR   : 0 SI PAS D'ERREUR RENCONTREE, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A&M University at QATAR    FEVRIER 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___ondeinit.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/donthe.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      CHARACTER*(*)     KNOMOB
C
      DOUBLE PRECISION  XYZNOEUD(3), ONDEC0(2)
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6),NOOBLA(12),NOOBPS(8)
      INTEGER           NUMIOB(4),MNDOEL(4)
C
C     INITIALISATION DU VECTEUR UG0(NBNOMA,2) DE L'ONDE NLSE INITIALE
C     AUX NOEUDS DU MAILLAGE DE L'OBJET ET A L'INSTANT TEMPS INITIAL
C     ===============================================================
C     LE VECTEUR GLOBAL DE L'ONDE NLSE INITIALE EST DECLARE
C     NTDL=2 NOMBRE DE NOEUDS
      NBNOMA = NTDL / 2
      IF( MNUG0 .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', NTDL, MNUG0 )
C        LE VECTEUR GLOBAL DE L'ONDE NLSE INITIALE EST MIS A ZERO
         CALL AZEROD( NTDL, MCN(MNUG0) )
      ENDIF
C
C     -------------------------------------------------------
C     TABLEAU DE L'ONDE COMPLEXE NLSE INITIALE
C     1-ERE COMPOSANTE = PARTIE REELLE     DE L'ONDE COMPLEXE
C     2-EME COMPOSANTE = PARTIE IMAGINAIRE DE L'ONDE COMPLEXE
C     -------------------------------------------------------
      NBVECT = 0
      CALL LXTSOU( NTLXOB, 'ONDEINIT', NTONIN, MNONIN )
      IF( MNONIN .LE. 0 ) GOTO 80
C
C     ICI LE TABLEAU DE L'ONDE NLSE INITIALE DE L'OBJET EXISTE
C     LE CODE DE CALCUL DE L'ONDE NLSE INITIALE
      LTOND0 = MCN( MNONIN + WTOND0 )
      IF( LTOND0 .EQ. 1 ) THEN
C
C        ONDE NLSE(NBNOMA,2) INITIALE CONSTANTE
C        --------------------------------------
         MN = ( MNUG0 - 1 ) / MOREE2
         DO K = 1, 2
            DO I=1,NBNOMA
               DMCN( MN  + I ) = RMCN( MNONIN + WAOND0 -1 + K )
            ENDDO
            MN = MN + NBNOMA
         ENDDO
         GOTO 80
C
      ELSE IF( LTOND0 .EQ. -1 ) THEN
C
C        L'ONDE NLSE(NBNOMA,2) EST OBTENUE PAR FONCTION
C        ----------------------------------------------
         CALL LXNMNO( NTOBJE, KNOMOB, NOOB, I )
         MN = (MNUG0-1) / MOREE2
         MM = MNXYZN + WYZNOE - 3
         DO 20 I=1,NBNOMA
            MM = MM + 3
            XYZNOEUD(1) = RMCN(MM  )
            XYZNOEUD(2) = RMCN(MM+1)
            XYZNOEUD(3) = RMCN(MM+2)
            CALL REONIN( 5, NOOB, XYZNOEUD, MNONIN,  ONDEC0 )
C           PARTIE REELLE
            DMCN( MN + I ) = ONDEC0(1)
C           PARTIE IMAGINAIRE
            DMCN( MN + I + NBNOMA ) = ONDEC0(2)
 20      CONTINUE
         GOTO 80
C
      ELSE IF( LTOND0 .EQ. 2 ) THEN
C
C        L'ONDE NLSE INITIALE EST STOCKEE DANS L'ACTUEL (LTOND0=2)
C        ~>OBJET>NomObjet>'VECTEUR"ONDENLSE'
C        ---------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"ONDENLSE', NTVECT, MNVECT )
         IF( MNVECT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PAS DE VECTEUR"ONDENLSE INITIAL'
            ELSE
               KERR(1) = 'NO VECTOR INITIAL ONDENLSE'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
         IF( NTDL .NE. MCN(MNVECT+WBCOVE) ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:8),'(I8)') NTDL
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOMBRE DE COMPOSANTES DU VECTEUR"ONDENLSE'
               KERR(2) = 'DIFFERENT DE ' // KERR(5)(1:8) // ' DL'
            ELSE
               KERR(1) = 'NUMBER of COMPONENTS of VECTOR"ONDENLSE'
               KERR(2) = 'NOT EQUAL TO ' // KERR(5)(1:8) // ' DoF'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
C        RECHERCHE DU TEMPS INITIAL TPSINI COMME LE PLUS PROCHE TEMPS
C        DE TPSINI PARMI LES TEMPS STOCKES DU CALCUL PRECEDENT
C        ------------------------------------------------------------
C        NOMBRE DE TEMPS STOCKES
         N      = MCN( MNVECT + WBCPIN )
         MNTIME = MNVECT + WECTEU + MOREE2*NTDL*MCN(MNVECT+WBVECT) - 1
C        NBVECT LE NUMERO DU VECTEUR STOCKE AU TEMPS LE PLUS PROCHE DE TPSINI
         NBVECT = 1
         R      = RINFO( 'GRAND' )
         DO 30 I=1,N
            RR = ABS( RMCN(MNTIME+I) - TPSINI )
            IF( RR .LT. R ) THEN
               NBVECT = I
               R      = RR
            ENDIF
 30      CONTINUE
C
         MNV = MNVECT + WECTEU + MOREE2 * NTDL * (NBVECT-1)
C        TRANSFERT DE L'ONDE NLSE INITIALE NBVECT DANS LE TABLEAU MNUG0
         CALL TRTATD( MCN(MNV), MCN(MNUG0), NTDL )
C
C        AFFICHAGE DU NOUVEAU TEMPS INITIAL RETENU
         TPSINI   = RMCN(MNTIME+NBVECT)
         TEMPS    = TPSINI
         TEMPSINI = TPSINI
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'ONDE INITIALE au TEMPS INITIAL=',TEMPS
         ELSE
            WRITE(IMPRIM,*) 'INITIAL WAVE at INITIAL TIME=',TEMPS
         ENDIF
C        INITIALISATION TERMINEE
         GOTO 9999
C
      ELSE IF( LTOND0 .EQ. 3 ) THEN
C
C        L'ONDE NLSE INITIALE EST STOCKEE DANS L'ACTUEL (LTOND0=3)
C        ~>OBJET>NomObjet>'VECTEUR"ONDENLSE'
C        ---------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"VALEURPROPRE', NTVECT, MNVECT )
         IF( MNVECT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PAS DE VECTEUR"VALEURPROPRE INITIAL'
            ELSE
               KERR(1) = 'NO VECTOR INITIAL EIGENVECTOR'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
         NBCOVE = MCN(MNVECT+WBCOVE)
         IF( NBNOMA .NE. NBCOVE ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:8),'(I8)') NBNOMA
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOMBRE DE COMPOSANTES DU VECTEUR"VALEURPROPRE'
               KERR(2) = 'DIFFERENT DE ' // KERR(5)(1:8) // ' NOEUDS'
            ELSE
               KERR(1) = 'NUMBER of COMPONENTS of VECTOR"VALEURPROPRE'
               KERR(2) = 'NOT EQUAL TO ' // KERR(5)(1:8) // ' NODES'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
C        RECHERCHE DU VECTEUR PROPRE DE NUMERO NUVPON
         NUVPON = MCN(MNONIN+WUVPON)
         IF( NUVPON .LE. 0 .OR. NUVPON .GT. MCN(MNVECT+WBVECT) ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:8),'(I8)') NUVPON
            WRITE(KERR(6)(1:8),'(I8)') MCN(MNVECT+WBVECT)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = KERR(5)(1:8) //
     %                 ' NUMERO INCORRECT DE VECTEUR PROPRE'
               KERR(2) = 'VALEUR CORRECTE ENTRE 0 et ' // KERR(6)(1:8)
            ELSE
               KERR(1) = KERR(5)(1:8) //
     %                 ' INCORRECT NUMBER of EIGENVECTOR'
               KERR(2) = 'CORRECT VALUE BETWEEN 0 and ' // KERR(6)(1:8)
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
C        NUVPON LE NUMERO DU VECTEUR PROPRE UTILISE COMME ONDE INITIALE
C        POUR SA PARTIE REELLE ET ZERO POUR LA PARTIE IMAGINAIRE
         NBCOVE = MCN(MNVECT+WBCOVE)
         MNV    = MNVECT + WECTEU + MOREE2 * NBCOVE * (NUVPON-1)
C
C        TRANSFERT DU VECTEUR PROPRE NUVPON COMME PARTIE REELLE DE L'ONDE
C        INITIALE DANS LE TABLEAU MNUG0
         CALL TRTATD( MCN(MNV), MCN(MNUG0), NBNOMA )
C
C        MISE A ZERO DE LA PARTIE IMAGINAIRE DE UG0
         CALL AZEROD( NBNOMA, MCN(MNUG0+MOREE2*NBNOMA) )
C
C        LA VALEUR PROPRE NUVPON
         MNTIME = MNVECT + WECTEU
     %          + MOREE2 * NBCOVE * MCN(MNVECT+WBVECT) - 1
         WRITE(IMPRIM,*)
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'VALEUR PROPRE',NUVPON,' =',RMCN(MNTIME+NUVPON)
         ELSE
            WRITE(IMPRIM,*) 'EIGENVALUE',NUVPON,' =',RMCN(MNTIME+NUVPON)
         ENDIF
C
C        INITIALISATION TERMINEE
C        UG0 N'EST PAS DANS LE TMS 'VECTEUR"ONDENLSE'
         NBVECT = 0
         GOTO 9999
C
      ENDIF
C
C     --------------------------------------------------------------
C     TRAITEMENT DE L'ONDE NLSE INITIALE PARTICULARISEE SUR DES PLSV
C     --------------------------------------------------------------
 80   IF( IEONIN .LE. 0 ) GOTO 9999
C     AU MOINS UN PLSV SUPPORTE UN TMS 'ONDEINIT'
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      MM0 = MNXYZN + WYZNOE - 3
      MN0 = (MNUG0-1) / MOREE2
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"No TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 95 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, SURFACES ET VOLUMES
            CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                  'ONDEINIT', MXDOTH, LPTEIN,
     %                   NOTYOB )
C
C           LE RECENSEMENT DE L'ONDE NLSE INITIALE
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO 90 J=1,NBNOE
C
C              L'ADRESSE MCN DU TABLEAU ONDEINIT
               MN1 = NOTYOB(3,J)
C
C              EXISTE-T-IL UNE "ONDEINIT" EN CE NOEUD ?
               IF( MN1 .GT. 0 ) THEN
C
C                 CALCUL DE L'ONDE NLSE INITIALE EN CE NOEUD
C                 LE TYPE OBJET DU NOEUD J DE L'EF
                  NTYOB = NOTYOB(1,J)
                  NOOB  = NOTYOB(2,J)
C                 LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
                  NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
                  MM    = MM0 + 3 * NONOE
                  XYZNOEUD(1) = RMCN(MM  )
                  XYZNOEUD(2) = RMCN(MM+1)
                  XYZNOEUD(3) = RMCN(MM+2)
                  CALL REONIN( NTYOB, NOOB, XYZNOEUD, MN1,  ONDEC0 )
C                 RANGEMENT PAR COMPOSANTES 1:PR puis 2:PI
C                 PARTIE REELLE
                  DMCN( MN0 + NONOE ) = ONDEC0(1)
C                 PARTIE IMAGINAIRE
                  DMCN( MN0 + NONOE + NBNOMA ) = ONDEC0(2)
C
               ENDIF
C
 90         CONTINUE
 95      CONTINUE
100   CONTINUE
C
9999  RETURN
ccc
ccc      print *,'TABLEAU de l'onde NLSE INITIALE='
ccc      mnn=MNXYZN+wyznoe
ccc      MN = (MNUG0-1) / MOREE2
ccc      do 200 k=1,ntdl
ccc         print *,'XYZ=',(RMCN(mnn+kk),kk=0,2),
ccc     %           ' UG01(',k,')=',DMCN(MN+k)
ccc         mnn = mnn + 3
ccc 200  continue
ccc
      END
