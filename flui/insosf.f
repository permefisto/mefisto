      SUBROUTINE INSOSF( MODESF, INTERP, NTLXOB,   MNDFOB,
     %                   NCAS0,  NCAS1,  SOLUTION, NMSOLU,
     %                   TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER L'INTEGRALE DES VECTEURS SOLUTION NCAS0 A NCAS1
C -----    (TEMPERATURE ou PRESSION ou ...) SUR LES SURFACES 3D
C           DE L'OBJET 3D
C
C ENTREES:
C --------
C MODESF : <=0 DEMANDE DU CALCUL SUR TOUTES LES SURFACES DE L'OBJET
C          >0  DEMANDE DU CALCUL SUR LES SURFACES NOMMEES PAR L'UTILISATEUR
C INTERP : NO D'INTERPOLATION A PRENDRE EN COMPTE
C          1 POLYNOME DE LAGRANGE DE DEGRE 1
C          2 POLYNOME DE LAGRANGE DE DEGRE 2
C          3 POLYNOME DE BREZZI-FORTIN
C
C NTLXOB : NO TMS DU LEXIQUE DE L'OBJET
C MNDFOB : ADRESSE MCN DU TABLEAU DE DEFINITION DE L'OBJET
C
C NBCOMP : NOMBRE DE COMPOSANTE D'UN VECTEUR SOLUTION
C NCAS0  : NUMERO DU PREMIER CAS A TRAITER PARMI LES NCAS0:NCAS1 VECTEURS
C NCAS1  : NUMERO DU DERNIER CAS A TRAITER PARMI LES NCAS0:NCAS1 VECTEURS
C solution: TABLEAU(ncas0:ncas1) de pointeurs sur le tableau solution(NBCOMP)
C NMSOLU : NOM DE LA SOLUTION (PRESSION, VITESSE, ...)
C TIMES  : INSTANT DU CALCUL DE CHACUN DES NCAS0:NCAS1 VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY SEPTEMBRE 2010
C23456---------------------------------------------------------------012
      PARAMETER   ( MXTYEL=7 )
      PARAMETER   ( LIGCON=0 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/donflu.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*90      KNOM
      CHARACTER*4       NOMELE(2)
      CHARACTER*(*)     NMSOLU

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(:), allocatable :: solution

      REAL              TIMES(NCAS0:NCAS1)
      INTEGER           NUMIOB(4),  NUMAOB(4), MNDOEL(4),  MXDOEL(4)
      INTEGER           NONOEF(10), NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NONOFK(8)
C
      MOREE2 = MOTVAR(6)
      MNNOSU = 0
      MNINSO = 0
      MNEXINT= 0
C
C     LECTURE DU NOMBRE ET DES NOMS DES SURFACES DE TRACE DE LA SOLUTION
C     ------------------------------------------------------------------
      IF( MODESF .LE. 0 ) THEN
         NBNOSU = -1
      ELSE
         NBNOSU = 0
      ENDIF
      CALL LINMSURF( MNDFOB,  NBNOSU, MNNOSU )
      IF( NBNOSU .LE. 0 ) RETURN
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF" DE l'OBJET
C     FIND THE TMS XYZSOMMET XYZNOEUD XYZPOINT NPEF"xxxx  of the OBJECT
C     Cf $MEFISTO/td/da/a___xyznoeud   a___npef
C     ===================================================================
      CALL MIMAOB(      1, NTLXOB, MXDOFL,
     %             NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES PLSV
C         NUMAOB          LES 4 NUMEROS MAXIMA DES PLSV
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX MIN MAX DES PLSV
C     MNDOEL THE 4 ADRESSES MCN of THE MIN MAX ARRAY MCN ADDRESSES OF PLSV
C     TABLEAUX DECRIVANT LA THERMIQUE DE L'OBJET COMPLET
C     SI EF BREZZI-FORTIN: XYZP EST REDUIT AUX SOMMETS
C                          XYZN CONTIENT LES SOMMETS ET BARYCENTRES
C
C     NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
C     BREZZI-FORTIN => SOMMETS + BARYCENTRES
      NBNOEU = MCN( MNXYZN + WNBNOE )
C
C     TABLEAU DE L'INTEGRALE DE LA SOLUTION SUR LES FACES DE CHAQUE SURFACE
C     ---------------------------------------------------------------------
C     POUR LES CAS NCAS0 A NCAS1
      NBSUNM = NUMAOB(3) - NUMIOB(3) + 1
      MOINSO = NBSUNM * (NCAS1 - NCAS0 + 1)
      IF( MOINSO .LE. 0 ) GOTO 9000
      CALL TNMCDC( 'REEL2', MOINSO, MNINSO )
      CALL AZEROD( MOINSO, MCN(MNINSO) )
C
C     EXISTENCE OU NON DU CALCUL
      CALL TNMCDC( 'ENTIER', NBSUNM, MNEXINT )
      CALL AZEROI( NBSUNM, MCN(MNEXINT) )
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ----------------------------------------
      NBFACT = 0
      DO 40 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM =  MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF DE CE TYPE
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
C           POINTS DIFFERENTS DES NOEUDS
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
C        --------------------------------------
         DO 30 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
            CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           PARCOURS DES FACES DE L'EF NUELEM
            DO NF = 1, NFACE
C
C              LA FACE NF EST ELLE SUR UNE DES SURFACES NOMMEES
               DO K = 1, NBNOSU
C
C                 NO DE LA SURFACE NOMMEE K
                  NOSURF = MCN(MNNOSU-1+K)
                  IF( NOSURF .EQ. NOOBSF(NF) ) THEN
C
C                    OUI: LA FACE EST SUR LA SURFACE NOSURF => CALCUL
C                    UNE FACE DE PLUS A TRAITER
                     NBFACT = NBFACT + 1
C                    TEMOIN D'EXISTENCE
                     MCN( MNEXINT + NOSURF - NUMIOB(3) ) = 1
C
                     IF( INTERP .EQ. 1 .OR. INTERP .EQ. 3 ) THEN
C
C                       INTERPOLATION DEGRE 1 AVEC LES SEULS SOMMETS
C                       LES SOMMETS DE LA FACE NF DE L'EF
                        NB = NBSOFA(NF)
                        DO N = 1, NB
                           NONOFK( N ) = NONOEF( NOSOFA(N,NF) )
                        ENDDO
C
                     ELSE
C
C                       INTERPOLATION DEGRE 2 AVEC LES NOEUDS
C                       RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE NF
                        CALL ELNOFA( NUTYEL, NF, NB, NONOFK )
C
C                       LES NBNOFK NOEUDS DE LA FACE NF DE L'EF
                        DO N = 1, NB
                           NONOFK( N ) = NONOEF( NONOFK(N) )
                        ENDDO
C
                     ENDIF
C
C                    CALCUL DE L'INTEGRALE DE LA SOLUTION SUR LA FACE
                     DO NCAS = NCAS0, NCAS1
                        MNICAS = MNINSO + MOREE2 * NBSUNM * (NCAS-NCAS0)
                        CALL INTSOLFA( NOSURF, NB, NONOFK,
     %                                 NBCOOR, MCN(MNXYZN+WYZNOE),
     %                                 NCAS,   solution,
     %                                 NUMIOB(3),NUMAOB(3), MCN(MNICAS))
                     ENDDO
C
                  ENDIF
               ENDDO
            ENDDO
 30      CONTINUE
 40   CONTINUE
C
C     AFFICHAGE DE L'INTEGRALE DE LA SOLUTION SUR LES SURFACES NOMMEES
C     POUR LES VECTEURS NCAS0 A NCAS1
C     ================================================================
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'NOMBRE DE FACES DES SURFACES NOMMEES=',NBFACT
      ELSE
         WRITE(IMPRIM,*) 'NUMBER of FACES of NAMED SURFACES=',NBFACT
      ENDIF
C
C     ADRESSE DU TABLEAU DU CAS NCAS
      LNM = NUDCNB( NMSOLU )
      MN  = ( MNINSO + 1 ) / MOREE2 - NUMIOB(3)
      DO NCAS = NCAS0, NCAS1
C
         WRITE(IMPRIM,*)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'Integrale du VECTEUR ',NMSOLU(1:LNM),NCAS,
     %                      ' au TEMPS', TIMES(NCAS)
         ELSE
            WRITE(IMPRIM,*) 'Integral of  ',NMSOLU(1:LNM),' VECTOR ',
     %                      NCAS,' at TIME', TIMES(NCAS)
         ENDIF
C
C        AFFICHAGE DE L'INTEGRALE DE LA SOLUTION DU CAS NCAS
C        SUR LES SURFACES NOMMEES SI ELLE EST NON NULLE
         N = 0
         DO K = NUMIOB(3), NUMAOB(3)
C           L'INTEGRALE DE LA SURFACE K A-T-ELLE ETE CALCULEE
            IF( MCN( MNEXINT + K - NUMIOB(3) ) .NE. 0 ) THEN
C              OUI: IMPRESSION DE L'INTEGRALE
               WRITE(KERR(MXLGER),'(G15.7)') DMCN(MN+K)
               N = N + 1
C              LE NOM DE LA SURFACE
               CALL NMOBNU( 'SURFACE', K, KNOM )
               NN = NUDCNB( KNOM )
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(N) = 'INTEGRALE ' // NMSOLU(1:LNM)
     %                    // KERR(MXLGER)(1:15)
     %                    // ' sur la SURFACE ' // KNOM(1:NN)
               ELSE
                  KERR(N) = 'INTEGRAL ' // NMSOLU(1:LNM)
     %                    // KERR(MXLGER)(1:15)
     %                    // ' on the SURFACE ' // KNOM(1:NN)
               ENDIF
            ENDIF
         ENDDO
         IF( N .GT. 0 ) THEN
            NBLGRC(NRERR) = N
            CALL LERESU
         ENDIF
C
C        PASSAGE AU CAS SUIVANT
         MN = MN + NBSUNM
      ENDDO
C     SORTIE DU TRACE DES ARETES ET FACES
C     ===================================
 9000 IF( MNINSO  .GT. 0 ) CALL TNMCDS( 'REEL2',  MOINSO, MNINSO )
      IF( MNNOSU  .GT. 0 ) CALL TNMCDS( 'ENTIER', NBNOSU, MNNOSU )
      IF( MNEXINT .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSUNM, MNEXINT )
      RETURN
      END
