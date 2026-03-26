      SUBROUTINE HAARSU( XYZSOM, NMSURF, NBMOTS, NTARSU, MNARSU, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CREER SI CE N'EST DEJA FAIT LE TABLEAU DES ARETES D'UNE SURFACE
C -----  CHAINER LES ARETES FRONTALIERES
C ENTREES:
C --------
C NMSURF : NOM DE LA SURFACE
C NBMOTS : NOMBRE DE MOTS D'INFORMATIONS A STOCKER EN PLUS PAR ARETE(>=4)
C
C SORTIES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE DE LA SURFACE
C NTARSU : NUMERO DU TMS DES ARETES DE LA SURFACE
C          0 SI NON CREE PAR CAUSE D'ERREUR
C MNARSU : ADRESSE MCN DU TABLEAU DES ARETES DE LA SURFACE
C          cf ~/td/d/a___arete
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR
C
C
C QUELQUES INFORMATIONS SUR LE TMS ARETE :
C MOARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU ARETES DES ARETES
C MXARET : NOMBRE D'ARETES DU TABLEAU ARETES
C
C LE TABLEAU ARETES CONTIENT :
C          ARETES(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          ARETES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          ARETES(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          SI SOMMET 1 < SOMMET 2    => ARETE   DIRECTE DANS LA FACE
C                      <             => ARETE INDIRECTE
C          ARETES(4,I)= NUMERO DE LA 1-ERE FACE CONTENANT CETTE ARETE
C                       >0 SI ARETE   DIRECTE DANS CETTE FACE
C                       <0 SI ARETE INDIRECTE DANS CETTE FACE
C          SI L'ARETE APPARTIENT A 2 FACES ALORS
C          ARETES(5,I)= NUMERO DE LA 2-EME FACE CONTENANT CETTE ARETE
C                       >0 SI ARETE   DIRECTE DANS CETTE FACE
C                       <0 SI ARETE INDIRECTE DANS CETTE FACE
C          ARETES(6,I)= NUMERO DE L'ARETE FRONTALIERE SUIVANTE
C                       0 SI C'EST LA DERNIERE
C                       0 SI ARETE NON FRONTALIERE
C          L1ARFB = NUMERO DE LA PREMIERE ARETE FRONTALIERE DANS ARETES
C          ARETES(7,I)= NUMERO DE L'ARETE A TG, 0 SINON
C
C          NUTGAR(1:2,1:NBARTG)= NUMERO DES 2 TANGENTES DES ARETES A TG
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___arete.inc"
      include"./incl/ntmnlt.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      REAL              XYZSOM(3,*)
      CHARACTER*(*)     NMSURF
      INTEGER           NOSOEL(12),NGS(2)
      EQUIVALENCE      (NS1,NGS(1)),(NS2,NGS(2))
      LOGICAL           AVANT
C
      IERR = 0
C
C     LA SURFACE INITIALE
C     ===================
C     LE TABLEAU LEXIQUE DE CETTE SURFACE
      CALL LXLXOU( NTSURF, NMSURF, NTLXSF, MN )
      IF( NTLXSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMSURF
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE INCONNUE'
         ELSE
            KERR(2) = 'UNKNOWN SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'NSEF', NTFASF, MNFASF )
      IF( NTFASF .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMSURF
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE SANS TMS NSEF'
         ELSE
            KERR(2) = 'SURFACE WITHOUT NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTSOSF, MNSOSF )
      IF( NTSOSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMSURF
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE SANS XYZSOMMET'
         ELSE
            KERR(2) = 'SURFACE WITHOUT XYZSOMMET TMS'
         ENDIF
          CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DU MAILLAGE DE LA SURFACE
      NBSOSF = MCN( MNSOSF + WNBSOM )
C     LE NOMBRE DE TANGENTES DU MAILLAGE DE LA SURFACE
      NBTGSF = MCN( MNSOSF + WNBTGS )
C
C     LE TABLEAU 'ARETE' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'ARETE', NTARSU, MNARSU )
      IF( NTARSU .GT. 0 ) THEN
C        LES ARETES SONT ELLES ANTERIEURES AUX SOMMETS ?
         IF( AVANT( MCN(MNSOSF), MCN(MNARSU) ) ) THEN
C           LES SOMMETS SONT ANTERIEURS AUX ARETES
C           LE TABLEAU ARETE EST ACCEPTE
            RETURN
         ENDIF
C        LE TABLEAU ARETE EST DETRUIT POUR ETRE ENSUITE REGENERE
         CALL LXTSDS( NTLXSF, 'ARETE' )
         NTARSU = 0
         MNARSU = 0
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNFASF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBFASU,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     NOMBRE D'ENTIERS POUR UNE ARETE
      MOARET = 3 + NBMOTS
C
C     MAJORATION DU NOMBRE DE ARETES : FORMULE A MODIFIER EVENTUELLEMENT
      IF( NBFASU .LT. 100 ) THEN
         MXARET = NBFASU * 4
      ELSE
         MXARET = NBFASU * 3 + 64
      ENDIF
C
C     LE NOMBRE MAXIMUM D'ARETES A TG (TRIANGLE OU QUADRANGLE)
      MXARTG = 4 * MCN( MNFASF + WBEFTG )
C
C     LE NOMBRE D'ARETES AVEC TG
      NBARTG = 0
C
C     ADRESSAGE ET OUVERTURE DU TABLEAU ARETES
C     ----------------------------------------
      CALL LXTNDC( NTLXSF, 'ARETE', 'ENTIER',
     %             W1LGFR + MOARET*MXARET + 2*MXARTG )
      CALL LXTSOU( NTLXSF, 'ARETE', NTARSU, MNARSU )
C
C     LE TABLEAU DES ARETES EST INITIALISE A ZERO
      MNARET = MNARSU + W1LGFR
      CALL AZEROI( MOARET*MXARET, MCN(MNARET) )
C
C     LE TABLEAU DES NUMEROS DES 2 TG DES FACES A TG
      MNTGAR = MNARET + MOARET*MXARET
C
C     LA 1-ERE ARETE LIBRE
      LIBREF = MXARET
C
C     LA BOUCLE SUR LES NO SOMMET DU MAILLAGE
C     ----------------------------------------
C     LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
      NBAR2F = 0
C
      DO 1000 N=1,NBFASU
C
C        LE NUMERO DES SOMMETS DE LA FACE N DE LA SURFACE
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNFASF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( IERR .NE. 0 ) GOTO 9990
C
C        BOUCLE SUR LES NCOGEL ARETES DE LA FACE
C        ---------------------------------------
         DO 200 NOA=1,NCOGEL
C           LE NUMERO LOCAL DES SOMMETS
            NS2 = NOA + 1
            IF( NS2 .GT. NCOGEL ) NS2 = 1
C           LE NUMERO GLOBAL DES SOMMETS DE L'ARETE
            NS1 = NOSOEL(NOA)
            NS2 = NOSOEL(NS2)
            IF( NS1 .GT. NS2 ) THEN
               L     = NS1
               NS1   = NS2
               NS2   = L
               ISENS =-1
            ELSE
               ISENS = 1
            ENDIF
C
C           1:NO LE PLUS FAIBLE DES SOMMETS DE L'ARETE
C           2:NO LE PLUS GRAND  DES SOMMETS DE L'ARETE
C           SI MIN = SOMMET 2  ARETE DIRECTE   DANS CETTE FACE
C                    SOMMET NS ARETE INDIRECTE DANS CETTE FACE
C
C           RECHERCHE OU ADJONCTION DE L'ARETE
C           ----------------------------------
            CALL HACHAG(2,NGS,MOARET,MXARET,MCN(MNARET),3,
     &                  LIBREF,NOARE)
            MN = MNARET + MOARET * ABS(NOARE) - MOARET
C
C           STOCKAGE ADRESSE DE L ELEMENT DANS ARETE(4,...;NOARE)
C           -----------------------------------------------------
C           PAR ARETE:st1,st2,lien,ef1,ef2,arete bord ou interne suivante,no are
            K = 2
 140        K = K + 1
            IF( K .LT. 5 ) THEN
C              RECHERCHE D'UN NUMERO NUL DE FACE
               IF( MCN(MN + K) .NE. 0) GOTO 140
C
C              NO ELEMENT > 0 SI L'ARETE EST   DIRECTE DANS L ELEMENT
C                         < 0 SI L'ARETE EST INDIRECTE DANS L ELEMENT
               MCN( MN + K ) = N * ISENS
C
            ELSE
C
C              ERREUR IL Y A PLUS DE 2 FACES CONTENANT CETTE ARETE
               IF( NBAR2F .EQ. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'TROP DE FACES CONTIENNENT CETTE ARETE'
                  ELSE
                     KERR(1) = 'TOO MANY FACES WITH THIS EDGE'
                  ENDIF
                  CALL LEREUR
                  IERR = 6
               ENDIF
               IF( NOSOEL(4) .EQ. 0 ) THEN
                  NT2 = 3
               ELSE
                  NT2 = 4
               ENDIF
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10200) N,(NOSOEL(L),L=1,NT2)
                  WRITE(IMPRIM,10201) NOA,(NGS(L),L=1,2),
     %                            ABS(MCN(MN+3)), ABS(MCN(MN+4))
                  DO 150 NT1=1,2
                     WRITE(IMPRIM,10150) NGS(NT1),
     %                    (XYZSOM(L,NGS(NT1)),L=1,3)
 150              CONTINUE
               ELSE
                  WRITE(IMPRIM,20200) N,(NOSOEL(L),L=1,NT2)
                  WRITE(IMPRIM,20201) NOA,(NGS(L),L=1,2),
     %                            ABS(MCN(MN+3)), ABS(MCN(MN+4))
                  DO 151 NT1=1,2
                     WRITE(IMPRIM,20150) NGS(NT1),
     %                    (XYZSOM(L,NGS(NT1)),L=1,3)
 151              CONTINUE
               ENDIF
               NBAR2F = NBAR2F + 1
            ENDIF
C
C           ARETE AVEC TANGENTES?
            IF( NUEFTG .GT. 0 ) THEN
C
C              LA FACE A DES TANGENTES
C              -----------------------
C              NUMERO DES 2 TANGENTES DE CETTE ARETE NOA SELON SON SENS
               NT1 = 2 * NOA - 1
               IF( NOA .NE. NCOGEL ) THEN
                  NT2 = NT1 + 3
               ELSE
                  NT2 = 2
               ENDIF
               IF( ISENS .LT. 0 ) THEN
C                 IL FAUT PERMUTER LES TANGENTES
                  L   = NT1
                  NT1 = NT2
                  NT2 = L
               ENDIF
               NT1 = NOSOEL( 4 + NT1 )
               NT2 = NOSOEL( 4 + NT2 )
C
C              L'ARETE A T ELLE DEJA DES TANGENTES ?
               NUARTG = MCN( MN + 6 )
               IF( NUARTG .EQ. 0 ) THEN
C
C                 NON. CELLES DE L'ARETE NF DE LA FACE
C                      DEFINISSENT CELLES DE L'ARETE
                  MNTF = MNTGAR + 2 * NBARTG
C                 UNE ARETE A TANGENTE DE PLUS
                  NBARTG = NBARTG + 1
                  MCN( MN   + 6 ) = NBARTG
                  MCN( MNTF     ) = NT1
                  MCN( MNTF + 1 ) = NT2
C
               ELSE
C
C                 OUI. CELLES DE L'ARETE NOA
C                      SONT REFONDUES AVEC CELLES PRECEDENTES
C                      AVEC PRIORITE A LA PRECEDENTE NON NULLE,
C                      PUIS, A LA NOUVELLE TG
                  MNTF = MNTGAR + 2 * NUARTG - 2
                  IF( NT1 .NE. 0 ) THEN
                     IF( MCN( MNTF ) .EQ. 0 ) THEN
C                       AJOUT DE LA TG SI PAS DE PRECEDENTE
                        MCN( MNTF ) = NT1
                     ENDIF
                  ENDIF
                  IF( NT2 .NE. 0 ) THEN
                     IF( MCN( MNTF + 1 ) .EQ. 0 ) THEN
C                       AJOUT DE LA TG SI PAS DE PRECEDENTE
                        MCN( MNTF + 1 ) = NT2
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
C
 200     CONTINUE
 1000 CONTINUE
C
      IF( NBAR2F .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)NBAR2F,' ARETES APPARTIENNENT A PLUS DE 2 FACES'
         ELSE
         WRITE(IMPRIM,*)NBAR2F,' EDGES IN MORE OF 2 FACES'
         ENDIF
      ENDIF
C
10200 FORMAT(' HAARSU: FACE:',I7,' DE SOMMETS:',T34,4I6)
10201 FORMAT(' L''ARETE',I2,' DE SOMMETS:',T34,2I7,
     %       ' APPARTIENT A PLUS DE 2 FACES => ERREUR'/
     %       ' NUMERO DES 2 FACES:',2I7/
     %       ' VERIFIER L''UNICITE DE CHAQUE SURFACE DANS LES UNIONS'/)
10150 FORMAT(' SOMMET',I7,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
C
20200 FORMAT(' HAARSU: FACE:',I7,' of VERTICES:',T34,4I6)
20201 FORMAT(' EDGE',I2,' with VERTICES:',T34,2I7,
     %       ' BELONGS TO MORE 2 FACES => ERROR'/
     %       ' NUMBER of 2 FACES:',2I7/
     %       ' VERIFY the UNICITY of EACH SURFACE in the UNIONS'/)
20150 FORMAT(' VERTEX',I7,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
C
C     CHAINAGE DES ARETES FRONTALIERES C-A-D APPARTENANT A UNE SEULE FACE
      MN     = MNARET
      L1ARFB = 0
      NBARFB = 0
      DO 2000 N=1,MXARET
C
C        ELIMINATION DES ARETES NON UTILISEES DANS MNARET
         IF( MCN( MN ) .EQ. 0 ) GOTO 1999
C
C        ELIMINATION DES ARETES APPARTENANT A 2 FACES
         IF( MCN( MN + 4 ) .NE. 0 ) GOTO 1999
C
C        L'ARETE EST FRONTALIERE. ELLE EST CHAINEE AVEC LA PRECEDENTE
         MCN( MN + 5 ) = L1ARFB
         L1ARFB        = N
         NBARFB        = NBARFB + 1
C
C        PASSAGE A L'ARETE SUIVANTE
 1999    MN = MN + MOARET
 2000 CONTINUE
C
C     LE TABLEAU 'ARETE' EST COMPLETE
C     LE NOMBRE D'ENTIERS PAR ARETE
      MCN( MNARSU + WOARET ) = MOARET
C     LA MAJORATION DU NOMBRE DE ARETES
      MCN( MNARSU + WXARET ) = MXARET
C     LE NOMBRE DE ARETES FRONTALIERES
      MCN( MNARSU + WBARFB ) = NBARFB
C     LE NUMERO DE LA PREMIERE ARETE FRONTALIERE
      MCN( MNARSU + W1ARFB ) = L1ARFB
C     LE NOMBRE DE ARETES INTERFACES . NON CALCUL'E
      MCN( MNARSU + WBARIN ) = 0
C     LE NUMERO DE LA PREMIERE ARETE INTERFACE . NON CALCUL'E
      MCN( MNARSU + W1ARIN ) = 0
C     LE NOMBRE D'ARETES AVEC TG
      MCN( MNARSU + WBARTG ) = NBARTG
C     LE NUMERO MINIMAL DES LIGNES FRONTIERES . NON CALCUL'E
      MCN( MNARSU + WUMILF ) = 0
C     LE NUMERO MAXIMAL DES LIGNES FRONTIERES . NON CALCUL'E
      MCN( MNARSU + WUMXLF ) = 0
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARSU) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARSU + MOTVAR(6) ) = NONMTD( '~>>>ARETE' )
C
C     REDUCTION EVENTUELLE
      IF( NBARTG .LT. MXARTG ) THEN
         CALL TAMSRA( NTARSU, W1LGFR + MOARET*MXARET + 2*NBARTG )
      ENDIF
      RETURN
C
C     ERREUR DETECTEE: LE TABLEAU ARETE INCOMPLET EST DETRUIT
 9990 CALL LXTSDS( NTLXSF, 'ARETE' )
      END
