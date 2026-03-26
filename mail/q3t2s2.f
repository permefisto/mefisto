      SUBROUTINE Q3T2S2( NT,     NTOP,   NSOP,   GRAND,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER 2 TETRAEDRES A PARTIR DE 3 TETRAEDRES ET UN SOMMET
C -----    DOUBLE OPPOSE A UN TETRAEDRE DE MAUVAISE QUALITE
C
C ENTREES :
C ---------
C NT     : NO DU TETRAEDRE EN COURS D'AMELIORATION DE SA QUALITE
C NTOP   : NUMERO DES 4 TETRAEDRES OPPOSES PAR LES 4 FACES
C          0 SI PAS DE TETRAEDRE OPPOSE
C NSOP   : NUMERO DES 4 SOMMETS OPPOSES DES TETRAEDRES OPPOSES
C GRAND  : LE PLUS GRAND REEL STOCKABLE
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C
C MODIFIES :
C ----------
C N1TEVI : NUMERO DU PREMIER TETRAEDRE VIDE DANS NSTETR
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C N1TESO : NUMERO DU PREMIER TETRAEDRE NON UTILISE DANS NOTESO
C          CHAINAGE DES PLACES VIDES DANS NOTESO(2,*)
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C MXTETA : NUMERO MAXIMAL DES TETRAEDRES ACTUELS
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C IERR   : = 0 SI PAS D'ERREUR ET PAS DE MODIFICATION
C          = 1 SI PAS DE SOMMET DOUBLE OPPOSE A 2 FACES DE NT
C          =-1 SI MODIFICATION EFFECTUEE
C          > 0 EN CAS DE SATURATION D'UN TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1992
C....................................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NTOP(4), NSOP(4)
      INTEGER           NSTETR(NBSOTE,*),
     %                  NUDMEF(*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*)
      REAL              XYZSOM(3,*),VOLUMT(*),QUALIT(*)
      REAL              ARMIN,ARMAX,SURFTR(4)
      INTEGER           NOSOFATE(3,4)
C                       NO 3 SOMMETS DES 4 FACES ORIENTEES POUR QUE
C                       NORMALE SOIT VERS L'INTERIEUR DU TETRAEDRE
      DATA              NOSOFATE / 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

C     RECHERCHE D'UN SOMMET OPPOSE DOUBLE
      DO 10 I=1,3
         NS0 = NSOP(I)
         IF( NS0 .EQ. 0 ) GOTO 10
         DO J=I+1,4
            IF( NS0 .EQ. NSOP(J) ) GOTO 20
         ENDDO
 10   ENDDO
      IERR = 1
      RETURN
C
C     IL EXISTE NS0 UN SOMMET DOUBLE APPARTENANT AUX 2 TETRAEDRES
C     ==========================================================
C     CALCUL DU NOMBRE DE TETRAEDRES ET SOMMETS OPPOSES
 20   NFR = 0
      DO K=1,4
         IF( NTOP(K) .GT. 0 ) NFR = NFR + 1
      ENDDO

      IF( NFR .EQ. 4 ) THEN
C        IL PEUT EXISTER UNE AUTRE POSSIBILITE DE TROUVER UN AUTRE
C        SOMMET DOUBLE
         K1 = I+1
         NS5 = NSOP(K1)
         DO K2=K1+1,4
            IF( NSOP(K2) .EQ. NS5 ) GOTO 30
         ENDDO
C        PAS DE SECOND SOMMET DOUBLE
         GOTO 40

C        ENTRE NS0 et NS5 LES 2 SOMMETS OPPOSES DOUBLES OU ENCORE
C        CONSERVER L'ARETE 13 OU L'ARETE 24
C        LEQUEL CHOISIR?
C        CELUI QUI DONNE LE MAX DES MIN DES QUALITES DES 2 CHOIX
C        DES 4 TETRAEDRES (4 SOMMETS DE NT + NS0 + NS5)
C        -------------------------------------------------------
 30      NS1 = NSTETR( NOSOFATE(1,I), NT )
         NS2 = NSTETR( NOSOFATE(3,I), NT )
         NS3 = NSTETR( NOSOFATE(2,I), NT )

C        4-EME SOMMET DU TETRAEDRE NT
         NS4 = NOSOFATE(1,I)
         IF( NS4 .EQ. 1 ) THEN
            NS4 = 4
         ELSE
            NS4 = NS4 - 1
         ENDIF
         NS4 = NSTETR( NS4, NT )

C        JEU AVEC L'ARETE 13 DE NT
C        TETRAEDRE 1230
         CALL QUATET( XYZSOM( 1, NS1 ),
     %                XYZSOM( 1, NS2 ),
     %                XYZSOM( 1, NS3 ),
     %                XYZSOM( 1, NS0 ),
     %                ARMIN, ARMAX, SURFTR, V1, Q1 )

C        TETRAEDRE 1340
         CALL QUATET( XYZSOM( 1, NSTETR( NOSOFATE(1,J), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(3,J), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(2,J), NT ) ),
     %                XYZSOM( 1, NS0 ),
     %                ARMIN, ARMAX, SURFTR, V2, Q2 )

C        TETRAEDRE 1325
         CALL QUATET( XYZSOM( 1, NS1 ),
     %                XYZSOM( 1, NS3 ),
     %                XYZSOM( 1, NS2 ),
     %                XYZSOM( 1, NS5 ),
     %                ARMIN, ARMAX, SURFTR, V3, Q3 )

C        TETRAEDRE 1435
         CALL QUATET( XYZSOM( 1, NSTETR( NOSOFATE(1,J), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(2,J), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(3,J), NT ) ),
     %                XYZSOM( 1, NS5 ),
     %                ARMIN, ARMAX, SURFTR, V4, Q4 )

         QUA13 = MIN( Q1, Q2, Q3, Q4 )
ccc         VOL13 = V1 + V2 + V3 + V4
         VOL13 = ABS(V1) + ABS(V2) + ABS(V3) + ABS(V4)

C        JEU AVEC L'ARETE 24 DE NT
C        TETRAEDRE 2340
         CALL QUATET( XYZSOM( 1, NSTETR( NOSOFATE(1,K1), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(2,K1), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(3,K1), NT ) ),
     %                XYZSOM( 1, NS0 ),
     %                ARMIN, ARMAX, SURFTR, V1, Q1 )

C        TETRAEDRE 2410
         CALL QUATET( XYZSOM( 1, NSTETR( NOSOFATE(1,K2), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(2,K2), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(3,K2), NT ) ),
     %                XYZSOM( 1, NS0 ),
     %                ARMIN, ARMAX, SURFTR, V2, Q2 )

C        TETRAEDRE 2435
         CALL QUATET( XYZSOM( 1, NSTETR( NOSOFATE(1,K1), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(3,K1), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(2,K1), NT ) ),
     %                XYZSOM( 1, NS5 ),
     %                ARMIN, ARMAX, SURFTR, V3, Q3 )
C        TETRAEDRE 2145
         CALL QUATET( XYZSOM( 1, NSTETR( NOSOFATE(1,K2), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(3,K2), NT ) ),
     %                XYZSOM( 1, NSTETR( NOSOFATE(2,K2), NT ) ),
     %                XYZSOM( 1, NS5 ),
     %                ARMIN, ARMAX, SURFTR, V4, Q4 )

         QUA24 = MIN( Q1, Q2, Q3, Q4 )
ccc         VOL24 = V1 + V2 + V3 + V4
         VOL24 = ABS(V1) + ABS(V2) + ABS(V3) + ABS(V4)

C        VERIFICATION A SUPPRIMER ENSUITE
         IF( ABS(VOL13-VOL24) .GT. 1E-4*(VOL13+VOL24) ) THEN
            PRINT *,'q3t2s2: VOLUMES DIFFERENTS V13=',VOL13,
     %      ' V24=',VOL24,' V1=',V1,' V2=',V2,' V3=',V3,' V4=',V4,
     %      ' A REVOIR'
         ENDIF

         IF( QUA13 .GT. QUA24 ) THEN
C           L'ARETE 13 EST CONSERVEE
C           L'ARETE 24 EST SUPPRIMEE
            PRINT *,'q3t2s2: REUSSI suppression arete 24:',NS4,NS2,
     %              ' Qualite 24=',QUA24, ' conservation arete 13:',
     %              NS1,NS3,' Qualite 13=',QUA13
            I = K1
            J = K2
            NS0 = NS5
         ENDIF

      ENDIF


C     LE SOMMET NS0 DOUBLE OPPOSE EST CHOISI
C     --------------------------------------
 40   NT1 = NTOP(I)
      NT2 = NTOP(J)
C     LA FACE COMMUNE DES 2 TETRAEDRES NT1 NT2
      CALL NOFC2TE( NT1, NT2, NBSOTE, NSTETR, NF1,NF2 )
C
      NUMATE = 0
      IF( NBDM .GT. 0 ) THEN
C        PAS DE MODIFICATION SI LES 3 TRIANGLES NE SONT PAS
C        DANS LE MEME VOLUME
         IF( NUDMEF(NT ) .NE. NUDMEF(NT1) .OR.
     %       NUDMEF(NT1) .NE. NUDMEF(NT2) ) RETURN
         NUMATE = NUDMEF(NT)
      ENDIF
C
      NS1 = NSTETR( NOSOFATE(1,NF1), NT1 )
      NS2 = NSTETR( NOSOFATE(2,NF1), NT1 )
      NS3 = NSTETR( NOSOFATE(3,NF1), NT1 )

C     ROTATION POUR METTRE NS0 EN POSITION NS3
      IF( NS0 .NE. NS3 ) THEN
         L   = NS1
         NS1 = NS2
         NS2 = NS3
         NS3 = L
      ENDIF
      IF( NS0 .NE. NS3 ) THEN
         L   = NS1
         NS1 = NS2
         NS2 = NS3
         NS3 = L
      ENDIF
ccc      PRINT *,'q3t2s2: disparition arete   :',NS1,NS2,
ccc     %        ' du tetredre',nt

      IF( NF1 .EQ. 1 ) THEN
         NS4 = 4
      ELSE
         NS4 = NF1 - 1
      ENDIF
C     LE SOMMET OPPOSE DANS NT1
      NS4 = NSTETR( NS4, NT1 )

      IF( NF2 .EQ. 1 ) THEN
         NS5 = 4
      ELSE
         NS5 = NF2 - 1
      ENDIF
C     LE SOMMET OPPOSE DANS NT2
      NS5 = NSTETR( NS5, NT2 )

C     COMPARAISON DES QUALITES AVANT ET APRES
      QUALI0 = MIN( QUALIT(NT), QUALIT(NT1), QUALIT(NT2) )
      CALL QUATET( XYZSOM(1,NS1), XYZSOM(1,NS4),
     %             XYZSOM(1,NS5), XYZSOM(1,NS3),
     %             ARMIN, ARMAX, SURFTR, V, QUALI1 )
      CALL QUATET( XYZSOM(1,NS5), XYZSOM(1,NS4),
     %             XYZSOM(1,NS2), XYZSOM(1,NS3),
     %             ARMIN, ARMAX, SURFTR, V, Q )
      QUALI1 = MIN( QUALI1, Q )

      IF( QUALI0 .GE. QUALI1 ) THEN
C        TRANSFORMATION INFERIEURE
         IERR = 0
         RETURN
      ENDIF

C     SUPPRESSION DU TETRAEDRE NT
      CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %             NO1TSO, N1TESO, NOTESO )
      QUALIT( NT ) = GRAND
      IF( NUMATE .GT. 0 ) NUDMEF( NT ) = 0

C     SUPPRESSION DU TETRAEDRE NT1
      CALL DS1TET( NT1,    NBSOTE, N1TEVI, NSTETR,
     %             NO1TSO, N1TESO, NOTESO )
      QUALIT( NT1 ) = GRAND
      IF( NUMATE .GT. 0 ) NUDMEF( NT1 ) = 0

C     SUPPRESSION DU TETRAEDRE NT2
      CALL DS1TET( NT2,    NBSOTE, N1TEVI, NSTETR,
     %             NO1TSO, N1TESO, NOTESO )
      QUALIT( NT2 ) = GRAND
      IF( NUMATE .GT. 0 ) NUDMEF( NT2 ) = 0
C
      IF( NBDM .GT. 0 ) THEN
         NUDMEF( NT  ) = 0
         NUDMEF( NT1 ) = 0
         NUDMEF( NT2 ) = 0
      ENDIF

C     CREATION DES TETRAEDRES 1453 5423
      CALL Q05S2T( NS1, NS2, NS3, NS4, NS5, NUMATE, NUDMEF,
     %             XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %             NO1TSO, N1TESO, NOTESO,
     %             MXTETA, VOLUMT, QUALIT, IERR )
      IERR = -1

      RETURN
      END
