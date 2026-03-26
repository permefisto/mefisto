      SUBROUTINE DFTOP3( MOELEM, NBELEM, MNELEM,
     %                   MOELTG, NBELTG, MNELTG, NUELTG,
     %                   LIENEL, NUDREL,
     %                   NBOBPR, NUOBPR, MNMNSO,
     %                   NDSCOU, NUSOUN,
     %                   NDTGOU, NUTGUN,
     %                   NDIM,   XYZ,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AJOUTER LES EF DES OBJETS PREMIERS AUX TABLEAUX DE HACHAGE DES
C -----   NOEUDSOMMETS ARETES TRIANGLES QUADRANGLES
C         TETRAEDRES PENTAEDRES HEXAEDRES 6-CUBES
C         DE L'OBJET COMPLET
C
C ENTREES:
C --------
C MOELEM : NOMBRE DE MOTS DE CHAQUE TYPE D'ELEMENT FINI
C NBELEM : NBELEM(1)=NOMBRE DE NOEUDSOMMETS
C          NBELEM(2)=SEGMENTS
C          NBELEM(3)=TRIANGLES
C          NBELEM(4)=QUADRANGLES
C          NBELEM(5)=TETRAEDRES
C          NBELEM(6)=PENTAEDRES
C          NBELEM(7)=HEXAEDRES
C          NBELEM(8)=6-CUBES
C          NBELEM(9)=PYRAMIDES
C MNELEM : ADRESSE MCN DES EVENTUELS 9 TABLEAUX DES ELEMENTS FINIS
C          0 SI PAS D'ELEMENTS DE CE TYPE
C MOELTG : NOMBRE DE TG PAR TYPE D'EF
C NBELTG : NOMBRE D'EF A TG PAR TYPE
C MNELTG : ADRESSE MCN DES EVENTUELS 8 TABLEAUX DES ELEMENTS FINIS A TG
C LIENEL : POSITION DU CHAINAGE POUR LE SP HACHAG
C NUDREL : POSITION DE LA 1-ERE COLONNE SUSCEPTIBLE D'ETRE LIBRE
C NBOBPR : NOMBRE D'OBJETS PREMIERS
C NUOBPR : LE TYPEOBJET DE CHAQUE OBJET PREMIER
C MNMNSO : ADRESSE MCN DU TABLEAU 'NSEF' DE CHAQUE OBJET PREMIER
C NDSCOU : NDSCOU(0)=0 NDSCOU(I)=POINTE DANS NUSOUN SUR LE DERNIER
C                                SOMMET DE L'OBJET PREMIER I
C NUSOUN : NUMERO DU SOMMET DANS L'UNION DES SOMMETS
C NDTGOU : NDTGOU(0)=0 NDTGOU(I)=POINTE DANS NUSOUN SUR LA DERNIERE
C                                TANGENTE DE L'OBJET PREMIER I
C NUSOUN : NUMERO DE LA TG DANS L'UNION DES TANGENTES
C NDIM   : DIMENSION (1 OU 2 OU 3 OU 6) DE L'ESPACE DES EF
C XYZ    : LES 3 COORDONNEES DES SOMMETS DES EF DE L'OBJET UNION
C
C SORTIES:
C --------
C NUELTG : NOMBRE D'EF A TG DECLARES DANS DFTOP3
C IERR   : 0 SI PAS D'ERREUR, NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1996
C MODIF  : ALAIN PERRONNET  TEXAS A & M UNIVERSITY            JULY  2005
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/a___nsef.inc"
      INTEGER           MOELEM(9), NBELEM(9), MNELEM(9),
     %                  MOELTG(9), NBELTG(9), MNELTG(9), NUELTG(9)
      INTEGER           LIENEL(9),
     %                  NUDREL(9), NUOBPR(2,NBOBPR), MNMNSO(NBOBPR),
     %                  NDSCOU(0:NBOBPR), NUSOUN(1:*),
     %                  NDTGOU(0:NBOBPR), NUTGUN(1:*)
      INTEGER           NOSOEL(64)
      REAL              XYZ(3,*)
      CHARACTER*8       KNOM
      CHARACTER*24      KNOM1,KNOM2
C
      IERR = 0
C
C     INITIALISATION DES HACHAGES
      DO 10 I=1,9
         NUDREL(I) = NBELEM(I)
         NUELTG(I) = 0
 10   CONTINUE
C
C     BOUCLE SUR LES OBJETS PREMIERS
C     ==============================
      NERR = 0
      DO 1000 I=1,NBOBPR
C
C        LE NUMERO DU TYPE DE L'OBJET
CCC         NUTYOB = NUOBPR(1,I)
C        LE NUMERO DE L'OBJET
         NUOBJT = NUOBPR(2,I)
C
C        LE TABLEAU NSEF DE L'OBJET
         MNSOOB = MNMNSO(I)
c
C        LES PARAMETRES DE RECUPERATION DES NO DES SOMMETS DU MAILLAGE
         CALL NSEFPA( MCN(MNSOOB),
     %                NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %                NX    , NY    , NZ    ,
     %                IERR  )
         IF( IERR .NE. 0 ) GOTO 1000
C
C        LA BOUCLE SUR LE TABLEAU NSEF DU MAILLAGE DE CET OBJET PREMIER
C        --------------------------------------------------------------
C        LE DEBUT DU TABLEAU DES NUMEROS DES NOEUDS DU MAILLAGE
         DO 900 N=1,NBEFOB
C
C           LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
            CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNSOOB, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C
C           LE NOMBRE DE SOMMETS DE CET ELEMENT FINI
            NBSO = NBSOME( NCOGEL )
C           TRANSFORMATION DU NUMERO DES SOMMETS DE CET OBJET PREMIER I
C           EN LES NUMEROS DE SOMMETS DE L'UNION
            DO 12 J=1,NBSO
C              LE NUMERO DU SOMMET DANS L'UNION APRES IDENTIFICATION
               NOSOEL(J) = NUSOUN( NDSCOU(I-1) + NOSOEL(J) )
 12         CONTINUE
C
            IF( NUEFTG .GT. 0 ) THEN
C              EF A TG : LE DECALAGE DANS NOSOEL PUIS LE NUMERO DES TG
               NU1TG = NBSOEF + 1
               DO 14 J = 1, MOELTG(NCOGEL)
C                 LE NUMERO INITIAL DE LA TG J DE L'EF A TG
                  NT = NOSOEL( NBSOEF + J )
                  IF( NT .EQ. 0 ) THEN
C                    PAS DE TANGENTE
                     GOTO 14
                  ELSEIF( NT .LT. 0 ) THEN
C                    TANGENTE OPPOSEE
                     LSIGNE = -1
                  ELSE
C                    TANGENTE TELLE QUE
                     LSIGNE = 1
                  ENDIF
C                 LE NUMERO DE LA TG APRES IDENTIFICATION DANS L'UNION
                  NT = NUTGUN( NDTGOU(I-1) + ABS(NT) )
                  IF( NT .LT. 0 ) LSIGNE = -LSIGNE
                  NOSOEL(NBSOEF+J) = LSIGNE * ABS(NT)
 14            CONTINUE
            ELSE
C              EF SANS TG
               NU1TG = 0
            ENDIF
C
            IF( NDIM .EQ. 2 ) THEN
               IF( NBSO .EQ. 3 .OR. NBSO .EQ. 4 ) THEN
C                 SURFACE ALGEBRIQUE DU TRIANGLE DES SOMMETS 1 2 3
                  S = SURTR2( XYZ(1,NOSOEL(1)), XYZ(1,NOSOEL(2)),
     %                        XYZ(1,NOSOEL(3)) )
                  IF( S .LT. 0 ) THEN
C                    REORIENTATION DU TRIANGLE OU QUADRANGLE
C                    PERMUTATION DES SOMMETS 2 ET NBSO
                     L            = NOSOEL(2)
                     NOSOEL(  2 ) = NOSOEL(NBSO)
                     NOSOEL(NBSO) = L
                     IF( NUEFTG .GT. 0 ) THEN
C                       PERMUTATION DES TANGENTES DE L'EF
                        IF( NBSO .EQ. 4 ) THEN
C                          QUADRANGLE
                           L          = NOSOEL(5)
                           NOSOEL(5)  = NOSOEL(6)
                           NOSOEL(6)  = L
                           L          = NOSOEL( 7)
                           NOSOEL( 7) = NOSOEL(12)
                           NOSOEL(12) = L
                           L          = NOSOEL( 9)
                           NOSOEL( 9) = NOSOEL(10)
                           NOSOEL(10) = L
                           L          = NOSOEL(11)
                           NOSOEL(11) = NOSOEL( 8)
                           NOSOEL( 8) = L
                        ELSE
C                          TRIANGLE
                           L          = NOSOEL(5)
                           NOSOEL(5)  = NOSOEL(6)
                           NOSOEL(6)  = L
                           L          = NOSOEL( 7)
                           NOSOEL( 7) = NOSOEL(10)
                           NOSOEL(10) = L
                           L          = NOSOEL(8)
                           NOSOEL(8)  = NOSOEL(9)
                           NOSOEL(9)  = L
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
C
C           NOEUDSOMMET  ARETE  TRIANGLE  QUADRANGLE
C           TETRAEDRE  PENTAEDRE  HEXAEDRE  PYRAMIDE
C           TRI CROISSANT DES SOMMETS SANS MODIFIER LA REGLE DE
C           VOLUME POSITIF DANS LA PRESENTATION DES SOMMETS D'UN ELEMENT FINI
C           ATTENTION: LA PYRAMIDE EST UNE EXCEPTION CAR 1234 est un QUADRANGLE
C                      et LE SOMMET 5 APPARTIENT A 4 TRIANGLES
            CALL HACREN( NCOGEL, NBSO, NU1TG, NOSOEL )
C
C           RECHERCHE ET AJOUT EVENTUEL SELON LE HACHAGE
            CALL HACHAG( NBSO, NOSOEL, MOELEM(NCOGEL),
     %                   NBELEM(NCOGEL),
     %                   MCN(MNELEM(NCOGEL)), LIENEL(NCOGEL),
     %                   NUDREL(NCOGEL), NUELEM )
C           EN SORTIE: NUELEM <0 SI EF AJOUTE
C                             >0 SI EF RETROUVE
C
            MNN = MNELEM(NCOGEL) + MOELEM(NCOGEL)*ABS(NUELEM) - 2
            IF( NUELEM .LT. 0 ) THEN
C
C              EF AJOUTE => LES TG SONT AJOUTEES DANS ELTG(NCOGEL)
               IF( NUEFTG .GT. 0 ) THEN
C                 EF A TG : STOCKAGE DU NUMERO DES TGS DE L'EF
                  NT = MNELTG(NCOGEL)+NUELTG(NCOGEL)*MOELTG(NCOGEL)-1
                  NUELTG( NCOGEL ) = NUELTG( NCOGEL ) + 1
                  MCN( MNN + 1 ) = NUELTG( NCOGEL )
C                 LES NUMEROS DES TANGENTES ELLES MEMES
                  DO 16 J=1,MOELTG(NCOGEL)
                     MCN( NT + J ) = NOSOEL( NBSOEF + J )
 16               CONTINUE
               ENDIF
C
            ELSE
C
C              EF RETROUVE
               IF( MCN(MNN) .GT. 0 .AND. MCN(MNN) .NE. I ) THEN
C                 ERREUR  : MEME EF RETROUVE DANS 2 OBJETS PREMIERS DIFFERENTS
C                 REMARQUE: SI UN EF APPARAIT 2 FOIS DANS UN MEME OBJET
C                           ALORS LE SECOND EST IGNORE
                  GOTO ( 21, 22, 23, 23, 24, 24, 24, 24, 24 ), NCOGEL
 21               KNOM = 'POINTS'
                  GOTO 25
 22               KNOM = 'LIGNES'
                  GOTO 25
 23               KNOM = 'SURFACES'
                  GOTO 25
 24               KNOM = 'VOLUMES'
C
 25               CALL NMOBNU( KNOM, NUOBJT, KNOM1 )
C                 LE NUMERO DE L'OBJET AUPARAVANT
                  CALL NMOBNU( KNOM, NUOBPR(2,MCN(MNN)), KNOM2 )
C
                  NBLGRC(NRERR) = 6
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10025) (NOSOEL(J),J=1,NBSO)
                     WRITE(IMPRIM,10026) KNOM,KNOM1,KNOM2
                     WRITE(IMPRIM,10027)
     %              (NOSOEL(J),(XYZ(K,NOSOEL(J)),K=1,3),J=1,NBSO)
C
      KERR(1) ='DFTOP3: AU MOINS UN ELEMENT FINI APPARTENANT A 2 OBJETS'
                  KERR(2) ='DE TYPE ' // KNOM
                  KERR(3) = KNOM1 // ' ET ' // KNOM2
                  KERR(4) ='IMPOSSIBLE DE CHOISIR LEQUEL EST LE BON'
            KERR(5) ='CHANGER la PRECISION d''IDENTIFICATION des POINTS'
                  KERR(6) ='ou MODIFIER le MAILLAGE'
                  ELSE
                     WRITE(IMPRIM,20025) (NOSOEL(J),J=1,NBSO)
                     WRITE(IMPRIM,20026) KNOM,KNOM1,KNOM2
                     WRITE(IMPRIM,20027)
     %              (NOSOEL(J),(XYZ(K,NOSOEL(J)),K=1,3),J=1,NBSO)
C
             KERR(1) ='DFTOP3: At LEAST ONE FINITE ELEMENT IN 2 OBJECTS'
                  KERR(2) ='of TYPE ' // KNOM
                  KERR(3) = KNOM1 // ' and ' // KNOM2
                  KERR(4) = 'IMPOSSIBLE to CHOOSE WHICH is the GOOD'
                  KERR(5) ='CHANGE PRECISION to IDENTIFY the POINTS'
                  KERR(6) ='or MODIFY the MESH'
                  ENDIF
10025 FORMAT(' DFTOP3: L''ELEMENT FINI DE SOMMETS',8I7)
10026 FORMAT(' APPARTIENT AUX ',A,2(1X,A)/)
10027 FORMAT(' SOMMET',I7,' X=',G15.7,' Y=',G15.7,' Z=',G15.7)
C
20025 FORMAT(' DFTOP3: FINITE ELEMENT with VERTICES',8I7)
20026 FORMAT(' BELONG to ',A,2(1X,A)/)
20027 FORMAT(' VERTEX',I7,' X=',G15.7,' Y=',G15.7,' Z=',G15.7)
                  CALL LEREUR
                  NERR = 1
                  GOTO 9000
               ENDIF
C
               IF( NUEFTG .GT. 0 ) THEN
C                 EF A TG : STOCKAGE DU NUMERO DES TGS DE L'EF
                  NUELT = MCN( MNN + 1 )
                  IF( NUELT .EQ. 0 ) THEN
C                    LE PRECEDENT EF N'AVAIT PAS DE TG => IL EST CREE
                     NT = MNELTG(NCOGEL)+NUELTG(NCOGEL)*MOELTG(NCOGEL)-1
                     NUELTG( NCOGEL ) = NUELTG( NCOGEL ) + 1
                     MCN( MNN + 1 ) = NUELTG( NCOGEL )
C                    LES NUMEROS DES TANGENTES ELLES MEMES
                     DO 30 J=1,MOELTG(NCOGEL)
                        MCN( NT + J ) = NOSOEL( NBSOEF + J )
 30                  CONTINUE
                  ELSE
C                    LE PRECEDENT EF AVAIT DES TG => LES TG SONT REFONDUES
                     NT = MNELTG(NCOGEL) + (NUELT-1)*MOELTG(NCOGEL)-1
                     DO 35 J=1,MOELTG(NCOGEL)
C                       SI L'ANCIENNE TG N'EXISTE PAS ALORS LA NOUVELLE PREND LA
                        IF( MCN(NT+J).EQ.0 ) MCN(NT+J)=NOSOEL(NBSOEF+J)
C                       SI LA TG ANCIENNE EXISTAIT ALORS ELLE EST CONSERVEE
 35                  CONTINUE
                  ENDIF
               ENDIF
            ENDIF
C
C           STOCKAGE DU DERNIER OBJET CONTENANT L'EF POUR AFFICHER TOUS LES CAS
C           STOCKAGE DU NO OBJET PREMIER DE L'EF
            MCN( MNN ) = I
 900     CONTINUE
 1000 CONTINUE
C
C     CODE D'ERREUR
 9000 IERR = MAX( IERR, NERR )
C
      DO 9010 I=1,9
         IF( NUELTG(I) .GT. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
             WRITE(IMPRIM,*) NUELTG(I),' EF de type',I,' avec TANGENTES'
             WRITE(IMPRIM,*) NBELTG(I),' Elements Finis a TANGENTES'
            ELSE
              WRITE(IMPRIM,*) NUELTG(I),' EF of type',I,' with TANGENTS'
              WRITE(IMPRIM,*) NBELTG(I),' Finite Element with TANGENTS'
            ENDIF
         ENDIF
 9010 CONTINUE
C
      RETURN
      END
