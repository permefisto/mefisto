      SUBROUTINE TRVMTR( MISTRE, NDIM,   KNOMOB, NTLXOB, MODECO,
     %                   NBTYEL, MNELEM, MNPOGE, NDPGST, NCAS,   IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LE CRITERE DE VON MISES OU TRESCA DE CET OBJET
C ----- Seuil de PLASTICITE = LIMITE de l'ELASTICITE du MATERIAU
C       En 3d:    si la i-eme contrainte principale
C       Von MISES = sqrt( (s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 ) / sqrt(2)
C       TRESCA    = MAX( abs(s1-s2), abs(s2-s3), abs(s3-s1) ) / 2
C       En 2d:
C       Von MISES = TRESCA = Abs( s1 - s2 )
C
C ENTREES:
C --------
C MISTRE : 1 POUR CRITERE DE VON MISES
C          2 POUR CRITERE DE TRESCA
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET A TRAITER
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS FINIS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NCAS   : NUMERO DU CAS A TRAITER
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR RENCONTREE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      PARAMETER     (LIGCON=0)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      REAL              CONTMN, CONTMX
      INTEGER           NBEF(4), NBST(4), MNCRIT(4), NOINTC(4)
      include"./incl/nomele.inc"
C
      CONTMN = 0
      CONTMX = 0
      MNCOMT = 0
      DO 5 I=1,4
         MNCRIT(I) = 0
 5    CONTINUE
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'EF DU MAILLAGE
C     ================================================
      DO 30 I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNELEM + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE DE TELS ELEMENTS FINIS
         NBELEM = MCN(MNELE + WBELEM )
         NBEF(I+1) = NBELEM
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        CALL ELNUNM( NUTYEL, NOMELE )
         KNOM = 'CONTRAINTE"' // NOMELE(2,NUTYEL)
C        OUVERTURE DU TABLEAU CONTRAINTE
         CALL LXTSOU( NTLXOB, KNOM, NTCONT, MNCONT )
         IF( NTCONT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'CONTRAINTES NON CALCULEES'
            ELSE
               KERR(1) = 'STRESSES NOT COMPUTED'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9000
         ENDIF
C
C        ADRESSE MCN DU TABLEAU VALEURS DES CONTRAINTES AUX POINTS
C        D'INTEGRATION
         LDCOPR = MCN(MNCONT+WDCOPR)
C
C        OUVERTURE DES TABLEAUX CONTENANT LES NUMEROS ET LES COORDONNEES
C        DES POINTS D'INTEGRATION UTILISES POUR L'EXTRAPOLATION
         CALL TNMCDC( 'ENTIER', 30,  MNNOPI )
         CALL TNMCDC( 'REEL2',  3*8, MNCOEX )
C
C        SUIVANT LE TYPE D'EF RECHERCHE DU NUMERO D'INTERPOLATION
C        AINSI QUE DES POINTS D'INTEGRATION UTILISES POUR L'EXTRAPOLATION
C        DU CRITERE AU NIVEAU DES NOEUDS
         CALL PIVMTR( NUTYEL,
     %                NOINTC(I+1), NBPIEX, MCN(MNNOPI), MCN(MNCOEX))
         NBST(I+1) = NBPIEX
C
C        OUVERTURE DU TABLEAU CRITERE DE VON MISES OU TRESCA
         MNCOMT = 0
         MOCOMT = NBPIEX * NBELEM
         CALL TNMCDC( 'REEL2', MOCOMT, MNCOMT )
         MNCRIT(I+1) = MNCOMT
C
C        OUVERTURE DU TABLEAU UTILISE POUR LA RESOLUTION DU
C        SYSTEME LINEAIRE QUI REALISE L'EXTRAPOLATION AUX SOMMETS
         CALL TNMCDC( 'REEL2', NBPIEX*(NBPIEX+NBELEM), MNSYST )
C
C        CALCUL DU CRITERE DE VON MISES OU DE TRESCA AUX SOMMETS DES EF
C        POUR LE VECTEUR DEPLACEMENT NCAS
         NBELFI = MCN(MNCONT+WBELFI)
         CALL NOVMTR( MISTRE, NBELFI, MCN(MNCONT+WBPIEF),
     %                MCN(MNCONT+WDIMES), MCN(MNCONT+WBJECA),
     %                MCN(MNCONT+LDCOPR), NCAS,
     %                NOINTC(I+1), NBPIEX, MCN(MNNOPI), MCN(MNCOEX),
     %                MCN(MNCOMT), MCN(MNSYST), CONTMN, CONTMX )
C        MCN(MNCOMT)=CRITERE(NBPIEX,NBELFI) A TRACER
C
C        DESTRUCTION DES TABLEAUX MC DEVENUS INUTILES
         IF( MNNOPI .GT. 0 ) CALL TNMCDS( 'ENTIER',  30, MNNOPI )
         IF( MNCOEX .GT. 0 ) CALL TNMCDS( 'REEL2',  3*8, MNCOEX )
         IF( MNSYST .GT. 0 ) CALL TNMCDS( 'REEL2',
     %                            NBPIEX*(NBPIEX+NBELEM), MNSYST )
C
 30   CONTINUE
C
C     TRACE DES TABLEAUX CRITERE(NBPIEX,NBELFI)
C     LECTURE DES DONNEES DU TRACE
 100  CALL LIMTCL( 'valzone', NMTCL )
      IF( NMTCL .LT.  0 ) GOTO 9000
      IF( NMTCL .EQ. 90 ) GOTO  200
      GOTO( 200, 100, 100, 100, 150, 160 ) , NMTCL
C
C     COULEUR des ARETES du MAILLAGE
C     ..............................
 150  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
C        COULEUR RESERVEE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 100
C
C     TYPE du TRAIT des ARETES du MAILLAGE
C     ....................................
 160  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAFR = I
      GOTO 100
C
C     LE TRACE EFFECTIF DU CRITERE
C     ============================
 200  IF( NDIM .EQ. 2 ) THEN
C        TRACE 2D
         CALL TRVMTR21( MISTRE, KNOMOB, NBTYEL, MNELEM, NDPGST,
     %                  MODECO, NCAS,   CONTMN, CONTMX,
     %                  NBST,   MNCRIT,
     %                  MCN(MNPOGE+WNBPOI), MCN(MNPOGE+WYZPOI) )
      ELSE
         CALL TRVMTR31( MISTRE, KNOMOB, NBTYEL, MNELEM,
     %                  MODECO, NCAS,   CONTMN, CONTMX,
     %                  NBST,   MNCRIT, NOINTC,
     %                  NDPGST, MNPOGE, MCN(MNPOGE+WNBPOI) )
      ENDIF
      GOTO 100
C
C     DESTRUCTION DES TABLEAUX DU CRITERE SUR LES DIFFERENTS TYPES D'EF
 9000 DO 9010 I=1, NBTYEL
         IF( MNCRIT(I) .GT. 0 )
     %       CALL TNMCDS( 'REEL2', NBST(I)*NBEF(I), MNCRIT(I) )
 9010 CONTINUE
C
      RETURN
      END
