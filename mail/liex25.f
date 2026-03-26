      SUBROUTINE LIEX25( NTLXLI, LADEFI,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LES ARETES ET VECTEURS TANGENTS A CES ARETES D'UN     +
C ----- POLYGONE REGULIER INSCRIT DANS L'INTERSECTION D'UN ELLIPSOIDE +
C       ET D'UN PLAN (LIGNE DE TYPE 25)                               +
C                                                                     +
C ENTREES :                                                           +
C ---------                                                           +
C NTLXLI : NUMERO DU TABLEAU TMS DU LEXIQUE DE LA LIGNE               +
C LADEFI : TABLEAU ENTIER DES DONNEES DE LA DEFINITION DE LA LIGNE    +
C          CF ~/td/d/a_ligne__definition                              +
C                                                                     +
C SORTIES :                                                           +
C ---------                                                           +
C NTNSEF: NUMERO DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE     +
C MNNSEF: ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE+
C         CF ~/td/d/a___nsef                                          +
C NTXYZS: NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE                  +
C MNXYZS: ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE                  +
C IERR  : 0 SI PAS D'ERREUR                                           +
C         1 SI NOMBRE D'ARETES DESIRE < 4                             +
C         2 SI LES TROIS POINTS DONNES NE SUFFISENT PAS A DEFINIR UN  +
C           PLAN                                                      +
C         3 SI LES AXES DE L'ELLIPSOIDE DONNES NE SONT PAS            +
C           PERPENDICULAIRES                                          +
C         4 SI L'INTERSECTION EST VIDE OU REDUITE A UN POINT          +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS       : JOSE URQUIZA  & ARNAUD MONTENAY  DECEMBRE 1994      +
C MODIFICATIONS : RODOLFO ARAYA & PASCAL MAILLOT   JANVIER  1996      +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "./incl/langue.inc"
      include "./incl/gsmenu.inc"
      include "./incl/a_ligne__definition.inc"
      include "./incl/a___xyzsommet.inc"
      include "./incl/a___nsef.inc"
      include "./incl/ntmnlt.inc"
C
C     LE SUPER-TABLEAU OU TOUS LES TMC ET TMS OUVERTS SONT STOCKES
C     ICI LE TABLEAU DMCN REEL DOUBLE PRECISION N'EST PAS UTILE
      include"./incl/pp.inc"
      COMMON        MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE ( MCN(1), RMCN(1) )
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
C     LE TABLEAU DES DONNEES DE LA DEFINITION DE LA LIGNE
C     VU SOUS FORME D'UN TABLEAU D'ENTIERS
      INTEGER LADEFI(0:*)
C
C     VARIABLES LOCALES AUXILIAIRES
C
C     XYZPL : MATRICE DES COORDONNEES DES 3 POINTS DU PLAN
C     X0    : COORDONNEES DU CENTRE DE L'ELLIPSOIDE
C     XYZEL : MATRICE DES QUATRE POINTS SUR LES AXES DE L'ELLIPSOIDE
C     R     : MATRICE DES COORDONNEES DES TROIS VECTEURS ORTHOGONAUX
C             DE L'ELLIPSOIDE
C     POINTS: TABLEAU (3*NOMBRE DE POINTS) DES COORDONNEES DE
C             L'INTERSECION
C
      REAL  X0(3), R(3,3)
      REAL  XYZPL(3,3), XYZEL(3,4)
C
C     RECUPERATION DU NOMBRE D'ARETES DE LA LIGNE
      NBARLI = LADEFI(WBARLI)
C
C     INITIALISATION DES ADRESSES DU SUPER-TABLEAU
      MNTHET = 0
      MNPOTS = 0
      MNTGS  = 0
C
C     RECUPERATION DU NOM DES 3 POINTS DE DEFINITION DU PLAN
      DO 2 I=1,3
C
C        LE NUMERO DU POINT DANS LE LEXIQUE DES POINTS
         NUPT = LADEFI( WUPTPL - 1 + I )
C
 1       IF( NUPT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'UN POINT de DEFINITION du PLAN est INCONNU'
            ELSE
               KERR(1) = 'AN UNKNOWN POINT to DEFINE the PLANE'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
C
C        RECUPERATION DES 3 COORDONNEES DES POINTS DU PLANL
         CALL LXNLOU( NTPOIN, NUPT, NTPOI, MN )
         CALL LXTSOU( NTPOI , 'XYZSOMMET', NTSOM, MNSOM )
         IF( NTSOM .LE. 0 ) THEN
            NUPT = 0
            GOTO 1
         ENDIF
         MN = MNSOM + WYZSOM
         XYZPL( 1 , I ) = RMCN( MN )
         XYZPL( 2 , I ) = RMCN( MN + 1 )
         XYZPL( 3 , I ) = RMCN( MN + 2 )
 2    ENDDO
C
C     RECUPERATION DU NOM DES 4 POINTS DE DEFINITION DE L'ELLIPSOIDE
      DO 4 I=1,4
C
C        LE NUMERO DU POINT DANS LE LEXIQUE DES POINTS
         NUPT = LADEFI( WUPTEL - 1 + I )
C
 3       IF( NUPT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'UN POINT de DEFINITION ELLPISOIDE est INCONNU'
            ELSE
               KERR(1) = 'AN UNKNOWN POINT to DEFINE the ELLIPSOID'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9999
         ENDIF
C
C        RECUPERATION DES 3 COORDONNEES DES POINTS DE L'ELLIPSOIDE
         CALL LXNLOU( NTPOIN, NUPT, NTPOI, MN )
         CALL LXTSOU( NTPOI , 'XYZSOMMET',  NTSOM, MNSOM )
         IF( NTSOM .LE. 0 ) THEN
            NUPT = 0
            GOTO 3
         ENDIF
         MN = MNSOM + WYZSOM
         XYZEL( 1 , I ) = RMCN( MN )
         XYZEL( 2 , I ) = RMCN( MN + 1 )
         XYZEL( 3 , I ) = RMCN( MN + 2 )
 4    ENDDO
C
C     RECUPERATION DES COORDONNEES DU CENTRE DE L'ELLIPSOIDE
      DO 15 I=1,3
          X0(I)=XYZEL(I,1)
 15   ENDDO
C
C     GENERATION DES VECTEURS D'AXE
      DO 40 J=1,3
          DO 35 I=1,3
              R(I,J)=XYZEL(I,J+1)-X0(I)
 35       ENDDO
 40   ENDDO
C
C     DECLARATION DES 3 TABLEAUX POINTS, TANGENTES ET THETA
      CALL TNMCDC( 'REEL', NBARLI * 3, MNPOTS )
      CALL TNMCDC( 'REEL', NBARLI * 3, MNTGS  )
      CALL TNMCDC( 'REEL', NBARLI + 1, MNTHET )
C
C     GENERATION DES POINTS DE L'INTERSECTION ET VECTEURS TANGENTS
C     A CES POINTS.
      CALL GENERE( XYZPL, X0, R, NBARLI,
     %             RMCN(MNPOTS), RMCN(MNTGS), RMCN(MNTHET), IERR )
C
C     VERIFICATION DU BON FONCTIONNEMENT
      IF(IERR .NE. 0) THEN
          GOTO 9999
      ENDIF
C
C     NOMBRE DE TANGENTES CALCULEES
      NBTGS = NBARLI
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE LA LIGNE
C     -----------------------------------------------
      CALL LXTNDC(NTLXLI, 'XYZSOMMET', 'ENTIER',
     %            WYZSOM+3*NBARLI+3*NBTGS)
      CALL LXTSOU(NTLXLI, 'XYZSOMMET', NTXYZS, MNXYZS)
C
C     NOMBRE DE SOMMETS DE LA LIGNE
      MCN( MNXYZS + WNBSOM )=NBARLI
C
C     NOMBRE DE TANGENTES DE LA LIGNE
      MCN( MNXYZS + WNBTGS )=NBTGS
      MCN( MNXYZS + WBCOOR )= 3
C
C     ADRESSE DU DEBUT DES COORDONNEES DU PREMIER POINT DE LA LIGNE
      MNS=MNXYZS+WYZSOM
C
C     ADRESSE DU DEBUT DES COORDONNEES DU PREMIER VECTEUR TANGENT DE LA LIGNE
      MNT=MNS+3*NBARLI
C
C     GENERATION DES COORDONNEES DES NBARLI POINTS ET NBTGS VECTEURS TANGENTS
      DO 75 I=1,NBARLI
          DO 72 J=1,3
              RMCN(MNS+(I-1)*3+J-1)=RMCN(MNPOTS-1+J+3*(I-1))
              RMCN(MNT+(I-1)*3+J-1)=RMCN(MNTGS -1+J+3*(I-1))
 72       ENDDO
 75   ENDDO
C
C     AJOUT DE LA DATE
      CALL ECDATE(MCN(MNXYZS))
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN(MNXYZS+MOTVAR(6))=NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' DE TYPE LIGNE STRUCTUREE
C     -------------------------------------------------------
      CALL LXTNDC(NTLXLI,'NSEF','ENTIER', WUSOEF + 6 * NBARLI )
      CALL LXTSOU(NTLXLI,'NSEF',NTNSEF,MNNSEF)
C
C     TYPE DE L'OBJET : ICI LIGNE
      MCN( MNNSEF + WUTYOB ) = 2
C
C     LIGNE FERMEE
      MCN( MNNSEF + WUTFMA ) = 1
C
C     NOMBRE DE SOMMETS PAR ARETE
      MCN( MNNSEF + WBSOEF ) = 2
C
C     NOMBRE DE TANGENTES PAR ARETE
      MCN( MNNSEF + WBTGEF ) = 2
C
C     NOMBRE D'ARETES DU SEGMENT STRUCTURE
      MCN( MNNSEF + WBEFOB ) = NBARLI
C
C     NOMBRE DE EF A POINTEUR ET TG
      MCN( MNNSEF + WBEFTG ) = NBARLI
      MCN( MNNSEF + WBEFAP ) = NBARLI
C
C     TYPE DU MAILLAGE : ICI SEGMENT NON STRUCTURE
      MCN( MNNSEF + WUTYMA ) = 0
C
C     NUMERO DES SOMMETS ET DES TANGENTES DES ARETES
      MN  = MNNSEF + WUSOEF
      MNP = MN  + 2 * NBARLI
      MNG = MNP + NBARLI
      MNT = MNG + NBARLI
      DO 200 I=1,NBARLI
C
C        NUMERO DES 2 SOMMETS DE L'ARETE I
         MCN( MN   ) = I
         MCN( MN+1 ) = I+1
         MN = MN + 2
C
C        LE POINTEUR
         MCN( MNP ) = I
         MNP = MNP + 1
C
C        LE CODE GEOMETRIQUE : ELLIPSE => 2
         MCN( MNG ) = 2
         MNG = MNG + 1
C
C        NUMERO DES 2 TANGENTES DE L'ARETE I
         MCN( MNT   ) =   I
         MCN( MNT+1 ) = -(I+1)
         MNT = MNT + 2
 200  ENDDO
C
C     LA LIGNE EST FERMEE
      MCN( MN -1 ) = 1
      MCN( MNT-1 ) =-1
C
C     AJOUT DE LA DATE
      CALL ECDATE(MCN(MNNSEF))
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF')
C
C     ERREUR
C     ======
C     DESTRUCTION DES 3 TABLEAUX POINTS TANGENTES THETA
 9999 IF( MNPOTS .GT. 0 ) CALL TNMCDS( 'REEL', NBARLI * 3, MNPOTS )
      IF( MNTGS  .GT. 0 ) CALL TNMCDS( 'REEL', NBARLI * 3, MNTGS  )
      IF( MNTHET .GT. 0 ) CALL TNMCDS( 'REEL', NBARLI + 1, MNTHET )
      END
C
C
      SUBROUTINE GENERE( XABC, X0, R, NBPOINTS,
     %                   POINTS, TANGENTES ,THETA, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LES ARETES ET LES VECTEURS TANGENTS A CES ARETES D'UN +
C ----- POLYGONE REGULIER INSCRIT DANS L'INTERSECTION D'UN ELLIPSOIDE +
C       ET D'UN PLAN (LIGNE DE TYPE 25)                               +
C                                                                     +
C ENTREES :                                                           +
C ---------                                                           +
C XABC     : TABLEAU DES TROIS POINTS DEFINISSANT LE PLAN             +
C X0       : COORDONNEES DU CENTRE DE L'ELLIPSOIDE                    +
C R        : TABLEAU REEL DES AXES DE L'ELLIPSOIDE                    +
C NBPOINTS : NOMBRE DE POINTS SUR L'INTERSECTION                      +
C THETA    : PARAMETRES DES POINTS SUR L'INTERSECTION                 +
C                                                                     +
C SORTIES :                                                           +
C ---------                                                           +
C POINTS   : TABLEAU DES COORDONNEES DES POINTS DE L'INTERSECTION     +
C TANGENTES: TABLEAU DES COORDONNEES DES VECTEURS TANGENTS            +
C IERR     : 0 SI PAS D'ERREUR                                        +
C            1 SI NOMBRE D'ARETES DESIRE < 4                          +
C            2 SI LES TROIS POINTS DONNES NE SUFFISENT PAS A DEFINIR  +
C              UN PLAN                                                +
C            3 SI LES AXES DE L'ELLIPSOIDES DONNES NE SONT PAS        +
C              PERPENDICULAIRES                                       +
C            4 SI L'INTERSECTION EST VIDE OU REDUITE A UN POINT       +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS       : JOSE URQUIZA & ARNAUD MONTENAY : DECEMBRE 1994      +
C MODIFICATIONS : RODOLFO ARAYA & PASCAL MAILLOT : JANVIER 1996       +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "./incl/langue.inc"
      include "./incl/gsmenu.inc"
C
      REAL XABC(3,3),X0(3),R(3,3),POINTS(3,NBPOINTS),THETA(0:NBPOINTS)
      REAL TANGENTES(3,NBPOINTS)
      INTEGER NBPOINTS, IERR
C
C     VARIABLES LOCALES AUXILIAIRES
C
C     R            : MATRICE DES VECTEURS UNITAIRES SUIVANT LES AXES DE L'ELLIPS
C     A, B, C      : DEMI-LONGUEURS DES AXES DE L'ELLIPSOIDE
C     R2           : COORDONNEES D'UN TRIEDRE DIRECT LIE AU PLAN
C     QUAD         : MATRICE DE LA FORME QUADRATIQUE ASSOSCIEE A L'ELLIPSOIDE
C                    DANS LE REPERE DE CENTRE LE CENTRE DE L'ELLIPSOIDE ET DE
C                    VECTEURS DE BASE, CEUX FORMES PAR R2
C     GAM0         : DISTANCE ALGEBRIQUE DU PLAN AU CENTRE DE L'ELLIPSOIDE
C                    SUIVANT LA NORMALE AU PLAN
C     CENTRE       : VECTEUR DES COORDONNEES DU CENTRE DE L'ELLIPSE DANS LE
C                    REPERE FORME PAR LE CENTRE DE L'ELLIPSOIDE ET R2
C     LAM1, LAM2   : DEMI-LONGUEURS DES AXES DE L'ELLIPSE
C     VECT1, VECT2 : VECTEURS DIRECTEURS NORMES DES AXES
C
      REAL A,B,C,R2(3,3),DET,L1,L2,L3
      REAL QUAD(3,3),GAM0,CH(3,3),CH1(3,3),CENTRE(3),VECT1(3),VECT2(3)
      REAL LAM1,LAM2,EPS1,EPS2,EPS3,MINAXE,MINEL
C
C     VALEURS SEUILS DES 3 ERREURS
C     EPS1 : LIMITE SUR LA DEGENERESCENCE DU TRIEDRE OABC
      EPS1=0.01
C     EPS2 : LIMITE SUR L'ECART A L'ORTHOGONALITE DES TROIS AXES.
      EPS2=0.01
C     EPS3 : TOLERANCE SUR L'EQUIDISTANCE ENTRE LES POINTS
      EPS3=0.001
C
C     TAILLE DE L'INTERVALLE DU PARAMETRE POUR UNE ARETE
      H = REAL( ATAN(1.D0)*8.D0 / NBPOINTS )
C
C     VERIFICATION DE LA COHERENCE DES DONNEES
      L1=SQRT((XABC(1,2)-XABC(1,1))**2+(XABC(2,2)-XABC(2,1))**2
     %+(XABC(3,2)-XABC(3,1))**2)
      L2=SQRT((XABC(1,3)-XABC(1,1))**2+(XABC(2,3)-XABC(2,1))**2
     %+(XABC(3,3)-XABC(3,1))**2)
      R2(1,3)=(XABC(2,2)-XABC(2,1))*(XABC(3,3)-XABC(3,1))
     %-(XABC(3,2)-XABC(3,1))*(XABC(2,3)-XABC(2,1))
      R2(2,3)=(XABC(3,2)-XABC(3,1))*(XABC(1,3)-XABC(1,1))
     %-(XABC(1,2)-XABC(1,1))*(XABC(3,3)-XABC(3,1))
      R2(3,3)=(XABC(1,2)-XABC(1,1))*(XABC(2,3)-XABC(2,1))
     %-(XABC(2,2)-XABC(2,1))*(XABC(1,3)-XABC(1,1))
      L3=SQRT(R2(1,3)**2+R2(2,3)**2+R2(3,3)**2)
      A=SQRT(R(1,1)**2+R(2,1)**2+R(3,1)**2)
      B=SQRT(R(1,2)**2+R(2,2)**2+R(3,2)**2)
      C=SQRT(R(1,3)**2+R(2,3)**2+R(3,3)**2)
      IF(NBPOINTS .LT. 4) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE INSUFFISANT d''ARETES de la LIGNE'
         ELSE
            KERR(1) = 'AUGMENT the EDGE NUMBER of the LINE'
         ENDIF
         CALL LEREUR
         IERR=1
         GOTO 9999
      ELSE IF(L3.LT.(EPS1*L1*L2)) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = '3 POINTS ALIGNES pour DEFINIR le PLAN'
         ELSE
            KERR(1)='3 POINTS ON A SAME STRAIGHT LINE to DEFINE a PLANE'
         ENDIF
         CALL LEREUR
         IERR=2
         GOTO 9999
      ELSE IF(ABS(R(1,1)*R(1,2)+R(2,1)*R(2,2)+R(3,1)*R(3,2)).GT.(EPS2*A
     %*B))THEN
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AXES 1 ET 2 de l''ELLIPSOIDE NON ORTHOGONAUX'
         ELSE
            KERR(1)='AXES 1 and 2 of the ELLIPSOIDE ARE NOT ORTHOGONAL'
         ENDIF
         CALL LEREUR
         IERR=3
         GOTO 9999
      ELSE IF(ABS(R(1,1)*R(1,3)+R(2,1)*R(2,3)+R(3,1)*R(3,3)).GT.(EPS2*A
     %*C))THEN
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AXES 1 ET 3 de l''ELLIPSOIDE NON ORTHOGONAUX'
         ELSE
            KERR(1)='AXES 1 and 3 of the ELLIPSOIDE ARE NOT ORTHOGONAL'
         ENDIF
         CALL LEREUR
         IERR=3
         GOTO 9999
      ELSE IF(ABS(R(1,3)*R(1,2)+R(2,3)*R(2,2)+R(3,3)*R(3,2)).GT.(EPS2*C
     %*B))THEN
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AXES 2 ET 3 de l''ELLIPSOIDE NON ORTHOGONAUX'
         ELSE
            KERR(1)='AXES 2 and 3 of the ELLIPSOIDE ARE NOT ORTHOGONAL'
         ENDIF
         CALL LEREUR
         IERR=3
         GOTO 9999
      ENDIF
C
C     CONSTRUCTION DU REPERE ORTHONORME DIRECT LIE AU PLAN I.E. LA MATRICE R2
C     1ER VECTEUR
      R2(1,1)=(XABC(1,2)-XABC(1,1))/L1
      R2(2,1)=(XABC(2,2)-XABC(2,1))/L1
      R2(3,1)=(XABC(3,2)-XABC(3,1))/L1
C     3EME VECTEUR
      DO 32 I=1,3
          R2(I,3)=R2(I,3)/L3
 32   ENDDO
C     2EME VECTEUR
      R2(1,2)=R2(2,3)*R2(3,1)-R2(3,3)*R2(2,1)
      R2(2,2)=R2(3,3)*R2(1,1)-R2(1,3)*R2(3,1)
      R2(3,2)=R2(1,3)*R2(2,1)-R2(2,3)*R2(1,1)
C
C     REORTHOGONALISATION DES TROIS AXES ET NORMALISATION DES
C     VECTEURS DIRECTEURS
      DO 35 I=1,3
          R(I,1)=R(I,1)/A
 35   ENDDO
      R(1,2)=R(2,3)*R(3,1)-R(2,1)*R(3,3)
      R(2,2)=R(3,3)*R(1,1)-R(3,1)*R(1,3)
      R(3,2)=R(1,3)*R(2,1)-R(2,3)*R(1,1)
      L2=SQRT(R(1,2)**2+R(2,2)**2+R(3,2)**2)
      DO 40 I=1,3
          R(I,2)=R(I,2)/L2
 40   ENDDO
      R(1,3)=R(2,1)*R(3,2)-R(2,2)*R(3,1)
      R(2,3)=R(3,1)*R(1,2)-R(3,2)*R(1,1)
      R(3,3)=R(1,1)*R(2,2)-R(2,1)*R(1,2)
C
C     CALCUL DE GAM0
      GAM0=(R2(1,3)*(XABC(1,1)-X0(1))+R2(2,3)*(XABC(2,1)-X0(2))
     %+R2(3,3)*(XABC(3,1)-X0(3)))
C
C     INITIALISATION DE LA FORME QUADRATIQUE QUAD
      DO 50 J=1,3
          DO 49 I=1,3
              CH(I,J)=0.0
              DO 48 K=1,3
                  CH(I,J)=CH(I,J)+R(K,I)*R2(K,J)
 48           ENDDO
 49       ENDDO
 50   ENDDO
      DO 55 J=1,3
          CH1(1,J)=CH(1,J)/(A**2)
          CH1(2,J)=CH(2,J)/(B**2)
          CH1(3,J)=CH(3,J)/(C**2)
 55   ENDDO
      DO 70 J=1,3
          DO 65 I=1,3
              QUAD(I,J)=0.0
              DO 60 K=1,3
                  QUAD(I,J)=QUAD(I,J)+CH(K,I)*CH1(K,J)
 60           ENDDO
 65       ENDDO
 70   ENDDO
C
C     CONSTRUCTION DES COORDONNEES DU CENTRE DE L'ELLIPSEE
      CENTRE(1)=GAM0*(-QUAD(2,2)*QUAD(1,3)+QUAD(1,2)*QUAD(2,3))/
     %(QUAD(1,1)*QUAD(2,2)-QUAD(1,2)**2)
      CENTRE(2)=GAM0*(QUAD(1,2)*QUAD(1,3)-QUAD(1,1)*QUAD(2,3))/
     %(QUAD(1,1)*QUAD(2,2)-QUAD(1,2)**2)
      CENTRE(3)=GAM0
C
C     PREMIERE VERIFICATION DE L'EXISTENCE DE L'INTERSECTION
      DET=1.0-QUAD(3,3)*GAM0**2+QUAD(1,1)*CENTRE(1)**2+2.0*QUAD(1,2)
     %*CENTRE(1)*CENTRE(2)+QUAD(2,2)*CENTRE(2)**2
      IF(DET .LE. 0.0) THEN
          IERR=4
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1)='INTERSECTION VIDE OU REDUITE A UN POINT'
          ELSE
             KERR(1)='EMPTY INTERSECTION or REDUCED to ONE POINT'
          ENDIF
          CALL LEREUR
          GOTO 9999
      ENDIF
C
C     CALCUL DES DEMI LONGUEURS D'AXE DE L'ELLIPSE ET DES VECTEURS
C     DIRECTEURS DES AXES DANS LE REPERE DE R2
      IF(QUAD(1,2) .EQ. 0.0) THEN
          LAM1=QUAD(1,1)/DET
          LAM2=QUAD(2,2)/DET
          VECT1(1)=1.0
          VECT1(2)=0.0
          VECT2(1)=0.0
          VECT2(2)=1.0
      ELSE
          LAM1=(QUAD(1,1)+QUAD(2,2)
     %    +SQRT((QUAD(1,1)-QUAD(2,2))**2+4*QUAD(1,2)**2))/2.0
          LAM2=(QUAD(1,1)+QUAD(2,2)
     %    -SQRT((QUAD(1,1)-QUAD(2,2))**2+4*QUAD(1,2)**2))/2.0
          VECT1(1)=1.0
          VECT1(2)=(LAM1-QUAD(1,1))/QUAD(1,2)
          LAM1=LAM1/DET
          LAM2=LAM2/DET
          DET=SQRT(VECT1(1)**2+VECT1(2)**2)
          VECT1(1)=VECT1(1)/DET
          VECT1(2)=VECT1(2)/DET
          VECT2(1)=-VECT1(2)
          VECT2(2)=VECT1(1)
      ENDIF
C
C
C     DEUXIEME VERIFICATION DE L'EXISTENCE DE L'INTERSECTION
      MINAXE=A
      IF(MINAXE.GT.B) MINAXE=B
      IF(MINAXE.GT.C) MINAXE=C
      MINEL=1.0/SQRT(LAM1)
      IF(MINEL.GT.(1.0/SQRT(LAM2))) MINEL=1.0/SQRT(LAM2)
      IF(MINEL.LE.0.01*MINAXE) THEN
          IERR=4
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1)='INTERSECTION VIDE OU REDUITE A UN POINT'
          ELSE
             KERR(1)='EMPTY INTERSECTION or REDUCED to ONE POINT'
          ENDIF
          CALL LEREUR
          GOTO 9999
      ENDIF
C
C     REECRITURE DES VECTEURS VECT1 ET VECT2 ET DES COORDOONEES DU
C     CENTRE DE L'ELLIPSE DANS LE REPERE D'ORIGINE
      L1=VECT1(1)
      L2=VECT1(2)
      VECT1(1)=R2(1,1)*L1+R2(1,2)*L2
      VECT1(2)=R2(2,1)*L1+R2(2,2)*L2
      VECT1(3)=R2(3,1)*L1+R2(3,2)*L2
      L1=VECT2(1)
      L2=VECT2(2)
      VECT2(1)=R2(1,1)*L1+R2(1,2)*L2
      VECT2(2)=R2(2,1)*L1+R2(2,2)*L2
      VECT2(3)=R2(3,1)*L1+R2(3,2)*L2
      L1=CENTRE(1)
      L2=CENTRE(2)
      L3=CENTRE(3)
      CENTRE(1)=X0(1)+R2(1,1)*L1+R2(1,2)*L2+R2(1,3)*L3
      CENTRE(2)=X0(2)+R2(2,1)*L1+R2(2,2)*L2+R2(2,3)*L3
      CENTRE(3)=X0(3)+R2(3,1)*L1+R2(3,2)*L2+R2(3,3)*L3
C
C     CALCUL DES PARAMETRES DU POLYGONE A NBARLI SOMMETS SUR L'ELLIPSE
      CALL ELLIPSE(1.0/SQRT(LAM1),1.0/SQRT(LAM2),NBPOINTS,EPS3,THETA)
C
C     CONSTRUCTION DES TABLEAUX POINTS ET TANGENTES DE LA LIGNE
C     ---------------------------------------------------------
C     GENERATION DES COORDONNEES DES NBPOINTS POINTS ET TANGENTES
      DO 75 I=1,NBPOINTS
          DO 72 J=1,3
              POINTS(J,I)=CENTRE(J)+VECT1(J)*COS(THETA(I))/SQRT(LAM1)+
     %        VECT2(J)*SIN(THETA(I))/SQRT(LAM2)
              TANGENTES(J,I)= -VECT1(J)*SIN(THETA(I))/SQRT(LAM1)+
     %        VECT2(J)*COS(THETA(I))/SQRT(LAM2)
              TANGENTES(J,I) = TANGENTES(J,I) * H
 72       ENDDO
 75   ENDDO
C     PAS D'ERREUR RENCONTREE
      IERR=0
C
C     ERREUR
C     ======
 9999 RETURN
      END
C
      SUBROUTINE ELLIPSE( A,B,NBPOINTS,EPSILON,THETAN )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LES PARAMETRES DE POINTS REPARTIS DE FACON            +
C ----- EQUIDISTANTE SUR UNE ELLIPSE DE DEMI AXES A ET B              +
C                                                                     +
C ENTREES :                                                           +
C ---------                                                           +
C A, B    : DEMI-AXES DE L'ELLIPSE                                    +
C NBPOINTS: NOMBRE DE COTES DU POLYGONE                               +
C EPSILON : TOLERANCE SUR L'EQUIDISTANCE ENTRE LES POINTS             +
C                                                                     +
C SORTIES :                                                           +
C ---------                                                           +
C THETAN  : TABLEAU DES PARAMETRES DES POINTS RESULTAT                +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS       : JOSE URQUIZA & ARNAUD MONTENAY : DECEMBRE 1994      +
C MODIFICATIONS : RODOLFO ARAYA & PASCAL MAILLOT : JANVIER 1996       +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     PARAMETRES DE L'ELLIPSE
      REAL A,B
C
C     NOMBRE DE POINTS VOULUS
      INTEGER NBPOINTS
C
C     PARAMETRES D'ECART PAR RAPPORT A L'EGALITE DE LONGUEUR
C     ENTRE LES SEGMENTS
      REAL EPSILON
C
C     PARAMETRES DES POINTS SUR L'ELLIPSE
      REAL THETAN (0:NBPOINTS)
C
C     VARIABLES LOCALES
C
C     TAMPONS POUR LES PARAMETRES
      REAL THETAB (0:100)
C
C     DISTANCE ENTRE LES POINTS SUCCESSIFS
      REAL    DISTAN, RESTE, DISTAB
      INTEGER I
C
C     INITIALISATIONS
      DISTAN=4*SQRT(A**2+B**2)/REAL(NBPOINTS)
      DISTAB=DISTAN
      DO 20 I=1,NBPOINTS
          THETAN(I)=0.0
          THETAB(I)=0.0
 20   ENDDO
 30   DO 40 I=2,NBPOINTS
          CALL DICHOTOMI(A,B,THETAB(I-1),DISTAB,0.1*EPSILON,
     %    NBPOINTS,THETAB(I))
          IF ( THETAB(I) .GT. 8.*ATAN(1.) ) THEN
            GOTO 60
          ENDIF
 40   ENDDO
      RESTE=SQRT((A*COS(THETAB(NBPOINTS))-A)**2+(B*SIN(
     %THETAB(NBPOINTS)))**2)
      IF (ABS(RESTE-DISTAB).LT.(EPSILON*DISTAB)) THEN
          DO 45 I=1,NBPOINTS
            THETAN(I)=THETAB(I)
 45       ENDDO
          DISTAN=DISTAB
          GOTO 99
      ELSEIF ((RESTE-DISTAB).GT.(EPSILON*DISTAB)) THEN
          DO 50 I=1,NBPOINTS
            THETAN(I)=THETAB(I)
 50       ENDDO
          DISTAN=DISTAB
          DISTAB=DISTAB+(RESTE-DISTAN)/REAL(NBPOINTS)
          GOTO 30
      ENDIF
C
 60   DO 55 I=1,NBPOINTS
         THETAB(I)=THETAN(I)
 55   ENDDO
      DISTAB=(DISTAN+DISTAB)/2.0
      GOTO 30
C
 99   RETURN
      END
C
      SUBROUTINE DICHOTOMI(A,B,T1,LONGUEUR,ERREUR,NBPOINTS,T2)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  GENERER LE PARAMETRE DU POINT SITUE SUR UNE ELLIPSE DE       +
C -----  DEMI-AXES A ET B, ET DISTANT DE 'LONGUEUR' DU POINT          +
C        DE PARAMETRS T1.                                             +
C                                                                     +
C ENTREES :                                                           +
C ---------                                                           +
C A,B      : DEMI-AXES DE L'ELLIPSE                                   +
C LONGUEUR : LONGUEUR DE LA DISTANCES INTERPOINT                      +
C NBPOINTS : NOMBRE DE POINTS                                         +
C T1       : PARAMETRE DU PREMIER POINT                               +
C ERREUR   : TOLERANCE SUR LA CONVERGENCE DE LA DICHOTOMIE            +
C                                                                     +
C SORTIES :                                                           +
C ---------                                                           +
C T2      : PARAMETRE DU SECOND POINT                                 +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS       : JOSE URQUIZA & ARNAUD MONTENAY : DECEMBRE 1994      +
C MODIFICATIONS : RODOLFO ARAYA & PASCAL MAILLOT : JANVIER 1996       +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     PARAMETRES D'ENTREE
      REAL A,B,T1,LONGUEUR,ERREUR
      INTEGER NBPOINTS
C
C     PARAMETRES DE SORTIE
      REAL T2
C
C     VARIABLES LOCALES
      REAL PAS,X1,X2,D,M

      PAS= 8.0 * ATAN(1.) /(REAL(NBPOINTS)*3.0)
      X1=0.0
      X2=T1+PAS
 5    D=SQRT((A*(COS(X2)-COS(T1)))**2+(B*(SIN(X2)-SIN(T1)))**2)
      IF(D.LT.LONGUEUR) THEN
          X1=X2
          X2=X2+PAS
          GOTO 5
      ENDIF
 10   M=(X1+X2)/2.0
      D=SQRT((A*(COS(M)-COS(T1)))**2+(B*(SIN(M)-SIN(T1)))**2)
      IF(D.GT.(LONGUEUR*(1.0+ERREUR))) THEN
          X2=M
      ELSE IF(D.LT.(LONGUEUR*(1.0-ERREUR))) THEN
          X1=M
      ELSE
          T2=M
          GOTO 15
      ENDIF
      GOTO 10
 15   RETURN
      END
