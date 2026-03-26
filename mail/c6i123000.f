      SUBROUTINE C6I123000( NOPROJ, NBEF6C, NSEF6C, NBST6C, XYZUVW6C,
     %                      NBVECT, TEMPE6C,
     %                      POLY6C, POLY3C, DPOLY3C,
     %                      MXEF3C, NBEF3C, XYZ3C, TEMPE3C, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER les 3 COORDONNEES  XYZ3C des NBEF3C 3Q1C
C -----   INTERSECTION DES 6CUBES AVEC U=V=W=0 et LA TEMPERATURE
C         EN LES SOMMETS DES 3CUBES.
C         PAS D'IDENTIFICATION DES MEMES SOMMETS DANS L'INTERSECTION
C         EN FAIT LA COORDONNEE 0 EST LA COORDONNEE MEDIANE ENTRE MIN ET MAX
C
C ENTREES:
C --------
C NOPROJ  : TYPE DE PROJECTION 0 CI-DESSOUS FIXE LA COORDONNEE A ZERO
C           0 : A DEFINIR LES TABLEAUX INDFIXES0(3) et COOFIXES(3)
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C NBEF6C : NOMBRE DE 6CUBES DU MAILLAGE
C NSEF6C : NSEF6C(NEF,NS) NO DU SOMMET NS DU 6CUBE NEF DU MAILLAGE
C NBST6C : NOMBRE DE SOMMETS DES 6CUBES DU MAILLAGE
C XYZUVW6C: 6 COORDONNEES DES SOMMETS DES 6CUBES
C TEMPE6C: TEMPERATURE AUX NBST6C SOMMETS DU MAILLAGE
C POLY6C : POLY6C(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY6C(I,J,K,L,M,N,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C                                               U**(L-1) V**(M-1) W**(N-1)
C POLY3C : POLY3C(L,Q) COEFFICIENT L DU Q-EME POLYNOME  I.E.
C          POLY3C(I,J,K,Q)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C DPOLY3C: DPOLY3C(L,Q,M) COEFFICIENT L DU Q-EME POLYNOME DERIVE DANS
C                       la DIRECTION M
C          DPOLY3C(I,J,K,Q,M)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C MXEF3C : MAXIMUM DE 3Q1C DECLARABLES
C
C SORTIES:
C --------
C NBEF3C : NOMBRE DE 3CUBES CREES DANS LE MAILLAGE A U=V=W=0
C XYZ3C  : 3 COORDONNEES DES SOMMETS DES 3CUBES A U=V=W=0
C TEMPE3C: TEMPERATURE AUX 8 * NBEF3C SOMMETS DU MAILLAGE 3CUBES
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF DEGENERE, 2 MXEF3C TROP PETIT...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L LIONS UPMC Paris Novembre2006
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/langue.inc"
      include"./incl/xyzext.inc"
      include"./incl/gsmenu.inc"
      include"./incl/dirpro.inc"
      REAL     XYZUVW6C(6,NBST6C), XYZ3C(3,8,MXEF3C)
      INTEGER  NSEF6C(NBEF6C,64)
      INTEGER  INDGARDES(6), INDFIXES0(3)
      REAL     COOFIXES(3)
C
      DOUBLE PRECISION  POLY6C(64,64), POLY3C(8,8), DPOLY3C(8,8,3)
      DOUBLE PRECISION  TEMPE6C(NBST6C,NBVECT), TEMPE3C(8,MXEF3C,NBVECT)
      DOUBLE PRECISION  V, XYZFe(3), A, B, C
      DOUBLE PRECISION  UVW(1:3,8), abcdef(6,8), def(3), XYZ(3)
      REAL              XMIN(1:3), XMAX(1:3), XK
      CHARACTER*1       LETTRE(6)
      DATA              LETTRE / 'X', 'Y', 'Z', 'U', 'V', 'W' /
      INTEGER           PUISS2(6)
      DATA              PUISS2 / 1, 2, 4, 8, 16, 32 /
C
      IERR = 0
      GOTO( 1, 2, 3, 4, 5 ), NOPROJ+1
C
C     0 : 'DEFINIR 3 DIRECTIONS et COORDONNEES'
 1    DO 12 K=1,3
C
C        LECTURE DU NUMERO DE LA DIRECTION FIXEE (1 a 6)
         CALL INVITE( 138 )
         N = 3+K
         NCVALS = 4
         CALL LIRENT( NCVALS , N )
         IF( NCVALS .LE. 0 ) THEN
            IERR = 1
            RETURN
         ENDIF
         INDFIXES0(K) = N
C
C        LECTURE DE LA VALEUR DE LA COORDONNE FIXEE
 11      CALL INVITE( 139 )
         CALL LIRRSP( NCVALS , R )
         IF( NCVALS .LE. 0 ) THEN
            IERR = 2
            RETURN
         ENDIF
         IF( R .LT. COOEXT(N,1) .OR. R .GT. COOEXT(N,2) ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'COORDONNEE TROP PETITE ou TROP GRANDE'
            ELSE
               KERR(1) = 'TOO SMALL or TOO LARGE COORDINATE'
            ENDIF
            CALL LEREUR
            GOTO 11
         ENDIF
         COOFIXES(K) = R
 12   CONTINUE
C
C     CONSTRUCTION DU TABLEAU INDGARDES(1:3)
      N = 0
      DO 15 K=1,6
         IF( K .NE. INDFIXES0(1) .AND.
     %       K .NE. INDFIXES0(2) .AND.
     %       K .NE. INDFIXES0(3) ) THEN
            N = N + 1
            INDGARDES(N) = K
         ENDIF
 15   CONTINUE
      IF( N .GE. 4 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'COORDONNEE FIXEE en DOUBLE'
            KERR(2) = 'REDONNER les 3 DIRECTIONS a FIXER parmi les 6'
         ELSE
            KERR(1) = 'TWICE FIXED COORDINATE'
            KERR(2) = 'GIVE AGAIN the 3 FIXED DIRECTIONS among 6'
         ENDIF
         CALL LEREUR
         GOTO 1
      ENDIF
      DO 16 K=1,3
         I = INDGARDES(K)
         KDIRPRO(I:I) = LETTRE(I)
         I = INDFIXES0(K)
         KDIRPRO(I:I) = '0'
 16   CONTINUE
      WRITE(IMPRIM,*) 'Projection XYZUVW => ',KDIRPRO
      GOTO 10
C
C     1 : 'X Y Z 0 0 0'
 2    INDGARDES(1) = 1
      INDGARDES(2) = 2
      INDGARDES(3) = 3
      INDFIXES0(1) = 4
      INDFIXES0(2) = 5
      INDFIXES0(3) = 6
      WRITE(IMPRIM,*) 'Projection XYZUVW => XYZ000'
      GOTO 6
C
C     2 : 'X Y 0 U 0 0'
 3    INDGARDES(1) = 1
      INDGARDES(2) = 2
      INDGARDES(3) = 4
      INDFIXES0(1) = 3
      INDFIXES0(2) = 5
      INDFIXES0(3) = 6
      WRITE(IMPRIM,*) 'Projection XYZUVW => XY0U00'
      GOTO 6
C
C     3 : 'X 0 0 U V 0'
 4    INDGARDES(1) = 1
      INDGARDES(2) = 4
      INDGARDES(3) = 5
      INDFIXES0(1) = 2
      INDFIXES0(2) = 3
      INDFIXES0(3) = 6
      WRITE(IMPRIM,*) 'Projection XYZUVW => X00UV0'
      GOTO 6
C
C     4 : '0 0 0 U V W'
 5    INDGARDES(1) = 4
      INDGARDES(2) = 5
      INDGARDES(3) = 6
      INDFIXES0(1) = 1
      INDFIXES0(2) = 2
      INDFIXES0(3) = 3
      WRITE(IMPRIM,*) 'Projection XYZUVW => 000UVW'
C
 6    DO 7 K=1,3
C        LE MILIEU ENTRE MIN ET MAX DE CETTE COORDONNEE IND
         IND = INDFIXES0(K)
         COOFIXES(K) = ( COOEXT(IND,1) + COOEXT(IND,2) ) / 2
 7    CONTINUE
C
 10   NBEF3C = 0
C
      DO 100 NEF6C=1,NBEF6C
C
C        LE 6CUBE NEF6C INTERSECTE T IL U=V=W=0?
C        RECHERCHE DU SOMMET MIN DES Xk>0 DU 6CUBE NEF6C k=4,5,6
C        RECHERCHE DU SOMMET MAX DES Xk<0 DU 6CUBE NEF6C k=4,5,6
         DO 18 K=1,3
            XMIN(K) = 1e27
            XMAX(K) =-1e27
 18      CONTINUE
C
C        RECHERCHE DES COORDONNEES MIN POSITIVES ET
C        MAX NEGATIVES DANS LES DIRECTIONS FIXEES
C        PARMI LES 64 SOMMETS DU 6CUBE
         DO 30 NS64=1,64
C           LE NO DU SOMMET NS64 DU 6CUBE NEF6C
            NS = NSEF6C(NEF6C,NS64)
            DO 20 K=1,3
               IND = INDFIXES0(K)
C              LE MILIEU ENTRE MIN ET MAX DE CETTE COORDONNEE IND
               XM = COOFIXES(K)
C              LA COORDONNEE IND DU SOMMET NS64
               XK = XYZUVW6C(IND,NS)
               IF( XK .GE. XM .AND. XK .LT. XMIN(K) ) XMIN(K) = XK
               IF( XK .LE. XM .AND. XK .GT. XMAX(K) ) XMAX(K) = XK
 20         CONTINUE
 30      CONTINUE
C
         IF( XMIN(1) .NE. 1e27 .AND. XMAX(1) .NE. -1e27 .AND.
     %       XMIN(2) .NE. 1e27 .AND. XMAX(2) .NE. -1e27 .AND.
     %       XMIN(3) .NE. 1e27 .AND. XMAX(3) .NE. -1e27 ) THEN
C
C            LE 6CUBE A UNE INTERSECTION AVEC XINDFIXES0(1:3)=COOFIXES
C            => 1 3CUBE DE PLUS
             IF( NBEF3C .GE. MXEF3C ) THEN
                IERR = 3
                GOTO 9000
             ENDIF
             NBEF3C = NBEF3C + 1
C
C            CONSTRUCTION DES 8 SOMMETS IJK DU 3CUBE
C            AVEC a=0 ou 1, b=0 ou 1, c=0 ou 1
             IJK = 0
             C = 0d0
             DO 70 K=0,1
                B = 0d0
                DO 60 J=0,1
                   A = 0d0
                   DO 50 I=0,1
C
C                     PARCOURS DES 8 SOMMETS AYANT les MEMES abc 0 ou 1
                      IJK = IJK + 1
                      LMN = 0
C
C                     NO DU PREMIER SOMMET X=A Y=B Z=C
                      N64ABC = 1 + I * PUISS2( INDGARDES(1) )
     %                           + J * PUISS2( INDGARDES(2) )
     %                           + K * PUISS2( INDGARDES(3) )
                      DO 39 N=0,1
                         DO 36 M=0,1
                            DO 33 L=0,1
C                              NO SOMMET L,M,N DE 1 A 64 AYANT MEME ABC
                               N64 = N64ABC
     %                             + L * PUISS2( INDFIXES0(1) )
     %                             + M * PUISS2( INDFIXES0(2) )
     %                             + N * PUISS2( INDFIXES0(3) )
C
C                              LE NO DU SOMMET AYANT MEME ABC
                               NS64 = NSEF6C(NEF6C,N64)
C
                               LMN = LMN + 1
                               UVW(1,LMN) = XYZUVW6C(INDFIXES0(1), NS64)
                               UVW(2,LMN) = XYZUVW6C(INDFIXES0(2), NS64)
                               UVW(3,LMN) = XYZUVW6C(INDFIXES0(3), NS64)
 33                         CONTINUE
 36                      CONTINUE
 39                   CONTINUE
C
C                     CALCUL de def DANS LE 3CUBE DE REFERENCE
C                     TEL QUE [3Q1C(def)]=[COOFIXES]
                      XYZ(1) = COOFIXES(1)
                      XYZ(2) = COOFIXES(2)
                      XYZ(3) = COOFIXES(3)
                      CALL F3Q1INV( POLY3C, DPOLY3C, UVW, XYZ,
     %                              def, IERR )
                      IF( IERR .NE. 0 ) RETURN
C
C                     LE POINT SUR LE 6CUBE DE REFERENCE
                      abcdef(INDGARDES(1),IJK) = A
                      abcdef(INDGARDES(2),IJK) = B
                      abcdef(INDGARDES(3),IJK) = C
                      abcdef(INDFIXES0(1),IJK) = def(1)
                      abcdef(INDFIXES0(2),IJK) = def(2)
                      abcdef(INDFIXES0(3),IJK) = def(3)
C
                      XYZFe(1) = 0d0
                      XYZFe(2) = 0d0
                      XYZFe(3) = 0d0
C
C                     CALCUL DE LA TEMPERATURE EN CE POINT Fe(abcdef(IJK))
                      DO 41 NBV=1,NBVECT
                         TEMPE3C(IJK,NBEF3C,NBV) = 0d0
 41                   CONTINUE
C
                      DO 45 L=1,64
C
C                        VALEUR DU POLYNOME L AU POINT abcdef
                         CALL PN6DVA( 2, POLY6C(1,L),
     %                                abcdef(1,IJK), abcdef(2,IJK),
     %                                abcdef(3,IJK), abcdef(4,IJK),
     %                                abcdef(5,IJK), abcdef(6,IJK),
     %                                V )
C
C                        Fe( abcdef(IJK) ) = XYZFe
                         NS = NSEF6C(NEF6C,L)
                         DO 42 KK=1,3
                            IND = INDGARDES(KK)
                            XYZFe(KK) = XYZFe(KK) + V * XYZUVW6C(IND,NS)
 42                      CONTINUE
C
C                        TEMPERATURES AU SOMMET L
                         DO 44 NBV=1,NBVECT
                            TEMPE3C(IJK,NBEF3C,NBV) =
     %                      TEMPE3C(IJK,NBEF3C,NBV)
     %                    + V * TEMPE6C(NSEF6C(NEF6C,L),NBV)
 44                      CONTINUE
C
 45                   CONTINUE
C
C                     COORDONNEES XYZ DU POINT IJK DU 3Q1C
                      XYZ3C(1,IJK,NBEF3C) = REAL( XYZFe(1) )
                      XYZ3C(2,IJK,NBEF3C) = REAL( XYZFe(2) )
                      XYZ3C(3,IJK,NBEF3C) = REAL( XYZFe(3) )
C
                      A = 1d0
 50                CONTINUE
                   B = 1d0
 60             CONTINUE
                C = 1d0
 70          CONTINUE
          ENDIF
 100   CONTINUE
C
 9000  WRITE(IMPRIM,*) 'NOPROJ=',NOPROJ,': EXTRACTION',NBEF3C,
     %' 3CUBES pour MXEF3C=',MXEF3C,' parmi ',NBEF6C,' 6CUBES'
C
       RETURN
       END
