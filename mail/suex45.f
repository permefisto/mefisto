      SUBROUTINE SUEX45( NTLXSU, LADEFI,
     %                   NTFASU, MNFASU, NTSOFA, MNSOFA, IERR)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTERSECTION, DE LA DIFFERENCE OU DE L'UNION DE
C -----    DEUX MAILLAGES DE SURFACES
C          EN 2D ==>PLANS<== C-A-D 2 SURFACES DE R2
C          EN 3D C-A-D 2 SURFACES FERMEES (D'UN VOLUME ou SOLIDE)
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TMS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU DES ENTIERS DE DEFINITION DE LA SURFACE
C          CF ~td/d/a_surface__definition
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C          CF ~td/d/a___nsef
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : BALDENSPERGER-MATHIEU DEA ANALYSE NUMERIQUE UPMC JANVIER 1999
C AJOUT3D: PERRONNET Alain LJLL UPMC& St Pierre du Perray SEPTEMBRE 2011
C2345X7..............................................................012
      INCLUDE "./incl/langue.inc"
      INCLUDE "./incl/gsmenu.inc"
      INCLUDE "./incl/a_surface__definition.inc"
      INCLUDE "./incl/a_ligne__definition.inc"
      INCLUDE "./incl/a___arete.inc"
      INCLUDE "./incl/a___xyzsommet.inc"
      INCLUDE "./incl/a___nsef.inc"
      INCLUDE "./incl/ntmnlt.inc"
C
      include"./incl/pp.inc"
      COMMON       MCN(MOTMCN)
      REAL        RMCN(1)
      EQUIVALENCE (MCN(1), RMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      INTEGER          NTFASU, MNFASU, NTSOFA, MNSOFA, IERR
      INTEGER          LADEFI(0:*)
C
      INTEGER          NUTRI1, NUTRI2, NUOPER
      INTEGER          TMSTR1, MCTRI1, TMSTR2, MCTRI2
      INTEGER          TMSEF1, MCNEF1, TMSEF2, MCNEF2
      INTEGER          TMSXY1, MCNXY1, TMSXY2, MCNXY2
      INTEGER          TAILL1, TAILL2, MCSEG1, MCSEG2
      INTEGER          MCEDG1, MCXED1, NBEDG1, MCCPC1, NBCPC1
      INTEGER          MCEDG2, MCXED2, NBEDG2, MCCPC2, NBCPC2
      INTEGER          MCNCOU, NBCOUR, NBPTS, PTDEP1, PTDEP2, ASUIVR
C
      INTEGER          EDGSZ1, EDGSZ2
      INTEGER          ISEG1, ISEG1S, ISEG2, ISEG2S
      DOUBLE PRECISION P11X, P11Y, P12X, P12Y
      DOUBLE PRECISION P21X, P21Y, p22X, P22Y, P13X, P13Y
      DOUBLE PRECISION P11P12(2), P11P21(2), P11P22(2), P21P22(2)
      DOUBLE PRECISION P12P21(2), DETERM(5), INTERS(2), INTEXY(2)
      DOUBLE PRECISION P11P13(2)
      INTEGER          SIGNE1, SIGNE2, SIGNE3
      INTEGER          LIEN1, LIEN2
      CHARACTER*24     NOMCOU, NUMERO
      INTEGER          LXCOUR, MCCOUR, NTNSEF, MCNSEF
      INTEGER          NTXYZS, MCXYZS
      INTEGER          I, J, K, INSIDE, TOKEEP
      INTEGER          MCLADS
      INTEGER          NUCEXT
      REAL             CXMIN
CCC   INTEGER          POS
C
C     NUMERO DES 2 SURFACES A TRAITER
      NUTRI1 = LADEFI(WUTRI1)
      NUTRI2 = LADEFI(WUTRI2)
C
C     NUMERO DE L'OPERATEUR LOGIQUE A TRAITER
C     ( 1: 'intersection' , 2: 'addition' , 3: 'soustraction' )
      NUOPER = LADEFI(WUOPER)
C
C     OUVERTURE DE LA SURFACE 1
      CALL LXNLOU (NTSURF, NUTRI1, TMSTR1, MCTRI1)
      IF (TMSTR1.LE.0) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LA SURFACE 1 EST INCONNUE'
         ELSE
            KERR(1) = 'SURFACE 1 IS UNKNOWN'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     OUVERTURE DE LA SURFACE 2
      CALL LXNLOU (NTSURF, NUTRI2, TMSTR2, MCTRI2)
      IF (TMSTR2.LE.0) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LA SURFACE 2 EST INCONNUE'
         ELSE
            KERR(1) = 'SURFACE 2 IS UNKNOWN'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     Ouverture des TMS contenant le maillage des deux surfaces
      CALL LXTSOU (TMSTR1, 'NSEF', TMSEF1, MCNEF1)
      IF (TMSEF1.LE.0) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE 1 SANS tms NSEF'
         ELSE
            KERR(1) = 'SURFACE 1 WITHOUT NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
      CALL LXTSOU (TMSTR2, 'NSEF', TMSEF2, MCNEF2)
      IF (TMSEF2.LE.0) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE 2 SANS tms NSEF'
         ELSE
            KERR(1) = 'SURFACE 2 WITHOUT NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
      IF ( MCN(MCNEF1+WUTYOB).NE.3 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE 1 N''EST PAS UNE SURFACE'
         ELSE
            KERR(1) = 'SURFACE 1 is NOT a SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 5
         RETURN
      ENDIF
C
      IF ( MCN(MCNEF2+WUTYOB).NE.3 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE 2 N''EST PAS UNE SURFACE'
         ELSE
            KERR(1) = 'SURFACE 2 is NOT a SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 6
         RETURN
      ENDIF
C
C     Ouverture du TMS contenant les coordonnees des sommets
      CALL LXTSOU (TMSTR1, 'XYZSOMMET', TMSXY1, MCNXY1)
      IF (TMSXY1 .LE. 0) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE 1 N''A PAS DE tms XYZSOMMET'
         ELSE
            KERR(1) = 'SURFACE 1 HAS NOT a tms XYZSOMMET'
         ENDIF
         CALL LEREUR
         IERR = 7
         RETURN
      ENDIF
      CALL DIMCOO( MCN(MCNXY1+WNBSOM), MCN(MCNXY1+WYZSOM), NDIM1 )
C
C     NDIM: DIMENSION DE L'ESPACE DE DEFINITION DES 2 SURFACES
C           =2 LES 2 SURFACES SONT DANS LE MEME PLAN
C           =3 CHACUNE DES 2 SURFACES EST FERMEE
C             (TOUTE ARETE APPARTIENT A 2 EF)
C
      CALL LXTSOU (TMSTR2, 'XYZSOMMET', TMSXY2, MCNXY2)
      IF (TMSXY2 .LE. 0) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE 2 N''A PAS DE tms XYZSOMMET'
         ELSE
            KERR(1) = 'SURFACE 2 HAS NOT a tms XYZSOMMET'
         ENDIF
         CALL LEREUR
         IERR = 8
         RETURN
      ENDIF
      CALL DIMCOO( MCN(MCNXY2+WNBSOM), MCN(MCNXY2+WYZSOM), NDIM )
C
      IF( NDIM1 .NE. NDIM ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = '2 SURFACES NON dans le MEME ESPACE R2 ou R3'
         ELSE
            KERR(1) = '2 SURFACES NOT in the SAME SPACE R2 or R3'
         ENDIF
         CALL LEREUR
         IERR = 9
         RETURN
      ENDIF
C
      IF( NDIM1 .EQ. 3 .AND. NDIM .EQ. 3 ) THEN
C
C        =======================================================
C        LES 2 SURFACES SONT DANS R3 C-A-D FRONTIERE D'UN VOLUME
C        =======================================================
         CALL OPLOR3( NUOPER, NTLXSU,
     %                MCNEF1, MCNXY1, MCNEF2, MCNXY2,
     %                NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
         GOTO 9999
C
      ENDIF
C
C     ===================================================
C     LES 2 SURFACES SONT DANS R2 C-A-D DANS UN MEME PLAN
C     ===================================================
C     Construction des contours de chaque surface :
C     Surface 1 : * construction de la hash-table
C                 * construction des contours
C                 * destruction de la hash-table
      CALL HASH45 (MCNEF1, MCNXY1, MCSEG1, TAILL1)
      CALL EDGE45 (MCSEG1, TAILL1, MCNXY1, MCEDG1, MCXED1, NBEDG1,
     &             MCCPC1, NBCPC1)
      CALL TNMCDS ('ENTIER', TAILL1 * 3, MCSEG1)
      EDGSZ1 = NBEDG1
C
C     Surface 2 : * construction de la hash-table
C                 * construction des contours
C                 * destruction de la hash-table
      CALL HASH45 (MCNEF2, MCNXY2, MCSEG2, TAILL2)
      CALL EDGE45 (MCSEG2, TAILL2, MCNXY2, MCEDG2, MCXED2, NBEDG2,
     &             MCCPC2, NBCPC2)
      CALL TNMCDS ('ENTIER', TAILL2 * 3, MCSEG2)
      EDGSZ2 = NBEDG2
CCCC
CCCC     Fermeture des LX et TMS des 2 SURFACES INITIALES
CCC      CALL LXNLFE (NTSURF, NUTRI1)
CCC      CALL LXNLFE (NTSURF, NUTRI2)
CCCC
CCC   LE CONTOUR EN ARETES DU MAILLAGE DE LA SURFACE 1
CCC   WRITE (IMPRIM,*)
CCC   WRITE (IMPRIM,*) 'NOMBRE ARETES CONTOUR SURFACE 1:',NBEDG1
CCC   DO I = 0, NBEDG1 - 1
CCC      POS = MCEDG1 + I * 4
CCC      WRITE (IMPRIM,*) 'POINT ', I+1, (MCN(POS+J),J=1,3)
CCC   ENDDO
CCC   WRITE (IMPRIM,*)
CCC   DO I = 0, NBCPC1 - 1
CCC   WRITE (IMPRIM,*) 'LE POINT ', MCN(MCCPC1 + I),
CCC  &                   ' EST SUR LA ',I + 1,'-EME COMPOSANTE CONNEXE'
CCC   ENDDO
CCC
CCC   LE CONTOUR EN ARETES DU MAILLAGE DE LA SURFACE 2
CCC   WRITE (IMPRIM,*)
CCC   WRITE (IMPRIM,*) 'NOMBRE ARETES CONTOUR SURFACE 2:',NBEDG2
CCC   DO I = 0, NBEDG2 - 1
CCC     POS = MCEDG2 + I * 4
CCC     WRITE (IMPRIM,*) 'POINT ', I+1, (MCN(POS+J),J=1,3)
CCC   ENDDO
CCC   WRITE (IMPRIM,*)
CCC   DO I = 0, NBCPC2 - 1
CCC     WRITE (IMPRIM,*) 'LE POINT ', MCN(MCCPC2 + I),
CCC  &                   ' EST SUR LA ',I + 1,'-EME COMPOSANTE CONNEXE'
CCC   ENDDO
C
C     AJOUT DES POINTS D'INTERSECTIONS ENTRE LES CONTOURS:
C     Boucle sur les segments de la surface 1
C     (parcourus dans l'ordre ou ils apparaissent dans la liste)
      ISEG1 = 1
 100  P11X = RMCN(MCXED1 + (ISEG1-1)*2)
      P11Y = RMCN(MCXED1 + (ISEG1-1)*2 + 1)
C     L'indice du segment suivant = ISEG1S
      ISEG1S = MCN(MCEDG1 + (ISEG1-1)*4 + 1)
      P12X = RMCN(MCXED1 + (ISEG1S-1)*2)
      P12Y = RMCN(MCXED1 + (ISEG1S-1)*2 + 1)
C     Le vecteur P11P12 :
      P11P12(1) = P12X - P11X
      P11P12(2) = P12Y - P11Y
C
C     Boucle sur les contours de la surface 2
      DO 410 ICONT = 0, NBCPC2 - 1
C
         ISEG2 = MCN(MCCPC2 + ICONT)
         P21X = RMCN(MCXED2 + (ISEG2-1)*2)
         P21Y = RMCN(MCXED2 + (ISEG2-1)*2 + 1)
C        Le vecteur P11P21
         P11P21(1) = P21X - P11X
         P11P21(2) = P21Y - P11Y
C        Le determinant P11P12 * P11P21
         DETERM(1) = P11P12(1) * P11P21(2) - P11P12(2) * P11P21(1)
         IF (DETERM(1) .LT. 0.0D0) THEN
            SIGNE1 = -1
         ELSE IF (DETERM(1) .GT. 0.0D0) THEN
            SIGNE1 = 1
         ELSE
            SIGNE1 = 0
         ENDIF
C
C        Boucle sur les segments du contour ICONT de la surface 2
C        (parcourus en suivant la courbe).
C
C           L'indice du segment suivant = ISEG1S
 200        ISEG2S = MCN(MCEDG2 + (ISEG2-1)*4 + 1)
            P22X = RMCN(MCXED2 + (ISEG2S-1)*2)
            P22Y = RMCN(MCXED2 + (ISEG2S-1)*2 + 1)
C           Le vecteur P11P22 :
            P11P22(1) = P22X - P11X
            P11P22(2) = P22Y - P11Y
C           Le determinant P11P12 * P11P22
            DETERM(2) = P11P12(1) * P11P22(2) - P11P12(2) * P11P22(1)
            IF (DETERM(2) .LT. 0.0D0) THEN
               SIGNE2 = -1
            ELSE IF (DETERM(2) .GT. 0.0D0) THEN
               SIGNE2 = 1
            ELSE
               SIGNE2 = 0
            ENDIF
C           Si les determinants sont non nuls et de meme signe,
C           il n'y a pas d'intersection on passe au segment suivant
C           (c'est le cas le plus frequent)
            IF (SIGNE1 * SIGNE2 .GT. 0) GOTO 400
C           Si les determinants sont nuls tous les deux, les
C           segments sont sur une meme droite.
C           On dit alors : pas d'intersection
            IF ((SIGNE1 .EQ. 0) .AND. (SIGNE2 .EQ. 0)) GOTO 400
C
C           Il y a une possibilite d'intersection
C           On teste les determinants croises
            P21P22(1) = P22X - P21X
            P21P22(2) = P22Y - P21Y
            DETERM(3) = - P21P22(1) * P11P21(2) + P21P22(2) * P11P21(1)
            P12P21(1) = P21X - P12X
            P12P21(2) = P21Y - P12Y
            DETERM(4) = - P21P22(1) * P12P21(2) + P21P22(2) * P12P21(1)
            IF (DETERM(3) .LT. 0.0D0) THEN
               SIGNE1 = -1
            ELSE IF (DETERM(3) .GT. 0.0D0) THEN
               SIGNE1 = 1
            ELSE
               SIGNE1 = 0
            ENDIF
            IF (DETERM(4) .LT. 0.0D0) THEN
               SIGNE3 = -1
            ELSE IF (DETERM(4) .GT. 0.0D0) THEN
               SIGNE3 = 1
            ELSE
               SIGNE3 = 0
            ENDIF
            IF (SIGNE1 * SIGNE3 .GT. 0) GOTO 400
            IF ((SIGNE1 .EQ. 0) .AND. (SIGNE3 .EQ. 0)) GOTO 400
C
C           On est sur que les segments se coupent, et qu'ils ne sont
C           pas paralleles : on resout le systeme 2x2
            DETERM(5) = - P11P12(1) * P21P22(2) + P11P12(2) * P21P22(1)
            DETERM(5) = 1.0D0 / DETERM(5)
            INTERS(1) = (-P21P22(2)*P11P21(1)+P21P22(1)*P11P21(2))
     &                  * DETERM(5)
            INTERS(2) = (-P11P12(2)*P11P21(1)+P11P12(1)*P11P21(2))
     &                  * DETERM(5)
C
C           Numeriquement, on peut etre a l'exterieur du segment.
C           Or, on est sur qu'il y a intersection, donc on
C           ramene l'intersection a l'interieur des segments.
CCC         WRITE(*,*) 'INTER',INTEXY(1),INTEXY(2), INTERS(1),INTERS(2)
            IF (INTERS(1) .LT. 0.0D0) THEN
               INTERS(1) = 0.0D0
            ELSE IF (INTERS(1) .GT. 1.0D0) THEN
               INTERS(1) = 1.0D0
            ENDIF
            IF (INTERS(2) .LT. 0.0D0) THEN
               INTERS(2) = 0.0D0
            ELSE IF (INTERS(2) .GT. 1.0D0) THEN
               INTERS(2) = 1.0D0
            ENDIF
C
C           Calcul des coordonnees de l'intersection
            INTEXY(1) = P11X + P11P12(1) * INTERS(1)
            INTEXY(2) = P11Y + P11P12(2) * INTERS(1)
CCC         WRITE (IMPRIM,*) 'INTER',INTEXY(1),INTEXY(2),
CCC  &                        INTERS(1),INTERS(2)
CCC         WRITE (IMPRIM,*) 'EXT1', P11X, P11Y, P12X, P12Y
CCC         WRITE (IMPRIM,*) 'EXT2', P21X, P21Y, P22X, P22Y
C
C           Si le point d'intersection est au milieu d'un segment,
C           on ajoute ce point dans la courbe (il faut l'inserer
C           dans la liste doublement chainee). On utilise une
C           tolerance de 5% de la longueur du segment, pour ne pas
C           creer de segment trop court qui produirait un mauvais
C           triangle.
            IF ((INTERS(1) .GT. 0.05) .AND. (INTERS(1) .LT. 0.95)) THEN
               IF (NBEDG1 .LE. EDGSZ1) THEN
C                 Il n'y a plus de place pour ajouter un nouveau point.
C                 On en fait !
C                 On augmente la taille des tableaux MCEDG1 et MCXED1
                  EDGSZ1 = EDGSZ1 + 50
                  CALL TNMCAU ('ENTIER',NBEDG1*4,EDGSZ1*4,NBEDG1*4,
     &                         MCEDG1)
                  CALL TNMCAU ('REEL', NBEDG1*2, EDGSZ1*2, NBEDG1*2,
     &                         MCXED1)
               ENDIF
               NBEDG1 = NBEDG1 + 1
               MCN(MCEDG1+(NBEDG1-1)*4) = 0
               MCN(MCEDG1+(NBEDG1-1)*4+3) = 0
C              Points suivant et precedent du nouveau point
               MCN(MCEDG1+(NBEDG1-1)*4+1) = ISEG1S
               MCN(MCEDG1+(NBEDG1-1)*4+2) = ISEG1
C              Point precedent du point suivant
               MCN(MCEDG1+(ISEG1S-1)*4+2) = NBEDG1
C              Point suivant du point precedent
               MCN(MCEDG1+(ISEG1-1)*4+1) = NBEDG1
C              Coordonnees
               RMCN(MCXED1+(NBEDG1-1)*2  ) = REAL( INTEXY(1) )
               RMCN(MCXED1+(NBEDG1-1)*2+1) = REAL( INTEXY(2) )
C              Numero du point ajoute pour le lien
               LIEN1 = NBEDG1
CCC            WRITE (IMPRIM,*) 'NEW1',LIEN1
            ELSE
C
C             Si l'intersection est sur une extremite, on n'ajoute
C             pas de point
               IF (INTERS(1) .LT. 0.5D0) THEN
CCC               WRITE (IMPRIM,*) 'OLD1',ISEG1
                  LIEN1 = ISEG1
               ELSE
CCC               WRITE (IMPRIM,*) 'OLD1',ISEG1S
                  LIEN1 = ISEG1S
               ENDIF
            ENDIF
C
C           Maintenant, on repete la meme chose avec le segment
C           de la deuxieme courbe.
            IF ((INTERS(2) .GT. 0.05) .AND. (INTERS(2) .LT. 0.95)) THEN
               IF (NBEDG2 .LE. EDGSZ2) THEN
C                 Il n'y a plus de place pour ajouter un nouveau point.
C                 On en fait !
C                 On augmente la taille des tableaux MCEDG1 et MCXED1
                  EDGSZ2 = EDGSZ2 + 50
                  CALL TNMCAU ('ENTIER',NBEDG2*4,EDGSZ2*4,NBEDG2*4,
     &                         MCEDG2)
                  CALL TNMCAU ('REEL', NBEDG2*2, EDGSZ2*2, NBEDG2*2,
     &                         MCXED2)
               ENDIF
               NBEDG2 = NBEDG2 + 1
               MCN(MCEDG2+(NBEDG2-1)*4) = 0
               MCN(MCEDG2+(NBEDG2-1)*4+3) = 0
C              Points suivant et precedent du nouveau point
               MCN(MCEDG2+(NBEDG2-1)*4+1) = ISEG2S
               MCN(MCEDG2+(NBEDG2-1)*4+2) = ISEG2
C              Point precedent du point suivant
               MCN(MCEDG2+(ISEG2S-1)*4+2) = NBEDG2
C              Point suivant du point precedent
               MCN(MCEDG2+(ISEG2-1)*4+1) = NBEDG2
C              Coordonnees
               RMCN(MCXED2+(NBEDG2-1)*2  ) = REAL( INTEXY(1) )
               RMCN(MCXED2+(NBEDG2-1)*2+1) = REAL( INTEXY(2) )
C              Numero du point ajoute pour le lien
               LIEN2 = NBEDG2
CCC            WRITE (IMPRIM,*) 'NEW2',LIEN2
            ELSE
C
C              Si l'intersection est sur une extremite, on n'ajoute
C              pas de point
               IF (INTERS(2) .LT. 0.5D0) THEN
                  LIEN2 = ISEG2
CCC               WRITE (IMPRIM,*) 'OLD2',ISEG2
               ELSE
                  LIEN2 = ISEG2S
CCC               WRITE (IMPRIM,*) 'OLD2',ISEG2S
               ENDIF
            ENDIF
C
C           Mise en place des liens entre les courbes
            IF (MCN(MCEDG1+(LIEN1-1)*4+3) .NE. LIEN2) THEN
               IF (MCN(MCEDG1+(LIEN1-1)*4+3) .EQ. 0) THEN
                  MCN(MCEDG1+(LIEN1-1)*4+3) = LIEN2
               ELSE
CCC               WRITE (IMPRIM,*) RMCN(MCXED1+(LIEN1-1)*2),
CCC  &                             RMCN(MCXED1+(LIEN1-1)*2+1)
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'IL Y A PLUSIEURS INTERSECTIONS CONFONDUES'
C                  CALL LEREUR
                  IERR = 9
               ENDIF
            ENDIF
            IF (MCN(MCEDG2+(LIEN2-1)*4+3) .NE. LIEN1) THEN
               IF (MCN(MCEDG2+(LIEN2-1)*4+3) .EQ. 0) THEN
                  MCN(MCEDG2+(LIEN2-1)*4+3) = LIEN1
               ELSE
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'IL Y A PLUSIEURS INTERSECTIONS CONFONDUES'
C                  CALL LEREUR
                  IERR = 10
               ENDIF
            ENDIF
C
C           Si on a insere un point sur un segment de la courbe 1, il faut
C           mettre a jour l'extremite du segment
            IF ((INTERS(1) .GT. 0.05) .AND. (INTERS(1) .LT. 0.95)) THEN
C              L'indice du segment suivant = ISEG1S
               ISEG1S = LIEN1
               P12X = RMCN(MCXED1 + (ISEG1S-1)*2)
               P12Y = RMCN(MCXED1 + (ISEG1S-1)*2 + 1)
C              Le vecteur P11P12 :
               P11P12(1) = P12X - P11X
               P11P12(2) = P12Y - P11Y
            ENDIF
C
C           On cherche le point suivant sur la courbe,
C           si c'est le point de depart, on s'arrete.
 400        P21X = P22X
            P21Y = P22Y
            P11P21(1) = P11P22(1)
            P11P21(2) = P11P22(2)
            SIGNE1 = SIGNE2
            ISEG2  = ISEG2S
         IF (ISEG2 .NE. MCN(MCCPC2 + ICONT)) GOTO 200
C
C        Fin de boucle sur les segments du contour ICONT
C
 410  CONTINUE
C
C     Fin de la boucle sur les segments de la surface 1
C
      ISEG1 = ISEG1 + 1
CCC   WRITE(IMPRIM,*) 'iseg1 = ',iseg1,'nbedg1=',nbedg1
      IF (ISEG1 .LE. NBEDG1) GOTO 100
C
C     --------------------------------------------------
C     PHASE 3 : contruction des courbes de l'objet final
C     On construit d'abord les courbes qui contiennent
C     un point d'intersection entre les deux surfaces.
C
C     Allocation d'un tableau temporaire pour stocker les courbes
C     en cours de construction. La taille du tableau est egale
C     au nombre total des points sur les contours des deux surfaces.
      CALL TNMCDC ('ENTIER', NBEDG1 + NBEDG2, MCNCOU)
C
C     La courbe exterieure est determinee en cherchant celle
C     qui porte le point dont la coordonnee X est minimale.
      NUCEXT = 0
      CXMIN = 1.0E30
C
C     Recherche d'un point d'intersection non encore utilise
      NBCOUR = 0
C     Boucle sur les courbes de l'objet final
 500     ISEG1 = 0
 510     ISEG1 = ISEG1 + 1
         IF (ISEG1 .GT. NBEDG1) GOTO 900
         IF (MCN(MCEDG1+(ISEG1-1)*4+3) .EQ. 0) GOTO 510
         IF (MCN(MCEDG1+(ISEG1-1)*4) .EQ. -1) GOTO 510
CCC      WRITE(*,*) 'point de depart',iseg1
         NBCOUR = NBCOUR + 1
         NBPTS = 0
         ISEG2 = MCN(MCEDG1+(ISEG1-1)*4+3)
C        On efface met -1 dans la premiere colonne
         MCN(MCEDG1+(ISEG1-1)*4) = -1
         MCN(MCEDG2+(ISEG2-1)*4) = -1
C        On memorise le point de la courbe d'ou on part.
C        On s'arretera en arrivant a nouveau sur ce point.
         PTDEP1 = ISEG1
         PTDEP2 = ISEG2
C
C        On est sur un point d'intersection. Selon la position relative
C        des courbes, et l'operation demandee (union, intersection, ...)
C        on determine laquelle des deux courbes on suit.
C        Pour l'intersection, on suite la courbe qui part vers
C        l'interieur de l'autre surface. Pour l'union, on suit la courbe
C        qui part vers l'exterieur.
C        La soustraction est traitee comme une intersection avec le
C        complementaire de la deuxieme surface. Ceci revient exactement
C        a inverser le sens de parcours des courbes (au lieu de prendre
C        le point suivant, on prend le point precedent).
C
C        Une fois determinee la courbe que l'on doit suivre, il suffit
C        de construire la nouvelle courbe en prenant les points dans
C        l'ordre. Si on rencontre un point d'intersection entre les
C        contours des deux surfaces, il faut changer de courbe et suivre
C        la nouvelle courbe dans le sens de parcours qui lui est associe.
C        Ainsi, on construit un contours de la surface demandee.
         P11X = RMCN(MCXED1+(ISEG1-1)*2)
         P11Y = RMCN(MCXED1+(ISEG1-1)*2+1)
C        Coordonnees du point suivant sur la surface 1
         ISEG1S = MCN(MCEDG1+(ISEG1-1)*4+1)
         P12X = RMCN(MCXED1+(ISEG1S-1)*2)
         P12Y = RMCN(MCXED1+(ISEG1S-1)*2+1)
C        Coordonnees du point precedent
         ISEG1P = MCN(MCEDG1+(ISEG1-1)*4+2)
         P13X = RMCN(MCXED1+(ISEG1P-1)*2)
         P13Y = RMCN(MCXED1+(ISEG1P-1)*2+1)
C        Coordonnees du point suivant sur l'autre surface
         ISEG2 = MCN(MCEDG1+(ISEG1-1)*4+3)
         ISEG2S = MCN(MCEDG2+(ISEG2-1)*4+1)
         P22X = RMCN(MCXED2+(ISEG2S-1)*2)
         P22Y = RMCN(MCXED2+(ISEG2S-1)*2+1)
C        Vecteurs dont on a besoin
         P11P12(1) = P12X - P11X
         P11P12(2) = P12Y - P11Y
         P11P13(1) = P13X - P11X
         P11P13(2) = P13Y - P11Y
         P11P22(1) = P22X - P11X
         P11P22(2) = P22Y - P11Y
C
C        Est-ce que le vecteur P11P22 est a l'interieur du secteur angulaire
C        defini par les vecteurs P11P12 et P11P13 ?
         DETERM(1) = P11P12(1) * P11P22(2) - P11P12(2) * P11P22(1)
         DETERM(2) = P11P13(1) * P11P22(2) - P11P13(2) * P11P22(1)
         DETERM(3) = P11P12(1) * P11P13(2) - P11P12(2) * P11P13(1)
         IF (DETERM(3) .GE. 0.0D0) THEN
            IF ((DETERM(1) .GT. 0.0D0) .AND. (DETERM(2) .LT. 0.0D0))THEN
               SIGNE1 = 1
            ELSE
               SIGNE1 = 0
            ENDIF
         ELSE
            IF ((DETERM(1) .GT. 0.0D0) .OR. (DETERM(2) .LT. 0.0D0)) THEN
               SIGNE1 = 1
            ELSE
               SIGNE1 = 0
            ENDIF
         ENDIF
C
C        Ici, on regarde quelle est l'operation demandee
         IF (NUOPER .EQ. 1) THEN
C           On demande une intersection : on choisit la courbe qui
C           part a l'interieur de l'autre surface.
            IF (SIGNE1 .EQ. 1) THEN
               ASUIVR = 2
            ELSE
               ASUIVR = 1
            ENDIF
         ELSE IF (NUOPER .EQ. 2) THEN
C           On demande une union : on choisit la courbe qui part
C           a l'exterieur de l'autre surface.
            IF (SIGNE1 .EQ. 1) THEN
               ASUIVR = 1
            ELSE
               ASUIVR = 2
            ENDIF
         ELSE
C           On demande une soustraction : on choisit la courbe qui part
C           a l'interieur de l'autre surface, l'interieur et l'exterieur
C           de la deuxieme surface etant echanges.
            IF (SIGNE1 .EQ. 1) THEN
               ASUIVR = 1
            ELSE
               ASUIVR = 2
            ENDIF
         ENDIF
C
C        BOUCLE DE PARCOURS DES POINTS D'UNE COURBE
C        A mesure qu'on trouve des intersections on met -1 dans la
C        premiere colonne pour montrer qu'ils ont ete parcourus.
 530     NBPTS = NBPTS + 1
         IF (ASUIVR .EQ. 1) THEN
            MCN(MCNCOU+(NBPTS-1)) = MCXED1+(ISEG1-1)*2
            MCN(MCEDG1+(ISEG1-1)*4) = -1
         ELSE
            MCN(MCNCOU+(NBPTS-1)) = MCXED2+(ISEG2-1)*2
            MCN(MCEDG2+(ISEG2-1)*4) = -1
         ENDIF
CCC      WRITE (IMPRIM,*) 'ajout point',asuivr,iseg1,iseg2
C
C        On cherche le point suivant sur la courbe qu'on est
C        en train de suivre. Si on fait une soustraction et qu'on suit
C        la courbe 2, on prend le point "precedent" puisqu'on inverse
C        le sens de parcours des courbes
         IF (ASUIVR .EQ. 1) THEN
           ISEG1 = MCN(MCEDG1+(ISEG1-1)*4+1)
C          Si il y a une intersection avec l'autre surface, on change
C          de courbe.
           IF (MCN(MCEDG1+(ISEG1-1)*4+3) .NE. 0) THEN
              MCN(MCEDG1+(ISEG1-1)*4) = -1
              ISEG2 = MCN(MCEDG1+(ISEG1-1)*4+3)
CCC           WRITE (IMPRIM,*) 'intersection : ',iseg1,iseg2
              ASUIVR = 2
           ENDIF
        ELSE
           IF (NUOPER .NE. 3) THEN
              ISEG2 = MCN(MCEDG2+(ISEG2-1)*4+1)
           ELSE
              ISEG2 = MCN(MCEDG2+(ISEG2-1)*4+2)
           ENDIF
C          Si il y a une intersection avec l'autre surface, on change
C          de courbe.
           IF (MCN(MCEDG2+(ISEG2-1)*4+3) .NE. 0) THEN
              MCN(MCEDG2+(ISEG2-1)*4) = -1
              ISEG1 = MCN(MCEDG2+(ISEG2-1)*4+3)
CCC           WRITE (IMPRIM,*) 'intersection : ',iseg1,iseg2
              ASUIVR = 1
           ENDIF
        ENDIF
C
C       Si on est a nouveau au point de depart, on s'arrete.
        IF ((ISEG1 .EQ. PTDEP1).AND.(ISEG2 .EQ. PTDEP2)) GOTO 550
        GOTO 530
C       FIN DE LA BOUCLE DE PARCOURS D'UNE COURBE
C
C       Ajout de la courbe en tant que LIGNE dans le TMS et le lexique
 550    CONTINUE
        WRITE (NUMERO,'(I4)') NBCOUR
        NBC = NBCHIF( NBCOUR )
        NOMCOU = 'LIGNE_POUR_SUEX45___' // NUMERO(5-NBC:4)
        WRITE (IMPRIM,*) 'nom de la courbe: ',NOMCOU
        CALL LXLXOU (NTLIGN, NOMCOU, LXCOUR, MCCOUR)
        IF( LXCOUR .GT. 0 ) THEN
C           LA LIGNE EXISTE => ELLE EST DETRUITE
           CALL LXLXDS( NTLIGN, NOMCOU )
        ENDIF
        CALL LXLXDC (NTLIGN, NOMCOU, 4, 24)
        CALL LXLXOU (NTLIGN, NOMCOU, LXCOUR, MCCOUR)
C
        CALL LXTSOU (LXCOUR, 'DEFINITION', NTNSEF, MCNSEF)
        IF( NTNSEF .GT. 0 ) THEN
           CALL LXTSDS( LXCOUR, 'DEFINITION' )
        ENDIF
        CALL LXTNDC (LXCOUR, 'DEFINITION', 'ENTIER', WUTYLI)
        CALL LXTSOU (LXCOUR, 'DEFINITION', NTNSEF, MCNSEF)
        MCN(MCNSEF+WTYTRL) = 1
        MCN(MCNSEF+WUTYLI) = 10
        CALL ECDATE (MCN(MCNSEF))
        MCN(MCNSEF+MOTVAR(6)) = NONMTD('~>LIGNE>>DEFINITION')
        I = WUSOEF + 2 * NBPTS
        CALL LXTSOU (LXCOUR, 'NSEF', NTNSEF, MCNSEF)
        IF( NTNSEF .GT. 0 ) THEN
           CALL LXTSDS( LXCOUR, 'NSEF' )
        ENDIF
        CALL LXTNDC (LXCOUR, 'NSEF', 'ENTIER', I)
        CALL LXTSOU (LXCOUR, 'NSEF', NTNSEF, MCNSEF)
        CALL ECDATE (MCN(MCNSEF))
        MCN(MCNSEF+MOTVAR(6)) = NONMTD('~>>>NSEF')
        MCN(MCNSEF+WUTYOB) = 2
        MCN(MCNSEF+WUTFMA) = 1
        MCN(MCNSEF+WBSOEF) = 2
        MCN(MCNSEF+WBTGEF) = 0
        MCN(MCNSEF+WBEFOB) = NBPTS
        MCN(MCNSEF+WBEFTG) = 0
        MCN(MCNSEF+WBEFAP) = 0
        MCN(MCNSEF+WUTYMA) = 0
        DO I = 0, NBPTS - 1
           MCN(MCNSEF+WUSOEF+I*2  ) = I + 1
           MCN(MCNSEF+WUSOEF+I*2+1) = I + 2
        ENDDO
        MCN(MCNSEF+WUSOEF+NBPTS*2-1) = 1
        CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
        IF( NTXYZS .GT. 0 ) THEN
           CALL LXTSDS( LXCOUR, 'XYZSOMMET' )
        ENDIF
        I = WUSOEF + NBPTS * 3
        CALL LXTNDC (LXCOUR, 'XYZSOMMET', 'MOTS', I)
        CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
        CALL ECDATE (MCN(MCXYZS))
        MCN(MCXYZS+MOTVAR(6)) = NONMTD('~>>>XYZSOMMET')
        MCN(MCXYZS+WNBSOM) = NBPTS
        MCN(MCXYZS+WNBTGS) = 0
        MCN(MCXYZS+WBCOOR) = 3
        DO I = 0, NBPTS - 1
           RMCN(MCXYZS+WYZSOM+I*3  ) = RMCN(MCN(MCNCOU+I))
           RMCN(MCXYZS+WYZSOM+I*3+1) = RMCN(MCN(MCNCOU+I)+1)
           RMCN(MCXYZS+WYZSOM+I*3+2) = 0.0E0
           IF (RMCN(MCXYZS+WYZSOM+I*3) .LT. CXMIN) THEN
              CXMIN = RMCN(MCXYZS+WYZSOM+I*3)
              NUCEXT = NBCOUR
           ENDIF
        ENDDO
        CALL LXLXFE (NTLIGN, NOMCOU)
C
CCC     WRITE (IMPRIM,*) 'Ajout de la courbe :',NBCOUR
CCC     WRITE (IMPRIM,*) 'Nb de points :',NBPTS
        DO I=0, NBPTS-1
           ISEG1 = MCN(MCNCOU+I)
CCC        WRITE (IMPRIM,*) RMCN(ISEG1),RMCN(ISEG1+1)
        ENDDO
C
C     OK, la courbe est dans le TMS, on passe a la courbe suivante
      GOTO 500
C
C     Fin du traitement des courbes contenant un point d'intersection
C     On passe maintenant aux courbes n'ayant aucun point d'intersection.
C     Ces courbes sont soit ignorees, soit copiees integralement dans
C     le nouvel objet. Si on fait une intersection, une courbe est copiee
C     si elle est a l'interieur de l'autre surface.
C     Si on fait une union, une courbe est copiee si elle est a l'exterieur
C     de l'autre surface.
C     Si on fait une soustraction, idem que l'intersection avec interieur
C     et exterieur echanges pour la deuxieme surface.
C
C     Boucle sur les courbes de la surface 1
 900  DO 950 K = 0, NBCPC1 - 1
C
C        Est-ce que la courbe contient un point d'intersection ?
         PTDEP1 = MCN(MCCPC1 + K)
         ISEG1 = PTDEP1
         NBPTS = 0
C        Si on trouve un lien, on passe a la courbe suivante
 910     IF (MCN(MCEDG1+(ISEG1-1)*4+3) .NE. 0) GOTO 950
         ISEG1 = MCN(MCEDG1+(ISEG1-1)*4+1)
         NBPTS = NBPTS + 1
         IF (ISEG1 .NE. PTDEP1) GOTO 910
C        On a trouve une courbe sans point d'intersection
C        Est ce qu'elle est a l'interieur de la surface 2 ?
         INSIDE = 0
C        On initialise une direction aleatoire
         P11P12(1) = 0.3736528942546315896124575D0
         P11P12(2) = 0.8747551321654278954651324D0
C        On prend un point sur la courbe 1
         P11X = RMCN(MCXED1+(PTDEP1-1)*2)
         P11Y = RMCN(MCXED1+(PTDEP1-1)*2+1)
C
C        On boucle sur les courbes de la surface 2.
         DO J = 0, NBCPC2 - 1
C
            PTDEP2 = MCN(MCCPC2 + J)
            ISEG2 = PTDEP2
            P21X = RMCN(MCXED2+(ISEG2-1)*2)
            P21Y = RMCN(MCXED2+(ISEG2-1)*2+1)
            P11P21(1) = P21X - P11X
            P11P21(2) = P21Y - P11Y
            DETERM(1) = P11P12(1) * P11P21(2) - P11P12(2) * P11P21(1)
            IF (DETERM(1) .LT. 0.0D0) THEN
               SIGNE1 = -1
            ELSE
               SIGNE1 = 1
            ENDIF
 920        ISEG2 = MCN(MCEDG2+(ISEG2-1)*4+1)
            P22X = RMCN(MCXED2+(ISEG2-1)*2)
            P22Y = RMCN(MCXED2+(ISEG2-1)*2+1)
            P11P22(1) = P22X - P11X
            P11P22(2) = P22Y - P11Y
            DETERM(2) = P11P12(1) * P11P22(2) - P11P12(2) * P11P22(1)
            IF (DETERM(2) .LT. 0.0D0) THEN
               SIGNE2 = -1
            ELSE
               SIGNE2 = 1
            ENDIF
C
C           Si le signe du determinant change, on a une chance d'avoir
C           une intersection. On teste la suite...
            IF (SIGNE1 .NE. SIGNE2) THEN
               P21P22(1) = P22X - P21X
               P21P22(2) = P22Y - P21Y
               DETERM(3) = P11P21(1) * P21P22(2) - P11P21(2) * P21P22(1)
C              Si le segment coupe la demi droite, on incremente le nombre
C              d'intersections.
               IF (DETERM(3) * SIGNE1 .LE. 0.0D0) THEN
                  INSIDE = INSIDE + 1
               ENDIF
            ENDIF
            P21X = P22X
            P21Y = P22Y
            P11P21(1) = P11P22(1)
            P11P21(2) = P11P22(2)
            SIGNE1 = SIGNE2
            IF (ISEG2 .NE. PTDEP2) GOTO 920
C
         ENDDO
C
C        Fin de la boucle sur les courbes de la surface 2.
C        Si le nombre d'intersections est impair, la courbe 1
C        est a l'interieur de la surface 2.
         IF (MOD(INSIDE,2) .NE. 0) THEN
            INSIDE = 1
         ELSE
            INSIDE = 0
         ENDIF
C
C        Intersection : on garde la courbe si elle est a l'interieur
C        Union : on garde la courbe si elle est a l'exterieur
C        Soustraction : on garde la courbe si elle est a l'exterieur
         IF (     ((NUOPER .EQ. 1) .AND. (INSIDE .EQ. 1))
     %       .OR. ((NUOPER .EQ. 2) .AND. (INSIDE .EQ. 0))
     %       .OR. ((NUOPER .EQ. 3) .AND. (INSIDE .EQ. 0))) THEN
            TOKEEP = 1
         ELSE
            TOKEEP = 0
         ENDIF
         IF (TOKEEP .EQ. 0) GOTO 950
C
C        Si on garde la courbe, on l'ajoute dans la liste
         NBCOUR = NBCOUR + 1
         WRITE (NUMERO,'(I4)') NBCOUR
         NBC = NBCHIF( NBCOUR )
         NOMCOU = 'LIGNE_POUR_SUEX45___' // NUMERO(5-NBC:4)
         WRITE (IMPRIM,*) 'nom de la courbe: ',NOMCOU
         CALL LXLXOU (NTLIGN, NOMCOU, LXCOUR, MCCOUR)
         IF( LXCOUR .GT. 0 ) THEN
            CALL LXLXDS( NTLIGN, NOMCOU )
         ENDIF
         CALL LXLXDC (NTLIGN, NOMCOU, 4, 24)
         CALL LXLXOU (NTLIGN, NOMCOU, LXCOUR, MCCOUR)
         CALL LXTSOU (LXCOUR, 'DEFINITION', NTNSEF, MCNSEF)
         IF( NTNSEF .GT. 0 ) THEN
            CALL LXTSDS( LXCOUR, 'DEFINITION' )
         ENDIF
         CALL LXTNDC (LXCOUR, 'DEFINITION', 'ENTIER', WUTYLI)
         CALL LXTSOU (LXCOUR, 'DEFINITION', NTNSEF, MCNSEF)
         MCN(MCNSEF+WTYTRL) = 1
         MCN(MCNSEF+WUTYLI) = 10
         CALL ECDATE (MCN(MCNSEF))
         MCN(MCNSEF+MOTVAR(6)) = NONMTD('~>LIGNE>>DEFINITION')
C
         CALL LXTSOU (LXCOUR, 'NSEF', NTNSEF, MCNSEF)
         IF( NTNSEF .GT. 0 ) THEN
            CALL LXTSDS( LXCOUR, 'NSEF' )
         ENDIF
         I = WUSOEF + 2 * NBPTS
         CALL LXTNDC (LXCOUR, 'NSEF', 'ENTIER', I)
         CALL LXTSOU (LXCOUR, 'NSEF', NTNSEF, MCNSEF)
         CALL ECDATE (MCN(MCNSEF))
         MCN(MCNSEF+MOTVAR(6)) = NONMTD('~>>>NSEF')
         MCN(MCNSEF+WUTYOB) = 2
         MCN(MCNSEF+WUTFMA) = 1
         MCN(MCNSEF+WBSOEF) = 2
         MCN(MCNSEF+WBTGEF) = 0
         MCN(MCNSEF+WBEFOB) = NBPTS
         MCN(MCNSEF+WBEFTG) = 0
         MCN(MCNSEF+WBEFAP) = 0
         MCN(MCNSEF+WUTYMA) = 0
         DO I = 0, NBPTS - 1
            MCN(MCNSEF+WUSOEF+I*2) = I + 1
            MCN(MCNSEF+WUSOEF+I*2+1) = I + 2
         ENDDO
         MCN(MCNSEF+WUSOEF+NBPTS*2-1) = 1
C
         CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
         IF( NTXYZS .GT. 0 ) THEN
            CALL LXTSDS( LXCOUR, 'XYZSOMMET' )
         ENDIF
         I = WUSOEF + NBPTS * 3
         CALL LXTNDC (LXCOUR, 'XYZSOMMET', 'MOTS', I)
         CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
         CALL ECDATE (MCN(MCXYZS))
         MCN(MCXYZS+MOTVAR(6)) = NONMTD('~>>>XYZSOMMET')
         MCN(MCXYZS+WNBSOM) = NBPTS
         MCN(MCXYZS+WNBTGS) = 0
         ISEG1 = PTDEP1
         DO I = 0, NBPTS - 1
            RMCN(MCXYZS+WYZSOM+I*3) = RMCN(MCXED1+(ISEG1-1)*2)
            RMCN(MCXYZS+WYZSOM+I*3+1) = RMCN(MCXED1+(ISEG1-1)*2+1)
            RMCN(MCXYZS+WYZSOM+I*3+2) = 0.0E0
            ISEG1 = MCN(MCEDG1+(ISEG1-1)*4+1)
            IF (RMCN(MCXYZS+WYZSOM+I*3) .LT. CXMIN) THEN
               CXMIN = RMCN(MCXYZS+WYZSOM+I*3)
               NUCEXT = NBCOUR
            ENDIF
         ENDDO
         CALL LXLXFE (NTLIGN, NOMCOU)
 950  CONTINUE
C
C     Boucle sur les courbes de la surface 2
      DO 990 K = 0, NBCPC2 - 1
C
C        Est-ce que la courbe contient un point d'intersection ?
         PTDEP2 = MCN(MCCPC2 + K)
         ISEG2 = PTDEP2
         NBPTS = 0
C        Si on trouve un lien, on passe a la courbe suivante
 970     IF (MCN(MCEDG2+(ISEG2-1)*4+3) .NE. 0) GOTO 990
         ISEG2 = MCN(MCEDG2+(ISEG2-1)*4+1)
         NBPTS = NBPTS + 1
         IF (ISEG2 .NE. PTDEP2) GOTO 970
C        On a trouve une courbe sans point d'intersection
C        Est ce qu'elle est a l'interieur de la surface 1 ?
         INSIDE = 0
C        On initialise une direction aleatoire
         P11P12(1) = 0.3736528942546315896124575D0
         P11P12(2) = 0.8747551321654278954651324D0
C        On prend un point sur la courbe 2
         P11X = RMCN(MCXED2+(PTDEP2-1)*2)
         P11Y = RMCN(MCXED2+(PTDEP2-1)*2+1)
C
C        On boucle sur les courbes de la surface 1.
         DO J = 0, NBCPC1 - 1
            PTDEP1 = MCN(MCCPC1 + J)
            ISEG1 = PTDEP1
            P21X = RMCN(MCXED1+(ISEG1-1)*2)
            P21Y = RMCN(MCXED1+(ISEG1-1)*2+1)
            P11P21(1) = P21X - P11X
            P11P21(2) = P21Y - P11Y
            DETERM(1) = P11P12(1) * P11P21(2) - P11P12(2) * P11P21(1)
            IF (DETERM(1) .LT. 0.0D0) THEN
               SIGNE1 = -1
            ELSE
               SIGNE1 = 1
            ENDIF
 980        ISEG1 = MCN(MCEDG1+(ISEG1-1)*4+1)
            P22X = RMCN(MCXED1+(ISEG1-1)*2)
            P22Y = RMCN(MCXED1+(ISEG1-1)*2+1)
            P11P22(1) = P22X - P11X
            P11P22(2) = P22Y - P11Y
            DETERM(2) = P11P12(1) * P11P22(2) - P11P12(2) * P11P22(1)
            IF (DETERM(2) .LT. 0.0D0) THEN
               SIGNE2 = -1
            ELSE
               SIGNE2 = 1
            ENDIF
C
C           Si le signe du determinant change, on a une chance d'avoir
C           une intersection. On teste la suite...
            IF (SIGNE1 .NE. SIGNE2) THEN
               P21P22(1) = P22X - P21X
               P21P22(2) = P22Y - P21Y
               DETERM(3) = P11P21(1) * P21P22(2) - P11P21(2) * P21P22(1)
C              Si le segment coupe la demi droite, on incremente le nombre
C              d'intersections.
               IF (DETERM(3) * SIGNE1 .LE. 0.0D0) THEN
                  INSIDE = INSIDE + 1
CCC               WRITE (IMPRIM,*) 'intersection !!!',inside
               ENDIF
            ENDIF
            P21X = P22X
            P21Y = P22Y
            P11P21(1) = P11P22(1)
            P11P21(2) = P11P22(2)
            SIGNE1 = SIGNE2
            IF (ISEG1 .NE. PTDEP1) GOTO 980
         ENDDO
C
C        Fin de la boucle sur les courbes de la surface 2.
C        Si le nombre d'intersections est impair, la courbe 1
C        est a l'interieur de la surface 2.
         IF (MOD(INSIDE,2) .NE. 0) THEN
            INSIDE = 1
         ELSE
            INSIDE = 0
         ENDIF
C
C        Intersection : on garde la courbe si elle est a l'interieur
C        Union : on garde la courbe si elle est a l'exterieur
C        Soustraction : on garde la courbe si elle est a l'interieur
         IF (     ((NUOPER .EQ. 1) .AND. (INSIDE .EQ. 1))
     %       .OR. ((NUOPER .EQ. 2) .AND. (INSIDE .EQ. 0))
     %       .OR. ((NUOPER .EQ. 3) .AND. (INSIDE .EQ. 1))) THEN
            TOKEEP = 1
         ELSE
            TOKEEP = 0
         ENDIF
         IF (TOKEEP .EQ. 0) GOTO 990
C
C        Si on garde la courbe, on l'ajoute dans la liste
         NBCOUR = NBCOUR + 1
         WRITE (NUMERO,'(I4)') NBCOUR
         NBC = NBCHIF( NBCOUR )
         NOMCOU = 'LIGNE_POUR_SUEX45___' // NUMERO(5-NBC:4)
         WRITE (IMPRIM,*) 'nom de la courbe: ',NOMCOU
         CALL LXLXOU (NTLIGN, NOMCOU, LXCOUR, MCCOUR)
         IF( LXCOUR .GT. 0 ) THEN
            CALL LXLXDS( NTLIGN, NOMCOU )
         ENDIF
         CALL LXLXDC (NTLIGN, NOMCOU, 4, 24)
         CALL LXLXOU (NTLIGN, NOMCOU, LXCOUR, MCCOUR)
C
         CALL LXTSOU (LXCOUR, 'DEFINITION', NTNSEF, MCNSEF)
         IF( NTNSEF .GT. 0 ) THEN
            CALL LXTSDS( LXCOUR, 'DEFINITION' )
         ENDIF
         CALL LXTNDC (LXCOUR, 'DEFINITION', 'ENTIER', WUTYLI)
         CALL LXTSOU (LXCOUR, 'DEFINITION', NTNSEF, MCNSEF)
         MCN(MCNSEF+WTYTRL) = 1
         MCN(MCNSEF+WUTYLI) = 10
         CALL ECDATE (MCN(MCNSEF))
         MCN(MCNSEF+MOTVAR(6)) = NONMTD('~>LIGNE>>DEFINITION')
C
         CALL LXTSOU (LXCOUR, 'NSEF', NTNSEF, MCNSEF)
         IF( NTNSEF .GT. 0 ) THEN
            CALL LXTSDS( LXCOUR, 'NSEF' )
         ENDIF
         I = WUSOEF + 2 * NBPTS
         CALL LXTNDC (LXCOUR, 'NSEF', 'ENTIER', I)
         CALL LXTSOU (LXCOUR, 'NSEF', NTNSEF, MCNSEF)
         CALL ECDATE (MCN(MCNSEF))
         MCN(MCNSEF+MOTVAR(6)) = NONMTD('~>>>NSEF')
         MCN(MCNSEF+WUTYOB) = 2
         MCN(MCNSEF+WUTFMA) = 1
         MCN(MCNSEF+WBSOEF) = 2
         MCN(MCNSEF+WBTGEF) = 0
         MCN(MCNSEF+WBEFOB) = NBPTS
         MCN(MCNSEF+WBEFTG) = 0
         MCN(MCNSEF+WBEFAP) = 0
         MCN(MCNSEF+WUTYMA) = 0
         DO I = 0, NBPTS - 1
            MCN(MCNSEF+WUSOEF+I*2) = I + 1
            MCN(MCNSEF+WUSOEF+I*2+1) = I + 2
         ENDDO
         MCN(MCNSEF+WUSOEF+NBPTS*2-1) = 1
C
         CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
         IF( NTXYZS .GT. 0 ) THEN
            CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
         ENDIF
         I = WUSOEF + NBPTS * 3
         CALL LXTNDC (LXCOUR, 'XYZSOMMET', 'MOTS', I)
         CALL LXTSOU (LXCOUR, 'XYZSOMMET', NTXYZS, MCXYZS)
         CALL ECDATE (MCN(MCXYZS))
         MCN(MCXYZS+MOTVAR(6)) = NONMTD('~>>>XYZSOMMET')
         MCN(MCXYZS+WNBSOM) = NBPTS
         MCN(MCXYZS+WNBTGS) = 0
         ISEG2 = PTDEP2
         DO I = 0, NBPTS - 1
            RMCN(MCXYZS+WYZSOM+I*3) = RMCN(MCXED2+(ISEG2-1)*2)
            RMCN(MCXYZS+WYZSOM+I*3+1) = RMCN(MCXED2+(ISEG2-1)*2+1)
            RMCN(MCXYZS+WYZSOM+I*3+2) = 0.0E0
            ISEG2 = MCN(MCEDG2+(ISEG2-1)*4+1)
            IF (RMCN(MCXYZS+WYZSOM+I*3) .LT. CXMIN) THEN
               CXMIN = RMCN(MCXYZS+WYZSOM+I*3)
               NUCEXT = NBCOUR
            ENDIF
         ENDDO
         CALL LXLXFE (NTLIGN, NOMCOU)
C
 990  CONTINUE
C
C     Liberation de la memoire temporaire utilisee
      CALL TNMCDS ('ENTIER', EDGSZ1 * 4, MCEDG1)
      CALL TNMCDS ('ENTIER', EDGSZ2 * 4, MCEDG2)
      CALL TNMCDS ('REEL',   EDGSZ1 * 2, MCXED1)
      CALL TNMCDS ('REEL',   EDGSZ2 * 2, MCXED2)
      CALL TNMCDS ('ENTIER', NBCPC1,     MCCPC1)
      CALL TNMCDS ('ENTIER', NBCPC2,     MCCPC2)
      CALL TNMCDS ('ENTIER', NBEDG1 + NBEDG2, MCNCOU)
C
C     APPEL DE SUEX09 POUR TRIANGULATION
C     ==================================
      CALL TNMCDC ('ENTIER', WULFTR + NBCOUR, MCLADS)
C     Taille ideale TRES GRAND => Choix de la longueur de l'arete maximale
      RMCN(MCLADS+WRETMX) = RINFO( 'GRAND' )
C     Nombre de lignes
      MCN(MCLADS+WBLFTR) = NBCOUR
C     Nombre de points internes futurs sommets
      MCN(MCLADS+WBPTIT) = 0
C
C     Numeros TMS des lignes fermees
      DO I = 1, NBCOUR
         WRITE (NUMERO,'(I4)') I
         NBC = NBCHIF( I )
         NOMCOU = 'LIGNE_POUR_SUEX45___' // NUMERO(5-NBC:4)
         CALL LXNMNO (NTLIGN, NOMCOU, LXCOUR, J)
C        LXCOUR EST LE NUMERO DE LA LIGNE DANS LE LEXIQUE DES LIGNES
         IF( LXCOUR .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = 'LIGNE VIDE: ' // NOMCOU
            KERR(2) = 'ABANDON DE LA TRIANGULATION 45'
            CALL LEREUR
            IERR = 11
            GOTO 1100
         ENDIF
         MCN(MCLADS+WULFTR+I-1) = LXCOUR
      ENDDO
C
C     On met la courbe exterieure en premier
      I = MCN(MCLADS+WULFTR+NUCEXT-1)
      MCN(MCLADS+WULFTR+NUCEXT-1) = MCN(MCLADS+WULFTR)
      MCN(MCLADS+WULFTR) = I
      WRITE (IMPRIM,*) 'NUMERO COURBE EXTERIEURE', NUCEXT
C
C     Appel du sous-programme SUEX09
      J = 9
      CALL SUEX09 ( J, NTLXSU, MCN(MCLADS), RMCN(MCLADS),
     &              NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'VERIFIER QUE LES 2 SURFACES'
         KERR(2) = 'SONT DANS UN MEME PLAN'
         CALL LEREUR
      ENDIF
C
C     On libere la memoire du TMC LADEFI
 1100 CALL TNMCDS ('ENTIER', WULFTR + NBCOUR, MCLADS)
C     On efface les TMS des lignes de contour
      DO I = 1, NBCOUR
         WRITE (NUMERO,'(I4)') I
         NBC    = NBCHIF( I )
         NOMCOU = 'LIGNE_POUR_SUEX45___' // NUMERO(5-NBC:4)
         CALL LXNMNO ( NTLIGN, NOMCOU, LXCOUR, J )
         IF( LXCOUR .GT. 0 ) THEN
            CALL LXLXDS ( NTLIGN, NOMCOU )
         ENDIF
      ENDDO
C
C     Et c'est fini !
 9999 RETURN
      END
