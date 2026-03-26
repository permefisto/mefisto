      SUBROUTINE TRPRTR( NBARPI, MXSOMM, NBSOMM, PXYD,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, MXARTR, N1ARTR, NOARTR, NOARST,
     %                   IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRIANGULER LES POINTS DE LA FRONTIERE ET INTERNES IMPOSES
C -----    PAR SUBDIVISION DU OU DES TRIANGLES CONTENANT LE POINT
C          ET MISE EN DELAUNAY DES TRIANGLES
C
C ENTREES:
C --------
C NBARPI : NOMBRE DE SOMMETS DE LA FRONTIERE + NOMBRE DE POINTS INTERNES
C          IMPOSES PAR L'UTILISATEUR
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS PXYD
C NBSOMM : NOMBRE DE SOMMETS UTILISES (NBARPI + 4)
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOSOAR
C MXSOAR : NOMBRE MAXIMAL D'ARETES STOCKABLES DANS LE TABLEAU NOSOAR
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C MXARTR : NOMBRE MAXIMAL DE TRIANGLES STOCKABLES DANS LE TABLEAU NOARTR
C
C MODIFIES:
C ---------
C N1SOAR : NUMERO DE LA PREMIERE ARETE VIDE DANS LE TABLEAU NOSOAR
C          UNE ARETE I DE NOSOAR EST VIDE  <=>  NOSOAR(1,I)=0
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          HACHAGE DES ARETES = NOSOAR(1)+NOSOAR(2)*2
C NOARST : NOARST(I) NUMERO DANS NOSOAR D'UNE ARETE DE SOMMET I
C
C SORTIES:
C --------
C N1ARTR : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOARTR
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOARTR(2,.)
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1 = 0 SI TRIANGLE VIDE => ARETE2 = TRIANGLE VIDE SUIVANT
C IERR   : =0 SI PAS D'ERREUR
C          =1 SI LE TABLEAU NOSOAR EST SATURE
C          =2 SI LE TABLEAU NOARTR EST SATURE
C          =3 SI AUCUN DES TRIANGLES NE CONTIENT L'UN DES POINTS
C          =4 SI POINT TROP PROCHE D'UN AUTRE POINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  JANVIER 2008
C2345X7..............................................................012
      DOUBLE PRECISION  EPSARETE
      PARAMETER        (EPSARETE=1D-15)
      PARAMETER        (LCHAIN=6)
C
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      DOUBLE PRECISION  PXYD(3,MXSOMM)
C
      INTEGER           NOSOAR(MOSOAR,MXSOAR),
     %                  NOARTR(MOARTR,MXARTR),
     %                  NOARST(MXSOMM)
C
      DOUBLE PRECISION  AIRETR, CB(3), XP,YP, X1,X2,X3, Y1,Y2,Y3
      INTEGER           NOSOTR(3), NUTR(4)
C
C     BOUCLE SUR LES POINTS A TRIANGULER
C     ==================================
      NBTRVU = 0
      NTRDER = 2
      DO 100 NPT=1, NBARPI
C
C        RECHERCHE D'UN TRIANGLE CONTENANT LE POINT NPT
C        METHODE BESTIALE MAIS SURE. A AMELIORER ENSUITE
C        -----------------------------------------------
         DO 90 NTR=NTRDER,1,-1
C
            IF( NOARTR(1,NTR) .EQ. 0 ) GOTO 90
            NBTRVU = NBTRVU + 1
C
C           LE TRIANGLE NTR EST ACTIF
C           LE NO DES 3 SOMMETS DU TRIANGLE NTR
            CALL NUSOTR( NTR, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR )
C
C           PERMUTATION CIRCULAIRE POUR AMENER EN DEBUT LE PLUS PETIT NO
C           POUR OBTENIR TOUJOURS LE MEME RESULTAT
            NSMIN = MIN( NOSOTR(1), NOSOTR(2), NOSOTR(3) )
            IF( NSMIN .EQ. NOSOTR(2) ) THEN
               NOSOTR(2) = NOSOTR(3)
               NOSOTR(3) = NOSOTR(1)
               NOSOTR(1) = NSMIN
               NSMIN            = NOARTR( 2, NTR )
               NOARTR( 2, NTR ) = NOARTR( 3, NTR )
               NOARTR( 3, NTR ) = NOARTR( 1, NTR )
               NOARTR( 1, NTR ) = NSMIN
            ELSE IF( NSMIN .EQ. NOSOTR(3) ) THEN
               NOSOTR(3) = NOSOTR(2)
               NOSOTR(2) = NOSOTR(1)
               NOSOTR(1) = NSMIN
               NSMIN            = NOARTR( 3, NTR )
               NOARTR( 3, NTR ) = NOARTR( 2, NTR )
               NOARTR( 2, NTR ) = NOARTR( 1, NTR )
               NOARTR( 1, NTR ) = NSMIN
            ENDIF
C
C           CHAINAGE DES ARETES A RENDRE DELAUNAY
            N1ADEL = 0
C
C           NOMBRE DE SOUS TRIANGLES GENERES
            NBSTR  = 0
C
C           AIRE ET COORDONNEES BARYCENTRIQUES DU POINT NPT
            XP = PXYD( 1, NPT )
            YP = PXYD( 2, NPT )
C
            N1 = NOSOTR( 1 )
            X1 = PXYD( 1, N1 )
            Y1 = PXYD( 2, N1 )
C
            N2 = NOSOTR( 2 )
            X2 = PXYD( 1, N2 )
            Y2 = PXYD( 2, N2 )
C
            N3 = NOSOTR( 3 )
            X3 = PXYD( 1, N3 )
            Y3 = PXYD( 2, N3 )
C
C           2 FOIS LA SURFACE DU TRIANGLE = DETERMINANT DE LA MATRICE
C           DE CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT P
            AIRETR=( X2 - X1 ) * ( Y3 - Y1 ) - ( X3 - X1 ) * ( Y2 - Y1 )
C
            IF( AIRETR .LE. 0D0 ) THEN
C
C              LE TRIANGLE NTR EST DEGENERE
C              ----------------------------
               print *,'trprtr: Aire=',AIRETR,' du triangle de sommets'
               print *,'sommet ',N1,' X1=',X1,' Y1=',Y1
               print *,'sommet ',N2,' X2=',X2,' Y2=',Y2
               print *,'sommet ',N3,' X3=',X3,' Y3=',Y3
               GOTO 90
C
            ELSE
C
C              LE TRIANGLE NTR EST NON DEGENERE
C              --------------------------------
C              CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU POINT NPT
               CB(1)=( ( X2-XP )*( Y3-YP )-( X3-XP )*( Y2-YP ) )/ AIRETR
               CB(2)=( ( X3-XP )*( Y1-YP )-( X1-XP )*( Y3-YP ) )/ AIRETR
               CB(3)=( ( X1-XP )*( Y2-YP )-( X2-XP )*( Y1-YP ) )/ AIRETR
C
ccc               print *,'Pt',NPT,' Triangle',NTR,':',NOSOTR,' CB=',CB
C
C              LE POINT NPT EST IL EXTERIEUR AU TRIANGLE NTR?
               IF( CB(1) .LT. 0D0 .OR.
     %             CB(2) .LT. 0D0 .OR.
     %             CB(3) .LT. 0D0 ) GOTO 90
C
C              NON. EST IL FRANCHEMENT A L'INTERIEUR?
               IF( CB(1) .GT. EPSARETE .AND.
     %             CB(2) .GT. EPSARETE .AND.
     %             CB(3) .GT. EPSARETE ) THEN
C
C                 OUI. POINT FRANCHEMENT A L'INTERIEUR DU TRIANGLE NTR
C                 LE TRIANGLE NTR EST SUBDIVISE EN 3 SOUS-TRIANGLES DE SOMMET NP
C                 --------------------------------------------------------------
C                 CHAINAGE DES 3 ARETES A RENDRE DELAUNAY
                  DO 30 I=1,3
                     N = ABS( NOARTR(I,NTR) )
                     NOSOAR(LCHAIN,N) = N1ADEL
                     N1ADEL = N
 30               CONTINUE
                  CALL TR3STR( NPT,    NTR,
     %                         MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                         MOARTR, MXARTR, N1ARTR, NOARTR,
     %                         NOARST,
     %                         NUTR,   IERR )
                  IF( IERR .NE. 0 ) RETURN
                  NBSTR = 3
C                 PASSAGE AU POINT NPT SUIVANT
                  GOTO 95
C
               ELSE
C
C                 NOMBRE DE COORDONNEES BARYCENTRIQUES < EPSARETE
                  NOARPT2 = 0
                  NBCB0   = 0
                  DO 40 I=1,3
                     IF( CB(I) .LE. EPSARETE ) THEN
                        NBCB0   = NBCB0 + 1
                        NOARPT2 = I
                     ENDIF
 40               CONTINUE
C
                  IF( NBCB0 .GE. 2 ) THEN
                     print *,'trprtr: point TROP PROCHE d''un ST'
                     print *,'point ',NPT,' XP=',XP,' YP=',YP
                     print *,'Triangle NTR=',NTR,' CB=',CB
                     IERR = 4
                     RETURN
                  ENDIF
C
C                 UNE SEULE COORDONNEE < EPSARETE
C                 => POINT SUR OU PROCHE DE L'ARETE DE NUMERO NOARPT
C                 LE TRIANGLE NTR ET OPPOSE SONT SUBDIVISES EN 4 SOUS-TRIANGLES
C                 AVEC LE SOMMET COMMUN LE POINT NPT
C                 -------------------------------------------------------------
                  NOARPT = MOD(NOARPT2,3) + 1
                  CALL TR4STR( NPT,    NTR,    NOARPT, N1ADEL,
     %                         MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                         MOARTR, MXARTR, N1ARTR, NOARTR,
     %                         NOARST,
     %                         NUTR,   IERR )
                  IF( IERR .EQ. 5 ) THEN
                     IERR = 0
                     GOTO 100
                  ELSE IF( IERR .NE. 0 ) THEN
                     RETURN
                  ENDIF
                  NBSTR = 4
C                 PASSAGE AU POINT NPT SUIVANT
                  GOTO 95
               ENDIF
            ENDIF
C
 90      CONTINUE
         print *,'trprtr: AUCUN TRIANGLE contient le point',NPT
         IERR = 3
         GOTO 100
C
C        NUMERO MAXIMAL DES TRIANGLES ACTUELS
 95      DO 97 I=1,NBSTR
            NTRDER = MAX( NTRDER, NUTR(I) )
 97      CONTINUE
C
         IF( TRATRI ) THEN
C     LES TRACES DES TRIAGLES AJOUTES SONT DEMANDES
CCC   CALL EFFACE
CCCC  LE CADRE OBJET GLOBAL EN UNITES UTILISATEUR
CCC   XX1 = MIN(PXYD(1,NOSOTR(1)),PXYD(1,NOSOTR(2)),PXYD(1,NOSOTR(3)))
CCC   XX2 = MAX(PXYD(1,NOSOTR(1)),PXYD(1,NOSOTR(2)),PXYD(1,NOSOTR(3)))
CCC   YY1 = MIN(PXYD(2,NOSOTR(1)),PXYD(2,NOSOTR(2)),PXYD(2,NOSOTR(3)))
CCC   YY2 = MAX(PXYD(2,NOSOTR(1)),PXYD(2,NOSOTR(2)),PXYD(2,NOSOTR(3)))
CCC   IF( XX1 .GE. XX2 ) XX2 = XX1 + (YY2-YY1)
CCC   IF( YY1 .GE. YY2 ) YY2 = YY1 + (XX2-XX1)*0.5
CCC   CALL ISOFENETRE( XX1-(XX2-XX1), XX2+(XX2-XX1),
CCC   %                YY1-(YY2-YY1), YY2+(YY2-YY1) )
            DO 96 I=1,NBSTR
C              TRACE DU TRIANGLE NUTR(I)
               CALL MTTRTR( PXYD, NUTR(I), MOARTR, NOARTR,
     %                      MOSOAR, NOSOAR, I, NCBLAN )
 96         CONTINUE
C
C           SAISIE D'UN POINT PAR CLIC DE LA SOURIS
C           OU ENTREE D'UN CARACTERE POUR VOIR LE TRACE
            CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
         ENDIF
C
C        MISE EN DELAUNAY DES ARETES PERIPHERIQUES DES SOUS TRIANGLES CREES
         CALL TEDELA( PXYD,   NOARST,
     %                MOSOAR, MXSOAR, N1SOAR, NOSOAR, N1ADEL,
     %                MOARTR, MXARTR, NOARTR, MODIFS )
C
cccC        affichage des sous triangles crees
ccc         DO 99 I=1,NBSTR
cccC           LE NO DES 3 SOMMETS DU TRIANGLE NUTR(I)
ccc            CALL NUSOTR( NUTR(I), MOSOAR, NOSOAR,
ccc     %                   MOARTR,  NOARTR, NOSOTR )
ccc            print *,'PT',NPT,' sous-triangle',NOSOTR
ccc 99      CONTINUE
 100  CONTINUE
C
ccc      DO 200 I=1,NTRDER
ccc         IF( noartr(1,i) .eq. 0 ) goto 200
ccc            CALL NUSOTR( i, MOSOAR, NOSOAR,
ccc     %                   MOARTR,  NOARTR, NOSOTR )
ccc         print *,'Triangle',i,' noartr',(noartr(k,i),k=1,3),
ccc     %           ' nosotr',nosotr
ccc         do 199 j=1,3
ccc            n = abs(noartr(j,i))
ccc            print *,'Arete',n,' : ',(nosoar(k,n),k=1,5)
ccc 199     continue
ccc 200  CONTINUE
ccc
      IF(  LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'TRPRTR: NOMBRE TRIANGLES VUS NBT=',NBTRVU,
     %                  ' POUR NBS=',NBSOMM,' SOMMETS'
      ELSE
       WRITE(IMPRIM,*)'TRPRTR: VISITED TRIANGLES NUMBER NBT=',NBTRVU,
     %                  ' FOR NBS=',NBSOMM,' VERTICES'
      ENDIF
      WRITE(IMPRIM,19000) DBLE(NBTRVU)/NBSOMM,
     %                    DBLE(NBTRVU)/DBLE(NBSOMM)**2
19000 FORMAT(' NBT/NBS=',G12.3,'  NBT/(NBS**2)=',G12.3)
C
      RETURN
      END
