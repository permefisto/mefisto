      SUBROUTINE QUALEF( NCOGEL,   NOSOEL, NBSOM, XYZSOM,
     %                   SURFVOLU, QUALIT, IERR   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA QUALITE D'UN ELEMENT FINI DE CODE GEOMETRIQUE NCOGEL
C ----- ET DE NUMEROS DE SOMMETS RANGES DANS NOSOEL
C
C ENTREES:
C --------
C NCOGEL : NUMERO DU CODE GEOMETRIQUE DU SOUS-OBJET NUELEM
C          1:POINT, 2:SEGMENT, 3:TRIANGLE, 4:QUADRANGLE,
C          5:TETRAEDRE, 6:PENTAEDRE, 7:HEXAEDRE, 8:6-CUBE, 9:PYRAMIDE
C NOSOEL : NUMERO DES NBSOEF SOMMETS DU MAILLAGE
C NBSOM  : NOMBRE TOTAL DE SOMMETS DU MAILLAGE
C XYZSOM : COORDONNEES DES NBSOM SOMMETS
C
C SORTIES:
C --------
C SURFVOLU : SURFACE ou VOLUME DE L'EF SELON SON TYPE
C QUALIT : VALEUR DE LA QUALITE DE L'ELEMENT FINI VALEUR DANS [0,1]
C          =1 QUALITE MAXIMALE, =0 ELEMENT FINI DEGENERE
C IERR   : = 0 SI PAS DE PROBLEME RENCONTRE
C          > 0 SI UNE ERREUR EST  RENCONTREE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS        MAI 1998
C MODIFS : ALAIN PERRONNET Laboratoire J.L. LIONS UPMC Paris   Mars 2007
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      REAL              XYZSOM(1:3,1:NBSOM)
      REAL              SURFTR(4)
      DOUBLE PRECISION  VOLTER
      INTEGER           NOSOEL(1:*)
      INTRINSIC         REAL

      IERR   = 0
      SURFVOLU = 0.

C     LES NUMEROS DES SOMMETS SONT ILS CORRECTS ?
C     ===========================================
      NBS  = NBSOME(NCOGEL)
      DO 10 I=1,NBS
         IF( NOSOEL(I) .LE. 0 .OR. NOSOEL(I) .GT. NBSOM ) THEN
C           ERREUR DETECTEE : NO DE SOMMET INCORRECT
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'qualef: Erreur: No SOMMET ',NOSOEL(I),
     %          '<=0 ou SUPERIEUR a NBSOM=',NBSOM,' dans EF de SOMMETS',
     %          (NOSOEL(IERR),IERR=1,NBS)
            ELSE
               PRINT*,'qualef: ERROR: VERTEX Number',NOSOEL(I),
     %              '<=0 or SUPERIOR to NBSOM=',NBSOM,' in FE VERTICES',
     %              (NOSOEL(IERR),IERR=1,NBS)
            ENDIF
C           NUMERO DE SOMMET IMPOSE A 0 POUR MONTRER L'ERREUR
            NOSOEL(I) = 0
C           QUALITE MARQUEE ERREUR!
            QUALIT = -5.
            IERR   = 1
            RETURN
         ENDIF
 10   CONTINUE

      GOTO ( 100, 200, 300, 400, 500, 600, 700, 800, 900 ), NCOGEL

C     ERREUR DETECTEE : NCOGEL INCORRECT
      NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:2),'(I2)') NCOGEL
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'qualef: CODE GEOMETRIQUE EF INCORRECT NCOGEL='
     %          // KERR(MXLGER)(1:2)
      ELSE
         KERR(1) = 'qualef: GEOMETRIC FE CODE is INCORRECT NCOGEL='
     %          // KERR(MXLGER)(1:2)
      ENDIF
      CALL LEREUR
      IERR = 2
      RETURN

C     POINT
C     =====
  100 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'QUALITE D''UN POINT n''a PAS de SENS'
      ELSE
         KERR(1) = 'QUALITY of a POINT is WITHOUT SENSE'
      ENDIF
      CALL LEREUR
      IERR = 1
      RETURN

C     LIGNE: QUALITE=0 MAIS LONGUEUR DE L'ARETE
C     ======
  200 NS1 = NOSOEL(1)
      NS2 = NOSOEL(2)
      SURFVOLU = SQRT( (XYZSOM(1,NS2)  - XYZSOM(1,NS1) ) **2
     %               + (XYZSOM(2,NS2)  - XYZSOM(2,NS1) ) **2
     %               + (XYZSOM(3,NS2)  - XYZSOM(3,NS1) ) **2 )
      QUALIT = 0.
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'QUALITE d''une ARETE=0  mais LONGUEUR=             '
         WRITE(KERR(1)(40:53),'(G14.6)') SURFVOLU
      ELSE
         KERR(1) = 'QUALITY of an EDGE=0 but LENGTH=              '
         WRITE(KERR(1)(33:46),'(G14.6)') SURFVOLU
      ENDIF
      CALL LEREUR
      IERR = 1
      RETURN


C     TRIANGLE
C     ========
 300  SURFVOLU = SURTRR( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(2)),
     %                   XYZSOM(1,NOSOEL(3)) )
      IF( SURFVOLU .LE. 0 ) THEN
         CALL EFSURFNE( XYZSOM, 3, NOSOEL, SURFVOLU )
         IF( SURFVOLU .LT. 0 ) THEN
            QUALIT = -1.0
         ELSE
            QUALIT = 0.0
         ENDIF
      ELSE
         CALL QUATRI( NOSOEL, XYZSOM,  QUALIT )
      ENDIF
      GOTO 9999

C     QUADRANGLE
C     ==========
  400 SURFVOLU = SURTRR( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(2)),
     %                   XYZSOM(1,NOSOEL(3)) )
     %         + SURTRR( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(3)),
     %                   XYZSOM(1,NOSOEL(4)) )
      IF( SURFVOLU .LE. 0 ) THEN
         CALL EFSURFNE( XYZSOM, 4, NOSOEL, SURFVOLU )
         IF( SURFVOLU .LT. 0 ) THEN
            QUALIT = -1.0
         ELSE
            QUALIT = 0.0
         ENDIF
      ELSE
         CALL QUAQUA( NOSOEL, XYZSOM,  QUALIT )
      ENDIF
      GOTO 9999

C     TETRAEDRE
C     =========
  500 CALL QUATET( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(2)),
     %             XYZSOM(1,NOSOEL(3)), XYZSOM(1,NOSOEL(4)),
     %             A, B, SURFTR, SURFVOLU, QUALIT )
      IF( SURFVOLU .LE. 0 ) THEN
         CALL EFVONE( XYZSOM, 4, NOSOEL, SURFVOLU )
         IF( SURFVOLU .LT. 0 ) THEN
            QUALIT = -1.0
         ELSE
            QUALIT = 0.0
         ENDIF
      ENDIF
      GOTO 9999

C     PENTAEDRE
C     =========
  600 SURFVOLU = REAL( VOLTER( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(2)),
     %                       XYZSOM(1,NOSOEL(3)), XYZSOM(1,NOSOEL(4)) ))
      IF( SURFVOLU .LE. 0 ) THEN
C        AFFICHER L'EF DE TETRAEDRE 1234 AVEC UN SURFVOLU NEGATIF
         CALL EFVONE( XYZSOM, 6, NOSOEL, SURFVOLU )
         IF( SURFVOLU .LT. 0 ) THEN
            QUALIT = -1.0
         ELSE
            QUALIT = 0.0
         ENDIF
      ELSE
C        CALCUL DE LA QUALITE D'UN PENTAEDRE COMME LA
C        QUALITE MINIMALE DE SES 5 FACES
         CALL QUAPEN( NOSOEL, XYZSOM,  QUALIT )
      ENDIF
      GOTO 9999

C     HEXAEDRE
C     ========
 700  SURFVOLU = REAL( VOLTER( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(2)),
     %                       XYZSOM(1,NOSOEL(3)), XYZSOM(1,NOSOEL(5)) ))
      IF( SURFVOLU .LE. 0 ) THEN
         CALL EFVONE( XYZSOM, 8, NOSOEL, SURFVOLU )
         IF( SURFVOLU .LT. 0 ) THEN
            QUALIT = -1.0
         ELSE
            QUALIT = 0.0
         ENDIF
      ELSE
         CALL QUAHEX( NOSOEL, XYZSOM,  QUALIT )
      ENDIF
      GOTO 9999

C     6-CUBE    QUALITE A DEFINIR
C     ======
  800 QUALIT=1.0
      GOTO 9999

C     PYRAMIDE
C     ========
 900  SURFVOLU = REAL( VOLTER( XYZSOM(1,NOSOEL(1)), XYZSOM(1,NOSOEL(2)),
     %                       XYZSOM(1,NOSOEL(3)), XYZSOM(1,NOSOEL(5)) ))
      IF( SURFVOLU .LE. 0 ) THEN
         CALL EFVONE( XYZSOM, 5, NOSOEL, SURFVOLU )
         IF( SURFVOLU .LT. 0 ) THEN
            QUALIT = -1.0
         ELSE
            QUALIT = 0.0
         ENDIF
      ELSE
         CALL QUAPYR( NOSOEL, XYZSOM,  QUALIT )
      ENDIF

 9999 IF( QUALIT .LE. 0 ) THEN
ccc         write(imprim,*)
         IF( LANGAG .EQ. 0 ) THEN
            write(imprim,*)'qualef: EF Sommets=',(nosoel(i),i=1,nbs),
     %                     ' a une QUALITE=',qualit,'<=0!'
            do 950 i=1,nbs
               ns = nosoel(i)
               write(imprim,10950) nosoel(i),(XYZSOM(k,ns),k=1,3)
 950        continue
            if( ncogel .GE. 5 ) then
               write(imprim,*)'qualef: EF VOLUME=',SURFVOLU
            else if( ncogel .GE. 3 ) then
               write(imprim,*)'qualef: EF SURFACE=',SURFVOLU
            endif
         ELSE
            write(imprim,*)'qualef: FE VERTICES=',(nosoel(i),i=1,nbs),
     %                     ' with a QUALITY=',qualit,'<=0!'
            do 951 i=1,nbs
               ns = nosoel(i)
               write(imprim,20950) nosoel(i),(XYZSOM(k,ns),k=1,3)
 951        continue
            if( ncogel .GE. 5 ) then
               write(imprim,*)'qualef: FE VOLUME=',SURFVOLU
            else if( ncogel .GE. 3 ) then
               write(imprim,*)'qualef: FE SURFACE=',SURFVOLU
            endif
         ENDIF
         write(imprim,*)
10950    format(' Sommet',i10,' XYZ=',3G18.9)
20950    format(' Vertex',i10,' XYZ=',3G19.9)
      ENDIF

      QUALIT = MAX( -1.0, QUALIT )
      RETURN
      END
