        SUBROUTINE VOEXP1(NUFACE,NSENS,NPSOM,NUTYPE,XYZ,IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LA SURFACE D'UN PENTAEDRE DEFINI PAR LES
C -----    SURFACES STRUCTUREES DE CHACUNE DE SES 5 FACES
C
C ENTREES:
C --------
C XYZ    : LES COORDONNEES DES SOMMETS DES 5 FACES
C
C SORTIES:
C --------
C NUFACE : RECLASSEMENT DES FACES DANS L'ORDRE USUEL
C NSENS  : ORIENTATION DES FACES
C NPSOM  : NUMERO LOCAL DU PREMIER SOMMET DE LA FACE
C NUTYPE : TYPE DE LA FACE : TRIANGLE OU QUADRILATERE
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS   MARS 1989
C.......................................................................
      include"./incl/gsmenu.inc"
      INTEGER           NUFACE(5),INVFA(5),NPSOM(5),NSENS(5),NUTYPE(5),
     %                  TSUIV(3),TPREC(3),QSUIV(4),QPREC(4)
      REAL              XYZ(3,4,5),XYZ1(3),XYZ2(3)
      DOUBLE PRECISION  VOLTER
      DATA   TSUIV/2,3,1/,TPREC/3,1,2/,QSUIV/2,3,4,1/,QPREC/4,1,2,3/
C
C     INITIALISATION
C     ==============
C
      IERR = 0
      NT   = 0
      DO 1 NF=1,5
         IF(NUTYPE(NF).EQ.3) THEN
            IF (NT.EQ.0) THEN
C              LE PREMIER  TRIANGLE EST PLACE EN FACE 1
               NT = NT + 1
               NUFACE(NF) = 1
               INVFA(1)   = NF
            ELSE
C              LE DEUXIEME TRIANGLE EST PLACE EN FACE 5
               NUFACE(NF) = 5
               INVFA(5)   = NF
            END IF
         ELSE
            NUFACE(NF) = 0
         END IF
 1    CONTINUE
      NSENS(1)  = 1
C
C     CONSTRUCTION DU PENTAEDRE
C     =========================
 59   DO 60 NF=1,5
         IF( NUTYPE(NF) .NE. 3 ) NUFACE(NF) = 0
 60   CONTINUE
C
C     1) LA FACE 1
C
      NUM = 1
      NF  = INVFA(NUM)
      NSP = 1
      NSQ = TSUIV(NSP)
      NPSOM(1)  = NSP
C
C     2) RECHERCHE DE LA FACE 2
C
      DO 3 NUM=1,5
         IF(NUFACE(NUM).EQ.0) THEN
            DO 30 NS=1,4
               DO 4 J=1,3
                  XYZ1(J)=XYZ(J,NSP,NF)
                  XYZ2(J)=XYZ(J,NS,NUM)
 4             CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
C              1) SENS DIRECT
                  DO 5 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,QPREC(NS),NUM)
 5                CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(2)   = NS
                     NSENS(2)   = 1
                     INVFA(2)   = NUM
                     NUFACE(NUM) = 2
                     GO TO 2
                  END IF
C              2) SENS RETROGRADE
                  DO 6 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,QSUIV(NS),NUM)
 6                CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(2)   = NS
                     NSENS(2)   = - 1
                     INVFA(2)   = NUM
                     NUFACE(NUM) = 2
                     GO TO 2
                  END IF
               END IF
 30         CONTINUE
         END IF
 3    CONTINUE
 2    CONTINUE
C
C     3) RECHERCHE DE LA FACE 3
C
      NUM = 1
      NF  = INVFA(NUM)
      NSP = 2
      NSQ = TSUIV(NSP)
      DO 8 NUM=1,5
         IF(NUFACE(NUM).EQ.0) THEN
            DO 80 NS=1,4
               DO 9 J=1,3
                  XYZ1(J)=XYZ(J,NSP,NF)
                  XYZ2(J)=XYZ(J,NS,NUM)
 9             CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
C              1) SENS DIRECT
                  DO 10 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,QPREC(NS),NUM)
 10               CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(3)   = NS
                     NSENS(3)   = 1
                     INVFA(3)   = NUM
                     NUFACE(NUM) = 3
                     GO TO 7
                  END IF
C              2) SENS RETROGRADE
                  DO 11 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,QSUIV(NS),NUM)
 11               CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(3)   = NS
                     NSENS(3)   = - 1
                     INVFA(3)   = NUM
                     NUFACE(NUM) = 3
                     GO TO 7
                  END IF
               END IF
 80         CONTINUE
         END IF
 8    CONTINUE
 7    CONTINUE
C
C     4) RECHERCHE DE LA FACE 4
C
      NUM = 1
      NF  = INVFA(NUM)
      NSP = 3
      NSQ = TSUIV(NSP)
      DO 12 NUM=1,5
         IF(NUFACE(NUM).EQ.0) THEN
            DO 120 NS=1,4
               DO 13 J=1,3
                  XYZ1(J)=XYZ(J,NSP,NF)
                  XYZ2(J)=XYZ(J,NS,NUM)
 13            CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
C              1) SENS DIRECT
                  DO 14 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,QPREC(NS),NUM)
 14               CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(4)   = NS
                     NSENS(4)   = 1
                     INVFA(4)   = NUM
                     NUFACE(NUM) = 4
                     GO TO 16
                  END IF
C              2) SENS RETROGRADE
                  DO 15 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,QSUIV(NS),NUM)
 15               CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(4)   = NS
                     NSENS(4)   = - 1
                     INVFA(4)   = NUM
                     NUFACE(NUM) = 4
                     GO TO 16
                  END IF
               END IF
 120        CONTINUE
         END IF
 12   CONTINUE
 16   CONTINUE
C
C     5) ORIENTATION DE LA FACE 5
C
      NUM = 2
      NF  = INVFA(NUM)
      NSO = NPSOM(NUM)
      IF( NSENS(2) .EQ. 1) THEN
         NSP = QSUIV(NSO)
         NSQ = QSUIV(NSP)
      ELSE
         NSP = QPREC(NSO)
         NSQ = QPREC(NSP)
      END IF
C
      NUM = 5
      NFA = INVFA(NUM)
      DO 21 NS=1,3
         DO 22 J=1,3
            XYZ1(J)=XYZ(J,NSP,NF)
            XYZ2(J)=XYZ(J,NS,NFA)
 22      CONTINUE
         CALL XYZIDE(XYZ1,XYZ2,IDENT)
         IF (IDENT.EQ.1) THEN
            NPSOM(NUM) = NS
C           VERIFICATION DU SENS DE PARCOURS
C           1) SENS DIRECT
            DO 23 J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,TSUIV(NS),NFA)
 23         CONTINUE
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = 1
               GO TO 20
            END IF
C           2) SENS RETROGRADE
            DO 24 J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,TPREC(NS),NFA)
 24         CONTINUE
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = - 1
               GO TO 20
            END IF
         END IF
 21   CONTINUE
C     LE SOMMET N'EST PAS RETROUVE
      NBLGRC(NRERR) = 2
      WRITE(KERR(MXLGER)(1:3),'(I3)') NF
      WRITE(KERR(MXLGER)(4:6),'(I3)') NFA
      KERR(1) = 'ERREUR VOEXP1 : LES FACES' // KERR(MXLGER)(1:3)
     %          // ' ET '  // KERR(MXLGER)(4:6)
      KERR(2) = ' N''ONT PAS D''ARETE COMMUNE'
      CALL LEREUR
      IERR = 1
 20   CONTINUE
      IF ( IERR .NE. 0 ) RETURN
C
C     DERNIERE VERIFICATION : ORIENTATION
      NT1 = INVFA(1)
      NT2 = INVFA(5)
      IF( VOLTER( XYZ(1,1,NT1), XYZ(1,2,NT1),
     %            XYZ(1,3,NT1), XYZ(1,1,NT2) ) .LE. 0 ) THEN
C         DANS LA FACE 1 LES SOMMETS 2 ET 3 SONT PERMUTES
          NSENS(1)=-1
C         PERMUTATION DES SOMMETS 2 ET 3 DE LA FACE 1
          DO 55 J=1,3
             R            = XYZ(J,2,NT1)
             XYZ(J,2,NT1) = XYZ(J,3,NT1)
             XYZ(J,3,NT1) = R
 55       CONTINUE
          GOTO 59
      ENDIF
C
C     MISE A JOUR DE NUFACE
      DO 50 NF=1,5
         NUFACE(NF) = INVFA(NF)
50    CONTINUE
C
      RETURN
      END
