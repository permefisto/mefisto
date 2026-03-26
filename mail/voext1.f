        SUBROUTINE VOEXT1(NUFACE,NSENS,NPSOM,XYZ,IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LA SURFACE D'UN TETRAEDRE DEFINI PAR LES
C -----    SURFACES STRUCTUREES DE CHACUNE DE SES 4 FACES
C
C ENTREES:
C --------
C XYZ    : LES COORDONNEES DES 3 SOMMETS DES 4 FACES
C
C SORTIES:
C --------
C NUFACE : RECLASSEMENT DES FACES DANS L'ORDRE USUEL
C NSENS  : ORIENTATION DES FACES
C NPSOM  : NUMERO LOCAL DU PREMIER SOMMET DE LA FACE
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS   MARS 1989
C.......................................................................
      include"./incl/gsmenu.inc"
      INTEGER           NUFACE(4),INVFA(4),NPSOM(4),NSENS(4),
     %                  SSUIV(3),SPREC(3)
      REAL              XYZ(3,3,4),XYZ1(3),XYZ2(3)
      DOUBLE PRECISION  VOLTER
      DATA              SSUIV/2,3,1/,SPREC/3,1,2/
C
C     INITIALISATION
C     ==============
      IERR = 0
      DO 1 NF=1,4
         NSENS(NF) = 0
 1    CONTINUE
      NSENS(1) = 1
C
 59   DO 60 NF=1,4
         NUFACE(NF) = 0
         INVFA(NF)  = 0
         NPSOM(NF)  = 0
 60   CONTINUE
C
C     CONSTRUCTION DU TETRAEDRE
C     =========================
C
C     1) LA FACE 1 : C'EST LA PREMIERE DE LA LISTE
C
      NUM = 1
      NF  = 1
      NSP = 1
      NSQ = SPREC(NSP)
      NPSOM(NF)  = NSP
      INVFA(NUM) = NF
      NUFACE(NF) = 1
C
C     2) RECHERCHE DE LA FACE 2
C
      DO 30 NUM=2,4
         IF(NUFACE(NUM).EQ.0) THEN
         DO 3 NS=1,3
            DO 4 J=1,3
               XYZ1(J)=XYZ(J,NSP,NF)
               XYZ2(J)=XYZ(J,NS,NUM)
 4          CONTINUE
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
C           1) SENS DIRECT
               DO 5 J=1,3
                  XYZ1(J)=XYZ(J,NSQ,NF)
                  XYZ2(J)=XYZ(J,SSUIV(NS),NUM)
 5             CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
                  NPSOM(2)   = NS
                  NSENS(2)   = 1
                  INVFA(2)   = NUM
                  NUFACE(NUM) = 2
                  GO TO 2
               END IF
C           2) SENS RETROGRADE
               DO 6 J=1,3
                  XYZ1(J)=XYZ(J,NSQ,NF)
                  XYZ2(J)=XYZ(J,SPREC(NS),NUM)
 6             CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
                  NPSOM(2)   = NS
                  NSENS(2)   = - 1
                  INVFA(2)   = NUM
                  NUFACE(NUM) = 2
                  GO TO 2
               END IF
            END IF
 3          CONTINUE
         END IF
 30      CONTINUE
 2    CONTINUE
C
C     3) RECHERCHE DE LA FACE 3
C
      NF  = 1
      NSP = 1
      NSQ = SSUIV(NSP)
      DO 80 NUM=2,4
         IF(NUFACE(NUM).EQ.0) THEN
         DO 8 NS=1,3
            DO 9 J=1,3
               XYZ1(J)=XYZ(J,NSP,NF)
               XYZ2(J)=XYZ(J,NS,NUM)
 9          CONTINUE
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
C           1) SENS DIRECT
               DO 10 J=1,3
                  XYZ1(J)=XYZ(J,NSQ,NF)
                  XYZ2(J)=XYZ(J,SPREC(NS),NUM)
 10            CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
                  NPSOM(3)   = NS
                  NSENS(3)   = 1
                  INVFA(3)   = NUM
                  NUFACE(NUM) = 3
                  GO TO 7
               END IF
C           2) SENS RETROGRADE
               DO 11 J=1,3
                  XYZ1(J)=XYZ(J,NSQ,NF)
                  XYZ2(J)=XYZ(J,SSUIV(NS),NUM)
 11            CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
                  NPSOM(3)   = NS
                  NSENS(3)   = - 1
                  INVFA(3)   = NUM
                  NUFACE(NUM) = 3
                  GO TO 7
               END IF
            END IF
 8          CONTINUE
            END IF
 80      CONTINUE
 7    CONTINUE
C
C     4) RECHERCHE DE LA FACE 4
C
      NUM = 1
      NF  = INVFA(NUM)
      NSO = NPSOM(NUM)
      NSP = SSUIV(NSO)
      NSQ = SPREC(NSO)
C
      DO 20 NUM=2,4
         IF(NUFACE(NUM).EQ.0) THEN
            DO 21 NS=1,3
               DO 22 J=1,3
                  XYZ1(J)=XYZ(J,NSP,NF)
                  XYZ2(J)=XYZ(J,NS,NUM)
 22            CONTINUE
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
C                 VERIFICATION DU SENS DE PARCOURS
C                 1) SENS DIRECT
                  DO 23 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,SSUIV(NS),NUM)
 23               CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(4) = NS
                     NSENS(4) = 1
                     INVFA(4) = NUM
                     NUFACE(NUM) = 4
                     GO TO 20
                  END IF
C                 2) SENS RETROGRADE
                  DO 24 J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,SPREC(NS),NUM)
 24               CONTINUE
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(4) = NS
                     NSENS(4) = - 1
                     INVFA(4) = NUM
                     NUFACE(NUM) = 4
                     GO TO 20
                  END IF
               END IF
 21         CONTINUE
C           LE SOMMET N'EST PAS RETROUVE
            NBLGRC(NRERR) = 1
            KERR(1) = 'VOEXT1:LES FACES NE FORMENT PAS UN TETRAEDRE'
            CALL LEREUR
            IERR = 1
         END IF
 20   CONTINUE
      IF ( IERR .NE. 0 ) RETURN
C
C     5) DERNIERE VERIFICATION :
C        LES FACES 2, 3 ET 4 ONT UN SOMMET COMMUN
C
      NUM = 2
      NF2 = INVFA(NUM)
      NSO = NPSOM(NUM)
      IF (NSENS(2).EQ.1) THEN
         NS2 = SPREC(NSO)
      ELSE
         NS2 = SSUIV(NSO)
      END IF
      NUM = 3
      NF3 = INVFA(NUM)
      NSO = NPSOM(NUM)
      IF (NSENS(3).EQ.1) THEN
         NS3 = SSUIV(NSO)
      ELSE
         NS3 = SPREC(NSO)
      END IF
      NUM = 4
      NF4 = INVFA(NUM)
      NSO = NPSOM(NUM)
      IF (NSENS(4).EQ.1) THEN
         NS4 = SPREC(NSO)
      ELSE
         NS4 = SSUIV(NSO)
      END IF
      DO 31 J=1,3
         XYZ1(J)=XYZ(J,NS2,NF2)
         XYZ2(J)=XYZ(J,NS3,NF3)
 31   CONTINUE
      CALL XYZIDE(XYZ1,XYZ2,IDENT1)
      DO 32 J=1,3
         XYZ1(J)=XYZ(J,NS2,NF2)
         XYZ2(J)=XYZ(J,NS4,NF4)
 32   CONTINUE
      CALL XYZIDE(XYZ1,XYZ2,IDENT2)
      IF (IDENT1*IDENT2.NE.1) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VOEXT1:LES FACES NE FORMENT PAS UN TETRAEDRE'
         CALL LEREUR
         IERR = 1
      END IF
      IF ( IERR .NE. 0 ) RETURN
C
C     MISE A JOUR DE NUFACE
      DO 50 NF=1,4
         NUFACE(NF) = INVFA(NF)
50    CONTINUE
C
C     DERNIERE VERIFICATION : ORIENTATION
C
      IF( VOLTER(XYZ(1,1,1), XYZ(1,2,1), XYZ(1,3,1), XYZ1) .LE. 0 ) THEN
C         DANS LA FACE 1 LES SOMMETS 2 ET 3 SONT PERMUTES
          NSENS(1)=-1
C         PERMUTATION DES SOMMETS 2 ET 3 DE LA FACE 1
          DO 55 J=1,3
             R          = XYZ(J,2,1)
             XYZ(J,2,1) = XYZ(J,3,1)
             XYZ(J,3,1) = R
 55       CONTINUE
          GOTO 59
      ENDIF
C
      RETURN
      END
