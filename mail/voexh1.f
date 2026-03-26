        SUBROUTINE VOEXH1(NUFACE,NSENS,NPSOM,XYZ,IERR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LA SURFACE D'UN HEXAEDRE DEFINI PAR LES
C -----    SURFACES STRUCTUREES DE CHACUNE DE SES 6 FACES
C
C ENTREES:
C --------
C XYZ    : LES COORDONNEES DES SOMMETS DES 6 FACES
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
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS   DECEMBRE  1988
C.......................................................................
      include"./incl/gsmenu.inc"
      INTEGER           NUFACE(6),INVFA(6),NPSOM(6),NSENS(6),
     %                  SSUIV(4),SPREC(4)
      REAL              XYZ(3,4,6),XYZ1(3),XYZ2(3)
      DOUBLE PRECISION  VOLTER
      DATA              SSUIV/2,3,4,1/,SPREC/4,1,2,3/
C
C     INITIALISATION
C     ==============
      IERR = 0
      DO NF=1,6
         NUFACE(NF) = 0
         INVFA(NF) = 0
         NSENS(NF) = 0
         NPSOM(NF) = 0
      ENDDO
C
C     CONSTRUCTION DE L'HEXAEDRE
C     ==========================
C
C     1) LA FACE 1 : C'EST LA PREMIERE DE LA LISTE
C
      NUM = 1
      NF  = 1
      NSP = 1
      NSQ = SSUIV(NSP)
      NSR = SSUIV(NSQ)
      NPSOM(NF)  = NSP
      NSENS(NF)  = 1
      INVFA(NUM) = NF
      NUFACE(NF)  = 1
C     VERIFICATION DE L'ORIENTATION DE LA FACE 1
C     RECHERCHE D'UN AUTRE SOMMET EN DEHORS DE LA FACE 1
      NSS = 0
      NG  = 0
      DO NUM=2,6
         DO 92 NSNUM=1,4
            DO NSNF=1,4
               CALL XYZIDE(XYZ(1,NSNUM,NUM),XYZ(1,NSNF,NF),IDENT)
               IF (IDENT.EQ.1) GOTO 92
            ENDDO
            NG  = NUM
            NSS = NSNUM
            GOTO 94
 92      ENDDO
      ENDDO
C
C     ERREUR DE DONNEES : LES 6 SOMMETS SONT COPLANAIRES
      WRITE(KERR(MXLGER)(1:3),'(I3)') NF
      NBLGRC(NRERR) = 1
      KERR(1) =  'ERREUR VOEXH1: TOUS LES SOMMETS SONT DANS'
     %                 // ' LA FACE' // KERR(MXLGER)(1:3)
      CALL LEREUR
      IERR = 1
C     CALCUL DU VOLUME DU TETRAEDRE
 94   IF( VOLTER(XYZ(1,NSP,NF),XYZ(1,NSQ,NF),
     %           XYZ(1,NSR,NF),XYZ(1,NSS,NG)) .LE. 0 ) THEN
C        REORIENTATION DE LA FACE 1
         NSENS(NF) = -1
C        ECHANGE DES SOMMETS 2 ET 4 DE LA FACE 1
         DO J=1,3
            XYZJ2       = XYZ(J,2,NF)
            XYZ(J,2,NF) = XYZ(J,4,NF)
            XYZ(J,4,NF) = XYZJ2
         ENDDO
      ENDIF
C
C     2) RECHERCHE DE LA FACE 2
C
      DO 60 NUM=2,6
         IF(NUFACE(NUM).EQ.0) THEN
            DO 3 NS=1,4
               DO J=1,3
                  XYZ1(J)=XYZ(J,NSP,NF)
                  XYZ2(J)=XYZ(J,NS,NUM)
               ENDDO
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
C              1) SENS DIRECT
                  DO J=1,3
                     XYZ1(J)= XYZ(J,NSQ,NF)
                     XYZ2(J)= XYZ(J,SSUIV(NS),NUM)
                  ENDDO
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(2)   = NS
                     NSENS(2)   = 1
                     INVFA(2)   = NUM
                     NUFACE(NUM) = 2
                     GO TO 2
                  END IF
C              2) SENS RETROGRADE
                  DO J=1,3
                     XYZ1(J)=XYZ(J,NSQ,NF)
                     XYZ2(J)=XYZ(J,SPREC(NS),NUM)
                  ENDDO
                  CALL XYZIDE(XYZ1,XYZ2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     NPSOM(2)   = NS
                     NSENS(2)   = - 1
                     INVFA(2)   = NUM
                     NUFACE(NUM) = 2
                     GO TO 2
                  END IF
               END IF
 3          ENDDO
         END IF
 60   ENDDO
C
C     3) RECHERCHE DE LA FACE 3
C
 2    NF  = 1
      NSP = 1
      NSQ = SPREC(NSP)
      DO 80 NUM=2,6
         IF(NUFACE(NUM).EQ.0) THEN
         DO 8 NS=1,4
            DO J=1,3
               XYZ1(J)=XYZ(J,NSP,NF)
               XYZ2(J)=XYZ(J,NS,NUM)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
C           1) SENS DIRECT
               DO J=1,3
                  XYZ1(J)=XYZ(J,NSQ,NF)
                  XYZ2(J)=XYZ(J,SSUIV(NS),NUM)
               ENDDO
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
                  NPSOM(3)   = NS
                  NSENS(3)   = 1
                  INVFA(3)   = NUM
                  NUFACE(NUM) = 3
                  GO TO 7
               END IF
C           2) SENS RETROGRADE
               DO J=1,3
                  XYZ1(J)=XYZ(J,NSQ,NF)
                  XYZ2(J)=XYZ(J,SPREC(NS),NUM)
               ENDDO
               CALL XYZIDE(XYZ1,XYZ2,IDENT)
               IF (IDENT.EQ.1) THEN
                  NPSOM(3)   = NS
                  NSENS(3)   = - 1
                  INVFA(3)   = NUM
                  NUFACE(NUM) = 3
                  GO TO 7
               END IF
            END IF
 8          ENDDO
         END IF
 80   ENDDO
C
C     4) COUPLAGE DES FACES
C
 7    DO NUM=1,3
         NF = INVFA(NUM)
C        RECHERCHE DE LA FACE OPPOSEE
         DO 13 NFA=1,6
            IF( NUFACE(NFA) .EQ. 0 ) THEN
               DO NS =1,4
                  DO NSA=1,4
                     DO J=1,3
                        XYZ1(J)=XYZ(J,NS,NF)
                        XYZ2(J)=XYZ(J,NSA,NFA)
                     ENDDO
                     CALL XYZIDE(XYZ1,XYZ2,IDENT)
                     IF (IDENT.EQ.1) GOTO 13
                  ENDDO
               ENDDO
C              LES FACES NF ET NFA N'ONT PAS DE SOMMET COMMUN
               NUFACE(NFA) = NUM + 3
               INVFA(NUM+3) = NFA
            END IF
 13      ENDDO
      ENDDO
C     VERIFICATION
      DO NF=1,6
         IF( NUFACE(NF) .EQ. 0 ) THEN
            WRITE(KERR(MXLGER)(1:3),'(I3)') NF
            NBLGRC(NRERR) = 1
            KERR(1) =  'ERREUR VOEXH1 : LA FACE' // KERR(MXLGER)(1:3)
     %                 // ' N''EST PAS ASSOCIEE'
            CALL LEREUR
            IERR = 1
         END IF
      ENDDO
      IF ( IERR .NE. 0 ) RETURN
C
C     5) ORIENTATION DE LA FACE 4
C
      NUM = 3
      NF  = INVFA(NUM)
      NSO = NPSOM(NUM)
      IF( NSENS(3) .EQ. 1) THEN
         NSP = SPREC(NSO)
         NSQ = SPREC(NSP)
      ELSE
         NSP = SSUIV(NSO)
         NSQ = SSUIV(NSP)
      END IF
C
      NUM = 4
      NFA = INVFA(NUM)
      DO 21 NS=1,4
         DO J=1,3
            XYZ1(J)=XYZ(J,NSP,NF)
            XYZ2(J)=XYZ(J,NS,NFA)
         ENDDO
         CALL XYZIDE(XYZ1,XYZ2,IDENT)
         IF (IDENT.EQ.1) THEN
            NPSOM(NUM) = NS
C           VERIFICATION DU SENS DE PARCOURS
C           1) SENS DIRECT
            DO J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,SPREC(NS),NFA)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = 1
               GO TO 20
            END IF
C           2) SENS RETROGRADE
            DO J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,SSUIV(NS),NFA)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = - 1
               GO TO 20
            END IF
         END IF
 21   ENDDO
C     LE SOMMET N'EST PAS RETROUVE
      NBLGRC(NRERR) = 2
      WRITE(KERR(MXLGER)(1:3),'(I3)') NF
      WRITE(KERR(MXLGER)(4:6),'(I3)') NFA
      KERR(1) = 'ERREUR VOEXH1 : LES FACES' // KERR(MXLGER)(1:3)
     %          // ' ET '  // KERR(MXLGER)(4:6)
      KERR(2) = ' N''ONT PAS D''ARETE COMMUNE'
      CALL LEREUR
      IERR = 1

 20   IF ( IERR .NE. 0 ) RETURN
C
C     6) ORIENTATION DE LA FACE 5
C
      NUM = 1
      NF  = INVFA(NUM)
      NSO = NPSOM(NUM)
      NSP = SPREC(NSO)
      NSQ = SPREC(NSP)
C
      NUM = 5
      NFA = INVFA(NUM)
      DO 31 NS=1,4
         DO J=1,3
            XYZ1(J)=XYZ(J,NSP,NF)
            XYZ2(J)=XYZ(J,NS,NFA)
         ENDDO
         CALL XYZIDE(XYZ1,XYZ2,IDENT)
         IF (IDENT.EQ.1) THEN
            NPSOM(NUM) = NS
C           VERIFICATION DU SENS DE PARCOURS
C           1) SENS DIRECT
            DO J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,SSUIV(NS),NFA)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = 1
               GO TO 30
            END IF
C           2) SENS RETROGRADE
            DO J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,SPREC(NS),NFA)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = - 1
               GO TO 30
            END IF
         END IF
 31   ENDDO
C     LE SOMMET N'EST PAS RETROUVE
      NBLGRC(NRERR) = 2
      WRITE(KERR(MXLGER)(1:3),'(I3)') NF
      WRITE(KERR(MXLGER)(4:6),'(I3)') NFA
      KERR(1) = 'ERREUR VOEXH1 : LES FACES' // KERR(MXLGER)(1:3)
     %          // ' ET '  // KERR(MXLGER)(4:6)
      KERR(2) = ' N''ONT PAS D''ARETE COMMUNE'
      CALL LEREUR
      IERR = 1

 30   IF ( IERR .NE. 0 ) RETURN
C
C     7) ORIENTATION DE LA FACE 6
C
      NUM = 1
      NF  = INVFA(NUM)
      NSO = NPSOM(NUM)
      NSP = SSUIV(NSO)
      NSQ = SSUIV(NSP)
C
      NUM = 6
      NFA = INVFA(NUM)
      DO 41 NS=1,4
         DO J=1,3
            XYZ1(J)=XYZ(J,NSP,NF)
            XYZ2(J)=XYZ(J,NS,NFA)
         ENDDO
         CALL XYZIDE(XYZ1,XYZ2,IDENT)
         IF (IDENT.EQ.1) THEN
            NPSOM(NUM) = NS
C           VERIFICATION DU SENS DE PARCOURS
C           1) SENS DIRECT
            DO J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,SSUIV(NS),NFA)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = 1
               GO TO 40
            END IF
C           4) SENS RETROGRADE
            DO J=1,3
               XYZ1(J)=XYZ(J,NSQ,NF)
               XYZ2(J)=XYZ(J,SPREC(NS),NFA)
            ENDDO
            CALL XYZIDE(XYZ1,XYZ2,IDENT)
            IF (IDENT.EQ.1) THEN
               NSENS(NUM) = - 1
               GO TO 40
            END IF
         END IF
 41   ENDDO
C     LE SOMMET N'EST PAS RETROUVE
      NBLGRC(NRERR) = 2
      WRITE(KERR(MXLGER)(1:3),'(I3)') NF
      WRITE(KERR(MXLGER)(4:6),'(I3)') NFA
      KERR(1) = 'ERREUR VOEXH1 : LES FACES' // KERR(MXLGER)(1:3)
     %          // ' ET '  // KERR(MXLGER)(4:6)
      KERR(2) = ' N''ONT PAS D''ARETE COMMUNE'
      CALL LEREUR
      IERR = 1

 40   IF ( IERR .NE. 0 ) RETURN
C
C     MISE A JOUR DE NUFACE
      DO NF=1,6
         NUFACE(NF) = INVFA(NF)
      ENDDO
C
      RETURN
      END
