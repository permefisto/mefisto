      SUBROUTINE DFTOP1( NBOBPR, MNNSEF,  NBELEM, NBELTG, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE NOMBRE D'EF DE CHAQUE PLSV DE L'OBJET
C -----
C ENTREES :
C ---------
C NBOBPR : NOMBRE D'OBJETS PREMIERS
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF' DE CHAQUE OBJET PREMIER
C
C SORTIES :
C ---------
C NBELEM : NBELEM(1)=NOMBRE DE NOEUDSOMMETS
C          NBELEM(2)=SEGMENTS
C          NBELEM(3)=TRIANGLES
C          NBELEM(4)=QUADRANGLES
C          NBELEM(5)=TETRAEDRES
C          NBELEM(6)=PENTAEDRES
C          NBELEM(7)=HEXAEDRES
C          NBELEM(8)=6-CUBES      DE L'OBJET COMPLET
C          NBELEM(9)=PYRAMIDES
C NBELTG : NOMBRE D'EF A TG DE CHACUN DES TYPES POUR L'OBJET COMPLET
C IERR   : 0 SI PAS D'ERREUR DETECTEE , NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1989
C MODIF  : ALAIN PERRONNET  TEXAS A & M UNIVERSITY            JULY  2005
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      INTEGER           NBELEM(9), MNNSEF(NBOBPR), NBELTG(9)
C
C     MISE A ZERO DU NOMBRE DES ELEMENTS FINIS PAR OBJET PREMIER
      CALL AZEROI( 9, NBELEM )
      CALL AZEROI( 9, NBELTG )
C
C     BOUCLE SUR LES OBJETS PREMIERS
C     ==============================
      DO 500 I=1,NBOBPR
C        L'ADRESSE MCN DU TABLEAU NSEF DE L'OBJET PREMIER I
         MNSS = MNNSEF(I)
C
C        LES PARAMETRES DU MAILLAGE DE L'OBJET PREMIER I
         CALL NSEFPA( MCN(MNSS),
     %                NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %                NX    , NY    , NZ    ,
     %                IERR   )
         IF( IERR .GT. 0 ) GOTO 500
C
C        LE NOMBRE D'EF A TG DE CET OBJET PREMIER
         NBEFTG = MCN(MNSS+WBEFTG)
C
C        LE MAILLAGE EST IL STRUCTURE ?
         IF( NUTYMA .EQ. 0 ) THEN
C
C           MAILLAGE NON STRUCTURE
C           ======================
            IF( NBSOEF .EQ. 1 .OR. NBSOEF .EQ. 2 ) THEN
C
C              LES NOEUDSOMMETS ET SEGMENTS
               NBELEM(NBSOEF) = NBELEM(NBSOEF) + NBEFOB
               NBELTG(NBSOEF) = NBELTG(NBSOEF) + NBEFTG
C
            ELSE IF( NBSOEF .EQ. 4 ) THEN
C
C              LES TRIANGLES ET QUADRANGLES
C              L'ADRESSE DU TABLEAU NUSOEF
               MN  = MNSS + WUSOEF + 3
               MNT = MNSS + LDAPEF
               DO 3 J=1,NBEFOB
C                 LE 4-EME SOMMET EST IL NUL ?
                  IF( MCN( MN ) .EQ. 0 ) THEN
C                    TRIANGLE
                     NCOGEL = 3
                  ELSE
C                    QUADRANGLE
                     NCOGEL = 4
                  ENDIF
C                 UN EF DE PLUS DE CE TYPE
                  NBELEM(NCOGEL) = NBELEM(NCOGEL) + 1
                  MN = MN + 4
C
C                 LE NUMERO DE L'EF A TG?
                  IF( NBEFTG .GT. 0 .AND. NBTGEF .GT. 0 ) THEN
                     NUEFTG = MCN( MNT )
                     MNT = MNT + 1
                  ELSE
                     NUEFTG = 0
                  ENDIF
                  IF( NUEFTG .GT. 0 ) THEN
C                    UN EF A TG DE PLUS
                     NBELTG(NCOGEL) = NBELTG(NCOGEL) + 1
                  ENDIF
 3             CONTINUE
C
            ELSE IF( NBSOEF .EQ. 8 ) THEN
C
C              LES TETRAEDRES PYRAMIDES PENTAEDRES HEXAEDRES
               MN  = MNSS + WUSOEF + 7
               MNT = MNSS + LDAPEF
               DO 5 J=1,NBEFOB
C                 LE 8-EME SOMMET EST IL NON NUL ?
                  IF( MCN( MN ) .NE. 0 ) THEN
C                    HEXAEDRE
                     NCOGEL = 7
                  ELSE IF( MCN( MN-2 ) .NE. 0 ) THEN
C                    PENTAEDRE
                     NCOGEL = 6
                  ELSE IF( MCN( MN-3 ) .NE. 0 ) THEN
C                    PYRAMIDE
                     NCOGEL = 9
                  ELSE IF( MCN( MN-4 ) .NE. 0 ) THEN
C                    TETRAEDRE
                     NCOGEL = 5
                  ELSE
C                    EF INCONNU
                     GOTO 99
                  ENDIF
C                 UN EF DE PLUS DE CE TYPE
                  NBELEM(NCOGEL) = NBELEM(NCOGEL) + 1
                  MN = MN + 8
C
C                 LE NUMERO DE L'EF A TG?
                  IF( NBEFTG .GT. 0 .AND. NBTGEF .GT. 0 ) THEN
                     NUEFTG = MCN( MNT )
                     MNT = MNT + 1
                  ELSE
                     NUEFTG = 0
                  ENDIF
                  IF( NUEFTG .GT. 0 ) THEN
C                    UN EF A TG DE PLUS
                     NBELTG(NCOGEL) = NBELTG(NCOGEL) + 1
                  ENDIF
 5             CONTINUE
C
            ELSE IF( NBSOEF .EQ. 64 ) THEN
C
C              LES 6-CUBES NON-STRUCTURES
               NCOGEL = 8
               NBELEM(NCOGEL) = NBELEM(NCOGEL) + NBEFOB
               NBELTG(NCOGEL) = NBELTG(NCOGEL) + NBEFTG
            ENDIF
C
         ELSE
C
C           MAILLAGE STRUCTURE
C           ==================
C           LE NOMBRE D'EF A TG DE CE TYPE
            NBELTG(NUTYMA) = NBELTG(NUTYMA) + NBEFTG
C
            GOTO ( 10, 20, 30, 40, 50, 60, 70, 80, 99 ), NUTYMA
C
C           LE NOMBRE DE NSEF DE L'OBJET PREMIER
C           NOEUDSOMMETS
 10         NBELEM(1) = NBELEM(1) + NX
            GOTO 500
C           SEGMENT STRUCTURE
 20         NBELEM(2) = NBELEM(2) + NX
            GOTO 500
C           TRIANGLE STRUCTURE
 30         NBELEM(3) = NBELEM(3) + NX ** 2
            GOTO 500
C           QUADRANGLE STRUCTURE
 40         NBELEM(4) = NBELEM(4) + NX * NY
            GOTO 500
C           TETRAEDRE STRUCTURE
 50         NBELEM(5) = NBELEM(5) + NX ** 3
            GOTO 500
C           PENTAEDRE STRUCTURE
 60         NBELEM(6) = NBELEM(6) +( NX ** 2 ) * NY
            GOTO 500
C           HEXAEDRE STRUCTURE
 70         NBELEM(7) = NBELEM(7) + NX * NY * NZ
            GOTO 500
C           6-CUBE STRUCTURE
 80         NBELEM(8) = NBELEM(8) + NX ** NY
         ENDIF
         GOTO 500
C
 99      IF( LANGAG .EQ. 0 ) THEN
           WRITE(IMPRIM,*) 'DFTOP1: ERREUR NUTYMA=',NUTYMA,' EF INCONNU'
         ELSE
            WRITE(IMPRIM,*) 'DFTOP1: ERROR NUTYMA=',NUTYMA,' UNKNOWN FE'
         ENDIF
         IERR = 1
C
 500  CONTINUE
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'NOMBRE de SOMMETS    ',NBELEM(1)
         WRITE(IMPRIM,*) 'NOMBRE de SEGMENTS   ',NBELEM(2)
         WRITE(IMPRIM,*) 'NOMBRE de TRIANGLES  ',NBELEM(3)
         WRITE(IMPRIM,*) 'NOMBRE de QUADRANGLES',NBELEM(4)
         WRITE(IMPRIM,*) 'NOMBRE de TETRAEDRES ',NBELEM(5)
         WRITE(IMPRIM,*) 'NOMBRE de PYRAMIDES  ',NBELEM(9)
         WRITE(IMPRIM,*) 'NOMBRE de PENTAEDRES ',NBELEM(6)
         WRITE(IMPRIM,*) 'NOMBRE de HEXAEDRES  ',NBELEM(7)
         WRITE(IMPRIM,*) 'NOMBRE de 6-CUBES    ',NBELEM(8)
      ELSE
         WRITE(IMPRIM,*) 'NUMBER of VERTICES   ',NBELEM(1)
         WRITE(IMPRIM,*) 'NUMBER of SEGMENTS   ',NBELEM(2)
         WRITE(IMPRIM,*) 'NUMBER of TRIANGLES  ',NBELEM(3)
         WRITE(IMPRIM,*) 'NUMBER of QUADRANGLES',NBELEM(4)
         WRITE(IMPRIM,*) 'NUMBER of TETRAHEDRA ',NBELEM(5)
         WRITE(IMPRIM,*) 'NUMBER of PYRAMIDS   ',NBELEM(9)
         WRITE(IMPRIM,*) 'NUMBER of PENTAHEDRA ',NBELEM(6)
         WRITE(IMPRIM,*) 'NUMBER of HEXAHEDRA  ',NBELEM(7)
         WRITE(IMPRIM,*) 'NUMBER of 6-CUBES    ',NBELEM(8)
       ENDIF
C
      RETURN
      END
