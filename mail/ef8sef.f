      SUBROUTINE EF8SEF( NBSOVI, XYZSTI, NBEFVI, NSEFVI,
     %           MXSOVF, NBSOVF, XYZSTF, NBEFVF, NSEFVF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES 8 SOUS EF DES EF D'UN VOLUME
C -----
C
C ENTREES:
C --------
C NBSOVI : NOMBRE DE SOMMETS DU VOLUME INITIAL
C XYZSTI : 3 COORDONNEES DES NBSOVI SOMMETS DU VOLUME INITIAL
C NBEFVI : NOMBRE D'EF DU VOLUME INITIAL
C NSEFVI : NO DES 8 SOMMETS DES EF DU VOLUME INITIAL
C MXSOVF : NOMBRE MAXIMAL DE SOMMETS DECLARABLES POUR LE MAILLAGE FINAL
C
C SORTIES:
C --------
C NBSOVF : NOMBRE DE SOMMETS DU VOLUME FINAL
C XYZSTF : 3 COORDONNEES DES NBSOVF SOMMETS DU VOLUME FINAL
C NBEFVF : NOMBRE D'EF DU VOLUME FINAL
C NSEFVF : NO DES 8 SOMMETS DES EF DU VOLUME FINAL
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C.......................................................................
      include"./incl/sotgar.inc"
      include"./incl/sotgfc.inc"
      REAL     XYZSTI(3,NBSOVI), XYZSTF(3,MXSOVF)
      INTEGER  NSEFVI(8,NBEFVI), NSEFVF(8,NBEFVF)
C
      INTEGER  NUSEF(27)
C
      INTEGER  NUSTET(4,8)
      INTEGER  NUSPYR(5,10)
      INTEGER  NUSPEN(6,8)
      INTEGER  NUSHEX(8,8)
      DATA     NUSTET / 1,5,7,8,  2,6,5,9,  3,7,6,10,  8,9,10,4,
     %                  5,8,9,10, 5,9,6,10, 5,6,7,10,  5,7,8,10 /
C
      DATA     NUSPYR / 1,6,14,9,10,   6,2,7,14,11,
     %                  9,14,8,4,13,   14,7,3,8,12,
     %                 10,11,12,13,5,  10,13,12,11,14,
     %                  6,14,10,11,0,   7,14,11,12,0,
     %                  8,14,12,13,0,   9,14,13,10,0  /
C
      DATA     NUSPEN / 1,7,9,10,17,16,    2,8,7,11,18,17,
     %                  7,8,9,17,18,16,    8,3,9,18,12,16,
     %                 10,17,16,4,13,15,  17,11,18,13,5,14,
     %                 17,18,16,13,14,15, 18,12,16,14,6,15  /
C
      DATA     NUSHEX / 1,9,21,12,13,23,27,22,   9,2,10,21,23,14,25,27,
     %                 12,21,11,4,22,27,26,16,  21,10,3,11,27,25,15,26,
     %                 13,23,27,22,5,17,24,20,  23,14,25,27,17,6,18,24,
     %                 22,27,26,16,20,24,19,8,  27,25,15,26,24,18,7,19 /
C
      IERR = 0
      DO 5 NCOGEL=5,9
C        LE NUMERO DES 2 SOMMETS DES ARETES DE L'EF DE TYPE NCOGEL
         CALL SOARFA( NCOGEL, NUSOAR(1,1,NCOGEL), NOSTFA(1,1,NCOGEL) )
 5    CONTINUE
C
C     COPIE DES XYZ DES SOMMETS INITIAUX EN DEBUT DES FINAUX
      DO 20 J=1,NBSOVI
         DO 10 I=1,3
            XYZSTF(I,J) = XYZSTI(I,J)
 10      CONTINUE
 20   CONTINUE
      NBSOVF = NBSOVI
C
C     SUBDIVISION DES TETRAEDRES PYRAMIDES PENTAEDRES HEXAEDRES
      NBEFVF = 0
      DO 100 N=1,NBEFVI
C
         IF( NBSOVF .GT. MXSOVF ) THEN
            IERR = 1
            RETURN
         ENDIF
         IF( NSEFVI(5,N) .EQ. 0 ) THEN
C
C           TETRAEDRE
C           CONSTRUCTION DES NUMEROS DES 4 SOMMETS INITIAUX
            DO 22 I=1,4
               NUSEF(I) = NSEFVI(I,N)
 22         CONTINUE
C
C           CONSTRUCTION DES 6 MILIEUX DES ARETES ET DE LEUR NUMERO
            DO 24 J=1,6
               NBSOVF = NBSOVF + 1
               NUSEF(4+J) = NBSOVF
               NSA1 = NUSEF( NUSOAR(1,J,5) )
               NSA2 = NUSEF( NUSOAR(2,J,5) )
               XYZSTF(1,NBSOVF) = (XYZSTF(1,NSA1) + XYZSTF(1,NSA2))* 0.5
               XYZSTF(2,NBSOVF) = (XYZSTF(2,NSA1) + XYZSTF(2,NSA2))* 0.5
               XYZSTF(3,NBSOVF) = (XYZSTF(3,NSA1) + XYZSTF(3,NSA2))* 0.5
 24         CONTINUE
C
C           DECOUPAGE DU TETRAEDRE EN 8 SOUS-TETRAEDRES
            DO 29 J=1,8
               NBEFVF = NBEFVF + 1
               DO 25 I=1,4
C                 LE NUMERO DES 4 SOMMETS DU SOUS EF J
                  NSEFVF(  I,NBEFVF) = NUSEF( NUSTET(I,J) )
                  NSEFVF(4+I,NBEFVF) = 0
 25            CONTINUE
 29         CONTINUE
C
         ELSE IF( NSEFVI(6,N) .EQ. 0 ) THEN
C
C           PYRAMIDE
C           CONSTRUCTION DES NUMEROS DES 5 SOMMETS INITIAUX
            DO 32 I=1,5
               NUSEF(I) = NSEFVI(I,N)
 32         CONTINUE
C
C           CONSTRUCTION DES 8 MILIEUX DES ARETES ET DE LEUR NUMERO
            DO 34 J=1,8
               NBSOVF = NBSOVF + 1
               NUSEF(5+J) = NBSOVF
               NSA1 = NUSEF( NUSOAR(1,J,9) )
               NSA2 = NUSEF( NUSOAR(2,J,9) )
               XYZSTF(1,NBSOVF) = (XYZSTF(1,NSA1) + XYZSTF(1,NSA2))* 0.5
               XYZSTF(2,NBSOVF) = (XYZSTF(2,NSA1) + XYZSTF(2,NSA2))* 0.5
               XYZSTF(3,NBSOVF) = (XYZSTF(3,NSA1) + XYZSTF(3,NSA2))* 0.5
 34         CONTINUE
C
C           CONSTRUCTION DU BARYCENTRE DE LA FACE QUADRANGULAIRE
C           FACE 1 DE LA PYRAMIDE
            NBSOVF = NBSOVF + 1
            NUSEF(14) = NBSOVF
            NSA1 = NUSEF( NOSTFA(1,1,9) )
            NSA2 = NUSEF( NOSTFA(2,1,9) )
            NSA3 = NUSEF( NOSTFA(3,1,9) )
            NSA4 = NUSEF( NOSTFA(4,1,9) )
            XYZSTF(1,NBSOVF) = ( XYZSTF(1,NSA1) + XYZSTF(1,NSA2)
     %                         + XYZSTF(1,NSA3) + XYZSTF(1,NSA4) )* 0.25
            XYZSTF(2,NBSOVF) = ( XYZSTF(2,NSA1) + XYZSTF(2,NSA2)
     %                         + XYZSTF(2,NSA3) + XYZSTF(2,NSA4) )* 0.25
            XYZSTF(3,NBSOVF) = ( XYZSTF(3,NSA1) + XYZSTF(3,NSA2)
     %                         + XYZSTF(3,NSA3) + XYZSTF(3,NSA4) )* 0.25
C
C           DECOUPAGE DE LA PYRAMIDE EN 6 SOUS-PYRAMIDES ET 4 TETRAEDRES
            DO 39 J=1,10
               NBEFVF = NBEFVF + 1
               DO 38 I=1,5
C                 LE NUMERO DES 5 SOMMETS DU SOUS EF J
                  NSA1 = NUSPYR(I,J)
                  IF( NSA1 .GT. 0 ) NSA1 = NUSEF( NSA1 )
                  NSEFVF(I,NBEFVF) = NSA1
 38            CONTINUE
               NSEFVF(6,NBEFVF) = 0
               NSEFVF(7,NBEFVF) = 0
               NSEFVF(8,NBEFVF) = 0
 39         CONTINUE
C
         ELSE IF( NSEFVI(7,N) .EQ. 0 ) THEN
C
C           PENTAEDRE
C           CONSTRUCTION DES NUMEROS DES 6 SOMMETS INITIAUX
            DO 42 I=1,6
               NUSEF(I) = NSEFVI(I,N)
 42         CONTINUE
C
C           CONSTRUCTION DES 9 MILIEUX DES ARETES ET DE LEUR NUMERO
            DO 44 J=1,9
               NBSOVF = NBSOVF + 1
               NUSEF(6+J) = NBSOVF
               NSA1 = NUSEF( NUSOAR(1,J,6) )
               NSA2 = NUSEF( NUSOAR(2,J,6) )
               XYZSTF(1,NBSOVF) = (XYZSTF(1,NSA1) + XYZSTF(1,NSA2))* 0.5
               XYZSTF(2,NBSOVF) = (XYZSTF(2,NSA1) + XYZSTF(2,NSA2))* 0.5
               XYZSTF(3,NBSOVF) = (XYZSTF(3,NSA1) + XYZSTF(3,NSA2))* 0.5
 44         CONTINUE
C
C           CONSTRUCTION DES 3 BARYCENTRES DES FACES QUADRANGULAIRES
C           FACE 2 DU PENTAEDRE
            NBSOVF = NBSOVF + 1
            NUSEF(16) = NBSOVF
            NSA1 = NUSEF( NOSTFA(1,2,6) )
            NSA2 = NUSEF( NOSTFA(2,2,6) )
            NSA3 = NUSEF( NOSTFA(3,2,6) )
            NSA4 = NUSEF( NOSTFA(4,2,6) )
            XYZSTF(1,NBSOVF) = ( XYZSTF(1,NSA1) + XYZSTF(1,NSA2)
     %                         + XYZSTF(1,NSA3) + XYZSTF(1,NSA4) )* 0.25
            XYZSTF(2,NBSOVF) = ( XYZSTF(2,NSA1) + XYZSTF(2,NSA2)
     %                         + XYZSTF(2,NSA3) + XYZSTF(2,NSA4) )* 0.25
            XYZSTF(3,NBSOVF) = ( XYZSTF(3,NSA1) + XYZSTF(3,NSA2)
     %                         + XYZSTF(3,NSA3) + XYZSTF(3,NSA4) )* 0.25
C           FACE 3 DU PENTAEDRE
            NBSOVF = NBSOVF + 1
            NUSEF(17) = NBSOVF
            NSA1 = NUSEF( NOSTFA(1,3,6) )
            NSA2 = NUSEF( NOSTFA(2,3,6) )
            NSA3 = NUSEF( NOSTFA(3,3,6) )
            NSA4 = NUSEF( NOSTFA(4,3,6) )
            XYZSTF(1,NBSOVF) = ( XYZSTF(1,NSA1) + XYZSTF(1,NSA2)
     %                         + XYZSTF(1,NSA3) + XYZSTF(1,NSA4) )* 0.25
            XYZSTF(2,NBSOVF) = ( XYZSTF(2,NSA1) + XYZSTF(2,NSA2)
     %                         + XYZSTF(2,NSA3) + XYZSTF(2,NSA4) )* 0.25
            XYZSTF(3,NBSOVF) = ( XYZSTF(3,NSA1) + XYZSTF(3,NSA2)
     %                         + XYZSTF(3,NSA3) + XYZSTF(3,NSA4) )* 0.25
C           FACE 5 DU PENTAEDRE
            NBSOVF = NBSOVF + 1
            NUSEF(18) = NBSOVF
            NSA1 = NUSEF( NOSTFA(1,5,6) )
            NSA2 = NUSEF( NOSTFA(2,5,6) )
            NSA3 = NUSEF( NOSTFA(3,5,6) )
            NSA4 = NUSEF( NOSTFA(4,5,6) )
            XYZSTF(1,NBSOVF) = ( XYZSTF(1,NSA1) + XYZSTF(1,NSA2)
     %                         + XYZSTF(1,NSA3) + XYZSTF(1,NSA4) )* 0.25
            XYZSTF(2,NBSOVF) = ( XYZSTF(2,NSA1) + XYZSTF(2,NSA2)
     %                         + XYZSTF(2,NSA3) + XYZSTF(2,NSA4) )* 0.25
            XYZSTF(3,NBSOVF) = ( XYZSTF(3,NSA1) + XYZSTF(3,NSA2)
     %                         + XYZSTF(3,NSA3) + XYZSTF(3,NSA4) )* 0.25
C
C           DECOUPAGE DU PENTAEDRE EN 8 SOUS-PENTAEDRES
            DO 49 J=1,8
               NBEFVF = NBEFVF + 1
               DO 48 I=1,6
C                 LE NUMERO DES 6 SOMMETS DU SOUS EF J
                  NSEFVF(I,NBEFVF) = NUSEF( NUSPEN(I,J) )
 48            CONTINUE
               NSEFVF(7,NBEFVF) = 0
               NSEFVF(8,NBEFVF) = 0
 49         CONTINUE
         ELSE
C
C           HEXAEDRE
C           CONSTRUCTION DES NUMEROS DES 8 SOMMETS INITIAUX
            DO 62 I=1,8
               NUSEF(I) = NSEFVI(I,N)
 62         CONTINUE
C
C           CONSTRUCTION DES 12 MILIEUX DES ARETES ET DE LEUR NUMERO
            DO 64 J=1,12
               NBSOVF = NBSOVF + 1
               NUSEF(8+J) = NBSOVF
               NSA1 = NUSEF( NUSOAR(1,J,7) )
               NSA2 = NUSEF( NUSOAR(2,J,7) )
               XYZSTF(1,NBSOVF) = (XYZSTF(1,NSA1) + XYZSTF(1,NSA2))* 0.5
               XYZSTF(2,NBSOVF) = (XYZSTF(2,NSA1) + XYZSTF(2,NSA2))* 0.5
               XYZSTF(3,NBSOVF) = (XYZSTF(3,NSA1) + XYZSTF(3,NSA2))* 0.5
 64         CONTINUE
C
C           CONSTRUCTION DES 3 BARYCENTRES DES 6 FACES QUADRANGULAIRES
            DO 65 J=1,6
               NBSOVF = NBSOVF + 1
               NUSEF(20+J) = NBSOVF
               NSA1 = NUSEF( NOSTFA(1,J,7) )
               NSA2 = NUSEF( NOSTFA(2,J,7) )
               NSA3 = NUSEF( NOSTFA(3,J,7) )
               NSA4 = NUSEF( NOSTFA(4,J,7) )
               XYZSTF(1,NBSOVF)=( XYZSTF(1,NSA1) + XYZSTF(1,NSA2)
     %                          + XYZSTF(1,NSA3) + XYZSTF(1,NSA4) )*0.25
               XYZSTF(2,NBSOVF)=( XYZSTF(2,NSA1) + XYZSTF(2,NSA2)
     %                          + XYZSTF(2,NSA3) + XYZSTF(2,NSA4) )*0.25
               XYZSTF(3,NBSOVF)=( XYZSTF(3,NSA1) + XYZSTF(3,NSA2)
     %                          + XYZSTF(3,NSA3) + XYZSTF(3,NSA4) )*0.25
 65         CONTINUE
C
C           BARYCENTRE DE L'HEXAEDRE
            NBSOVF = NBSOVF + 1
            NUSEF(27) = NBSOVF
            XYZSTF(1,NBSOVF) = ( XYZSTF(1,NUSEF(1)) + XYZSTF(1,NUSEF(2))
     %                         + XYZSTF(1,NUSEF(3)) + XYZSTF(1,NUSEF(4))
     %                         + XYZSTF(1,NUSEF(5)) + XYZSTF(1,NUSEF(6))
     %                         + XYZSTF(1,NUSEF(7)) + XYZSTF(1,NUSEF(8))
     %                         ) * 0.125
            XYZSTF(2,NBSOVF) = ( XYZSTF(2,NUSEF(1)) + XYZSTF(2,NUSEF(2))
     %                         + XYZSTF(2,NUSEF(3)) + XYZSTF(2,NUSEF(4))
     %                         + XYZSTF(2,NUSEF(5)) + XYZSTF(2,NUSEF(6))
     %                         + XYZSTF(2,NUSEF(7)) + XYZSTF(2,NUSEF(8))
     %                         ) * 0.125
            XYZSTF(3,NBSOVF) = ( XYZSTF(3,NUSEF(1)) + XYZSTF(3,NUSEF(2))
     %                         + XYZSTF(3,NUSEF(3)) + XYZSTF(3,NUSEF(4))
     %                         + XYZSTF(3,NUSEF(5)) + XYZSTF(3,NUSEF(6))
     %                         + XYZSTF(3,NUSEF(7)) + XYZSTF(3,NUSEF(8))
     %                         ) * 0.125
C
C           DECOUPAGE DE L'HEXAEDRE EN 8 SOUS-HEXAEDRES
            DO 69 J=1,8
               NBEFVF = NBEFVF + 1
               DO 68 I=1,8
C                 LE NUMERO DES 8 SOMMETS DU SOUS EF J
                  NSEFVF(I,NBEFVF) = NUSEF( NUSHEX(I,J) )
 68            CONTINUE
 69         CONTINUE
C
         ENDIF
C
 100  CONTINUE
C
      RETURN
      END
