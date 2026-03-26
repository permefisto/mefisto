        SUBROUTINE JOTQ2L( NBSOM1, XYZ1,   NBTGL1, XYZTG1,
     %                     NUTYM1, NUTFM1, NBARE1, NOSOA1,
     %                     NUEFA1, NBEFT1, NUTGL1,
     %                     NBSOM2, XYZ2,   NBTGL2, XYZTG2,
     %                     NUTYM2, NUTFM2, NBARE2, NOSOA2,
     %                     NUEFA2, NBEFT2, NUTGL2,
     %                     NUDEST, NUSEN1, NBSOSU, XYZS  ,
     %                     NUDETG, NBTGSU, XYZT,
     %                     NUDEEF, NUTYMS, NBFASU, NOSOFA,
     %                     NUEFAP, NBEFTG, NUTGFA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LA TRIANGULATION-QUADRANGULATION D'UNE SURFACE
C -----    PAR JONCTION DES ARETES DE LIGNES 2 A 2
C ENTREES:
C --------
C NBSOM1 : NOMBRE DE SOMMETS DE LA LIGNE 1
C XYZ1   : 3 COORDONNEES DES NBSOM1 SOMMETS DE LA LIGNE 1
C NBTGL1 : NOMBRE DE TANGENTES STOCKEES DE LA LIGNE 1
C XYZTG1 : 3 COMPOSANTES DES NBTGL1 TANGENTES STOCKEES DE LA LIGNE 1
C NUTYM1 : CODE DE MAILLAGE DE LA LIGNE 1 (0 NON STRUCTURE, 2 STRUCTURE)
C NUTFM1 : Ligne Fermee ou NON-Ferme
C          -1: inconnu, 0: non-fermee , 1: ligne fermee
C NBARE1 : NOMBRE D'ARETES DE LA LIGNE 1
C NOSOA1 : NUMERO DES 2 SOMMETS DES NBARE1 ARETES DE LA LIGNE 1
C NUEFA1 : TABLEAU POINTEUR SUR LES ARETES AVEC TANGENTES DE LA LIGNE 1
C NBEFT1 : NOMBRE D'ARETES AVEC TANGENTES DE LA LIGNE 1
C NUTGL1 : NUMERO DES 2 TANGENTES DES ARETES AVEC TANGENTES DE LA LIGNE 1
C
C NBSOM2 : NOMBRE DE SOMMETS DE LA LIGNE 2
C XYZ2   : 3 COORDONNEES DES NBSOM2 SOMMETS DE LA LIGNE 2
C NBTGL2 : NOMBRE DE TANGENTES STOCKEES DE LA LIGNE 2
C XYZTG2 : 3 COMPOSANTES DES NBTGL2 TANGENTES STOCKEES DE LA LIGNE 2
C NUTYM2 : CODE DE MAILLAGE DE LA LIGNE 2 (0 NON STRUCTURE, 2 STRUCTURE)
C NUTFM2 : Ligne Fermee ou NON-Ferme
C          -1: inconnu, 0: non-fermee , 1: ligne fermee
C NBARE2 : NOMBRE D'ARETES DE LA LIGNE 2
C NOSOA2 : NUMERO DES 2 SOMMETS DES NBARE2 ARETES DE LA LIGNE 2
C NUEFA2 : TABLEAU POINTEUR SUR LES ARETES AVEC TANGENTES DE LA LIGNE 2
C NBEFT2 : NOMBRE D'ARETES AVEC TANGENTES DE LA LIGNE 2
C NUTGL2 : NUMERO DES 2 TANGENTES DES ARETES AVEC TANGENTES DE LA LIGNE 2

C SORTIES:
C --------
C NUDEST : NUMERO DU DERNIER SOMMET AVANT LA LIGNE 1 DANS LA SURFACE
C NUSEN1 : SENS DE PARCOURS DE LA LIGNE 1
C NBSOSU : NOMBRE DE SOMMETS DU MAILLAGE DE LA SURFACE
C XYZS   : 3 COORDONNEES DES NBSOSU SOMMETS DE LA SURFACE
C NUDETG : NUMERO DE LA DERNIERE TANGENTE AVANT CETTE COUCHE
C NBTGSU : NOMBRE TOTAL DE TANGENTES DE LA SURFACE
C XYZT   : 3 COMPOSANTES DES NBTGSU TANGENTES DE LA SURFACE
C NUDEEF : NUMERO DU DERNIER ELEMENT FINI AVANT CETTE COUCHE A AJOUTER
C NUTYMS : NUMERO DU TYPE DU MAILLAGE (4: STRUCTURE, 0:NON STRUCTURE)
C NBFASU : NOMBRE DE FACES DE LA SURFACE
C NOSOFA : 4 NUMEROS DES SOMMETS DES NBFASU FACES DE LA SURFACE
C NUEFAP : POINTEUR SUR LES EF A TG (EXISTE SI NBTGSU>0)
C NBEFTG : NOMBRE D'EF A TANGENTES DE LA SURFACE
C NUTGFA : NUMERO DES 8 TANGENTES DES EF DE LA SURFACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1996
C23456...............................................................012
      REAL     XYZ1(3,NBSOM1),XYZ2(3,NBSOM2),
     %         XYZS(3,NBSOSU),XYZT(3,NBTGSU),
     %         XYZTG1(3,NBTGL1),XYZTG2(3,NBTGL2)
      INTEGER  NOSOA1(2,NBARE1),NOSOA2(2,NBARE2),NOSOFA(4,NBFASU)
      INTEGER  NUEFA1(NBARE1),NUTGL1(2,NBEFT1),
     %         NUEFA2(NBARE2),NUTGL2(2,NBEFT2),
     %         NUEFAP(NBEFTG),NUTGFA(8,NBEFTG)

C     RECHERCHE D'INVERSION OU NON DE LA SECONDE LIGNE
C     NUSEN* EST LE SENS DE PARCOURS DE LA LIGNE
      IF( NUTFM1 .EQ. 1  .OR. NUTFM2 .EQ. 1 ) THEN
         NUSEN2 = NUSEN1
      ELSE
         D11 = DIST2P( XYZ1(1,1), XYZ2(1,1) )
         D1N = DIST2P( XYZ1(1,1), XYZ2(1,NBSOM2) )
         DN1 = DIST2P( XYZ1(1,NBSOM1), XYZ2(1,1) )
         DNN = DIST2P( XYZ1(1,NBSOM1), XYZ2(1,NBSOM2) )
         D   = MIN( D11, D1N, DN1, DNN )
         IF( D .EQ. D11 .OR. D .EQ. DNN ) THEN
            NUSEN2 = NUSEN1
         ELSE
            NUSEN2 = -NUSEN1
         ENDIF
      ENDIF

C     LES COORDONNEES DES SOMMETS ET DES TANGENTES
      IF( NUDEST .EQ. 0 ) THEN
C        LES SOMMETS DE LA 1ERE LIGNE SONT STOCKES
         CALL TRTATA( XYZ1,   XYZS, 3*NBSOM1 )
         IF( NBTGL1 .GT. 0 ) THEN
C           LES TANGENTES DE LA 1ERE LIGNE SONT STOCKES
            CALL TRTATA( XYZTG1, XYZT, 3*NBTGL1 )
         ENDIF
C        PAS D'EF A TG RECENSES
         NBEFTG = 0
C        LE DERNIER EF A TG AVANT CETTE COUCHE D'EF
         NUDETG = 0
      ENDIF
      DO 20 J=1,NBSOM2
         IF( NUSEN2 .LT. 0 ) THEN
C           LE NUMERO DU SOMMET EN RECULANT
            KL = NBSOM2 + 1 - J
         ELSE
            KL = J
         ENDIF
         DO 10 I=1,3
            XYZS(I,NUDEST+NBSOM1+J) = XYZ2(I,KL)
 10      CONTINUE
 20   CONTINUE
      CALL TRTATA( XYZTG2, XYZT(1,NUDETG+NBTGL1+1), 3*NBTGL2 )
C
C     MAILLAGE NON STRUCTURE DE LA SURFACE
C     ====================================
      NS1 = 0
      NS2 = 0
      NS3 = 0
      NS4 = 0
      NA1 = 0
      NA2 = 0
      NBA = MAX ( NBARE1, NBARE2 )
      DO 100 J=1,NBA
C
         IF( NUTYMS .EQ. 0 ) THEN
C
C           MAILLAGE NON STRUCTURE DE LA SURFACE
C           LE NUMERO NS1 NS2 DES 2 SOMMETS DE CETTE ARETE SUR LA LIGNE 1
            IF( J .LE. NBARE1 ) THEN
               IF( NUTYM1 .EQ. 2 ) THEN
C                 SEGMENT STRUCTURE LES XYZ ONT ETE REMIS DANS LE BON SENS
                  NS1 = J
                  NS2 = J + 1
                  IF( NUSEN1 .LT. 0 ) THEN
C                    INVERSION DU SENS DE PARCOURS
                     NA1 = NBARE1 + 1 - J
                  ENDIF
               ELSE
C                 LIGNE NON STRUCTUREE
                  IF( NUSEN1 .GT. 0 ) THEN
C                    SENS DE PARCOURS CLASSIQUE
                     NS1 = NOSOA1(1,J)
                     NS2 = NOSOA1(2,J)
                  ELSE
C                    INVERSION DU SENS DE PARCOURS
                     NA1 = NBARE1 + 1 - J
                     NS1 = NBSOM1 + 1 - NOSOA1(2,NA1)
                     NS2 = NBSOM1 + 1 - NOSOA1(1,NA1)
                  ENDIF
               ENDIF
            ENDIF
C
C           LE NUMERO NS3 NS4 DES 2 SOMMETS DE CETTE ARETE SUR LA LIGNE 2
            IF( J .LE. NBARE2 ) THEN
               IF( NUTYM2 .EQ. 2 ) THEN
C                 SEGMENT STRUCTURE LES XYZ ONT ETE REMIS DANS LE BON SENS
                  NS3 = J
                  NS4 = J + 1
                  IF( NUSEN2 .LT. 0 ) THEN
C                    INVERSION DU SENS DE PARCOURS
                     NA2 = NBARE2 + 1 - J
                  ENDIF
               ELSE
C                 LIGNE NON STRUCTUREE
                  IF( NUSEN2 .GT. 0 ) THEN
C                    SENS DE PARCOURS CLASSIQUE
                     NS3 = NOSOA2(1,J)
                     NS4 = NOSOA2(2,J)
                  ELSE
C                    INVERSION DU SENS DE PARCOURS
                     NA2 = NBARE2 + 1 - J
                     NS3 = NBSOM2 + 1 - NOSOA2(2,NA2)
                     NS4 = NBSOM2 + 1 - NOSOA2(1,NA2)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
C        LE NUMERO NT1 NT2 DES 2 TANGENTES DE L'ARETE DE LA LIGNE 1
         NT1 = 0
         NT2 = 0
         IF( NBEFT1 .GT. 0 ) THEN
            IF( NUSEN1 .GT. 0 ) THEN
C              SENS DE PARCOURS CLASSIQUE
               NUEFTG = NUEFA1(J)
               IF( NUEFTG .GT. 0 ) THEN
                  NT1 = NUTGL1(1,NUEFTG)
                  NT2 = NUTGL1(2,NUEFTG)
               ENDIF
            ELSE
C              INVERSION DU SENS DE PARCOURS
               NUEFTG = NUEFA1(NA1)
               IF( NUEFTG .GT. 0 ) THEN
                  NT1 = NUTGL1(2,NUEFTG)
                  NT2 = NUTGL1(1,NUEFTG)
               ENDIF
            ENDIF
         ENDIF
C
C        LE NUMERO NT3 NT4 DES 2 TANGENTES DE L'ARETE DE LA LIGNE 2
         NT3 = 0
         NT4 = 0
         IF( NBEFT2 .GT. 0 ) THEN
            IF( NUSEN2 .GT. 0 ) THEN
C              SENS DE PARCOURS CLASSIQUE
               NUEFTG = NUEFA2(J)
               IF( NUEFTG .GT. 0 ) THEN
                  NT3 = NUTGL2(1,NUEFTG)
                  NT4 = NUTGL2(2,NUEFTG)
               ENDIF
            ELSE
C              INVERSION DU SENS DE PARCOURS
               NUEFTG = NUEFA2(NA2)
               IF( NUEFTG .GT. 0 ) THEN
                  NT3 = NUTGL2(2,NUEFTG)
                  NT4 = NUTGL2(1,NUEFTG)
               ENDIF
            ENDIF
         ENDIF
C
         IF( NBA .EQ. NBARE1 .AND. NBA .EQ. NBARE2 ) THEN
C
C           LE QUADRANGLE NS 1 2 3 4
C           ------------------------
            NUDEEF = NUDEEF + 1
            IF( NUTYMS .EQ. 0 ) THEN
C              LE QUADRANGLE NS 1 2 3 4
               NOSOFA(1,NUDEEF) = NS1 + NUDEST
               NOSOFA(2,NUDEEF) = NS2 + NUDEST
               NOSOFA(3,NUDEEF) = NS4 + NUDEST + NBSOM1
               NOSOFA(4,NUDEEF) = NS3 + NUDEST + NBSOM1
            ENDIF
C           LES EVENTUELLES 8 TANGENTES DU QUADRANGLE
            IF( NBTGSU .GT. 0 ) THEN
               IF( NT1.NE.0 .OR. NT2.NE.0 .OR.
     %             NT3.NE.0 .OR. NT4.NE.0  ) THEN
C                 IL EXISTE AU MOINS UNE TANGENTE POUR CET EF
                  NBEFTG = NBEFTG + 1
                  NUEFAP(NUDEEF) = NBEFTG
                  IF( NT1 .GE. 0 ) THEN
                     NUTGFA(1,NBEFTG) = NT1 + NUDETG
                  ELSE
                     NUTGFA(1,NBEFTG) = NT1 - NUDETG
                  ENDIF
                  NUTGFA(2,NBEFTG) = 0
                  NUTGFA(3,NBEFTG) = 0
                  IF( NT2 .GE. 0 ) THEN
                     NUTGFA(4,NBEFTG) = NT2 + NUDETG
                  ELSE
                     NUTGFA(4,NBEFTG) = NT2 - NUDETG
                  ENDIF
                  IF( NT4 .GE. 0 ) THEN
                     NUTGFA(5,NBEFTG) = NT4 + NUDETG + NBTGL1
                  ELSE
                     NUTGFA(5,NBEFTG) = NT4 - NUDETG - NBTGL1
                  ENDIF
                  NUTGFA(6,NBEFTG) = 0
                  NUTGFA(7,NBEFTG) = 0
                  IF( NT3 .GE. 0 ) THEN
                     NUTGFA(8,NBEFTG) = NT3 + NUDETG + NBTGL1
                  ELSE
                     NUTGFA(8,NBEFTG) = NT3 - NUDETG - NBTGL1
                  ENDIF
               ELSE
C                 EF SANS TG
                  NUEFAP(NUDEEF) = 0
               ENDIF
            ENDIF
            GOTO 100
         ENDIF
C
         IF( NBARE1 .GT. NBARE2 ) THEN
C
C           LES TRIANGLES NS 1 2 3 ET NS 3 2 4
C           ----------------------------------
            NUDEEF = NUDEEF + 1
            IF( NUTYMS .EQ. 0 ) THEN
C              LE TRIANGLE NS 1 2 3
               IF( J .EQ. NBA ) NS3 = NS4
               NOSOFA(1,NUDEEF) = NS1 + NUDEST
               NOSOFA(2,NUDEEF) = NS2 + NUDEST
               NOSOFA(3,NUDEEF) = NS3 + NUDEST + NBSOM1
               NOSOFA(4,NUDEEF) = 0
            ENDIF
C           LES EVENTUELLES 8 TANGENTES DU TRIANGLE NS 1 2 3
            IF( NBTGSU .GT. 0 ) THEN
               IF( NT1.NE.0 .OR. NT2.NE.0 ) THEN
C                 IL EXISTE AU MOINS UNE TANGENTE POUR CET EF
                  NBEFTG = NBEFTG + 1
                  NUEFAP(NUDEEF) = NBEFTG
                  IF( NT1 .GE. 0 ) THEN
                     NUTGFA(1,NBEFTG) = NT1 + NUDETG
                  ELSE
                     NUTGFA(1,NBEFTG) = NT1 - NUDETG
                  ENDIF
                  NUTGFA(2,NBEFTG) = 0
                  NUTGFA(3,NBEFTG) = 0
                  IF( NT2 .GE. 0 ) THEN
                     NUTGFA(4,NBEFTG) = NT2 + NUDETG
                  ELSE
                     NUTGFA(4,NBEFTG) = NT2 - NUDETG
                  ENDIF
                  NUTGFA(5,NBEFTG) = 0
                  NUTGFA(6,NBEFTG) = 0
                  NUTGFA(7,NBEFTG) = 0
                  NUTGFA(8,NBEFTG) = 0
               ELSE
C                 EF SANS TG
                  NUEFAP(NUDEEF) = 0
               ENDIF
            ENDIF
C
            IF( J .EQ. NBA ) GOTO 100
C
C           LE TRIANGLE NS 3 2 4
            NUDEEF = NUDEEF + 1
            IF( NUTYMS .EQ. 0 ) THEN
               NOSOFA(1,NUDEEF) = NS3 + NUDEST + NBSOM1
               NOSOFA(2,NUDEEF) = NS2 + NUDEST
               NOSOFA(3,NUDEEF) = NS4 + NUDEST + NBSOM1
               NOSOFA(4,NUDEEF) = 0
            ENDIF
C           LES EVENTUELLES 8 TANGENTES DU TRIANGLE NS 3 2 4
            IF( NBTGSU .GT. 0 ) THEN
               IF( NT3.NE.0 .OR. NT4.NE.0 ) THEN
C                 IL EXISTE AU MOINS UNE TANGENTE POUR CET EF
                  NBEFTG = NBEFTG + 1
                  NUEFAP(NUDEEF) = NBEFTG
                  NUTGFA(1,NBEFTG) = 0
                  IF( NT3 .GE. 0 ) THEN
                     NUTGFA(2,NBEFTG) = NT3 + NUDETG + NBTGL1
                  ELSE
                     NUTGFA(2,NBEFTG) = NT3 - NUDETG - NBTGL1
                  ENDIF
                  NUTGFA(3,NBEFTG) = 0
                  NUTGFA(4,NBEFTG) = 0
                  IF( NT4 .GE. 0 ) THEN
                     NUTGFA(5,NBEFTG) = NT4 + NUDETG + NBTGL1
                  ELSE
                     NUTGFA(5,NBEFTG) = NT4 - NUDETG - NBTGL1
                  ENDIF
                  NUTGFA(6,NBEFTG) = 0
                  NUTGFA(7,NBEFTG) = 0
                  NUTGFA(8,NBEFTG) = 0
               ELSE
C                 EF SANS TG
                  NUEFAP(NUDEEF) = 0
               ENDIF
            ENDIF
C
         ELSE
C
C           LES TRIANGLES NS 3 1 4 ET NS 4 1 2
C           ----------------------------------
            NUDEEF = NUDEEF + 1
            IF( NUTYMS .EQ. 0 ) THEN
C              LE TRIANGLE NS 3 1 4
               IF( J .EQ. NBA ) NS1 = NS2
               NOSOFA(1,NUDEEF) = NS3 + NUDEST + NBSOM1
               NOSOFA(2,NUDEEF) = NS1 + NUDEST
               NOSOFA(3,NUDEEF) = NS4 + NUDEST + NBSOM1
               NOSOFA(4,NUDEEF) = 0
            ENDIF
C           LES EVENTUELLES 8 TANGENTES DU TRIANGLE NS 3 1 4
            IF( NBTGSU .GT. 0 ) THEN
               IF( NT3.NE.0 .OR. NT4.NE.0 ) THEN
C                 IL EXISTE AU MOINS UNE TANGENTE POUR CET EF
                  NBEFTG = NBEFTG + 1
                  NUEFAP(NUDEEF) = NBEFTG
                  NUTGFA(1,NBEFTG) = 0
                  IF( NT3 .GE. 0 ) THEN
                     NUTGFA(2,NBEFTG) = NT3 + NUDETG + NBTGL1
                  ELSE
                     NUTGFA(2,NBEFTG) = NT3 - NUDETG - NBTGL1
                  ENDIF
                  NUTGFA(3,NBEFTG) = 0
                  NUTGFA(4,NBEFTG) = 0
                  IF( NT4 .GE. 0 ) THEN
                     NUTGFA(5,NBEFTG) = NT4 + NUDETG + NBTGL1
                  ELSE
                     NUTGFA(5,NBEFTG) = NT4 - NUDETG - NBTGL1
                  ENDIF
                  NUTGFA(6,NBEFTG) = 0
                  NUTGFA(7,NBEFTG) = 0
                  NUTGFA(8,NBEFTG) = 0
               ELSE
C                 EF SANS TG
                  NUEFAP(NUDEEF) = 0
               ENDIF
            ENDIF
C
            IF( J .EQ. NBA ) GOTO 100
C
C           LE TRIANGLE NS 4 1 2
            NUDEEF = NUDEEF + 1
            IF( NUTYMS .EQ. 0 ) THEN
               NOSOFA(1,NUDEEF) = NS4 + NUDEST + NBSOM1
               NOSOFA(2,NUDEEF) = NS1 + NUDEST
               NOSOFA(3,NUDEEF) = NS2 + NUDEST
               NOSOFA(4,NUDEEF) = 0
            ENDIF
C           LES EVENTUELLES 8 TANGENTES DU TRIANGLE NS 4 1 2
            IF( NBTGSU .GT. 0 ) THEN
               IF( NT1.NE.0 .OR. NT2.NE.0 ) THEN
C                 IL EXISTE AU MOINS UNE TANGENTE POUR CET EF
                  NBEFTG = NBEFTG + 1
                  NUEFAP(NUDEEF) = NBEFTG
                  NUTGFA(1,NBEFTG) = 0
                  NUTGFA(2,NBEFTG) = 0
                  IF( NT1 .GE. 0 ) THEN
                     NUTGFA(3,NBEFTG) = NT1 + NUDETG
                  ELSE
                     NUTGFA(3,NBEFTG) = NT1 - NUDETG
                  ENDIF
                  NUTGFA(4,NBEFTG) = 0
                  NUTGFA(5,NBEFTG) = 0
                  IF( NT2 .GE. 0 ) THEN
                     NUTGFA(6,NBEFTG) = NT2 + NUDETG
                  ELSE
                     NUTGFA(6,NBEFTG) = NT2 - NUDETG
                  ENDIF
                  NUTGFA(7,NBEFTG) = 0
                  NUTGFA(8,NBEFTG) = 0
               ELSE
C                 EF SANS TG
                  NUEFAP(NUDEEF) = 0
               ENDIF
            ENDIF
         ENDIF
  100 CONTINUE

C     MISE A JOUR DES VARIABLES
      NUDEST = NUDEST + NBSOM1
      NUSEN1 = NUSEN2
      NUDETG = NUDETG + NBTGL1

      RETURN
      END
