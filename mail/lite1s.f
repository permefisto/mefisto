      SUBROUTINE LITE1S( NS,     N1TETS, NOTETR,
     %                   NBTET0, NBTET1, MXTEST, NOTEST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LISTER LE NUMERO NOTETR DES TETRAEDRES DE SOMMET NS
C -----
C
C ENTREES:
C --------
C NS     : NUMERO DU SOMMET COMMUN A TOUTES LES FACES
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NBTET0 : NOMBRE INITIAL DE TETRAEDRES STOCKES DANS LE TABLEAU NOTEST

C SORTIES:
C --------
C NBTET1 : NOMBRE FINAL DE TETRAEDRES STOCKES DANS LE TABLEAU NOTEST
C          =0 EN CAS DE SATURATION DE TABLEAU
C MXTEST : NOMBRE MAXIMAL DE NUMEROS  STOCKES DANS LE TABLEAU NOTEST
C NOTEST : NUMERO NOTETR DES TETRAEDRES DE SOMMET NS
C          DERRIERE LES NUMEROS DE TETRAEDRES INITIAUX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  OCTOBRE 2014
C....................................................................012
      PARAMETER        (MXTEFA=1024)
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*),
     %                  NOTEST(MXTEST)
      INTEGER           NOTEFA(2,MXTEFA)
C
C     INITIALISATION EN CAS D'ERREUR ET SORTIE
      NBTET1 = NBTET0
C
C     UN TETRAEDRE DE SOMMET NS
      NT = N1TETS( NS )
      IF( NT .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'LITE1S: SOMMET ',NS,' DANS AUCUN TETRAEDRE'
         NBTET1 = 0
         RETURN
      ENDIF
C
C     RECHERCHE DU NUMERO DE SOMMET DE NS DANS NT
      DO I=1,4
         IF( NS .EQ. NOTETR(I,NT) ) GOTO 10
      ENDDO
      WRITE(IMPRIM,*) 'LITE1S: SOMMET ',NS,' NON DANS LE TETRAEDRE ',NT
      NBTET1 = 0
      RETURN

C     STOCKAGE DES 3 FACES DU TETRAEDRE NT DE SOMMET NS
C     =================================================
C     NS EST LE SOMMET I DU TETRAEDRE NT
C     IL N'APPARTIENT PAS A LA FACE I+1 DU TETRAEDRE NT
 10   NBF = 0
      I   = MOD(I,4) + 1
C     I LA FACE I NE CONTENANT PAS LE SOMMET NS
      DO J=1,4
         IF( J .NE. I ) THEN
C           EXISTE T IL UN TETRAEDRE OPPOSE ?
            CALL NOFAOP( J, NT, NOTETR, K, NT1 )
            IF( K .LE. 0 ) THEN
C              PAS DE TETRAEDRE AU DELA DE LA FACE J DE NT
C              STOCKAGE DE CETTE FACE
               NT1 = NT
               K   = J
            ENDIF
C           LA FACE K DU TETRAEDRE NT1 EST STOCKEE
C           IL Y A DE LA PLACE POUR LES 4 PREMIERS
C           ADRESSE DEBUT DE LA FACE
            NBF = NBF + 1
C           LE TETRAEDRE OPPOSE
            NOTEFA(1,NBF) = NT1
C           LE NUMERO LOCAL DE LA FACE DANS NT1
            NOTEFA(2,NBF) = K
         ENDIF
      ENDDO
C
C     PARCOURS DES TETRAEDRES PAR LES FACES DE SOMMET NS
      NBFT = 1
C
 45   IF( NBFT .LE. NBF ) THEN

C        LA FACE EXISTE C'EST LA FACE K DU TETRAEDRE NT1
         NT1 = NOTEFA(1,NBFT)
         K   = NOTEFA(2,NBFT)
C        LA FACE EST TRAITEE
         NBFT = NBFT + 1

C        LE NUMERO J DU SOMMET NS DANS LE TETRAEDRE NT1
         DO J=1,4
            IF( NOTETR(J,NT1) .EQ. NS ) GOTO 60
         ENDDO
C        J LE NUMERO DE LA FACE DE NT1 NE CONTENANT PAS NS
 60      J = MOD(J,4)+1
C
C        LES 2 AUTRES EVENTUELLES FACES DE SOMMET NS SONT AJOUTEES
C        SI ELLES N'APPARTIENNENT PAS AUX FACES RECENSEES
C        ET SI ELLES EXISTENT DANS LA TETRAEDRISATION
         DO 90 I=1,4
            IF( I .EQ. K ) GOTO 90
C           LA FACE N'A PAS ETE TRAITEE
            IF( I .EQ. J ) GOTO 90
C           LA FACE I DE NT1 CONTIENT NS
C           LE TETRAEDRE OPPOSE
            CALL NOFAOP( I, NT1, NOTETR, I2, NT2 )
            IF( I2 .LE. 0 ) THEN
C              PAS DE TETRAEDRE AU DELA DE LA FACE I DE NT1
C              STOCKAGE DE CETTE FACE
               NT2 = NT1
               I2  = I
            ENDIF
            IF( I2 .GT. 0 ) THEN
C              CETTE FACE A T ELLE ETE STOCKEE ?
               DO L=1,NBF
                  NT = NOTEFA(1,L)
                  II = NOTEFA(2,L)
                  IF( NT .EQ. NT1 ) THEN
                     IF( II .EQ. I ) GOTO 90
                  ENDIF
                  IF( NT .EQ. NT2 ) THEN
                     IF( II .EQ. I2 ) GOTO 90
                  ENDIF
               ENDDO
C
C              FACE NON RETROUVEE. ELLE EST EMPILEE
               IF( NBF .GE. MXTEFA ) THEN
                  WRITE(IMPRIM,*) 'LITE1S: AUGMENTER MXTEFA=',MXTEFA
                  NBTET1 = 0
                  RETURN
               ENDIF
               NBF = NBF + 1
               NOTEFA(1,NBF) = NT2
               NOTEFA(2,NBF) = I2
            ENDIF
 90      CONTINUE
         GOTO 45
      ENDIF

C     LISTAGE DES TETRAEDRES DE SOMMET NS SANS REPETITION
      DO 100 I=1,NBF
         NT = NOTEFA(1,I)
         DO J=1,NBTET1
            IF( NOTEST(J) .EQ. NT ) GOTO 100
         ENDDO
C        LE TETRAEDRE NT DE SOMMET NS EST AJOUTE A NOTEST
         IF( NBTET1 .GE. MXTEST ) THEN
            WRITE(IMPRIM,*) 'LITE1S: AUGMENTER MXTEST=',MXTEST
            NBTET1 = 0
            RETURN
         ENDIF
         NBTET1 = NBTET1 + 1
         NOTEST( NBTET1 ) = NT
 100  CONTINUE

      RETURN
      END
