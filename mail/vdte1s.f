      SUBROUTINE VDTE1S( NS, N1TETS, NOTETR,  NBFANS, MXTEFA, MNTEFA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO DE TETRAEDRE ET LE NUMERO LOCAL DE TOUTES
C -----    LES FACES DES TETRAEDRES DE SOMMET NS

C ENTREES:
C --------
C NS     : NUMERO DU SOMMET COMMUN A TOUTES LES FACES
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C MODIFIES:
C ---------
C MXTEFA : NOMBRE D'ENTIERS DECLARES DU TABLEAU TEFA
C          ATTENTION CETTE VALEUR PEUT ETRE CHANGEE PAR LA SUBROUTINE

C SORTIES:
C --------
C NBFANS : NOMBRE DE FACES AYANT NS COMME SOMMET
C MNTEFA : ADRESSE MCN DU TABLEAU TEFA
C          TEFA(1,I) = NUMERO DU TETRAEDRE DE SOMMET NS
C          TEFA(2,I) = NUMERO DE 1 A 4 DE LA FACE DE CE TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JUILLET 1992
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*)

C     INITIALISATION EN CAS D'ERREUR ET SORTIE
      NBFANS = 0

C     UN TETRAEDRE DE SOMMET NS
      NT = N1TETS( NS )

      IF( NT .LE. 0 ) THEN
         WRITE(IMPRIM,*) 'vdte1s: SOMMET ',NS,' DANS AUCUN TETRAEDRE'
      ENDIF

      IF( NOTETR(1,NT) .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'vdte1s: SOMMET ',NS,
     %                   ' dans un TETRAEDRE ABANDONNE',NT
         WRITE(IMPRIM,*) 'NOTETR(',NT,')=',(NOTETR(I,NT),I=1,8)
         RETURN
      ENDIF

C     RECHERCHE DU NUMERO DE SOMMET DE NS DANS NT
      DO I=1,4
         IF( NS .EQ. NOTETR(I,NT) ) GOTO 20
      ENDDO
      WRITE(IMPRIM,*) 'vdte1s: SOMMET ',NS,' NON DANS LE TETRAEDRE ',NT
      WRITE(IMPRIM,*) 'NOTETR(',NT,')=',(NOTETR(I,NT),I=1,8)
      RETURN

C     STOCKAGE DES 3 FACES DU TETRAEDRE NT DE SOMMET NS
C     =================================================
C     NS EST LE SOMMET I DU TETRAEDRE NT
C     IL N'APPARTIENT PAS A LA FACE I+1 DU TETRAEDRE NT
 20   NBFANS = 0
      I   = MOD(I,4) + 1
C     LA FACE I NE CONTIENT PAS LE SOMMET NS
      DO 40 J=1,4
         IF( J .NE. I ) THEN
C           EXISTE T IL UN TETRAEDRE OPPOSE ?
            CALL NOFAOP( J, NT, NOTETR, K, NT1 )
            IF( K .LE. 0 ) THEN
C              PAS DE TETRAEDRE AU DELA DE LA FACE J DE NT
C              STOCKAGE DE CETTE FACE
               NT1 = NT
               K   = J
            ENDIF
C
C           LA FACE K DU TETRAEDRE NT1 EST STOCKEE
C           IL Y A DE LA PLACE POUR LES 4 PREMIERS
C           ADRESSE DEBUT DE LA FACE
            MNT = MNTEFA + 2 * NBFANS
            NBFANS = NBFANS + 1
C           LE TETRAEDRE OPPOSE
            MCN( MNT     ) = NT1
C           LE NUMERO LOCAL DE LA FACE DANS NT1
            MCN( MNT + 1 ) = K
         ENDIF
 40   CONTINUE
C
C     PARCOURS DES TETRAEDRES PAR LES FACES DE SOMMET NS
      NBFANST = 1
C
 45   IF( NBFANST .LE. NBFANS ) THEN
C        LA FACE EXISTE C'EST LA FACE K DU TETRAEDRE NT1
         MNT = MNTEFA + 2*NBFANST - 2
         NT1 = MCN(MNT)
         K   = MCN(MNT+1)
C        LA FACE EST TRAITEE
         NBFANST = NBFANST + 1
C
C        LE NUMERO J DU SOMMET NS DANS LE TETRAEDRE NT1
         DO 50 J=1,4
            IF( NOTETR(J,NT1) .EQ. NS ) GOTO 60
 50      CONTINUE
C        LE NUMERO DE LA FACE DE NT1 NE CONTENANT PAS NS
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
               MNTT = MNTEFA
               DO 80 L=1,NBFANS
                  IF( MCN(MNTT) .EQ. NT1 ) THEN
                     IF( MCN(MNTT+1) .EQ. I ) GOTO 90
                  ENDIF
                  IF( MCN(MNTT) .EQ. NT2 ) THEN
                     IF( MCN(MNTT+1) .EQ. I2 ) GOTO 90
                  ENDIF
                  MNTT = MNTT + 2
 80            CONTINUE
C
C              FACE NON RETROUVEE. ELLE EST EMPILEE
               IF( NBFANS*2+2 .GE. MXTEFA ) THEN
                  print *,'vdte1s:',nbf,' faces du sommet ns=',ns,
     %                    ' mxtefa=',mxtefa,' mntefa avant=',mntefa
                CALL TNMCAU('ENTIER',MXTEFA,MXTEFA+4096,2*NBFANS,MNTEFA)
                  MXTEFA = MXTEFA + 4096
                  print *,'vdte1s:',nbf,' faces du sommet ns=',ns,
     %                    ' mxtefa=',mxtefa,' mntefa apres=',mntefa
               ENDIF

               MNT = MNTEFA + 2 * NBFANS
               NBFANS = NBFANS + 1
               MCN(MNT  ) = NT2
               MCN(MNT+1) = I2
            ENDIF
 90      CONTINUE
         GOTO 45
      ENDIF

      RETURN
      END
