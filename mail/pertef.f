      SUBROUTINE PERTEF( NCOGEX, NBELEM, MOELEM, MNELEM, XYZ, NBEFPE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES COORDONNEES DES SOMMETS DES EF NON IDENTIFIES
C -----    COMME TRACES D'EF DE DIMENSION SUPERIEURE
C ENTREES:
C --------
C NCOGEX : CODE GEOMETRIQUE MAXIMAL A TRAITER
C          7:HEXAEDRE, ... , 4:QUADRANGLE, 3:TRIANGLE, 2:ARETE, 1:SOMMET
C NBELEM : NOMBRE INITIAL D'EF DE CE TYPE GEOMETRIQUE
C MOELEM : NOMBRE DE MOTS DE CHAQUE TYPE D'ELEMENT
C MNELEM : ADRESSE MCN DES EVENTUELS TABLEAUX DES 9 TYPES D'EF
C          0 SI PAS DE D'ELEMENTS FINIS DE CE TYPE
C          MCN(MNELEM(NCOGEL))=LISTE DE HACHAGE DES EF DE CODE NCOGEL
C          CETTE LISTE A ETE CALCULEE PAR UN CALL DFTOP3
C XYZ    : 3 COORDONNEES DES SOMMETS DU MAILLAGE DE L'OBJET
C
C SORTIES:
C --------
C NBEFPE : NOMBRE EFFECTIF D'EF PERDUS APRES PARCOURS DES TABLES DE HACHAGE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS        MAI 1996
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           MOELEM(9),MNELEM(9),NBELEM(9)
      REAL              XYZ(3,*),X(4),Y(4),Z(4)
C
C     LE NOMBRE D'EF PERDU
      NBEFPE = 0
      NOSOM  = 0
C
      DO 1000 NCOGEL=NCOGEX,1,-1
C
C        LE NOMBRE ACTUEL D'ELEMENTS AVEC CE CODE GEOMETRIQUE
         NBELT = NBELEM(NCOGEL)
         IF( NBELT .LE. 0 ) GOTO 1000
C
C        LE NOMBRE DE SOMMETS DE CE TYPE D'EF
         NBSO = NBSOME( NCOGEL )
C
C        BOUCLE SUR LES ELEMENTS FINIS DE TYPE NCOGEL
         MOEL = MOELEM( NCOGEL )
         MNEL = MNELEM( NCOGEL ) - 1 - MOEL
C
         DO 100 I=1,NBELT
C
C           TOUT ELEMENT FINI DEJA VU OU INEXISTANT EST SAUTE
            MNEL = MNEL + MOEL
            IF( MCN( MNEL + 1 )    .LE. 0 ) GOTO 100
            IF( MCN( MNEL + MOEL ) .LE. 0 ) GOTO 100
C
C           UN ELEMENT FINI EST PERDU
            IF( NBEFPE .EQ. 0 ) THEN
C
C              INITIALISATION DU TRACE
C              AUCUN ITEM SUR L'ECRAN
               CALL EFFACE
               CALL ITEMS0
C              LES PARAMETRES DU CADRE MAXIMAL
               CALL VISEE0
C
            ENDIF
C
C           UN EF PERDU DE PLUS
            NBEFPE = NBEFPE + 1
            IF( NBSO .EQ. 4 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10004) NBEFPE
               ELSE
                  WRITE(IMPRIM,20004) NBEFPE
               ENDIF
10004 FORMAT(' EF PERDU ',I5,': QUADRANGLE de SOMMETS' )
20004 FORMAT(' LOST FE ',I5,': QUADRANGLE of VERTICES' )
            ELSE IF( NBSO .EQ. 3 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10003) NBEFPE
               ELSE
                  WRITE(IMPRIM,20003) NBEFPE
               ENDIF
10003 FORMAT(' EF PERDU ',I5,': TRIANGLE de SOMMETS' )
20003 FORMAT(' LOST FE ',I5,': TRIANGLE of VERTICES' )
            ELSE IF( NBSO .EQ. 2 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10002) NBEFPE
               ELSE
                  WRITE(IMPRIM,20002) NBEFPE
               ENDIF
10002 FORMAT(' EF PERDU ',I5,': ARETE de SOMMETS' )
20002 FORMAT(' LOST FE ',I5,': EDGE of VERTICES' )
            ELSE
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10001) NBEFPE
               ELSE
                  WRITE(IMPRIM,20001) NBEFPE
               ENDIF
10001 FORMAT(' EF PERDU ',I5,': POINT de SOMMET' )
20001 FORMAT(' LOST FE ',I5,': POINT of VERTEX' )
            ENDIF
C
            DO 10 J=1,NBSO
C              LE NUMERO DU SOMMET ET SES XYZ SONT AFFICHES
               NOSOM = MCN(MNEL+J)
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10010) (XYZ(K,NOSOM),K=1,3)
               ELSE
                  WRITE(IMPRIM,20010) (XYZ(K,NOSOM),K=1,3)
               ENDIF
10010 FORMAT(' SOMMET X=',G15.6,' Y=',G15.6,' Z=',G15.6)
20010 FORMAT(' VERTEX X=',G15.6,' Y=',G15.6,' Z=',G15.6)
               X(J) = XYZ(1,NOSOM)
               Y(J) = XYZ(2,NOSOM)
               Z(J) = XYZ(3,NOSOM)
10          CONTINUE
C
C           TRACE DE L'EF PERDU
            IF( NDIMLI .EQ. 3 ) THEN
C
C              TRACE EN DIMENSION 3
               GOTO( 15, 20, 30, 30 ),NBSO
C
C              UN POINT=SOMMET
 15            IF( LANGAG .EQ. 0 ) THEN
                  CALL SYMBOLE3D( NCROUG, XYZ(1,NOSOM), '+ POINT PERDU')
               ELSE
                  CALL SYMBOLE3D( NCROUG, XYZ(1,NOSOM), '+ LOST POINT' )
               ENDIF
               GOTO 90
C
C              UNE ARETE
 20            CALL  TRAIT3D( NCROUG, XYZ(1,MCN(MNEL+1)),
     %                                XYZ(1,MCN(MNEL+2)) )
               GOTO 90
C
C              UN TRIANGLE OU QUADRANGLE
 30            CALL FAP13D( NCROUG, NCNOIR, 0, NBSO, X, Y, Z )
C
            ELSE
C
C              TRACE EN DIMENSION 2
               GOTO( 50, 60, 70, 70 ),NBSO
C
C              UN POINT=SOMMET
 50            IF( LANGAG .EQ. 0 ) THEN
                  CALL SYMBOLE2D( NCROUG, X(1), Y(1), '+ POINT PERDU' )
               ELSE
                  CALL SYMBOLE2D( NCROUG, X(1), Y(1), '+ LOST POINT' )
               ENDIF
               GOTO 90
C
C              UNE ARETE
 60            CALL  TRAIT2D( NCROUG, X(1), Y(1),
     %                                X(2), Y(2) )
               GOTO 90
C
C              UN TRIANGLE OU QUADRANGLE
 70            CALL FACE2D( NCROUG, NCNOIR, NBSO, X, Y )
C
            ENDIF
C
C           POUR NE PLUS REVOIR CET EF
 90         MCN( MNEL + MOEL ) = -MCN( MNEL + MOEL )
C
100      CONTINUE
1000  CONTINUE
C
      IF( NBEFPE .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            CALL TRFINS( 'TRACE des ELEMENTS FINIS NON RETROUVES' )
         ELSE
            CALL TRFINS( 'DRAWING of NOT RECOVERED FINITE ELEMENTS' )
         ENDIF
         CALL CLICSO
      ENDIF
C
      RETURN
      END
