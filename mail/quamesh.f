      SUBROUTINE QUAMESH( TEXT,   ITER,   GRAND,  XYZSOM,
     %                    NBSOTE, MXTETR, NOSOTE,
     %                    VOLUMT, QUALIT, QUAMIN0,
     %                    NBTU,   VOLTOT, QUAMOY, QUAMIN, NBPB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LA QUALITE D'UNE TETRAEDRISATION
C -----
C
C ENTREES:
C --------
C TEXT   : A AFFICHER
C ITER   : NO ITERATION A AFFICHER
C GRAND  : VALEUR TEMOIN DE TETRAEDRE NON UTILISE
C XYZSOM : 3 COORDONNEES DES SOMMETS
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NOSOTE(>3)
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARES
C NOSOTE : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NOSOTE DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C
C SORTIES:
C --------
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C NBTU   : NOMBRE DE TETRAEDRES UTILISES
C VOLTOT : VOLUME DU MAILLAGE
C QUAMOY : QUALITE MOYENNE  DU MAILLAGE
C QUAMIN : QUALITE MINIMALE DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Septembre 2012
C2345X7..............................................................012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     TEXT
      REAL              XYZSOM(3,*), VOLUMT(MXTETR), QUALIT(MXTETR)
      INTEGER           NOSOTE(NBSOTE,*)
      REAL              ARMIN, ARMAX, SURFTR(4)
C
      NBTQ0  = 0
      VOLTOT = 0
      QUAMIN = GRAND
      QUAMOY = 0
      NBTU   = 0
C
      DO NT=1,MXTETR
C
         IF( NOSOTE(1,NT) .GT. 0 ) THEN
C
C           UN TETRAEDRE UTILISE DE PLUS
            NBTU = NBTU + 1
C
C           QUALITE DU TETRAEDRE NT
            CALL QUATET( XYZSOM(1,NOSOTE(1,NT)),
     %                   XYZSOM(1,NOSOTE(2,NT)),
     %                   XYZSOM(1,NOSOTE(3,NT)),
     %                   XYZSOM(1,NOSOTE(4,NT)),
     %                   ARMIN, ARMAX, SURFTR, VOLUMT(NT),QUALIT(NT))
C
C           VOLUME DES TETRAEDRES DU MAILLAGE
            VOLTOT = VOLTOT + VOLUMT(NT)
C
C           QUALITE MOYENNE DES TETRAEDRES DU MAILLAGE
            QUAMOY = QUAMOY + QUALIT(NT)
            IF( QUALIT(NT) .LT. QUAMIN ) QUAMIN = QUALIT(NT)
C
            IF( QUALIT(NT) .LE. 0 ) THEN
               NBTQ0 = NBTQ0 + 1
               print*,'quamesh',ITER,': QUALIT(',NT,')=',QUALIT(NT),
     %                ' VOLUMT(',NT,')=',VOLUMT(NT),' NS=',
     %         (nosote(kk,nt),kk=1,4)
            ENDIF
C
         ELSE
C
C           TETRAEDRE VIDE
            VOLUMT(NT) = 0.0
            QUALIT(NT) = GRAND
C
         ENDIF
      ENDDO
C
C     QUALITE MOYENNE DU MAILLAGE
      QUAMOY = QUAMOY / NBTU
C
C     AFFICHAGE DE LA QUALITE MINIMALE DE L'ITERATION
      NBPB = 0
      IF( QUAMIN0 .GT. QUAMIN ) THEN
         NBPB = 1
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,11030) TEXT,ITER, QUAMIN0, QUAMIN,
     %                       NBTU, VOLTOT, QUAMOY, QUAMIN,
     %                       NBTQ0
         ELSE
         WRITE(IMPRIM,21030) TEXT,ITER, QUAMIN0, QUAMIN,
     %                       NBTU, VOLTOT, QUAMOY, QUAMIN,
     %                       NBTQ0
         ENDIF
      ENDIF
C
11030 FORMAT( A,' tetra',I8,' Regression Qualite de ',G14.6,' a ',G14.6/
     %        ' NOMBRE TETRAEDRES=',I8,
     %        ' VOLUME MAILLAGE=',G15.6,
     %        ' QUALITE MOYENNE=',G15.6,
     %        ' QUALITE MINIMALE=',G15.6,
     %        ' NB TETRAEDRES de QUALITE NULLE=',I8)
     %       
21030 FORMAT(A,'tetra',I8,' Quality Regression from',G14.6,' to ',G14.6/
     %        ' TETRAHEDRA NUMBER=',I8,
     %        ' MESH VOLUME=',G15.6,
     %        ' MEAN QUALITY=',G15.6,
     %        ' MINIMUM QUALITY=',G15.6,
     %        ' NB of NULL QUALITY TETRAHEDRA=',I8)
C
      RETURN
      END
