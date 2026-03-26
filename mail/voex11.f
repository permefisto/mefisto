      SUBROUTINE VOEX11( NTLXVO, LADEFI,
     %                   NTCUVO, MNCUVO, NTSOCU, MNSOCU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN PARALLELIPIPEDE RECTANGLE DE COTES
C -----    NBINHX NBINHY NBINHZ

C ENTREES:
C --------
C NTLXVO : NUMERO DU TABLEAU TS DU LEXIQUE DU HEXAEDRE
C LADEFI : TABLEAU DE DEFINITION DU VOLUME PARTITIONNEE
C          CF '~td/d/a_volume__definition'

C SORTIES:
C --------
C NTCUVO : NUMERO      DU TMS 'NSEF' DES NUMEROS DES CUBES
C MNCUVO : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES CUBES
C          cf '~td/d/a___nsef'
C NTSOCU : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOCU : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF '~td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS   MARS 1989
C.......................................................................
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)

      IERR = 0

C     NOMBRE D'ARETES EN X ET EN Y
C     ============================
      NBINHX = LADEFI(WBINHX)
      NBINHY = LADEFI(WBINHY)
      NBINHZ = LADEFI(WBINHZ)
      IF( NBINHX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBINHX
         KERR(1)='NOMBRE D''ARETES EN X INCORRECT ='//KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBINHY .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBINHY
         KERR(1)='NOMBRE D''ARETES EN Y INCORRECT ='//KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBINHY .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBINHZ
         KERR(1)='NOMBRE D''ARETES EN Z INCORRECT ='//KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

C     GENERATION DES SOMMETS
C     ======================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      NBSOM = ( NBINHX + 1 ) * ( NBINHY + 1 ) * ( NBINHZ + 1 )
      CALL LXTNDC( NTLXVO, 'XYZSOMMET', 'MOTS',  WYZSOM + 3 * NBSOM )
      CALL LXTSOU( NTLXVO, 'XYZSOMMET',  NTSOCU, MNSOCU )
C     LE NOMBRE DE SOMMETS
      MCN( MNSOCU + WNBSOM) = NBSOM
C     CALCUL DES COORDONNEES DES SOMMETS
      MN = MNSOCU + WYZSOM
      DO K=0,NBINHZ
         DO J=0,NBINHY,1
            DO I=0,NBINHX,1
               RMCN( MN     ) = I
               RMCN( MN + 1 ) = J
               RMCN( MN + 2 ) = K
               MN  = MN + 3
            ENDDO
         ENDDO
      ENDDO
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOCU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOCU + WNBTGS ) = 0
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOCU + WBCOOR ) = 3
      MCN( MNSOCU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C     GENERATION DES NSEF
C     ==========================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTNDC( NTLXVO,'NSEF','ENTIER',WBARZH+1)
      CALL LXTSOU( NTLXVO, 'NSEF',  NTCUVO, MNCUVO )
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE VOLUME
C     TYPE DE L'OBJET : VOLUME
      MCN( MNCUVO + WUTYOB ) = 4
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNCUVO + WUTFMA ) = -1
C     NUMERO DU TYPE DU MAILLAGE : HEXAEDRE STRUCTURE
      MCN( MNCUVO + WUTYMA ) = 7
C     NBARXQ NOMBRE D'ARETES EN X
      MCN( MNCUVO + WBARXH ) = NBINHX
C     NBARYQ NOMBRE D'ARETES EN Y
      MCN( MNCUVO + WBARYH ) = NBINHY
C     NBARYQ NOMBRE D'ARETES EN Z
      MCN( MNCUVO + WBARZH ) = NBINHZ
C     LE NOMBRE DE SOMMETS PAR CUBE
      MCN( MNCUVO + WBSOEF ) = 8
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : VOLUME C0
      MCN( MNCUVO + WBTGEF ) = 0
      MCN( MNCUVO + WBEFAP ) = 0
      MCN( MNCUVO + WBEFTG ) = 0
C     LE NOMBRE D'EF DE LA SURFACE
      MCN( MNCUVO + WBEFOB ) = NBINHX * NBINHY * NBINHZ
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNCUVO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCUVO + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

      RETURN
      END
