        SUBROUTINE VOEX13( NTLXVO , LADEFI ,
     %                     NTCUVO , MNCUVO , NTSOCU , MNSOCU , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN PENTAEDRE DROIT
C -----
C
C ENTREES:
C --------
C NTLXVO : NUMERO DU TABLEAU TS DU LEXIQUE DU HEXAEDRE
C LADEFI : TABLEAU DE DEFINITION DU VOLUME PARTITIONNEE
C          CF '~td/d/a_volume__definition'
C
C SORTIES:
C --------
C NTCUVO : NUMERO      DU TMS 'NSEF' DES NUMEROS DES CUBES
C MNCUVO : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES CUBES
C          CF '~td/d/a___nsef'
C NTSOCU : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOCU : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF '~td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : P.JOLY LABORATOIRE ANALYSE NUMERIQUE PARIS   MARS 1989
C.......................................................................
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
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
C
      IERR = 0
C
C     NOMBRE D'ARETES SUR UNE ARETE DE TRIANGLE ET EN HAUTEUR
C     =======================================================
      NBSA = LADEFI(WBACTP)
      NBSH = LADEFI(WBARHP)
      IF( NBSA .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NBSA
         KERR(1) = 'NOMBRE INCORRECT D''ARETES DES FACES TRIANGULAIRES='
     %             //KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE IF( NBSH .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NBSH
         KERR(1) = 'NOMBRE INCORRECT D''ARETES DANS LA HAUTEUR='
     %            //KERR(MXLGER)(1:10)
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS
C     ======================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      NBA = NBSA + 1
      NBH = NBSH + 1
      NBSOM = NBA * ( NBA + 1 ) * NBH / 2
      CALL LXTNDC( NTLXVO , 'XYZSOMMET' , 'MOTS' ,  WYZSOM + 3 * NBSOM )
      CALL LXTSOU( NTLXVO , 'XYZSOMMET' ,  NTSOCU , MNSOCU )
C     LE NOMBRE DE SOMMETS
      MCN( MNSOCU + WNBSOM) = NBSOM
C     CALCUL DES COORDONNEES DES SOMMETS
      MN = MNSOCU + WYZSOM
      DO 30 K=1,NBH
         DO 20 I=1,NBA
            DO 10 J=1,I
               RMCN( MN     ) = I - J
               RMCN( MN + 1 ) = J - 1
               RMCN( MN + 2 ) = K - 1
               MN  = MN + 3
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOCU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOCU + WNBTGS ) = 0
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOCU + WBCOOR ) = 3
      MCN( MNSOCU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     GENERATION DES NSEF
C     ==========================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
      NBCUVO = NBSA * NBSA *  NBSH
      CALL LXTNDC( NTLXVO,'NSEF','ENTIER',WUSOEF+8*NBCUVO)
      CALL LXTSOU( NTLXVO , 'NSEF' ,  NTCUVO , MNCUVO )
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE VOLUME
C     TYPE DE L'OBJET : VOLUME
      MCN( MNCUVO + WUTYOB ) = 4
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNCUVO + WUTFMA ) = -1
C     LE NOMBRE DE SOMMETS PAR CUBE
      MCN( MNCUVO + WBSOEF ) = 8
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : VOLUME C0
      MCN( MNCUVO + WBTGEF ) = 0
      MCN( MNCUVO + WBEFAP ) = 0
      MCN( MNCUVO + WBEFTG ) = 0
C     LE NOMBRE D'EF DU VOLUME
      MCN( MNCUVO + WBEFOB ) = NBSH * ( NBSA ** 2)
C     NUMERO DU TYPE DU MAILLAGE : PENTAEDRE STRUCTURE
      MCN( MNCUVO + WUTYMA ) = 6
C     NBACTP NOMBRE D'ARETES SUR LE COTE DU TRIANGLE
      MCN( MNCUVO + WBARTP ) = NBSA
C     NBACTP NOMBRE D'ARETES SUR LE COTE DU RECTANGLE
      MCN( MNCUVO + WBARZP ) = NBSH
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNCUVO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCUVO + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      END
