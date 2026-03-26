      SUBROUTINE AFDEPL( NBNOMX, NUMCAS, MNXYZN,
     %                   NBCODE, NBNOEU, NDSM,   DEPLAC,
     %                   DECMAX, NOFOTI, DEXMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES DEPLACEMENTS CALCULES DU CAS NUMCAS PARMI NDSM
C -----    LES DEPLACEMENTS EXACTS ET LES ERREURS SI LA FONCTION
C          UTILISATEUR  DEPLACEMENT_EXACT(TEMPS,X,Y,Z,NC) EXISTE
C
C ENTREES:
C --------
C NBNOMX : NUMERO MAXIMAL DES NOEUDS DE DEPLACEMENT A AFFICHER
C NUMCAS : NUMERO DE LA CARTE DES DEPLACEMENTS A AFFICHER
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NBCODE : NOMBRE DE COMPOSANTES DU DEPLACEMENT 1 2 OU 3
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE DE L'OBJET
C NDSM   : NOMBRE DE CARTES DE DEPLACEMENTS
C DEPLAC : LES NDSM CARTES DES DEPLACEMENTS AUX NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C DECMAX : NORME DU DEPLACEMENT MAXIMAL CALCULE POUR CETTE CARTE NUMCAS
C NOFOTI : 0 SI PAS DE FONCTION  DEPLACEMENT_EXACT(TEMPS,X,Y,Z,NC)
C          NO DE LA FONCTION DANS SON LEXIQUE SINON
C DEXMAX : NORME DU DEPLACEMENT MAXIMAL EXACT POUR CETTE CARTE NUMCAS
C          0D0 SI NOFOTI=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1998
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/ctemps.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      DOUBLE PRECISION  DEPLAC(NBCODE,NBNOEU,NDSM)
      DOUBLE PRECISION  DPARAF(5),DEPEXA(3),DIFF(3),
     %                  DEXMAX,DECMAX,PCERMX,ERRMAX,
     %                  DEX,DEC,DIF,PROSCD
C
      DEXMAX = 0D0
      DECMAX = 0D0
      PCERMX = 0D0
      ERRMAX = 0D0
C
C     L'ADRESSE MCN DES COORDONNEES DES NOEUDS DU MAILLAGE DE L'OBJET
      MN = MNXYZN + WYZNOE - 3
C
C     EXISTENCE OU NON DE LA FONCTION 'REGION'
      CALL LXNMNO( NTFONC, 'REGION', NOFORE, I )
C     NOFORE>0 SI CETTE FONCTION EXISTE
C
C     EXISTENCE OU NON DE LA FONCTION 'DEPLACEMENT_EXACT'
      NOFOTI = NOFODEEX()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
C
      IF( NOFOTI .GT. 0 ) THEN
C
C        AFFICHAGE DE LA NORME DU VECTEUR DEPLACEMENT AVEC LES ERREURS
C        =============================================================
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) NBNOEU, NUMCAS
         ELSE
            WRITE(IMPRIM,20000) NBNOEU, NUMCAS
         ENDIF
C
10000 FORMAT(/'Le DEPLACEMENT des',I6,' PREMIERS NOEUDS du VECTEUR DEPLA
     %CEMENT',I6/
     %'LES NORMES DU DEPLACEMENT CALCULE, DEPLACEMENT EXACT ET ERREURS'/
     %100('='))
20000 FORMAT(/'The DISPLACEMENT of',I6,
     %' FIRST NODES of the DISPLACEMENT VECTOR',I6/
     %'The NORMS of COMPUTED DISPLACEMENT, EXACT DISPLACEMENT and ERRORS
     %'/100('='))
C
         DO 10 I=1,NBNOEU
C
C           LES 5 PARAMETRES D'APPEL DE LA FONCTION 'DEPLACEMENT_EXACT'
C           LE TEMPS EN 1-ER PARAMETRE
            DPARAF(1) = TEMPS
            MN        = MN  + 3
            DPARAF(2) = RMCN(MN)
            DPARAF(3) = RMCN(MN+1)
            DPARAF(4) = RMCN(MN+2)
C
            DO 5 NC=1,NBCODE
C              LA COMPOSANTE DU DEPLACEMENT EN 5-EME PARAMETRE
               DPARAF(5) = NC
C              FONCTION DEPLACEMENT_EXACT(TEMPS,X,Y,Z,NC)
               CALL FONVAL( NOFOTI, 5, DPARAF, NCODEV, DEPEXA(NC) )
               IF( NCODEV .EQ. 0 ) GOTO 10
               DIFF(NC) = DEPEXA(NC) - DEPLAC(NC,I,NUMCAS)
 5          CONTINUE
C
C           NORME DU DEPLACEMENT EXACT CORRECTEMENT INITIALISE
            DEX = SQRT( PROSCD( DEPEXA, DEPEXA, NBCODE ) )
            DEC = SQRT( PROSCD( DEPLAC(1,I,NUMCAS),
     %                          DEPLAC(1,I,NUMCAS), NBCODE ) )
            DIF = SQRT( PROSCD( DIFF, DIFF, NBCODE ) )
C
            DEXMAX = MAX( DEXMAX, DEX )
            DECMAX = MAX( DECMAX, DEC )
            ERRMAX = MAX( ERRMAX, DIF )
C
            IF( ABS(DEX) .LT. 1D-10 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10010) I, (RMCN(MN+K),K=0,2),
     %                                DEX, DEC, DIF
               ELSE
                  WRITE(IMPRIM,20010) I, (RMCN(MN+K),K=0,2),
     %                                DEX, DEC, DIF
               ENDIF
            ELSE
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10011) I, (RMCN(MN+K),K=0,2),
     %                                DEX, DEC, DIF, DIF/DEX*100
               ELSE
                  WRITE(IMPRIM,20011) I, (RMCN(MN+K),K=0,2),
     %                                DEX, DEC, DIF, DIF/DEX*100
               ENDIF
               PCERMX = MAX( PCERMX, DIF/DEX )
            ENDIF
 10      CONTINUE
C
10010 FORMAT( 'NOEUD',I6,':X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' ||DEP EXACT||=',G15.7,
     %        ' ||DEP CALCUL||=',G15.7,
     %        ' ||Dif||=',G11.3,' = % /0')
20010 FORMAT( 'NODE',I6,':X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' ||EXACT DISP||=',G15.7,
     %        ' ||COMPU DISP||=',G15.7,
     %        ' ||Dif||=',G11.3,' = % /0')
C
10011 FORMAT( 'NOEUD',I6,':X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' ||DEP EXACT||=',G15.7,
     %        ' ||DEP CALCUL||=',G15.7,
     %        ' ||Dif||=',G11.3,' = %',G10.2)
20011 FORMAT( 'NODE',I6,':X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' ||EXACT DISP||=',G15.7,
     %        ' ||COMPU DISP||=',G15.7,
     %        ' ||Dif||=',G11.3,' = %',G10.2)
C
         IF( DECMAX .NE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10012) DEXMAX, DECMAX, ERRMAX,
     %                             100*ERRMAX/DECMAX, 100*PCERMX,NBNOEU
            ELSE
               WRITE(IMPRIM,20012) DEXMAX, DECMAX, ERRMAX,
     %                             100*ERRMAX/DECMAX, 100*PCERMX,NBNOEU
            ENDIF
         ENDIF
C
10012 FORMAT(/
     %'MAX ||DEPLACEMENT(Noeud)||     EXACT         =',G14.6/
     %'MAX ||DEPLACEMENT(Noeud)||     CALCULE       =',G14.6/
     %'MAX ||DEPL EXACT(Noeud) - DEPL CALC(Noeud)|| =',G14.6//
     %'MAX(||DEPL EXACT(N)-DEPL CALC(N)||) / MAX(||DEPL EXACT(N)||) = %'
     %,G10.2/
     %'MAX(||DEPL EXACT(N)-DEPL CALC(N)||  /     ||DEPL EXACT(N)||) = %'
     %,G10.2//
     %'MAX designe le MAXIMUM sur les',i7,' Noeuds du MAILLAGE')
20012 FORMAT(/
     %'MAX ||DISPLACEMENT(Node)||     EXACT       =',G14.6/
     %'MAX ||DISPLACEMENT(Node)||     COMPUTED    =',G14.6/
     %'MAX ||EXACT DISP(Node) - CALC DISP(Node)|| =',G14.6//
     %'MAX(||EXACT DISP(N)-CALC DISP(N)||) / MAX(||EXACT DISP(N)||) = %'
     %,G10.2/
     %'MAX(||EXACT DISP(N)-CALC DISP(N)||  /     ||EXACT DISP(N)||) = %'
     %,G10.2//
     %'MAX is the MAXIMUM on',i7,' mesh Nodes')
CC
      ELSE
C
C        AFFICHAGE DES COMPOSANTES DES DEPLACEMENTS SANS LES ERREURS
C        ===========================================================
         IF( NOFORE .GT. 0 ) THEN
C
C           LIMITATION PAR LA FONCTION 'REGION(t,x,y,z)'
C           --------------------------------------------
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10020) NUMCAS
            ELSE
               WRITE(IMPRIM,20020) NUMCAS
            ENDIF
C
            DO 30 I=1,NBNOEU
C
C              LES 4 PARAMETRES D'APPEL DE LA FONCTION 'CHOIX'
C              LE TEMPS EN 1-ER PARAMETRE
               DPARAF(1) = TEMPS
C              PUIS LES 3 COORDONNEES X Y Z DU NOEUD
               MN        = MN  + 3
               DPARAF(2) = RMCN(MN)
               DPARAF(3) = RMCN(MN+1)
               DPARAF(4) = RMCN(MN+2)
C              FONCTION CHOIX(TEMPS,X,Y,Z)
               CALL FONVAL( NOFORE, 4, DPARAF, NCODEV, DIF )
               IF( NCODEV .NE. 0 .AND. NINT(DIF) .NE. 0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10050) I, (RMCN(MN+K),K=0,2),
     %                    (DEPLAC(K,I,NUMCAS),K=1,NBCODE)
                  ELSE
                     WRITE(IMPRIM,20050) I, (RMCN(MN+K),K=0,2),
     %                    (DEPLAC(K,I,NUMCAS),K=1,NBCODE)
                  ENDIF
               ENDIF
 30         CONTINUE
C
         ELSE
C           LIMITATION AUX NBNOMX PREMIERS NOEUDS
C           -------------------------------------
            NC = MIN( NBNOEU, NBNOMX )
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10040) NUMCAS, NC
            ELSE
               WRITE(IMPRIM,20040) NUMCAS, NC
            ENDIF
            DO 50 I=1,NC
               MN = MN  + 3
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10050) I, (RMCN(MN+K),K=0,2),
     %                 (DEPLAC(K,I,NUMCAS),K=1,NBCODE)
               ELSE
                  WRITE(IMPRIM,20050) I, (RMCN(MN+K),K=0,2),
     %                 (DEPLAC(K,I,NUMCAS),K=1,NBCODE)
               ENDIF
 50         CONTINUE
         ENDIF
C
      ENDIF
C
10020 FORMAT(/'CAS ',I6,
     %' Le DEPLACEMENT dans la REGION DEFINIE par la FONCTION REGION'/
     %80(1H=))
20020 FORMAT(/'CAS ',I6,
     %' The DISPLACEMENT in the REGION DEFINED by the FUNCTION REGION'/
     %80(1H=))
C
10040 FORMAT(/'CAS ',I6,' LE DEPLACEMENT DES',I6,
     %       ' PREMIERS NOEUDS :'/ 80(1H=))
20040 FORMAT(/'CASE ',I6,' The DISPLACEMENT of',I6,
     %       ' FIRST NODES :'/ 80(1H=))
C
10050 FORMAT('NOEUD',I6,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' DEPLACEMENT=',3G15.7)
20050 FORMAT('NODE',I6,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' DISPLACEMENT=',3G15.7)
C
C     RECHERCHE DU POINT DE DEPLACEMENT MAXIMAL
C     =========================================
C     L'ADRESSE MCN DES COORDONNEES DES NOEUDS DU MAILLAGE DE L'OBJET
      MN     = MNXYZN + WYZNOE - 3
      DECMAX = -1D0
      MNMAX  = 0
      DO 120 I=1,NBNOEU
C
         MN = MN  + 3
C
C        NORME DU DEPLACEMENT CALCULE AU NOEUD I
         DEC = PROSCD( DEPLAC(1,I,NUMCAS),
     %                 DEPLAC(1,I,NUMCAS), NBCODE )
         IF( DEC .GT. DECMAX ) THEN
            DECMAX = DEC
            NOEMAX = I
            MNMAX  = MN
         ENDIF
C
 120  CONTINUE
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10120) TEMPS,NUMCAS,SQRT(DECMAX)
         WRITE(IMPRIM,10050) NOEMAX, (RMCN(MNMAX+K),K=0,2),
     %                      (DEPLAC(K,NOEMAX,NUMCAS),K=1,NBCODE)
      ELSE
         WRITE(IMPRIM,20120) TEMPS,NUMCAS,SQRT(DECMAX)
         WRITE(IMPRIM,20050) NOEMAX, (RMCN(MNMAX+K),K=0,2),
     %                      (DEPLAC(K,NOEMAX,NUMCAS),K=1,NBCODE)
      ENDIF
C
10120 FORMAT(/'Au TEMPS',G15.7,' le DEPLACEMENT MAXIMAL du CAS',
     %I6,' VAUT ',G15.7 )
20120 FORMAT(/'At TIME',G15.7,' the MAXIMUM DISPLACEMENT of CASE',
     %I6,' = ',G15.7 )
C
      RETURN
      END
