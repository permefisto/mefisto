      SUBROUTINE AFPRES( NTDLMX, NUMCAS, MNXYZN,
     %                   NTDL,   NDSM,   TEMP,
     %                   TECMIN, TECMAX, NOFOTI, TEXMIN, TEXMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES PRESSIONS CALCULEES DU CAS NUMCAS PARMI NDSM
C -----    LES PRESSIONS EXACTES ET LES ERREURS SI LA FONCTION
C          UTILISATEUR  PRESSION_EXACTE(t,x,y,z) EXISTE
C
C ENTREES:
C --------
C NTDLMX : NOMBRE MAXIMAL DE PRESSIONS A AFFICHER POUR CETTE CARTE NUMCAS
C NUMCAS : NUMERO DE LA CARTE DE PRESSION A AFFICHER
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NTDL   : NOMBRE TOTAL DE NOEUDS DU MAILLAGE DE L'OBJET
C NDSM   : NOMBRE DE CARTES DE PRESSIONS
C TEMP   : LES NDSM CARTES DE PRESSIONS AUX NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C TECMIN : PRESSION  CALCULEE  MINIMALE POUR CETTE CARTE NUMCAS
C TECMAX : PRESSION  CALCULEE  MAXIMALE POUR CETTE CARTE NUMCAS
C
C NOFOTI : 0 SI PAS DE FONCTION PRESSION_EXACTE(t,x,y,z)
C          NO DE LA FONCTION DANS SON LEXIQUE SINON
C TEXMIN : PRESSION  EXACTE  MINIMALE POUR CETTE CARTE NUMCAS
C          0D0 SI NOFOTI=0
C TEXMAX : PRESSION  EXACTE  MAXIMALE POUR CETTE CARTE NUMCAS
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
      DOUBLE PRECISION  TEMP(NTDL,NDSM)
      DOUBLE PRECISION  DPARAF(7),TEMEXA,TECMIN,TECMAX,TEXMIN,TEXMAX,
     %                  ERRMAX,ERP,DINFO,ERRL1,SOLL1
C
      TECMAX = TEMP(1,NUMCAS)
      TECMIN = TECMAX
C     NBCOOR : NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     EXISTENCE OU NON DE LA FONCTION 'REGION'
C     ========================================
      CALL LXNMNO( NTFONC, 'REGION', NOFORE, I )
C     NOFORE>0 SI CETTE FONCTION EXISTE
C
C     EXISTENCE OU NON DE LA FONCTION 'PRESSION_EXACTE'
C     ====================================================
      CALL LXNMNO( NTFONC, 'PRESSION_EXACTE', NOFOTI, I )
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NOFOTI .GT. 0 ) THEN
C
C        =======================================================
C        AFFICHAGE DE LA PRESSION AUX NOEUDS AVEC LES ERREURS
C        =======================================================
         TEXMIN =  DINFO( 'GRAND' )
         TEXMAX = -TEXMIN
         SOLL1  = 0D0
         ERRL1  = 0D0
         ERRMAX = 0D0
         ERP    = 0D0
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) NTDL
         ELSE
            WRITE(IMPRIM,20000) NTDL
         ENDIF
10000 FORMAT(/'PRESSION CALCULEE, PRESSION EXACTE, ERREUR aux',
     %         I8,' NOEUDS'/100('='))
20000 FORMAT(/'COMPUTED PRESSION, EXACT PRESSION, ERROR at',
     %         I8,' NODES'/100('='))
C
         DO 10 I=1,NTDL
C
C           LES 4 PARAMETRES D'APPEL DE LA FONCTION 'PRESSION_EXACTE'
C           LE TEMPS EN 1-ER PARAMETRE PUIS X Y Z
            DPARAF(1) = TEMPS
            MN        = MNXYZN + WYZNOE + NBCOOR * ( I - 1 )
            N         = 1
            DO 5 K=0,NBCOOR-1
               N = N + 1
               DPARAF(N) = RMCN(MN+K)
 5          CONTINUE
C
C           TEMEXA = PRESSION_EXACTE(TEMPS,X,Y,Z)
            CALL FONVAL( NOFOTI, N, DPARAF, NCODEV, TEMEXA )
            IF( NCODEV .GT. 0 ) THEN
C
C              PRESSION EXACTE CORRECTEMENT INITIALISEE
               DPARAF(1) = TEMEXA
               TEXMAX = MAX( TEXMAX, TEMEXA )
               TEXMIN = MIN( TEXMIN, TEMEXA )
               SOLL1  = SOLL1 + ABS( TEMEXA )
C
C              PRESSION CALCULEE
               DPARAF(2) = TEMP(I,NUMCAS)
               TECMAX = MAX( TECMAX, DPARAF(2) )
               TECMIN = MIN( TECMIN, DPARAF(2) )
C
C              ERREUR AU NOEUD: PRESSION EXACTE - TEMP CALCULEE
               DPARAF(3) = TEMEXA - TEMP(I,NUMCAS)
               ERRMAX = MAX( ERRMAX, ABS(DPARAF(3)) )
               ERRL1  = ERRL1 + ABS(DPARAF(3))
C
               IF( LANGAG .EQ. 0 ) THEN
C                 AFFICHAGE EN FRANCAIS
                  IF( ABS(TEMEXA) .LT. 1D-10 ) THEN
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,10010) I,(RMCN(MN+K),K=0,2),
     %                                        (DPARAF(K),K=1,3)
                     ELSE
                        IF( ABS(DPARAF(3)) .GT. 1D-10 )
     %                     WRITE(IMPRIM,10016) I,(RMCN(MN+K),K=0,5),
     %                                           (DPARAF(K),K=1,3)
                     ENDIF
                  ELSE
                     DPARAF(4) = ABS( DPARAF(3) / TEMEXA * 100 )
                     ERP = MAX( ERP, ABS(DPARAF(4)) )
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,10011)I,(RMCN(MN+K),K=0,2),
     %                                       (DPARAF(K),K=1,4)
                     ELSE
                        IF( ABS(DPARAF(3)) .GT. 1D-10 .AND.
     %                      ABS(DPARAF(4)) .GT. 1D-4 ) THEN
                           WRITE(IMPRIM,10017)I,(RMCN(MN+K),K=0,5),
     %                                          (DPARAF(K),K=1,4)
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
C                 AFFICHAGE EN ANGLAIS
                  IF( ABS(TEMEXA) .LT. 1D-10 ) THEN
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,20010) I,(RMCN(MN+K),K=0,2),
     %                                        (DPARAF(K),K=1,3)
                     ELSE
                        IF( ABS(DPARAF(3)) .GT. 1D-10 )
     %                     WRITE(IMPRIM,20016) I,(RMCN(MN+K),K=0,5),
     %                                           (DPARAF(K),K=1,3)
                     ENDIF
                  ELSE
                     DPARAF(4) = ABS( DPARAF(3) / TEMEXA * 100 )
                     ERP = MAX( ERP, ABS(DPARAF(4)) )
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,20011)I,(RMCN(MN+K),K=0,2),
     %                                       (DPARAF(K),K=1,4)
                     ELSE
                        IF( ABS(DPARAF(3)) .GT. 1D-10 .AND.
     %                      ABS(DPARAF(4)) .GT. 1D-4 ) THEN
                           WRITE(IMPRIM,20017)I,(RMCN(MN+K),K=0,5),
     %                                          (DPARAF(K),K=1,4)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
 10     CONTINUE
C
10010 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,' % = /0 ')
10011 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,' % =',G10.2)
10016 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,' % = /0 ')
10017 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,' % =',G10.2)
C
20010 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,' % = /0 ')
20011 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,' % =',G10.2)
20016 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,' % = /0 ')
20017 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,' % =',G10.2)
C
C        BILAN FINAL
         IF( TEXMAX-TEXMIN .NE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10012) NUMCAS, TEXMIN, TECMIN,
     %                             TEXMAX, TECMAX,
     %                             TEXMAX-TEXMIN, TECMAX-TECMIN,
     %                            ERRMAX,100*ERRMAX/(TEXMAX-TEXMIN),ERP,
     %                             ERRL1/SOLL1*100
            ELSE
               WRITE(IMPRIM,20012) NUMCAS, TEXMIN, TECMIN,
     %                             TEXMAX, TECMAX,
     %                             TEXMAX-TEXMIN, TECMAX-TECMIN,
     %                            ERRMAX,100*ERRMAX/(TEXMAX-TEXMIN),ERP,
     %                             ERRL1/SOLL1*100
            ENDIF
         ENDIF
C
10012 FORMAT(/'VECTEUR PRESSION NUMERO',I5/
     %'PRESSION  Minimale  EXACTE ET CALCULEE =',2G14.6/
     %'PRESSION  MAXIMALE  EXACTE ET CALCULEE =',2G14.6/
     %'PRESSION MAX - Min  EXACTE ET CALCULEE =',2G14.6//
     %'MAX |T EXACTE(Noeud) - T CALCULEE(Noeud)|  =',G14.6/
     %'MAX |T EXACTE-T CALCULEE|(Noeud)/(MAX T(Noeud)-MIN T(Noeud)) =',
     % G10.2,' %'/
     %'MAX(|T EXACTE(Noeud)-T CALCULEE(Noeud)| / |T EXACTE(Noeud)|) =',
     % G10.2,' %'/
     %'MAX indique ici un MAXIMUM sur les valeurs aux NOEUDS'//
     %'SOM |(T EXACTE-T CALCULEE)(Noeud)| / SOM |T EXACTE(Noeud)| =',
     % G10.2,' %'/
     %'SOM indique ici la SOMME des valeurs en tous les NOEUDS'/)
C
20012 FORMAT(/'PRESSURE VECTOR NUMBER',I5/
     %'Minimum    PRESSURE EXACT and COMPUTED =',2G14.6/
     %'MAXIMUM    PRESSURE EXACT and COMPUTED =',2G14.6/
     %'MAX - Min  PRESSURE EXACT and COMPUTED =',2G14.6//
     %'MAX |T EXACT(Node) - T COMPUTED(Node)|    =',G14.6/
     %'MAX |T EXACT-T COMPUTED|(Node)/(MAX T(Node)-MIN T(Node)) =',
     % G10.2,' %'/
     %'MAX(|T EXACT(Node) -T COMPUTED(Node)| / |T EXACT(Node)|) =',
     % G10.2,' %'/
     %'here MAX means a MAXIMUM on the values at NODES'//
     %'SUM |(T EXACT - T COMPUTED)(Node)| / SUM |T EXACT(Node)| =',
     % G10.2,' %'/
     %'here SUM means the sum of values at ALL NODES'/)
C
      ELSE
C
C        =======================================================
C        AFFICHAGE DE LA PRESSION AUX NOEUDS SANS LES ERREURS
C        =======================================================
         TEXMIN = 0D0
         TEXMAX = 0D0
         MN     = MNXYZN + WYZNOE - NBCOOR
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
            DO 30 I=1,NTDL
C
C              LES 4 PARAMETRES D'APPEL DE LA FONCTION 'REGION'
C              LE TEMPS EN 1-ER PARAMETRE
               DPARAF(1) = TEMPS
C              PUIS LES 3 COORDONNEES X Y Z DU NOEUD
               MN        = MN + NBCOOR
               N         = 1
               DO 28 K=0,NBCOOR-1
                  N = N + 1
                  DPARAF(N) = RMCN(MN+K)
 28            CONTINUE
C              FONCTION REGION(TEMPS,X,Y,Z)
               CALL FONVAL( NOFORE, N, DPARAF, NCODEV, ERP )
               IF( NCODEV .NE. 0 .AND. NINT(ERP) .NE. 0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10050) I, (RMCN(MN+K),K=0,2),
     %                                  (TEMP(I,K),K=1,NDSM)
                  ELSE
                     WRITE(IMPRIM,20050) I, (RMCN(MN+K),K=0,2),
     %                                  (TEMP(I,K),K=1,NDSM)
                  ENDIF
               ENDIF
 30         CONTINUE
C
         ELSE
C
C           LIMITATION AUX NTDLI PREMIERS DL
C           --------------------------------
            NTDLI = MIN( NTDL, NTDLMX )
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10049) NUMCAS, NTDLI
            ELSE
               WRITE(IMPRIM,20049) NUMCAS, NTDLI
            ENDIF
            DO 50 I=1,NTDLI
               MN = MN + NBCOOR
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10050) I,(RMCN(MN+K),K=0,2),
     %                               (TEMP(I,K),K=1,NDSM)
               ELSE
                  WRITE(IMPRIM,20050) I,(RMCN(MN+K),K=0,2),
     %                               (TEMP(I,K),K=1,NDSM)
               ENDIF
 50         CONTINUE
C
         ENDIF
C
      ENDIF
C
10020 FORMAT(/'CAS',I5,
     %' PRESSIONS dans la REGION DEFINIE par la FONCTION ''REGION'''/
     %80(1H=))
10049 FORMAT(/'CAS',I5,' PRESSION DES',I8,
     %        ' PREMIERS NOEUDS :'/80('='))
10050 FORMAT('DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %      (T58,': PRESSION =',4G15.7))
C
20020 FORMAT(/'CASE',I5,
     %' PRESSIONS in the REGION DEFINED by the ''REGION'' FUNCTION'/
     %80(1H=))
20049 FORMAT(/'CASE',I5,' PRESSION of',I8,
     %        ' FIRST NODES :'/80('='))
20050 FORMAT('DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %      (T58,': PRESSION =',4G15.7))
C
C     RECHERCHE DU POINT DE PRESSION MINIMALE ET MAXIMALE
C     ======================================================
      TECMAX = TEMP(1,NUMCAS)
      TECMIN = TECMAX
      IMIN   = 1
      IMAX   = 1
      DO 130 I=2,NTDL
         IF( TEMP(I,NUMCAS) .GT. TECMAX ) THEN
             TECMAX = TEMP(I,NUMCAS)
             IMAX = I
         ENDIF
         IF( TEMP(I,NUMCAS) .LT. TECMIN ) THEN
             TECMIN = TEMP(I,NUMCAS)
             IMIN = I
         ENDIF
 130  CONTINUE
C
      MN = MNXYZN + WYZNOE - NBCOOR
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10130) TEMPS,NUMCAS
         IF( NBCOOR .LE. 3 ) THEN
         WRITE(IMPRIM,10131) TECMIN,IMIN,(RMCN(MN+NBCOOR*IMIN+K),K=0,2)
         WRITE(IMPRIM,10132) TECMAX,IMAX,(RMCN(MN+NBCOOR*IMAX+K),K=0,2)
         ELSE
         WRITE(IMPRIM,10136) TECMIN,IMIN,(RMCN(MN+NBCOOR*IMIN+K),K=0,5)
         WRITE(IMPRIM,10137) TECMAX,IMAX,(RMCN(MN+NBCOOR*IMAX+K),K=0,5)
         ENDIF
      ELSE
         WRITE(IMPRIM,20130) TEMPS,NUMCAS
         IF( NBCOOR .LE. 3 ) THEN
         WRITE(IMPRIM,20131) TECMIN,IMIN,(RMCN(MN+NBCOOR*IMIN+K),K=0,2)
         WRITE(IMPRIM,20132) TECMAX,IMAX,(RMCN(MN+NBCOOR*IMAX+K),K=0,2)
         ELSE
         WRITE(IMPRIM,20136) TECMIN,IMIN,(RMCN(MN+NBCOOR*IMIN+K),K=0,5)
         WRITE(IMPRIM,20137) TECMAX,IMAX,(RMCN(MN+NBCOOR*IMAX+K),K=0,5)
         ENDIF
      ENDIF
C
10130 FORMAT(/'Au TEMPS',G15.7,' VECTEUR"PRESSION NUMERO =',I5)
10131 FORMAT('PRESSION CALCULEE MINIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
10132 FORMAT('PRESSION CALCULEE MAXIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
C
10136 FORMAT('PRESSION CALCULEE MINIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
10137 FORMAT('PRESSION CALCULEE MAXIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
C
20130 FORMAT(/'At TIME',G15.7,' VECTOR"PRESSURE NUMBER =',I5)
20131 FORMAT('MINIMUM COMPUTED PRESSURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
20132 FORMAT('MAXIMUM COMPUTED PRESSURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
C
20136 FORMAT('MINIMUM COMPUTED PRESSURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
20137 FORMAT('MAXIMUM COMPUTED PRESSURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
C
      RETURN
      END
