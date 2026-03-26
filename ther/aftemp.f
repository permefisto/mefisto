      SUBROUTINE AFTEMP( NTDLMX, NUMCAS, MNXYZN, NTDL,   NDSM,   TEMP,
     %                   TECMOY, TECMIN, IDLMIN, TECMAX, IDLMAX,
     %                   NOFOTI, TEXMIN, TEXMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES TEMPERATURES CALCULEES DU CAS NUMCAS PARMI NDSM
C -----    LES TEMPERATURES EXACTES ET LES ERREURS SI LA FONCTION
C          UTILISATEUR  TEMPERATURE_EXACTE(t,x,y,z) EXISTE

C ENTREES:
C --------
C NTDLMX : NOMBRE MAXIMAL DE TEMPERATURES A AFFICHER POUR CETTE CARTE NUMCAS
C NUMCAS : NUMERO DE LA CARTE DE TEMPERATURE A AFFICHER
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NTDL   : NOMBRE TOTAL DE NOEUDS DU MAILLAGE DE L'OBJET
C NDSM   : NOMBRE DE CARTES DE TEMPERATURES
C TEMP   : LES NDSM CARTES DE TEMPERATURES AUX NOEUDS DU MAILLAGE

C SORTIES:
C --------
C TECMOY : TEMPERATURE  CALCULEE  MOYENNE  POUR CETTE CARTE NUMCAS
C TECMIN : TEMPERATURE  CALCULEE  MINIMALE POUR CETTE CARTE NUMCAS
C IDLMIN : NUMERO DU DL DE LA TEMPERATURE MINIMALE
C TECMAX : TEMPERATURE  CALCULEE  MAXIMALE POUR CETTE CARTE NUMCAS
C IDLMAX : NUMERO DU DL DE LA TEMPERATURE MAXIMALE
C NOFOTI : =0 SI PAS DE FONCTION TEMPERATURE_EXACTE(t,x,y,z) ou
C          >0 NO DE LA FONCTION DANS SON LEXIQUE
C TEXMIN : TEMPERATURE  EXACTE  MINIMALE POUR CETTE CARTE NUMCAS
C          0D0 SI NOFOTI>0
C TEXMAX : TEMPERATURE  EXACTE  MAXIMALE POUR CETTE CARTE NUMCAS
C          0D0 SI NOFOTI>0
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

      DOUBLE PRECISION  TEMP(NTDL,NDSM)
      DOUBLE PRECISION  DPARAF(7),TEMEXA,TECMIN,TECMAX,TEXMIN,TEXMAX,
     %                  TECMOY,ERRMAX,ERP,DINFO,ERRL1,SOLL1,T

      TECMAX = TEMP(1,NUMCAS)
      TECMIN = TECMAX
C     NBCOOR : NOMBRE (3 ou 6) DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )

C     EXISTENCE OU NON DE LA FONCTION 'REGION'
C     ========================================
      CALL LXNMNO( NTFONC, 'REGION', NOFORE, I )
C     NOFORE>0 SI CETTE FONCTION EXISTE

C     EXISTENCE OU NON DE LA FONCTION 'TEMPERATURE_EXACTE'
C     ====================================================
      NOFOTI = NOFOTEEX()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NOFOTI .GT. 0 ) THEN

C        =======================================================
C        AFFICHAGE DE LA TEMPERATURE AUX NOEUDS AVEC LES ERREURS
C        =======================================================
         TEXMIN =  DINFO( 'GRAND' )
         TEXMAX = -TEXMIN
         SOLL1  = 0D0
         ERRL1  = 0D0
         ERRMAX = 0D0
         ERP    = 0D0

         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) NTDL
         ELSE
            WRITE(IMPRIM,20000) NTDL
         ENDIF
10000 FORMAT(/'TEMPERATURE CALCULEE, TEMPERATURE EXACTE, ERREUR aux',
     %         I8,' NOEUDS')
20000 FORMAT(/'COMPUTED TEMPERATURE, EXACT TEMPERATURE, ERROR at',
     %         I8,' NODES')

         DO 10 I=1,NTDL

C           LES 4 PARAMETRES D'APPEL DE LA FONCTION 'TEMPERATURE_EXACTE'
C           LE TEMPS EN 1-ER PARAMETRE PUIS X Y Z
            DPARAF(1) = TEMPS
            MN        = MNXYZN + WYZNOE + NBCOOR * ( I - 1 )
            N         = 1
            DO K=0,NBCOOR-1
               N = N + 1
               DPARAF(N) = RMCN(MN+K)
            ENDDO

C           TEMEXA = TEMPERATURE_EXACTE(TEMPS,X,Y,Z)
            CALL FONVAL( NOFOTI, N, DPARAF, NCODEV, TEMEXA )
            IF( NCODEV .GT. 0 ) THEN

C              TEMPERATURE EXACTE CORRECTEMENT INITIALISEE
               DPARAF(1) = TEMEXA
               TEXMAX = MAX( TEXMAX, TEMEXA )
               TEXMIN = MIN( TEXMIN, TEMEXA )
               SOLL1  = SOLL1 + ABS( TEMEXA )

C              TEMPERATURE CALCULEE
               DPARAF(2) = TEMP(I,NUMCAS)
               TECMAX = MAX( TECMAX, DPARAF(2) )
               TECMIN = MIN( TECMIN, DPARAF(2) )

C              ERREUR AU NOEUD: TEMPERATURE EXACTE - TEMP CALCULEE
               DPARAF(3) = TEMEXA - TEMP(I,NUMCAS)
               ERRMAX = MAX( ERRMAX, ABS(DPARAF(3)) )
               ERRL1  = ERRL1 + ABS(DPARAF(3))

               IF( LANGAG .EQ. 0 ) THEN
C                 AFFICHAGE EN FRANCAIS
                  IF( ABS(TEMEXA) .LT. 1D-10 ) THEN
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,10010) I,(RMCN(MN+K),K=0,2),
     %                                        (DPARAF(K),K=1,3)
                     ELSE
ccc                        IF( ABS(DPARAF(3)) .GT. 1D-10 )
                        WRITE(IMPRIM,10016) I,(RMCN(MN+K),K=0,5),
     %                                        (DPARAF(K),K=1,3)
                     ENDIF
                  ELSE
                     DPARAF(4) = ABS( DPARAF(3) / TEMEXA * 100 )
                     ERP = MAX( ERP, ABS(DPARAF(4)) )
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,10011)I,(RMCN(MN+K),K=0,2),
     %                                       (DPARAF(K),K=1,4)
                     ELSE
ccc                        IF( ABS(DPARAF(3)) .GT. 1D-10 .AND.
ccc     %                      ABS(DPARAF(4)) .GT. 1D-4 ) THEN
                           WRITE(IMPRIM,10017)I,(RMCN(MN+K),K=0,5),
     %                                          (DPARAF(K),K=1,4)
ccc                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
C                 AFFICHAGE EN ANGLAIS
                  IF( ABS(TEMEXA) .LT. 1D-10 ) THEN
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,20010) I,(RMCN(MN+K),K=0,2),
     %                                        (DPARAF(K),K=1,3)
                     ELSE
ccc                        IF( ABS(DPARAF(3)) .GT. 1D-10 )
                        WRITE(IMPRIM,20016) I,(RMCN(MN+K),K=0,5),
     %                                           (DPARAF(K),K=1,3)
                     ENDIF
                  ELSE
                     DPARAF(4) = ABS( DPARAF(3) / TEMEXA * 100 )
                     ERP = MAX( ERP, ABS(DPARAF(4)) )
                     IF( NBCOOR .LE. 3 ) THEN
                        WRITE(IMPRIM,20011)I,(RMCN(MN+K),K=0,2),
     %                                       (DPARAF(K),K=1,4)
                     ELSE
ccc                        IF( ABS(DPARAF(3)) .GT. 1D-10 .AND.
ccc     %                      ABS(DPARAF(4)) .GT. 1D-4 ) THEN
                           WRITE(IMPRIM,20017)I,(RMCN(MN+K),K=0,5),
     %                                          (DPARAF(K),K=1,4)
ccc                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
 10     ENDDO

10010 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,'  % = /0 ')
10011 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,'  % =',G10.2)
10016 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,'  % = /0 ')
10017 FORMAT( 'DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACTE=',G15.7,
     %        ' T CALCUL=',G15.7,
     %        '  Dif=',G11.3,'  % =',G10.2)

20010 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,'  % = /0 ')
20011 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,'  % =',G10.2)
20016 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,'  % = /0 ')
20017 FORMAT( 'DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %        ' U=',G13.5,' V=',G13.5,' W=',G13.5,
     %        ' T EXACT=',G15.7,
     %        ' T COMPUT=',G15.7,
     %        '  Dif=',G11.3,'  % =',G10.2)

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

10012 FORMAT(/'VECTEUR TEMPERATURE NUMERO',I5/
     %'TEMPERATURE  Minimale  EXACTE ET CALCULEE =',2G14.6/
     %'TEMPERATURE  MAXIMALE  EXACTE ET CALCULEE =',2G14.6/
     %'TEMPERATURE MAX - Min  EXACTE ET CALCULEE =',2G14.6//
     %'MAX |T EXACTE(Noeud) - T CALCULEE(Noeud)|  =',G14.6/
     %'MAX |T EXACTE-T CALCULEE|(Noeud)/(MAX T(Noeud)-MIN T(Noeud)) =',
     % G10.2,' %'/
     %'MAX(|T EXACTE(Noeud)-T CALCULEE(Noeud)| / |T EXACTE(Noeud)|) =',
     % G10.2,' %'/
     %'MAX indique ici un MAXIMUM sur les valeurs aux NOEUDS'//
     %'SOM |(T EXACTE-T CALCULEE)(Noeud)| / SOM |T EXACTE(Noeud)| =',
     % G10.2,' %'/
     %'SOM indique ici la SOMME des valeurs en tous les NOEUDS'/)

20012 FORMAT(/'TEMPERATURE VECTOR NUMBER',I5/
     %'Minimum    TEMPERATURE EXACT and COMPUTED =',2G14.6/
     %'MAXIMUM    TEMPERATURE EXACT and COMPUTED =',2G14.6/
     %'MAX - Min  TEMPERATURE EXACT and COMPUTED =',2G14.6//
     %'MAX |T EXACT(Node) - T COMPUTED(Node)|    =',G14.6/
     %'MAX |T EXACT-T COMPUTED|(Node)/(MAX T(Node)-MIN T(Node)) =',
     % G10.2,' %'/
     %'MAX(|T EXACT(Node) -T COMPUTED(Node)| / |T EXACT(Node)|) =',
     % G10.2,' %'/
     %'here MAX means a MAXIMUM on the values at NODES'//
     %'SUM |(T EXACT - T COMPUTED)(Node)| / SUM |T EXACT(Node)| =',
     % G10.2,' %'/
     %'here SUM means the sum of values at ALL NODES'/)

      ELSE

C        =======================================================
C        AFFICHAGE DE LA TEMPERATURE AUX NOEUDS SANS LES ERREURS
C        =======================================================
         TEXMIN = 0D0
         TEXMAX = 0D0

         IF( NOFORE .GT. 0 ) THEN

C           LIMITATION PAR LA FONCTION 'REGION(t,x,y,z)'
C           --------------------------------------------
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10020) NUMCAS
            ELSE
               WRITE(IMPRIM,20020) NUMCAS
            ENDIF

            MN = MNXYZN + WYZNOE - NBCOOR
            DO 30 I=1,NTDL

C              LES 4 PARAMETRES D'APPEL DE LA FONCTION 'CHOIX'
C              LE TEMPS EN 1-ER PARAMETRE
               DPARAF(1) = TEMPS
C              PUIS LES 3 COORDONNEES X Y Z DU NOEUD
               MN = MN + NBCOOR
               N  = 1
               DO K=0,NBCOOR-1
                  N = N + 1
                  DPARAF(N) = RMCN(MN+K)
               ENDDO
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
 30         ENDDO

         ELSE

C           LIMITATION AUX NTDLMX DL AU DELA DE NTDL/2
C           ------------------------------------------
C           NUMEROS EXTREMES DES DL TEMPERATURE A AFFICHER
            NDL1 = NTDL/2 + 1
            NDL2 = MIN( NTDL/2+NTDLMX, NTDL )

            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10049) NUMCAS, NDL1, NDL2
            ELSE
               WRITE(IMPRIM,20049) NUMCAS, NDL1, NDL2
            ENDIF
            DO I=NDL1,NDL2
               MN = MNXYZN + WYZNOE + (I-1) * NBCOOR
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10050) I,(RMCN(MN+K),K=0,2),
     %                               (TEMP(I,K),K=1,NDSM)
               ELSE
                  WRITE(IMPRIM,20050) I,(RMCN(MN+K),K=0,2),
     %                               (TEMP(I,K),K=1,NDSM)
               ENDIF
            ENDDO

         ENDIF

      ENDIF

10020 FORMAT('CAS',I5,
     %' TEMPERATURES dans la REGION DEFINIE par la FONCTION ''REGION''')
10049 FORMAT('CAS',I5,' TEMPERATURE DES NOEUDS',I8,' A',I8)
10050 FORMAT('DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %      (T58,': TEMPERATURE=',5G15.7))

20020 FORMAT('CASE',I5,
     %' TEMPERATURES in the REGION DEFINED by the ''REGION'' FUNCTION')
20049 FORMAT('CASE',I5,' TEMPERATURE NODES',I8,' to',I8)
20050 FORMAT('DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %      (T58,': TEMPERATURE=',5G15.7))

C     RECHERCHE DU POINT DE TEMPERATURE MINIMALE ET MAXIMALE
C     ======================================================
      TECMOY = 0D0
      TECMAX = TEMP(1,NUMCAS)
      TECMIN = TECMAX
      IDLMIN = 1
      IDLMAX = 1
      DO I=1,NTDL
         T = TEMP(I,NUMCAS)
         TECMOY = TECMOY + T
         IF( T .GT. TECMAX ) THEN
             TECMAX = T
             IDLMAX = I
         ENDIF
         IF( T .LT. TECMIN ) THEN
             TECMIN = T
             IDLMIN = I
         ENDIF
      ENDDO
      TECMOY = TECMOY / NTDL

      MN = MNXYZN + WYZNOE - NBCOOR
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10130) TEMPS,NUMCAS
         IF( NBCOOR .LE. 3 ) THEN
      WRITE(IMPRIM,10131) TECMIN,IDLMIN,(RMCN(MN+NBCOOR*IDLMIN+K),K=0,2)
      WRITE(IMPRIM,10132) TECMOY
      WRITE(IMPRIM,10133) TECMAX,IDLMAX,(RMCN(MN+NBCOOR*IDLMAX+K),K=0,2)
         ELSE
      WRITE(IMPRIM,10136) TECMIN,IDLMIN,(RMCN(MN+NBCOOR*IDLMIN+K),K=0,5)
      WRITE(IMPRIM,10132) TECMOY
      WRITE(IMPRIM,10137) TECMAX,IDLMAX,(RMCN(MN+NBCOOR*IDLMAX+K),K=0,5)
         ENDIF
      ELSE
         WRITE(IMPRIM,20130) TEMPS,NUMCAS
         IF( NBCOOR .LE. 3 ) THEN
      WRITE(IMPRIM,20131) TECMIN,IDLMIN,(RMCN(MN+NBCOOR*IDLMIN+K),K=0,2)
      WRITE(IMPRIM,20132) TECMOY
      WRITE(IMPRIM,20133) TECMAX,IDLMAX,(RMCN(MN+NBCOOR*IDLMAX+K),K=0,2)
         ELSE
      WRITE(IMPRIM,20136) TECMIN,IDLMIN,(RMCN(MN+NBCOOR*IDLMIN+K),K=0,5)
      WRITE(IMPRIM,20132) TECMOY
      WRITE(IMPRIM,20137) TECMAX,IDLMAX,(RMCN(MN+NBCOOR*IDLMAX+K),K=0,5)
         ENDIF
      ENDIF

10130 FORMAT('Au TEMPS',G15.7,' VECTEUR"TEMPERATURE NUMERO =',I5)
10131 FORMAT('TEMPERATURE CALCULEE MINIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
10132 FORMAT('TEMPERATURE CALCULEE MOYENNE  =',G15.7)
10133 FORMAT('TEMPERATURE CALCULEE MAXIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)

10136 FORMAT('TEMPERATURE CALCULEE MINIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
10137 FORMAT('TEMPERATURE CALCULEE MAXIMALE =',G15.7,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)

20130 FORMAT('At TIME',G15.7,' VECTOR"TEMPERATURE NUMBER =',I5)
20131 FORMAT('MINIMUM COMPUTED TEMPERATURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
20132 FORMAT('MEAN    COMPUTED TEMPERATURE =',G15.7)
20133 FORMAT('MAXIMUM COMPUTED TEMPERATURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)

20136 FORMAT('MINIMUM COMPUTED TEMPERATURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
20137 FORMAT('MAXIMUM COMPUTED TEMPERATURE =',G15.7,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %' U=',G13.5,' V=',G13.5,' W=',G13.5)
C
      RETURN
      END
