      SUBROUTINE AFNLSE( NBNOAF, NUMCAS, MNXYZN, NBNOMA, NDSM,   WAVE,
     %                   WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LA PARTIE REELLE ET IMAGINAIRE D'UNE ONDE NLSE
C -----    DU CAS NUMCAS PARMI LES NDSM CAS STOKES
C          LA PARTIE REELLE ET IMAGINAIRE EXACTE ET LES ERREURS
C          SI LES FONCTIONS UTILISATEUR
C          PARTIE_REELLE_EXACTE(t,x,y,z) et PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
C          ou Exact_Real_Part(t,x,y,z) et Exact_Imaginary_Part(t,x,y,z)
C          EXISTENT
C
C ENTREES:
C --------
C NBNOAF : NOMBRE MAXIMAL DE DL A AFFICHER POUR CETTE CARTE NUMCAS
C NUMCAS : NUMERO DE LA CARTE DU VECTEUR SOLUTION COMPLEXE A AFFICHER
C
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NBNOMA : NOMBRE TOTAL DE NOEUDS DU MAILLAGE DE L'OBJET
C NDSM   : NOMBRE DE CARTES DE WAVES
C WAVE   : NDSM CARTES DE WAVES COMPLEXES (PR et PI) AUX NOEUDS DU MAILLAGE
C          WAVE(NBNOMA,2,NDSM)
C
C SORTIES:
C --------
C WCAMIN : WAVE PR et PI CALCULEE MINIMALE POUR CETTE CARTE NUMCAS
C WCAMAX : WAVE PR et PI CALCULEE MAXIMALE POUR CETTE CARTE NUMCAS
C
C NOFOWE(1) =0 SI PAS DE FONCTION
C          Partie_Reelle_Exacte(t,x,y,z) ou Exact_Real_Part(t,x,y,z)
C          >0 NO DE LA FONCTION DANS SON LEXIQUE
C NOFOWE(2) =0 SI PAS DE FONCTION
C          Partie_Imaginaire_Exacte(t,x,y,z) ou Exact_Imaginary_Part(t,x,y,z)
C          >0 NO DE LA FONCTION DANS SON LEXIQUE
C WEXMIN : WAVE PR et PI EXACTE MINIMALE POUR CETTE CARTE NUMCAS
C          0D0 SI NOFOWE>0 POUR PR et PI
C WEXMAX : WAVE PR et PI EXACTE MAXIMALE POUR CETTE CARTE NUMCAS
C          0D0 SI NOFOWE>0 POUR PR et PI
C ERRMAX : ERREUR MAX POUR WAVE PR et PI (EXACTE - CALCULEE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
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
      DOUBLE PRECISION  WAVE(NBNOMA,2,NDSM)
      DOUBLE PRECISION  DPARAF(7,2), WAVEXA(2), SOMODU, SOMODIF, MODUMX
      DOUBLE PRECISION  WCAMIN(2), WCAMAX(2),WEXMIN(2), WEXMAX(2),
     %                  ERRMAX(2), ERP(2), DINFO, ERRL1(2), SOLL1(2)
      INTEGER           NOFOWE(2), IMIN(2), IMAX(2)
      INTRINSIC         MAX
C
C     LIMITATION DE L'AFFICHAGE A NBNOF NOEUDS
C     ----------------------------------------
      NBNOF = MIN( NBNOMA, NBNOAF )
      IF( NBNOF .LT. NBNOMA/2 ) THEN
         NBNOM2 = NBNOMA / 2
      ELSE
         NBNOM2 = 1
      ENDIF
C
C     INITIALISATION DES MIN ET MAX DES PR et PI
      DO M=1,2
         WCAMAX(M) = WAVE(1,M,NUMCAS)
         WCAMIN(M) = WCAMAX(M)
      ENDDO
C
C     NBCOOR : NOMBRE 3 DES COORDONNEES DES NOEUDS
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     EXISTENCE OU NON DE LA FONCTION 'REGION'
C     ========================================
      CALL LXNMNO( NTFONC, 'REGION', NOFORE, I )
C     NOFORE>0 SI CETTE FONCTION EXISTE
C
C     EXISTENCE OU NON DES 2 FONCTIONS
C     PARTIE_REELLE_EXACTE(t,x,y,z)     ou EXACT_REAL_PART(t,x,y,z)
C     PARTIE_IMAGINAIRE_EXACTE(t,x,y,z) ou EXACT_IMAGINARY_PART(t,x,y,z)
C     ==================================================================
      MODUMX    = 0D0
      SOMODU    = 0D0
      SOMODIF   = 0D0
      ERRMAX(1) = 0D0
      ERRMAX(2) = 0D0
      WEXMIN(1) = 0D0
      WEXMIN(2) = 0D0
      WEXMAX(1) = 0D0
      WEXMAX(2) = 0D0
      NOFOWE(1) = NOFOPREX()
      NOFOWE(2) = NOFOPIEX()
C     NOFOWE>0 SI CETTE FONCTION EXISTE
      IF( NOFOWE(1) .GT. 0 .AND. NOFOWE(2) .GT. 0 ) THEN
C
C        ================================================
C        AFFICHAGE DE LA WAVE AUX NOEUDS AVEC LES ERREURS
C        ================================================
         DO M=1,2
            WEXMIN(M) =  DINFO( 'GRAND' )
            WEXMAX(M) = -WEXMIN(M)
            SOLL1(M)  = 0D0
            ERRL1(M)  = 0D0
            ERRMAX(M) = 0D0
            ERP(M)    = 0D0
         ENDDO
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) TEMPS, NBNOF, NBNOMA
         ELSE
            WRITE(IMPRIM,20000) TEMPS, NBNOF, NBNOMA
         ENDIF
10000 FORMAT('Au TEMPS',G14.6,' ONDE CALCULEE, ONDE EXACTE, ERREUR aux'
     %        ,I8,' NOEUDS parmi ',I9,'NOEUDS')
C            /140('='))
20000 FORMAT('At TIME',G14.6,' COMPUTED WAVE, EXACT WAVE, ERROR at',
     %         I8,' NODES among 'I9,' NODES')
C            /140('='))
C
         MNXYZ = MNXYZN + WYZNOE - NBCOOR
         DO 10 I=1,NBNOMA
C
C           LES 4 PARAMETRES D'APPEL DE LA FONCTION 'PARTIE_REELLE_EXACTE...'
C           LE TEMPS EN 1-ER PARAMETRE PUIS X Y Z
            DPARAF(1,1)= TEMPS
            MN         = MNXYZ + NBCOOR * I
            N          = 1
            DO K=0,NBCOOR-1
               N = N + 1
               DPARAF(N,1) = RMCN(MN+K)
            ENDDO
C
C           WAVEXA(1) = PARTIE_REELLE_EXACTE(t,x,y,z)
C           WAVEXA(2) = PARTIE_IMAGINAIRE_EXACTE(t,x,y,z)
            DO M=1,2
               CALL FONVAL( NOFOWE(M), N, DPARAF, NCODEV, WAVEXA(M) )
            ENDDO
            IF( NCODEV .GT. 0 ) THEN
C
               DO M=1,2
C                 WAVE EXACTE CORRECTEMENT INITIALISEE
                  DPARAF(1,M) = WAVEXA(M)
                  WEXMAX(M) = MAX( WEXMAX(M), WAVEXA(M) )
                  WEXMIN(M) = MIN( WEXMIN(M), WAVEXA(M) )
                  SOLL1(M)  = SOLL1(M) + ABS( WAVEXA(M) )
C
C                 WAVE CALCULEE
                  DPARAF(2,M) = WAVE(I,M,NUMCAS)
                  WCAMAX(M) = MAX( WCAMAX(M), DPARAF(2,M) )
                  WCAMIN(M) = MIN( WCAMIN(M), DPARAF(2,M) )
C
C                 ERREUR AU NOEUD: WAVE EXACTE - WAVE CALCULEE
                  DPARAF(3,M) = WAVEXA(M) - WAVE(I,M,NUMCAS)
                  ERRMAX(M) = MAX( ERRMAX(M), ABS(DPARAF(3,M)) )
                  ERRL1(M)  = ERRL1(M) + ABS(DPARAF(3,M))
               ENDDO
C
C              ERREUR SUR LE MODULE: MODULE EXACT - MODULE CALCULE
               DPARAF(5,1) = WAVEXA(1)**2 + WAVEXA(2)**2
               DPARAF(5,1) = SQRT(DPARAF(5,1))
               DPARAF(6,1) = WAVE(I,1,NUMCAS)**2 + WAVE(I,2,NUMCAS)**2
               DPARAF(6,1) = SQRT(DPARAF(6,1))
               DPARAF(7,1) = ABS( DPARAF(5,1) - DPARAF(6,1) )
               IF( DPARAF(5,1) .LT. 1D-10 ) THEN
                  DPARAF(7,2) = 0D0
               ELSE
                  DPARAF(7,2) = DPARAF(7,1) / DPARAF(5,1) * 100D0
               ENDIF
C
C              SOMME DU MODULE AUX NOEUDS et SOMME DE LA DIFFERENCE DES MODULES
               SOMODU  = SOMODU  + DPARAF(5,1)
               SOMODIF = SOMODIF + DPARAF(7,1)
C
ccc               IF( I .LE. NBNOF ) THEN
               IF( I .GE. NBNOM2 .AND. I .LE. NBNOM2+NBNOF-1 ) THEN
C
C                 CALCUL DU POURCENTAGE SI C'EST POSSIBLE
                  DO M=1,2
                     IF( ABS(WAVEXA(M)) .LT. 1D-10 ) THEN
                        DPARAF(4,M) = 0D0
                     ELSE
                        DPARAF(4,M) = ABS( DPARAF(3,M)/WAVEXA(M) )*100D0
                        ERP(M) = MAX( ERP(M), ABS(DPARAF(4,M)) )
                     ENDIF
                  ENDDO
                  IF( LANGAG .EQ. 0 ) THEN
C                    AFFICHAGE EN FRANCAIS
                     WRITE(IMPRIM,10011)I,(RMCN(MN+K),K=0,2),
     %               (DPARAF(1,M),M=1,2),DPARAF(5,1),
     %               (DPARAF(2,M),M=1,2),DPARAF(6,1),
     %                DPARAF(3,1),DPARAF(3,2),DPARAF(7,1),
     %                DPARAF(4,1),DPARAF(4,2),DPARAF(7,2)
                  ELSE
C                    AFFICHAGE EN ANGLAIS
                     WRITE(IMPRIM,20011)I,(RMCN(MN+K),K=0,2),
     %               (DPARAF(1,M),M=1,2),DPARAF(5,1),
     %               (DPARAF(2,M),M=1,2),DPARAF(6,1),
     %                DPARAF(3,1),DPARAF(3,2),DPARAF(7,1),
     %                DPARAF(4,1),DPARAF(4,2),DPARAF(7,2)
                  ENDIF
               ENDIF
            ENDIF
 10      CONTINUE
C
10011 FORMAT('DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %   T60,' PR  Exact=',G14.6,
     %      '  PI  Exact=',G14.6,
     %      '  |U| Exact=',G14.6/
     %   T60,' PR Calcul=',G14.6,
     %      '  PI Calcul=',G14.6,
     %      '  |U|Calcul=',G14.6/
     %   T60,' PR  Diff =',G14.6,
     %      '  PI  Diff =',G14.6,
     %      '  |U| Diff =',G14.6/
     %   T60,' PR  Diff =',G11.3,' % ',
     %      '  Pi  Diff =',G11.3,' % ',
     %      '  |U| Diff =',G11.3,' % ')
C
20011 FORMAT('DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %   T60,' Exact RP =',G14.6,
     %      '  Exact IP =',G14.6,
     %      '  Exact |U|=',G14.6/
     %   T60,' Compu RP =',G14.6,
     %      '  Compu IP =',G14.6,
     %      '  Compu |U|=',G14.6/
     %   T60,' Diff  RP =',G14.6,
     %      '  Diff  IP =',G14.6,
     %      '  Diff  |U|=',G14.6/
     %   T60,' Diff  RP =',G11.3,' % ',
     %      '  Diff  IP =',G11.3,' % ',
     %      '  Diff  |U|=',G11.3,' % ')
C
C        BILAN FINAL
C        VALEUR MOYENNE DU MODULE ET DE L'ERREUR SUR LE MODULE
         SOMODU  = SOMODU  / NBNOMA
         SOMODIF = SOMODIF / NBNOMA
         IF( LANGAG .EQ. 0 ) THEN
C           PARTIE REELLE EN FRANCAIS
            IF( WEXMAX(1)-WEXMIN(1) .NE. 0 ) THEN
               WRITE(IMPRIM,10012) TEMPS, NUMCAS, WEXMIN(1), WCAMIN(1),
     %         WEXMAX(1), WCAMAX(1),
     %         WEXMAX(1)-WEXMIN(1), WCAMAX(1)-WCAMIN(1),
     %       ERRMAX(1),100*ERRMAX(1)/(WEXMAX(1)-WEXMIN(1)),
     %         ERP(1),
     %         ERRL1(1)/SOLL1(1)*100
            ENDIF
C           PARTIE IMAGINAIRE EN FRANCAIS
            IF( WEXMAX(2)-WEXMIN(2) .NE. 0 ) THEN
               WRITE(IMPRIM,10013) TEMPS, NUMCAS, WEXMIN(2), WCAMIN(2),
     %         WEXMAX(2), WCAMAX(2),
     %         WEXMAX(2)-WEXMIN(2), WCAMAX(2)-WCAMIN(2),
     %       ERRMAX(2),100*ERRMAX(2)/(WEXMAX(2)-WEXMIN(2)),
     %         ERP(2),
     %         ERRL1(2)/SOLL1(2)*100
            ENDIF
C           MODULE EN FRANCAIS
            WRITE(IMPRIM,10014) SOMODU, SOMODIF, SOMODIF/SOMODU*100
         ELSE
C           PARTIE REELLE EN ANGLAIS
            IF( WEXMAX(1)-WEXMIN(1) .NE. 0 ) THEN
               WRITE(IMPRIM,20012) TEMPS, NUMCAS, WEXMIN(1), WCAMIN(1),
     %         WEXMAX(1), WCAMAX(1),
     %         WEXMAX(1)-WEXMIN(1), WCAMAX(1)-WCAMIN(1),
     %       ERRMAX(1),100*ERRMAX(1)/(WEXMAX(1)-WEXMIN(1)),
     %         ERP(1),
     %         ERRL1(1)/SOLL1(1)*100
            ENDIF
C           PARTIE IMAGINAIRE EN ANGLAIS
            IF( WEXMAX(2)-WEXMIN(2) .NE. 0 ) THEN
               WRITE(IMPRIM,20013) TEMPS, NUMCAS, WEXMIN(2), WCAMIN(2),
     %         WEXMAX(2), WCAMAX(2),
     %         WEXMAX(2)-WEXMIN(2), WCAMAX(2)-WCAMIN(2),
     %       ERRMAX(2),100*ERRMAX(2)/(WEXMAX(2)-WEXMIN(2)),
     %         ERP(2),
     %         ERRL1(2)/SOLL1(2)*100
            ENDIF
C           MODULE EN ANGLAIS
            WRITE(IMPRIM,20014) SOMODU, SOMODIF, SOMODIF/SOMODU*100
         ENDIF
C
10012 FORMAT('Au temps',G14.6,' Vecteur OndeNLSE Numero',I5,
     %' ONDE Partie REELLE:'/
     %'Minimum Exact',G14.6,' Calc=',G14.6,
     %' Maximum Exact',G14.6,' Calc=',G14.6,
     %' Max-Min Exact',G14.6,' Calc=',G14.6/
     %'MAX |PR Exact(N)-PR Calc(N)|=',G14.6,
     %T63,'MAX |PR Exact-PR Calc|(N)/(MAX PR(N)-Min PR(N))=',G10.2,' %'/
     %'MAX(|PR Exact(N)-PR Calc(N)|/|PR Exact(N)|)=',G10.2,
     %T63,'Som |(PR Exact-PR Calc)(N)| / Som |PR Exact(N)|=',
     %G10.2,' %'/
     %'MAX indique un MAXIMUM sur les valeurs aux NOEUDS N et ',
     %' Som indique la SOMME des VALEURS en TOUS les NOEUDS N')
C
10013 FORMAT('Au temps',G14.6,' Vecteur OndeNLSE Numero',I5,
     %' Partie IMAGINAIRE:'/
     %'Minimum ONDE PI Exact',G14.6,' Calc =',G14.6,
     %' Maximum ONDE PI Exact',G14.6,' Calc =',G14.6,
     %' MAX-Min ONDE PI Exact',G14.6,' Calc =',G14.6/
     %'MAX |PI Exact(N)-PI Calc(N)|=',G14.6,
     %T63,'MAX |PI Exact-PI Calc|(N)/(MAX PI(N)-MIN PI(N))=',G10.2,' %'/
     %'MAX(|PI Exact(N)-PI Calc(N)|/|PI Exact(N)|)=', G10.2,
     %T63,'Som|(PI Exact-PI Calc)(N)| / Som |PI Exact(N)| =',
     %G10.2,' %')
C
10014 FORMAT(
     %'Som |Uex|(Noeud)/NbNoeuds=',G14.6,
     %' Som |Uex|-|Uca|(Noeud)/NbNoeuds=',G14.6,
     %' Som |Uex|-|Uca|(Noeud)/Som |Uex|(Noeud)=', G10.2,' %')
C
20012 FORMAT('At time',G14.6,' WAVE Vector Number',I5,
     %' WAVE REAL PART:'/
     %'Minimum RP Exact',G14.6,' Computed =',G14.6,
     %' Maximum RP Exact',G14.6,' Computed =',G14.6,
     %' MAX-Min RP Exact',G14.6,' Computed =',G14.6/
     %'MAX |RP Exact(N) - RP Computed(N)| =',G14.6,
     %T63,'MAX |RP Exact-RP Computed|(N) / (MAX RP(N)-MIN RP(N))=',
     %G10.2,' %'/
     %'MAX(|RP Exact(N) - RP Computed(N)| / |RP Exact(N)|)=',G10.2,
     %T63,'Sum |(RP Exact - RP Computed)(N)| / Sum |RP Exact(N)|=',
     %G10.2,' %'/
     %'Here MAX means a MAXIMUM on the values at NODES N and ',
     %'Sum means the SUM of VALUES at ALL NODES N')
C
20013 FORMAT('At time',G14.6,' WAVE Vector Number',I5,
     %' WAVE IMAGINARY PART:'/
     %'Minimum IP Exact',G14.6,' Computed =',G14.6,
     %' Maximum IP Exact',G14.6,' Computed =',G14.6,
     %' MAX-Min IP Exact',G14.6,' Computed =',G14.6/
     %'MAX|IP Exact(N) - IP Computed(N)| =',G14.6,
     %T63,'MAX|IP Exact-IP Computed|(N) / (MAX IP(N)-MIN IP(N))=',
     % G10.2,' %'/
     %'MAX(|IP Exact(N) - IP Computed(N)| /|IP Exact(N)|)=',
     % G10.2,
     %T63,'Sum|(IP Exact - IP Computed)(N)| / Sum |IP Exact(N)|=',
     % G10.2,' %')
C
20014 FORMAT(
     %'Sum |Uex|(Node)/NodeNb=',G14.6,
     %' Sum |Uex|-|Uca|(Node)/NodeNb=',G14.6,
     %' Sum |Uex|-|Uca|(Node)/Sum |Uex|(Node)=',G10.2,' %')
C
      ELSE
C
C        =======================================================
C        AFFICHAGE DE LA WAVE AUX NOEUDS SANS LES ERREURS
C        =======================================================
         MN = MNXYZN + WYZNOE - NBCOOR
         IF( NOFORE .GT. 0 ) THEN
C
C           LIMITATION PAR LA FONCTION 'REGION(t,x,y,z)'
C           --------------------------------------------
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10020) TEMPS, NUMCAS
            ELSE
               WRITE(IMPRIM,20020) TEMPS, NUMCAS
            ENDIF
C
            DO 30 I=1,NBNOMA
C
C              LES 4 PARAMETRES D'APPEL DE LA FONCTION 'CHOIX'
C              LE TEMPS EN 1-ER PARAMETRE
               DPARAF(1,1) = TEMPS
C              PUIS LES 3 COORDONNEES X Y Z DU NOEUD
               MN        = MN + NBCOOR
               N         = 1
               DO 28 K=0,NBCOOR-1
                  N = N + 1
                  DPARAF(N,1) = RMCN(MN+K)
 28            CONTINUE
C              FONCTION REGION(TEMPS,X,Y,Z)
               CALL FONVAL( NOFORE, N, DPARAF, NCODEV, ERP(1) )
               IF( NCODEV .NE. 0 .AND. NINT(ERP(1)) .NE. 0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10050) I,(RMCN(MN+K),K=0,2),
     %                                 ((WAVE(I,M,K),M=1,2),K=1,NDSM)
                  ELSE
                     WRITE(IMPRIM,20050) I,(RMCN(MN+K),K=0,2),
     %                                 ((WAVE(I,M,K),M=1,2),K=1,NDSM)
                  ENDIF
               ENDIF
 30         CONTINUE
C
         ELSE
C
C           LIMITATION A NBNOF NOEUDS
C           -------------------------
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10049) TEMPS, NUMCAS, NBNOF
            ELSE
               WRITE(IMPRIM,20049) TEMPS, NUMCAS, NBNOF
            ENDIF
            DO 50 I = NBNOM2, NBNOM2+NBNOF-1
               MN = MNXYZN + WYZNOE + NBCOOR * (I-1)
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10050) I,(RMCN(MN+K),K=0,2),
     %                              ((WAVE(I,M,K),M=1,2),K=1,NDSM)
               ELSE
                  WRITE(IMPRIM,20050) I,(RMCN(MN+K),K=0,2),
     %                              ((WAVE(I,M,K),M=1,2),K=1,NDSM)
               ENDIF
 50         CONTINUE
C
         ENDIF
C
      ENDIF
C
10020 FORMAT('Au TEMPS ',G14.6,' CAS',I5,
     %' VALEUR de l''ONDE dans la REGION DEFINIE par la FONCTION ''REGIO
     %N''')
C     %/80(1H=))
10049 FORMAT('Au TEMPS ',G14.6,' CAS',I6,' Valeur de l''ONDE en',I8,
     %        ' NOEUDS')
C   /97('='))
10050 FORMAT('DL',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %       (T58,': ONDE=',G14.6,' + i ',G14.6) )
C
20020 FORMAT('At TIME ',G14.6,' CASE',I5,
     %' WAVE VALUES in the REGION DEFINED by the ''REGION'' FUNCTION')
C    %/80(1H=))
20049 FORMAT('At TIME ',G14.6,' CASE',I6,' WAVE Value at',I8,' NODES')
C     %/97('='))
20050 FORMAT('DoF',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %       (T58,': WAVE=',G14.6,' + i ',G14.6) )
C
C     RECHERCHE DU POINT DE WAVE MINIMALE ET MAXIMALE
C     ===============================================
      DO M=1,2
         WCAMAX(M) = WAVE(1,M,NUMCAS)
         WCAMIN(M) = WCAMAX(M)
         IMIN(M)   = 1
         IMAX(M)   = 1
         DO I=2,NBNOMA
            IF( WAVE(I,M,NUMCAS) .GT. WCAMAX(M) ) THEN
               WCAMAX(M) = WAVE(I,M,NUMCAS)
               IMAX(M) = I
            ENDIF
            IF( WAVE(I,M,NUMCAS) .LT. WCAMIN(M) ) THEN
               WCAMIN(M) = WAVE(I,M,NUMCAS)
               IMIN(M) = I
            ENDIF
         ENDDO
      ENDDO
C
C     MAXIMUM DU MODULE
      MODUMX = 0D0
      DO I=1,NBNOMA
         MODUMX = MAX( MODUMX, WAVE(I,1,NUMCAS)**2+WAVE(I,2,NUMCAS)**2 )
      ENDDO
      MODUMX = SQRT( MODUMX )
C
      MN = MNXYZN + WYZNOE - NBCOOR
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10130) TEMPS, NUMCAS, MODUMX
         WRITE(IMPRIM,10131) 'PReel ONDE ',WCAMIN(1),IMIN(1),
     %                       (RMCN(MN+NBCOOR*IMIN(1)+K),K=0,2)
         WRITE(IMPRIM,10132) 'PReel ONDE ',WCAMAX(1),IMAX(1),
     %                       (RMCN(MN+NBCOOR*IMAX(1)+K),K=0,2)
         WRITE(IMPRIM,10131) 'PImag ONDE ',WCAMIN(2),IMIN(2),
     %                       (RMCN(MN+NBCOOR*IMIN(2)+K),K=0,2)
         WRITE(IMPRIM,10132) 'PImag ONDE ',WCAMAX(2),IMAX(2),
     %                       (RMCN(MN+NBCOOR*IMAX(2)+K),K=0,2)
      ELSE
         WRITE(IMPRIM,20130) TEMPS, NUMCAS, MODUMX
         WRITE(IMPRIM,20131) 'Wave Real Part ',WCAMIN(1),IMIN(1),
     %                       (RMCN(MN+NBCOOR*IMIN(1)+K),K=0,2)
         WRITE(IMPRIM,20132) 'Wave Real Part ',WCAMAX(1),IMAX(1),
     %                       (RMCN(MN+NBCOOR*IMAX(1)+K),K=0,2)
         WRITE(IMPRIM,20131) 'Wave Imag Part ',WCAMIN(2),IMIN(2),
     %                       (RMCN(MN+NBCOOR*IMIN(2)+K),K=0,2)
         WRITE(IMPRIM,20132) 'Wave Imag Part ',WCAMAX(2),IMAX(2),
     %                       (RMCN(MN+NBCOOR*IMAX(2)+K),K=0,2)
      ENDIF
C
10130 FORMAT('Au TEMPS',G14.6,' VECTEUR"ONDE NUMERO =',I5,
     %        ' MAX |ONDE(Noeud)| =',G14.6)
10131 FORMAT(A,'CALCULEE MINIMALE =',G14.6,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
10132 FORMAT(A,'CALCULEE MAXIMALE =',G14.6,
     %' au NOEUD',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
C
20130 FORMAT('At TIME',G14.6,' VECTOR"WAVE NUMBER =',I5,
     %        ' MAX |WAVE(Node)| =',G14.6)
20131 FORMAT(A,'MINIMUM COMPUTED =',G14.6,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
20132 FORMAT(A,'MAXIMUM COMPUTED =',G14.6,
     %' at NODE',I8,' : X=',G13.5,' Y=',G13.5,' Z=',G13.5)
C
      RETURN
      END
