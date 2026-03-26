      SUBROUTINE TRMXVIPR( KNOMOB, IERR, DCPU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER A CHAQUE PAS de TEMPS du CALCUL DES VITESSE-PRESSION de
C -----  MOYVIT = Somme  |VITESSE|(NoeudVitesse) / Nombre NoeudsVitesse
C                NoeudsVitesse
C        MAXVIT = Max  |VITESSE|(NoeudVitesse)
C                NoeudsVitesse

C        MXMIPR = Max-Min Pression(NoeudPression)
C                NoeudsPression

C        FLUX- de la VITESSE PRESSION(Temps)
C        FLUX+ de la VITESSE PRESSION(Temps)

C        No du PAS DT de TEMPS de l'integration du Pb INSTATIONNAIRE)

C        TEMPERATURE Moyenne (Temps)
C        TEMPERATURE Minimale(Temps)
C        TEMPERATURE Maximale(Temps)

C        POUR DES ELEMENTS FINIS DE TAYLOR-HOOD ou BREZZI-FORTIN

C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A CALCULER

C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C DCPU   : SECONDES DE CPU DE L'EXECUTION DE CE SOUS PROGRAMME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     Mars 2012
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY           Fevrier 2022
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY          Decembre 2022
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/a___vecteur.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  DINFO, DCPU
      CHARACTER*(*)     KNOMOB

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/100('=')/
     %'Trace des MAXIMA VITESSE-PRESSION de l''OBJET ',A/100('='))
20000 FORMAT(/100('=')/
     %'Drawings of VELOCITY-PRESSURE MAXIMA of the OBJECT ',A/100('='))
C
C     INITIALISATION DU TEMPS CALCUL INITIAL
      DCPU = DINFO( 'CPU' )
      IERR = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)

C     AFFICHAGE ET VERIFICATION DU NOM_DE_L'OBJET
C     ===========================================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF

C     RECHERCHE DU TMS  VECTEUR"MAXVITESSEPRESS  DE L'OBJET
C     =====================================================
      CALL LXTSOU(NTLXOB, 'VECTEUR"MAXVITESSEPRESS', NTVPMIMX, MNVPMIMX)
      IF( NTVPMIMX .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'MAXIMA VITESSE PRESSION NON CALCULEES'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'VELOCITY PRESSURE MAXIMA NOT COMPUTED'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF

C     N1VPMIMX = NOMBRE DE COMPOSANTES DE CHAQUE VECTEUR
C     1:Temps, 2:|VITESSE|Moyen, 3:|VITESSE|Max,
C     4:PressionMoyen, 5:Max-MinPression, 6:FluxNegatif, 7:FluxPositif
C     8:No du PAS DT de TEMPS de l'integration du Pb INSTATIONNAIRE)
      N1VPMIMX = MCN(MNVPMIMX+WBCOVE)
C     si N1VPMIMX=11 alors en plus
C     9:TemperatureMoyen 10: Temp Min 11: Temp Max
C     NOMBRE DE VECTEURS MAX VITESSE+PRESSION CALCULES (0:NBPAST)
      N2VPMIMX = MCN(MNVPMIMX+WBVECT)
C     NOMBRE DE PAS DE TEMPS CALCULES = NOMBRE DE VECTEURS VP -1
      NBPAST = N2VPMIMX-1

C     ADRESSE MCN DES VECTEURS
      MNMAXVIT = MNVPMIMX + WECTEU
      MNRESU = ( MNMAXVIT + 1 ) / MOREE2
      TEMPS  = REAL( DMCN( MNRESU ) )

C     AFFICHAGE DU TABLEAU DES MIN MAX DES VECTEURS VITESSE+PRESSION
C     --------------------------------------------------------------
      CALL AFVPMIMX( N1VPMIMX, NBPAST, MCN(MNVPMIMX+WECTEU) )

C     TEMPS DES CALCULS D => R
      MNTIME = 0
      CALL TNMCDC( 'REEL', N2VPMIMX, MNTIME )
      IF( MNTIME .LE. 0 ) GOTO 9999
C     TRANSFORMATION DU TEMPS REEL DOUBLE PRECISION EN SIMPLE PRECISION
      MN = MNRESU
      DO NOCOMP=1,N2VPMIMX
         RMCN(MNTIME-1+NOCOMP) = REAL( DMCN(MN) )
         MN = MN + N1VPMIMX
      ENDDO

C     TRACER la COURBE |VITESSE| Moyen = f(Temps)
C     =================================================================
      NOCOMP = 2
      CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'VelocityMean', NCBLEU,
     %              'Time', '|Velocity| Mean',
     %              '|Velocity| Mean(Time)',
     %              ' ',
     %              ' ' )

C     TRACER la COURBE |VITESSE| Maximum = f(Temps)
C     =================================================================
      NOCOMP = 3
      CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'VelocityMax', NCNOIR,
     %              'Time', '|Velocity| Max',
     %              '|Velocity| Maximum(Time)',
     %              ' ',
     %              ' ' )
C
C     TRACER la COURBE PRESSION Maximum-Minimum = f(Temps)
C     =================================================================
      NOCOMP = 5
      CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'PressureMax-Min', NCORAN,
     %              'Time', 'Pressure Max-Min',
     %              'Pressure Maximum-Minimum(Time)',
     %              ' ',
     %              ' ' )

C     TRACER la COURBE FLUX NEGATIF DE LA VITESSE PRESSION = f(Temps)
C     =================================================================
      NOCOMP = 6
      CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'NegativeVelFlux', NCMAGE,
     %              'Time', 'Negative Velocity Flux',
     %              'Negative Velocity Flux(Time)',
     %              ' ',
     %              ' ' )

C     TRACER la COURBE FLUX POSITIF DE LA VITESSE PRESSION = f(Temps)
C     =================================================================
      NOCOMP = 7
      CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'PositiveVelFlux', NCMAGE,
     %              'Time', 'Positive Velocity Flux',
     %              'Positive Velocity Flux(Time)',
     %              ' ',
     %              ' ' )


      IF( N1VPMIMX .EQ. 11 ) THEN

C        TRACER la COURBE TEMPERATURE Moyenne = f(Temps)
C        ==============================================================
         NOCOMP = 9
         CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'TemperatureMean', NCROUG,
     %              'Time', 'Temperature Mean',
     %              'Temperature Mean(Time)',
     %              ' ',
     %              ' ' )

C        TRACER la COURBE TEMPERATURE Min = f(Temps)
C        ==============================================================
         NOCOMP = 10
         CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'TemperatureMin', NCROUG,
     %              'Time', 'Temperature Min',
     %              'Temperature Min(Time)',
     %              ' ',
     %              ' ' )

C        TRACER la COURBE TEMPERATURE Max = f(Temps)
C        ==============================================================
         NOCOMP = 11
         CALL TRTABLE( N2VPMIMX, RMCN(MNTIME),
     %              NOCOMP,  N1VPMIMX, MCN(MNMAXVIT),
     %              'TemperatureMax', NCROUG,
     %              'Time', 'Temperature Max',
     %              'Temperature Max(Time)',
     %              ' ',
     %              ' ' )

      ENDIF

C     AFFICHAGE DU TEMPS CALCUL QUI PEUT ETRE IMPORTANT PAR L'EXECUTION
C     DES FONCTIONS EXACTES DE L'UTILISATEUR VITESSE_EXACTE, PRESSIONEXACTE
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'TEMPS Trace |VITESSE|Moy, |VITESSE|Max, PRESSI
     %ON Max-Min, FLUX NEGATIF Vitesse, FLUX POSITIF Vitesse, TEMPERATUR
     %ES=',DCPU, ' secondes'
      ELSE
         WRITE(IMPRIM,*)'CPU TIME to DRAW |Velocity|Mean, |Velocity|Max
     %PressureMax-Min, NEGATIVE Velocity FLUX, POSITIVE Velocity FLUX, T
     %EMPERATURES=',DCPU,' seconds'
      ENDIF

C     DESTRUCTION DES TEMPS DE CALCUL
      CALL TNMCDS( 'REEL', N2VPMIMX, MNTIME )

 9999 RETURN
      END
