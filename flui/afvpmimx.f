      SUBROUTINE AFVPMIMX( N1VPMIMX, NBPAST, VPMIMX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES MIN MAX AU NBPAST PAS DE TEMPS CALCULE
C -----    DU VECTEUR VITESSES+PRESSIONS+TEMPERATURES D'UN FLUIDE

C ENTREES:
C --------
C N1VPMIMX: 8 ou 11 NOMBRE DE VALEURS PAR PAS DE TEMPS
C NBPAST  : NOMBRE de PAS de TEMPS de CALCUL des VECTEURS VITESSE+PRESSION

C si N1VPMIMX=8 alors
C VPMIMX : VPMIMX(1,i)=Temps           au PAS de TEMPS i
C          VPMIMX(2,i)=|VITESSE|Moyenne
C          VPMIMX(3,i)=|VITESSE|Max
C          VPMIMX(4,i)=PressionMoyenne
C          VPMIMX(5,i)=Max-MinPression,
C          VPMIMX(6,i)=FluxVitesseNegatif,
C          VPMIMX(7,i)=FluxVitessePositif
C          VPMIMX(8,i)=Numero du PAS de TEMPS du CALCUL du
C                      VECTEUR VITESSE+PRESSION de MAX MIN A AFFICHER
C si N1VPMIMX=11 alors
C          VPMIMX( 9,i)=TEMPERATURE Moyenne  au PAS de TEMPS i
C          VPMIMX(10,i)=TEMPERATURE Minimale
C          VPMIMX(11,i)=TEMPERATURE Maximale
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Fevrier 2013
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray          Decembre 2022
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      DOUBLE PRECISION  VPMIMX(N1VPMIMX,0:NBPAST)
      REAL              TPSINI, TPSFIN
      DOUBLE PRECISION  VMAXMAX, VMOYMOY, PMXIMXI, PMOYMOY,
     %                  FLUMMOY, FLUPMOY, FLUPMMOY,
     %                  TMOYMOY, TMINMIN, TMAXMAX

      PRINT *

C     TEMPS INITIAL et FINAL
      TPSINI = VPMIMX(1,0)
      TPSFIN = VPMIMX(1,NBPAST)

      PRINT 19991
19991 FORMAT(254('='))

      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,
     %'Resultats a chaque pas de temps dans l''intervalle de temps (',
     %TPSINI, ',', TPSFIN,' )'
      ELSE
         PRINT *,
     %'RESULTS at every TIME STEP during the INTERVAL of TIME ('
     %,TPSINI, ',', TPSFIN,' )'
      ENDIF
      PRINT 19991

      VMOYMOY = 0D0
      VMAXMAX = VPMIMX(3,1)

      PMOYMOY = 0D0
      PMXIMXI = VPMIMX(5,1)

      FLUMMOY = 0D0
      FLUPMOY = 0D0
      FLUPMMOY= 0D0

C     PRESENCE DES TEMPERATURES MOY MIN MAX A CHAQUE PAS DE TEMPS?
      TMOYMOY = 0D0
      IF( N1VPMIMX .EQ. 11 ) THEN
         TMINMIN = VPMIMX(10,1)
         TMAXMAX = VPMIMX(11,1)
      ELSE
         TMINMIN = 0D0
         TMAXMAX = 0D0
      ENDIF

      DO K=0,NBPAST

C        NUMERO DU PAS DE TEMPS CALCULE
         NBPASDT = NINT( VPMIMX(8,K) )

         IF( N1VPMIMX .EQ. 11 ) THEN
C           PRESENCE DES TEMPERATURES MOY MIN MAX A CHAQUE PAS DE TEMPS
            IF( LANGAG .EQ. 0 ) THEN
               PRINT 19990,NBPASDT,(VPMIMX(I,K),I=1,7),
     %                     VPMIMX(6,K)+VPMIMX(7,K),
     %                     VPMIMX(10,K),VPMIMX(9,K),VPMIMX(11,K)
            ELSE
               PRINT 29990,NBPASDT,(VPMIMX(I,K),I=1,7),
     %                     VPMIMX(6,K)+VPMIMX(7,K),
     %                     VPMIMX(10,K),VPMIMX(9,K),VPMIMX(11,K)
            ENDIF
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
               PRINT 19993,NBPASDT,(VPMIMX(I,K),I=1,7),
     %                     VPMIMX(6,K)+VPMIMX(7,K)
            ELSE
               PRINT 29993,NBPASDT,(VPMIMX(I,K),I=1,7),
     %                     VPMIMX(6,K)+VPMIMX(7,K)
            ENDIF
         ENDIF

19990 FORMAT(' Pas Dt',I9,' au Temps=',g13.6,
     %' VitMoy=', g13.6,' VitMax=',g13.6,
     %'  PreMoy=',g13.6,' PreMax-Min=',g13.6,
     %'  FluV-=', g13.6,' FluV+=',g13.6,' FluV-+=',g13.6,
     %'  TemMin=',g13.6,' TemMoy=',g13.6,' TemMax=',g13.6)

29990 FORMAT(' Time STEP',I9,' at Time=',g13.6,
     %' VelMean=', g13.6,' VelMax=',g13.6,
     %'  PreMean=',g13.6,' PreMax-Min=',g13.6,
     %'  FluV-=',  g13.6,' FluV+=',g13.6,' FluV-+=',g13.6,
     %'  TemMin=', g13.6,' TemMoy=',g13.6,' TemMax=',g13.6)

19993 FORMAT(' Pas Dt',I9,' au Temps=',g13.6,
     %' VitMoy=', g13.6,' VitMax=',g13.6,
     %'  PreMoy=',g13.6,' PreMax-Min=',g13.6,
     %'  FluV-=', g13.6,' FluV+=',g13.6,' FluV-+=',g13.6)

29993 FORMAT(' Time STEP',I9,' at Time=',g13.6,
     %' VelMean=', g13.6,' VelMax=',g13.6,
     %'  PreMean=',g13.6,' PreMax-Min=',g13.6,
     %'  FluV-=',  g13.6,' FluV+=',g13.6,' FluV-+=',g13.6)

C        MOYENNE DES VITESSES MOYENNES
         VMOYMOY = VMOYMOY + VPMIMX(2,K)

C        MAXIMUM DE LA VITESSE MAXIMALE DE CE PAS DE TEMPS
         IF( VPMIMX(3,K) .GT. VMAXMAX ) VMAXMAX = VPMIMX(3,K)

C        MOYENNE DES PRESSIONS MOYENNES
         PMOYMOY = PMOYMOY + VPMIMX(4,K)

C        MAXIMUM DE MAX-MIN PRESSION DE CE PAS DE TEMPS
         IF( VPMIMX(5,K) .GT. PMXIMXI ) PMXIMXI = VPMIMX(5,K)

C        MOYENNE DES FLUX NEGATIFS
         FLUMMOY = FLUMMOY + VPMIMX(6,K)
C        MOYENNE DES FLUX POSITIFS
         FLUPMOY = FLUPMOY + VPMIMX(7,K)
C        MOYENNE DES FLUX ( NEGATIFS + POSITIFS )
         FLUPMMOY = FLUPMMOY + VPMIMX(6,K)+VPMIMX(7,K)

         IF( N1VPMIMX .EQ. 11 ) THEN
C           PRESENCE DES TEMPERATURES MOY MIN MAX A CHAQUE PAS DE TEMPS
C           MOYENNE DES TEMPERATURES MOYENNES
            TMOYMOY = TMOYMOY + VPMIMX(9,K)
C           MINIMUM DES TEMPERATURES MINIMALES
            IF( VPMIMX(10,K) .LT. TMINMIN ) TMINMIN = VPMIMX(10,K)
C           MAXIMUM DES TEMPERATURES MAXIMALES
            IF( VPMIMX(11,K) .GT. TMAXMAX ) TMAXMAX = VPMIMX(11,K)
         ENDIF

      ENDDO

C     LES VALEURS MOYENNES SUR TOUS LES PAS DE TEMPS DE TPSINI A TPSFIN
      VMOYMOY = VMOYMOY / (1+NBPAST)
      PMOYMOY = PMOYMOY / (1+NBPAST)
      FLUMMOY = FLUMMOY / (1+NBPAST)
      FLUPMOY = FLUPMOY / (1+NBPAST)
      FLUPMMOY= FLUPMMOY/ (1+NBPAST)

      IF( N1VPMIMX .EQ. 11 ) THEN
C        PRESENCE DES TEMPERATURES MOY MIN MAX A CHAQUE PAS DE TEMPS
C        MOYENNE DES TEMPERATURES MOYENNES
         TMOYMOY = TMOYMOY / (1+NBPAST)
      ENDIF

      PRINT 19991

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 19990, NBPASDT, VPMIMX(1,NBPAST), VMOYMOY, VMAXMAX,
     %         PMOYMOY, PMXIMXI, FLUMMOY, FLUPMOY, FLUPMMOY,
     %         TMINMIN, TMOYMOY, TMAXMAX
         PRINT *,
     %' Minima, Moyennes, Maxima des VITESSES PRESSIONs TEMPERATURES sur
     % l''INTERVALLE de TEMPS [',TPSINI,',', TPSFIN,' ]'
      ELSE
         PRINT 29990, NBPASDT, VPMIMX(1,NBPAST), VMOYMOY, VMAXMAX,
     %         PMOYMOY, PMXIMXI, FLUMMOY, FLUPMOY, FLUPMMOY,
     %         TMINMIN, TMOYMOY, TMAXMAX
         PRINT *,
     %'Minima, Means, Maxima of Velocity Pressure Temperature DURING the
     % INTERVAL of TIME (',TPSINI,',', TPSFIN,' )'
       ENDIF

      PRINT 19991
      PRINT *

      RETURN
      END
