      SUBROUTINE LINMFIVPT( KNOMOB,  MXVPFILE,
     %                      MNTIMES, NMVPFILE, NBVPFILE,
     %                      NCAST0,  NCAST1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RECUPERER DANS LE REPERTOIRE DU PROJET LE NOM DE FICHIER
C -----     DES NBVPFILE VECTEURS VITESSE+PRESSION+TEMPERATURE STOCKES
C           CONSTRUIRE LE TABLEAU TIMES DES NBVPFILE TEMPS
C           DE SAUVEGARDE DE CHACUN D'EUX

C ENTREES:
C --------
C MXVPFILE: NOMBRE MAXIMAL DE FICHIERS VITESSE+PRESSION ECRITS
C MNTIMES : ADRESSE MCN DU TABLEAU DES TEMPS

C SORTIES:
C --------
C MNTIMES : ADRESSE MCN DU TABLEAU DES NBVPFILE TEMPS =>
C           RMCN(MNTIMES:MNTIMES-1+NBVPFILE) = TIMES(NCAST0:NCAST1)
C NMVPFILE: NOM DES NBVPFILE FICHIERS VECTEURS VITESSE+PRESSION
C NBVPFILE: NOMBRE DE FICHIERS VECTEURS VITESSE+PRESSION RETROUVES
C           NOMBRE DE TEMPS DU TABLEAU MCN D'ADRESSE MNTIMES(1:NBVPFILE)
C NCAST0  : NUMERO DU CAS DU PREMIER FICHIER VITESSE+PRESSION+TEMPERATURE LU
C               OU DU TEMPS DU PREMIER FICHIER LU
C NCAST1  : NUMERO DU CAS DU DERNIER FICHIER VITESSE+PRESSION+TEMPERATURE LU
C               OU DU TEMPS DU DERNIER FICHIER LU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Veulettes & St Pierre du Perray     Aout 2020
C MODIFS : ALAIN PERRONNET Veulettes & St Pierre du Perray  Janvier 2022
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ctemps.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      (RMCN(1),MCN(1))

      CHARACTER*(*)     KNOMOB
      CHARACTER*128     KNMVPFI
      CHARACTER*128     NMVPFILE(MXVPFILE)
      CHARACTER*16      KTEMPS
      LOGICAL           LEXIST, LOPEN

C     RECUPERATION EVENTUELLE des TEMPS a partir du nom des FICHIERS des
C     VECTEURS VITESSE+PRESSION SAUVEGARDES SUR DISQUE dans le
C     REPERTOIRE du PROJET
C     ------------------------------------------------------------------
C     DECLARATION DU TABLEAU des TEMPS des VECTEURS VITESSE+PRESSION
      IF( MNTIMES .GT. 0 ) CALL TNMCDS( 'REEL', MXVPFILE, MNTIMES )
      CALL TNMCDC( 'REEL', MXVPFILE, MNTIMES )
      CALL AZEROR( MXVPFILE, RMCN(MNTIMES) )

      NBVPFILE = 0
      NCAST0 = 0
      NCAST1 = 0

C     RECUPERATION DANS LE SYSTEME DE LA LISTE DES FICHIERS DE NOM
C     DEBUTANT PAR LE NOM DE L'OBJET
      NBK = NUDCNB( KNOMOB )
      CALL SYSTEM( 'ls -l ' // KNOMOB(1:NBK) //'*_T=* >./NmFiVeVP' )
ccc      CALL SYSTEM( 'cat ./NmFiVeVP' )

C     LECTURE DU FICHIER NmFiVeVP => NOM DES FICHIERS DES
C     Vecteurs Vitesse+Pression aux differents Temps
      INQUIRE( FILE='NmFiVeVP', EXIST=LEXIST, OPENED=LOPEN )
      IF( .NOT. LEXIST ) THEN
         GOTO 9999
      ENDIF

C     OUVERTURE DU FICHIER 'NmFiVeVP'
      CALL TRUNIT( NOFILE )
      OPEN( FILE='NmFiVeVP', UNIT=NOFILE, STATUS='OLD', IOSTAT=NERR )
      IF( NERR .GT. 0 ) THEN
         GOTO 9999
      ENDIF

      DO NVP = 1, MXVPFILE

C        LECTURE DU NOM DU FICHIER NVP
         READ( UNIT=NOFILE, FMT='(A128)', END=10, IOSTAT=NERR ) KNMVPFI

C        SUPPRESSION des BLANCS et NOMBRE DE CARACTERES
         CALL SANSBL( KNMVPFI, NC1 )

C        NUMERO DU PREMIER CARACTERE DU NOM DE L'OBJET
         NC0 = INDEX( KNMVPFI, KNOMOB(1:NBK) )

         IF( NC0 .GT. 0 ) THEN
C           MISE A JOUR DU NOMBRE DE VECTEURS VITESSE+PRESSION RETROUVE
C           SUR LES FICHIERS DU REPERTOIRE DU PROJET
            NBVPFILE = NBVPFILE + 1
            NMVPFILE(NBVPFILE) = KNMVPFI(NC0:NC1) // '               '

ccc            IF( LANGAG .EQ. 0 ) THEN
ccc               PRINT*,'linmfivpt: Existence du fichier',NBVPFILE,': ',
ccc     %                 NMVPFILE(NBVPFILE)
ccc            ELSE
ccc               PRINT*,'linmfivpt: Existence of file',NBVPFILE,': ',
ccc     %                 NMVPFILE(NBVPFILE)
ccc            ENDIF

         ENDIF

      ENDDO

C     FERMETURE DU FICHIER 'NmFiVeVP'
 10   CLOSE( NOFILE, STATUS='DELETE' )

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'linmfivpt a RETROUVE',NBVPFILE,
     %          ' FICHIERS VECTEUR VITESSE+PRESSION+TEMPERATURE'
      ELSE
         PRINT*,'linmfivpt:',NBVPFILE,
     %          ' VELOCITY+PRESSURE+TEMPERATURE VECTOR FILES FOUND'
      ENDIF

      IF( NBVPFILE .LE. 0 ) THEN
         GOTO 9999
      ENDIF

C     NBVPFILE VECTEURS VITESSE+PRESSION RETROUVES
C     DANS LE REPERTOIRE DU PROJET
C     SUPPRESSION DES DROITS D'ACCES... PAR TRANSLATION A GAUCHE
      NBPB = 0
      NCASTMX = 0
      DO NVP = 1, NBVPFILE

C        POSITION du TEMPS du VECTEUR dans le NOM du FICHIER
         NCT  = INDEX( NMVPFILE(NVP), '_T=' )

C        RECUPERATION du TEMPS du VECTEUR NVP
         IF( NCT .GT. 0 ) THEN

            KTEMPS = NMVPFILE(NVP)(NCT+3:NCT+15)
C           TRADUCTION CARACTERES KTEMPS -> T REEL
            READ( KTEMPS,'(E13.7)' ) T

ccc            IF( LANGAG .EQ. 0 ) THEN
ccc               PRINT*,'linmfivpt: Existence du fichier ',
ccc     %                 NMVPFILE(NVP)(1:NCT+15),' du TEMPS=', T
ccc            ELSE
ccc               PRINT*,'linmfivpt: Existence of FILE ',
ccc     %                 NMVPFILE(NVP)(1:NCT+15),' du TEMPS=', T
ccc            ENDIF

         ELSE

            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'linmfivpt: NOM du FICHIER',NVP,
     %                ' NON CONFORME SANS _T=',NMVPFILE(NVP)
            ELSE
               PRINT*,'linmfivpt: NAME of FILE',NVP,
     %                ' INCORRECT WITHOUT _T=',NMVPFILE(NVP)
            ENDIF

            T = 0

         ENDIF

C        POSITION du NUMERO de CAS du VECTEUR dans le NOM du FICHIER
         NC = INDEX( NMVPFILE(NVP), '_N=' )

         IF( NC .GT. 0 ) THEN
C           TRADUCTION CARACTERES -> ENTIER
            READ( NMVPFILE(NVP)(NC+4:NC+7), '(I4)' ) NCAST

            IF( NVP .NE. NCAST ) THEN
               NBPB = NBPB + 1
               IF( LANGAG .EQ. 0 ) THEN
                  PRINT 10030, T, NCAST, NVP
10030             FORMAT('linmfivpt: Temps,',E15.7,' No du temps',I6,
     %                   ' et celui du FICHIER',I6, ' sont DIFFERENTS' )
               ELSE
                  PRINT 20030, T, NCAST, NVP
20030             FORMAT('linmfivpt: Time,',E15.7,' Time Number',I6,
     %                   ' and this of FILE',I6, ' are DIFFERENTS' )
               ENDIF
            ENDIF

C           LE NUMERO DU TEMPS
            RMCN( MNTIMES -1 + NCAST ) = T

            IF( NCAST .GT. NCASTMX ) THEN
               NCASTMX = NCAST
            ENDIF
         ENDIF

      ENDDO

C     POSITION du NUMERO de CAS du VECTEUR dans le NOM du 1-ER FICHIER
      NCT = INDEX( NMVPFILE(1), '_N=' )
C     RECUPERATION du NUMERO du CAS DU PREMIER VECTEUR LU
      IF( NCT .GT. 0 ) THEN
         KTEMPS = NMVPFILE(1)(NCT+4:NCT+7)
C        TRADUCTION CARACTERES -> ENTIER
         READ( KTEMPS,'(I4)' ) NCAST0
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'linmfivpt: NOM du FICHIER NON CONFORME SANS _N='
         ELSE
            PRINT*,'linmfivpt: FILE NAME INCORRECT WITHOUT _N='
         ENDIF
         PRINT*, NMVPFILE(1)
         NCAST0 = 1
      ENDIF

C     POSITION du NUMERO de CAS du VECTEUR dans le NOM du DERNIER FICHIER
      NCAST1 = NCASTMX

C     SI PB, AFFICHAGE DES TEMPS DES NBVPFILE VECTEURS VITESSE+PRESSION
      IF( NBPB .GT. 0 ) THEN
         PRINT*
         DO NVP=1,NBVPFILE
            IF( LANGAG .EQ. 0 ) THEN
               PRINT 19000, NVP, RMCN(MNTIMES-1+NVP), NMVPFILE(NVP)
19000          FORMAT('linmfivpt: Cas',I6,' au Temps=',E15.7,
     %              ' Nom du Fichier Vitesse+Pression+Temperature: ',A )
            ELSE
               PRINT 29000,NVP, RMCN(MNTIMES-1+NVP), NMVPFILE(NVP)
29000          FORMAT('linmfivpt: Case',I6,' at Time=',E15.7,
     %                ' Velocity+Pressure+TEMPERATURE File Name: ', A )
            ENDIF
         ENDDO
         PRINT*
      ENDIF

 9999 RETURN
      END
