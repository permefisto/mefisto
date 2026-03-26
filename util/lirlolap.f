      SUBROUTINE LIRLOLAP( NDIM,   MNXYZP, MNNPEF, TCAS0, TCAS1,
     %                     NBPART, MNPART, NCVALS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  LECTURE DU MODULE DES VITESSES INITIALES DES PARTICULES
C -----  DU RAYON DES BOULES PARTICULES
C        DES LONGITUDES LATITUDES MIN et MAX et 
C        DU NOMBRE DE SUBDIVISIONS DE LEURS ANGLES
C        POUR GENERER LES COORDONNEES XYZ et VITESSE INITIALE
C        DES PARTICULES DE PARCOURS A CALCULER ENSUITE
C        CONNAISSANT LA VITESSE DU FLUIDE EN TOUT POINT ET TEMPS
C        ET LE TEMPS DU DEPART DE LA PARTICULE

C        UNE PARTICULE EST REFUSEE SI AUCUN TETRAEDRE NE LA CONTIENT

C ENTREES:
C --------
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET 2 ou 3
C MNXYZP : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE L'OBJET
C TCAS0  : TEMPS DU PREMIER VECTEUR VITESSE+PRESSION a TRAITER
C TCAS1  : TEMPS DU DERNIER VECTEUR VITESSE+PRESSION a TRAITER
C          POUR CONTROLER LE TEMPS DE DEPART DES PARTICULES

C MODIFIES:
C ---------
C NBPART : NOMBRE DE PARTICULES DEJA DECLAREES (ET MNPART>0)
C MNPART : ADRESSE MCN DU TABLEAU DES COORDONNEES XYZ + XYZVIT + RAYON + T
C          de la BOULE des NBPART PARTICULES

C SORTIE :
C --------
C NCVALS : > 0 NBPART PARTICULES SONT DEFINIES DANS MCN(MNPART)
C          <=0 ABANDON DEMANDE LORS DES DONNEES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY         Novembre 2020
C23456---------------------------------------------------------------012
      PARAMETER         (MXPART=10000, MXRAYPAR=64, MXVITMOD=32)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      ( MCN(1), RMCN(1) )
      DOUBLE PRECISION  CB1, CB2, CB3, CB4, XP0, YP0
      REAL              XYZP(3), LONGIMIN, LONGIMAX, LATIMIN, LATIMAX,
     %                  LONGI, LATI, DEGRAD, TCAS0, TCAS1,
     %                  RAYPAR(MXRAYPAR), VITMOD(MXVITMOD)

C     DEGRE -> RADIANS
      DEGRAD = ATAN( 1. ) / 45.

C     NBPART est le NOMBRE DE PARTICULES DU PRECEDENT CALCUL
C     DESTRUCTION DU TABLEAU DES PARTICULES PRECEDENTES
 1    IF( NBPART .GT. 0 .AND. MNPART .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU DES NBPART XYZ PARTICULES
         CALL TNMCDS( 'REEL', 8*NBPART, MNPART )
      ENDIF

      NBPART  = 0
      NBPART1 = 0

C     LECTURE DU NOMBRE DE POINTS DE DEPART DES PARTICULES
C     ----------------------------------------------------
      CALL INVITE( 176 )
      NCVALS = 0
      CALL LIRENT( NCVALS, NBPTDEPAR )
      IF( NCVALS    .LE. 0 ) GOTO 9999
      IF( NBPTDEPAR .LE. 0 ) GOTO 9999

C     RESERVATION DU TABLEAU XYZ VitXYZ DES PARTICULES
C     ------------------------------------------------
      CALL TNMCDC( 'REEL', 8*MXPART, MNPART )

C     =======================================================
C     BOUCLE DE LECTURE DES DONNEES EN CHAQUE POINT DE DEPART
C     =======================================================
      DO NUPTDE=1, NBPTDEPAR

C        LECTURE DES XYZ DU POINT DE DEPART DES PARTICULES
C        -------------------------------------------------
         CALL INVITE( 142 )
         NCVALS = 0
         CALL LIRXYZ( NCVALS, XYZP )
         IF( NCVALS .LE. 0 ) GOTO 9999

         PRINT*
         PRINT*,'lirlolap: Point de DEPART',NUPTDE,';',XYZP

         IF( NDIM .EQ. 3 ) THEN
C           RECHERCHE EXHAUSTIVE DU TETRAEDRE NEF CONTENANT LE POINT XYZP
            CALL RETETRXYZ0( XYZP, MNNPEF, MNXYZP,
     %                       NEF, CB1, CB2, CB3, CB4,IERR )
         ELSE

C           PASSAGE DE REEL EN DOUBLE PRECISION
            XP0 = XYZP(1)
            YP0 = XYZP(2)

C           RECHERCHE EXHAUSTIVE DU TRIANGLE NEF CONTENANT LE POINT XYP0
            CALL RETRIAXY0( XP0, YP0, MNNPEF, MNXYZP,
     %                      NEF, CB1, CB2, CB3, IERR )
            CB4 = 0D0

         ENDIF

C        IERR=0 PAS D'ERREUR, NEF>0 EST LE NUMERO DE L'EF CONTENANT XYZP0
C            =1 PAS DE TETRAEDRE ou TRIANGLE CONTENANT LE POINT XYZP0
C            =2 EXISTENCE D'UN TETRAEDRE D'AIRE <=0
         IF( IERR .EQ. 1 ) THEN
C           PAS DE POINT DE DEPART DES PARTICULES INTERIEUR AU FLUIDE
C           -> SORTIE
            GOTO 9999
         ENDIF

C     LECTURE DU RAYON DE LA BOULE PARTICULE
C     -------------------------------------------------------------------------
C     Pour FAIRE un CHOIX de RAYON de la GOUTTELETTE avec COVID19
C     Revue "Pour la Science" du 31 mars 2020:
C     Jean-Michel COURTY, Edouard KIERLIK,
C     Laboratoire de Physique Sorbonne Universite
C    "Le Covid19 est transmis par l'expiration de gouttelettes d'eau
C     de diametre compris entre 1 micrometre (1E-6 m) et 100 micrometres
C     qui s'evaporent rapidement et liberent dans l'air des bacteries
C     (0,5 a 5 micrometres) et virus (0,02 a 0,3 micrometre,
C      0,1 micrometre (1E-7 m) pour le virus Covid19)"

C     Pour le CHOIX de la PLUS GROSSE des GOUTELETTES de RAYON 50 MICRO-METRES
C     RMCN( MN + 7 ) = RAYPAR = 50 * 1E-6
C     ------------------------------------------------------------------------
 6       CALL INVITE( 167 )
         NCVALS = 0
         CALL LIRENT( NCVALS, NBRAYPAR )
         IF( NCVALS   .LE. 0 ) GOTO 1
         IF( NBRAYPAR .LE. 0 ) GOTO 6
         IF( NBRAYPAR .GT. MXRAYPAR ) GOTO 6

         DO N = 1, NBRAYPAR
 7          CALL INVITE( 168 )
            NCVALS = 0
            CALL LIRRSP( NCVALS, RAYPAR(N) )
            IF( NCVALS .LE. 0 ) GOTO 1
            IF( RAYPAR(N) .LE. 0 ) GOTO 7
         ENDDO

C        LECTURE DES MODULES DE LA VITESSE DE DEPART DES PARTICULES
C        --------------------------------------------------------
 8       CALL INVITE( 169 )
         NCVALS = 0
         CALL LIRENT( NCVALS, NBVITMOD )
         IF( NCVALS   .LE. 0 ) GOTO 1
         IF( NBVITMOD .LE. 0 ) GOTO 8

         DO N = 1, NBVITMOD
 9          CALL INVITE( 160 )
            NCVALS = 0
            CALL LIRRSP( NCVALS, VITMOD(N) )
            IF( NCVALS .LE. 0 ) GOTO 1
            IF( VITMOD(N) .LE. 0 ) GOTO 9
         ENDDO

C        LECTURE DE LA LONGITUDE MIN ET MAX DES VITESSES DES PARTICULES
C        --------------------------------------------------------------
         CALL INVITE( 161 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, LONGIMIN )
         IF( NCVALS .LE. 0 ) GOTO 1
 10      IF( LONGIMIN .LT. -360. ) THEN
            LONGIMIN = LONGIMIN + 360.
            GOTO 10
         ENDIF
         IF( LONGIMIN .GT. 360. ) THEN
            LONGIMIN = LONGIMIN - 360.
            GOTO 10
         ENDIF

         CALL INVITE( 162 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, LONGIMAX )
         IF( NCVALS .LE. 0 ) GOTO 1
 11      IF( LONGIMAX .LT. -360. ) THEN
            LONGIMAX = LONGIMAX + 360.
            GOTO 11
         ENDIF
         IF( LONGIMAX .GT. 360. ) THEN
            LONGIMAX = LONGIMAX - 360.
            GOTO 11
         ENDIF

C        ICI -180 <= LONGIMIN LONGIMAX <= 180
         IF( LONGIMIN .GT. 0 .AND. LONGIMAX .LT. 0 ) THEN
C           DESIGNE L'INTERVALLE [0->2Pi] dans le SENS DIRECT
            LONGIMAX = LONGIMAX + 360
         ENDIF

C        LECTURE DU NOMBRE DE SUDIVISIONS DE LA LONGITUDE
C        ------------------------------------------------
 14      CALL INVITE( 163 )
         NCVALS = 0
         CALL LIRENT( NCVALS, NBLONGI )
         IF( NCVALS  .LE. 0 ) GOTO 1
         IF( NBLONGI .LE. 0 ) GOTO 14

C        LECTURE DE LA LATITUDE MIN ET MAX DES VITESSES DES PARTICULES
C     -------------------------------------------------------------
 12      CALL INVITE( 164 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, LATIMIN )
         IF( NCVALS .LE. 0 ) GOTO 1
         IF( LATIMIN .LT. -90. ) THEN
            GOTO 12
         ENDIF
         IF( LATIMIN .GT. 90. ) THEN
            GOTO 12
         ENDIF

 13      CALL INVITE( 165 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, LATIMAX )
         IF( NCVALS .LE. 0 ) GOTO 1
         IF( LATIMAX .LT. -90. ) THEN
            GOTO 13
         ENDIF
         IF( LATIMAX .GT. 90. ) THEN
            GOTO 13
         ENDIF

         IF( LATIMAX .LT. LATIMIN ) THEN
            R       = LATIMIN
            LATIMIN = LATIMAX
            LATIMAX = R
         ENDIF

C        LECTURE DU NOMBRE DE SUDIVISIONS DE LA LATITUDE
C        ------------------------------------------------
 15      CALL INVITE( 166 )
         NCVALS = 0
         CALL LIRENT( NCVALS, NBLATI )
         IF( NCVALS .LE. 0 ) GOTO 1
         IF( NBLATI .LE. 0 ) GOTO 15

C        NOMBRE DE PARTICULES POUR CE POINT DE DEPART DE PARCOURS A CALCULER
C        -------------------------------------------------------------------
         NBPA1RAY = (1+NBLONGI) * (1+NBLATI)
         NBPART1  = NBVITMOD * NBRAYPAR * NBPA1RAY
         IF( NBPART+NBPART1 .GT. MXPART ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'lirlolap:',NBPART+NBPART1,' >',MXPART,
     %                ' NOMBRE MAXIMAL de PARTICULES a TRAITER'
            ELSE
               PRINT*,'lirlolap:',NBPART+NBPART1,' >',MXPART,
     %                ' MAXIMUM PARTICLE NUMBER TO BE TREAT'
            ENDIF
            GOTO 9999
         ENDIF

C        LECTURE DU TEMPS DU DEPART DES NBPART1 PARTICULES DE CE POINT
C        -------------------------------------------------------------
         CALL INVITE( 172 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, TDEPAR0 )
         IF( NCVALS .LE. 0 ) GOTO 1
         IF( TDEPAR0 .LT. TCAS0 .OR. TDEPAR0 .GE. TCAS1 ) THEN
            TDEPAR0 = TCAS0
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='TEMPS du DEPART de la PREMIERE PARTICULE NON COM
     %PRIS ENTRE                                 '
               WRITE( KERR(1)(60:73),'(G14.7)' ) TCAS0
               WRITE( KERR(1)(75:88),'(G14.7)' ) TCAS1
               KERR(2)='TEMPS du DEPART FIXE a                '
               WRITE( KERR(2)(24:37),'(G14.7)' ) TDEPAR0
            ELSE
               KERR(1)='START TIME of FIRST PARTICULE NOT BETWEEN       
     %                               '
               WRITE( KERR(1)(43:56),'(G14.7)' ) TCAS0
               WRITE( KERR(1)(58:71),'(G14.7)' ) TCAS1
               KERR(2)='START TIME FIXED to                '
               WRITE( KERR(2)(21:34),'(G14.7)' ) TDEPAR0
            ENDIF
            CALL LERESU
         ENDIF

         CALL INVITE( 173 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, TDEPAR1 )
         IF( NCVALS .LE. 0 ) GOTO 1
         IF( TDEPAR1 .LT. TCAS0 .OR. TDEPAR1 .GE. TCAS1 ) THEN
            TDEPAR1 = (TCAS0+TCAS1)/2
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='TEMPS du DEPART de la DERNIERE PARTICULE NON COM
     %PRIS ENTRE                                 '
               WRITE( KERR(1)(60:73),'(G14.7)' ) TCAS0
               WRITE( KERR(1)(75:88),'(G14.7)' ) TCAS1
               KERR(2)='TEMPS du DEPART FIXE a                '
               WRITE( KERR(2)(24:37),'(G14.7)' ) TDEPAR1
            ELSE
               KERR(1)='START TIME of LAST PARTICULE NOT BETWEEN        
     %                            '
               WRITE( KERR(1)(43:56),'(G14.7)' ) TCAS0
               WRITE( KERR(1)(58:71),'(G14.7)' ) TCAS1
               KERR(2)='START TIME FIXED to                  '
               WRITE( KERR(2)(21:34),'(G14.7)' ) TDEPAR1
            ENDIF
            CALL LERESU
 
         ENDIF

C        CALCUL DE LA VITESSE DES PARTICULES
C        -----------------------------------
         ANGLATI  = (LATIMAX - LATIMIN) / NBLATI
         ANGLONGI = (LONGIMAX-LONGIMIN) / NBLONGI

         NOPART   = 0
         MN       = MNPART + 8*NBPART - 1

         DO NoVITMOD = 1, NBVITMOD

C        VITMOD(NoVITMOD) MODULE DE LA VITESSE INITIALE DES PARTICULES
C        (IDENTIQUE POUR TOUS LES RAYONS et LES LONGITUDES et LATITUDES)
         VIMO = VITMOD( NoVITMOD )

         DO NoRAYPAR = 1, NBRAYPAR

C           RAYPAR(NoRAYPAR) EST LE RAYON DES NBPA1RAY BOULES-PARTICULES
C           (IDENTIQUE POUR TOUTES LONGITUDES et LATITUDES)
            RAPA = RAYPAR( NoRAYPAR )

C           VALEURS SELON LA LONGITUDE et LATITUDE INITIALE DE LA PARTICULE
            DO J = 0, NBLATI

C              LATITUDE DU POINT I J DE LA SPHERE DE RAYON VIMO
               LATI = ( LATIMIN + ANGLATI *
     %                          ( J + (NoVITMOD-1)/NBVITMOD) ) * DEGRAD

               COSLATI = COS( LATI )
               SINLATI = SIN( LATI )
               DO I = 0, NBLONGI

C                 XYZ DU DEPART DE LA PARTICULE
                  RMCN( MN + 1 ) = XYZP( 1 )
                  RMCN( MN + 2 ) = XYZP( 2 )
                  RMCN( MN + 3 ) = XYZP( 3 )

C                 LONGITUDE DU POINT I J DE LA SPHERE DE RAYON VIMO
                  LONGI = LONGIMIN + ANGLONGI *
     %                             ( I + (NoRAYPAR-1) / NBRAYPAR )

                  IF( LONGI .GT. 180 ) LONGI = LONGI - 360
                  IF( LONGI .LT.-180 ) LONGI = LONGI + 360
                  LONGI = LONGI * DEGRAD

C                 VITESSE XYZ DU DEPART DE LA PARTICULE
                  RMCN( MN + 4 ) = VIMO * COSLATI * COS( LONGI )
                  RMCN( MN + 5 ) = VIMO * COSLATI * SIN( LONGI )
                  RMCN( MN + 6 ) = VIMO * SINLATI

C                 RAYON DE LA BOULE-PARTICULE I
C                 POUR CETTE LONGITUDE et LATITUDE
                  RMCN( MN + 7 ) = RAPA

C                 LE NUMERO DE LA PARTICULE
                  NOPART = NOPART + 1

C                 TEMPS DE DEPART DE LA PARTICULE NOPART
                  MN = MN + 8
                  TDEP = TDEPAR0 +NOPART*(TDEPAR1-TDEPAR0) / NBPART1
                  RMCN( MN )= TDEP

C                 REMARQUE: SI TDEPAR0=TDEPAR1 alors LE TEMPS DE DEPART
C                                                    EST CONSTANT
               ENDDO

            ENDDO

            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'lirlolap:',NOPART,
     %                ' No de la DERNIERE PARTICULE de RAYON',RAPA,
     %                ' VITESSE INITIALE',VIMO,' TEMPS DEPART=',TDEP
            ELSE
               PRINT*,'lirlolap:',NOPART,
     %                ' LAST PARTICLE of RADIUS',RAPA,
     %                ' INITIAL VELOCITY',VIMO,' START TIME=',TDEP
            ENDIF

         ENDDO

         ENDDO

C        NOMBRE DE PARTICULES INITIALISEES
         NBPART = NBPART + NBPART1

C     ========================================================
C     FIN DES DONNEES DES PARTICULES DU POINT DE DEPART NUPTDE
C     ========================================================
      ENDDO

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'lirlolap:',NBPART,
     %          ' XYZ+V0XYZ+RAYON+TEMPS0 des PARTICULES sont STOCKEES'
      ELSE
         PRINT*,'lirlolap:',NBPART,
     %          ' XYZ+V0XYZ+RADIUS+TIME0 PARTICLES DATA are STORED'
      ENDIF

      IF( NBPART .LT. MXPART ) THEN
C        REDUCTION DU TABLEAU DES XYZ POUR NBPART1 PARTICULES
         CALL TNMCDC( 'REEL', 8*NBPART, MNPART1 )
         CALL TRTATA( RMCN(MNPART), RMCN(MNPART1), 8*NBPART )
         CALL TNMCDS( 'REEL', 8*MXPART, MNPART )
         MNPART = MNPART1
      ENDIF

      RETURN
      END
