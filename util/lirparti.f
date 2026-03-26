      SUBROUTINE LIRPARTI( NDIM,   MNXYZP, MNNPEF, TCAS0, TCAS1,
     %                     NBPART, MNPART, NCVALS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  LECTURE DES COORDONNEES XYZ, VITESSE INITIALE VXYZ,
C -----  Rayon DES BOULES PARTICULES et Temps du depart DES PARTICULES
C        DE PARCOURS A VISUALISER
C        UNE PARTICULE EST REFUSEE SI AUCUN ELEMENT FINI NE LA CONTIENT

C ENTREE :
C --------
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET 2 ou 3
C MNXYZP : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE

C TCAS0  : TEMPS DU PREMIER VECTEUR VITESSE+PRESSION a TRAITER
C TCAS1  : TEMPS DU DERNIER VECTEUR VITESSE+PRESSION a TRAITER
C          POUR CONTROLER LE TEMPS DE DEPART DES PARTICULES
C NBPART : NOMBRE DE PARTICULES DEJA DECLAREES (ET MNPART>0)

C SORTIES:
C --------
C NBPART : NOMBRE DE PARTICULES (REMPLACENT CELLES EXISTANTES)
C MNPART : ADRESSE MCN DU TABLEAU DES
C          X0 Y0 Z0 COORDONNEES INITIALES + VX0 VY0 VZ0 VITESSE INITIALES
C          + RAYON + TEMPS DE DEPART DES NBPART PARTICULES
C          RMCN(MNPART...)=PARTICU(8,NBPART)
C          REMARQUE:
C          EN 2D, BIEN QUE DEMANDES, Z0, VZ0, RAYON, NE SONT PAS UTILISES
C NCVALS : <=0 SI ABANDON DEMANDE
C          >0  SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY NOVEMBRE 2010
C MODIFS : ALAIN PERRONNET Veulettes sur mer                Juillet 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      ( MCN(1), RMCN(1) )
      DOUBLE PRECISION  CB1, CB2, CB3, CB4, XP0, YP0

 10   IF( NBPART .GT. 0 .AND. MNPART .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU DES NBPART (XYZ VXYZ R T) PARTICULES
         CALL TNMCDS( 'REEL', 8*NBPART, MNPART )
         NBPART = 0
      ENDIF

C     LECTURE DU NOMBRE NBPART DE PARTICULES
C     --------------------------------------
 15   CALL INVITE( 141 )
      NCVALS = 0
      CALL LIRENT( NCVALS, NBPART0 )
      IF( NCVALS  .LE. 0 ) GOTO 9999
      IF( NBPART0 .LE. 0 ) GOTO 15

C     DECLARATION DU TABLEAU DES 8 COORDONNEES DES NBPART0 PARTICULES
      CALL TNMCDC( 'REEL', 8*NBPART0, MNPART )
      CALL AZEROR( 8*NBPART0, RMCN(MNPART) )

      PRINT*
      NBPART1 = 0
      MNP     = MNPART
      DO 50 K = 1, NBPART0

C        LECTURE DES XYZ DU POINT DE DEPART DE LA PARTICULE K
C        ----------------------------------------------------
         CALL INVITE( 142 )
         NCVALS = 0
         CALL LIRXYZ( NCVALS, RMCN(MNP) )
         IF( NCVALS .LE. 0 ) GOTO 10

         IF( NDIM .EQ. 3 ) THEN

C           RECHERCHE EXHAUSTIVE DU TETRAEDRE NEF CONTENANT LE POINT K
            CALL RETETRXYZ0( RMCN(MNP), MNNPEF, MNXYZP,
     %                       NEF, CB1, CB2, CB3, CB4, IERR )
         ELSE

C           PASSAGE DE REEL EN DOUBLE PRECISION
            XP0 = RMCN(MNP )
            YP0 = RMCN(MNP+1)

C           RECHERCHE EXHAUSTIVE DU TRIANGLE NEF CONTENANT LE POINT XYP0
            CALL RETRIAXY0( XP0, YP0, MNNPEF, MNXYZP,
     %                      NEF, CB1, CB2, CB3, IERR )
            CB4 = 0D0

         ENDIF

C        IERR=0 PAS D'ERREUR, NEF>0 EST LE NUMERO DE L'EF CONTENANT XYZP0
C            =1 PAS DE TETRAEDRE ou TRIANGLE CONTENANT LE POINT XYZP0
C            =2 EXISTENCE D'UN TETRAEDRE D'AIRE <=0
         IF( IERR .EQ. 1 ) THEN
C           LA PARTICULE K N'EST PAS STOCKEE
C           PASSAGE A LA DEMANDE DES XYZ DE LA PARTICULE K
            GOTO 50
         ENDIF

C        LECTURE DES XYZ DE LA VITESSE INITIALE DE LA PARTICULE K
C        --------------------------------------------------------
         CALL INVITE( 157 )
         NCVALS = 0
         CALL LIRXYZ( NCVALS, RMCN(MNP+3) )
         IF( NCVALS .LE. 0 ) GOTO 10

C     LECTURE DU RAYON DE LA BOULE PARTICULE K
C     -----------------------------------------------------------------
C     Pour FAIRE un CHOIX de RAYON de la GOUTTELETTE avec COVID19
C     Revue "Pour la Science" du 31 mars 2020:
C     Jean-Michel COURTY, Edouard KIERLIK,
C     Laboratoire de Physique Sorbonne Universite
C    "Le Covid19 est transmis par l'expiration de gouttelettes d'eau
C     de diametre compris entre 1 micrometre (1E-6 m) et 100 micrometres
C     qui s'evaporent rapidement et liberent dans l'air des bacteries
C     (0,5 a 5 micrometres) et virus (0,02 a 0,3 micrometre,
C      0,1 micrometre (1E-7 m) pour le virus Covid19)"

C     SI CHOIX de la PLUS GROSSE des GOUTELETTES de RAYON 50 MICRO-METRES
C     RMCN(MN+7) = RAYPAR = 50 * 1E-6
C     -------------------------------------------------------------------
 20      CALL INVITE( 168 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, RAYPAR )
         IF( NCVALS .LE. 0 ) GOTO 10
         IF( RAYPAR .LT. 0 ) GOTO 20
         RMCN(MNP+6) = RAYPAR

C        LE TEMPS DU DEPART DE LA PARTICULE
C        ----------------------------------
 30      CALL INVITE( 171 )
         NCVALS = 0
         CALL LIRRSP( NCVALS, TDEPAR )
         IF( NCVALS .LE. 0 ) GOTO 10
         IF( TDEPAR .LT. TCAS0 .OR. TDEPAR .GT. TCAS1 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='TEMPS du DEPART de la PREMIERE PARTICULE NON COM
     %PRIS ENTRE                                      ' 
               WRITE( KERR(1)(60:73),'(G14.7)' ) TCAS0
               WRITE( KERR(1)(75:88),'(G14.7)' ) TCAS1
            ELSE
               KERR(1)='START TIME of FIRST PARTICULE NOT BETWEEN       
     %                                '     
               WRITE( KERR(1)(43:56),'(G14.7)' ) TCAS0
               WRITE( KERR(1)(58:71),'(G14.7)' ) TCAS1
            ENDIF
            CALL LERESU
            GOTO 30
         ENDIF
         RMCN( MNP+7 )= TDEPAR

C        UNE PARTICULE DE PLUS
         NBPART1 = NBPART1 + 1
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'lirparti: Particule',NBPART1,
     %         ' en X=',RMCN(MNP  ),' Y=',RMCN(MNP+1),' Z=',RMCN(MNP+2),
     %         ' Vitesse V0X=',RMCN(MNP+3),' V0Y=',RMCN(MNP+4),
     %         ' V0Z=',RMCN(MNP+5),' Rayon Boule=',RMCN(MNP+6),
     %         ' Temps0=',RMCN(MNP+7)
         ELSE
            PRINT*,'lirparti: Particle',NBPART1,
     %         ' at X=',RMCN(MNP  ),' Y=',RMCN(MNP+1),' Z=',RMCN(MNP+2),
     %         ' Velocity V0X=',RMCN(MNP+3),' V0Y=',RMCN(MNP+4),
     %         ' V0Z=',RMCN(MNP+5),' Bool Radius=',RMCN(MNP+6),
     %         ' Time0=',RMCN(MNP+7)
         ENDIF

         MNP = MNP + 8
C        PASSAGE A LA LECTURE DES XYZ VXYZ R T0 DE LA PARTICULE SUIVANTE

 50   ENDDO

C     NOMBRE FINAL DE PARTICULES INTERNES AU MAILLAGE
      NBPART = NBPART1
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'lirparti:',NBPART,' XYZ+VXYZ+R+T0 PARTICULES STOCKEES'
      ELSE
         PRINT*,'lirparti:',NBPART,' XYZ+VXYZ+R+T0 STORED PARTICLES'
      ENDIF

      IF( NBPART .LE. 0 ) THEN
C        AUCUNE PARTICULE DANS LE MAILLAGE
         NBPART = 1
         GOTO 10
      ENDIF

      IF( NBPART .LT. NBPART0 ) THEN
C        REDUCTION DU TABLEAU DES XYZ POUR NBPART1 PARTICULES
         CALL TNMCDC( 'REEL', 8*NBPART, MNPART1 )
         CALL TRTATA( RMCN(MNPART), RMCN(MNPART1), 8*NBPART )
         CALL TNMCDS( 'REEL', 8*NBPART0, MNPART )
         MNPART = MNPART1
      ENDIF

 9999 RETURN
      END
