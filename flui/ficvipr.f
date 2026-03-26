      SUBROUTINE FICVIPR( KNOMOB, NBVPFILE, NDIM,   NUTYEL, MNXYZN,
     %                    NBNOVI, NODDL,    NTDLVP, NCAS,   VXYZPN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECRIRE SUR UN FICHIER LES VITESSES ET LES PRESSIONS CALCULEES
C -----    AUX NOEUDS DU CAS NCAS PARMI LES NCAS CALCULES AU TEMPS

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NBVPFILE:NOMBRE D'APPELS ficvipr au DEBUT c-a-d
C          NOMBRE DE VECTEURS VITESSES-PRESSIONS ECRITS SUR UN FICHIER
C TEMPS  : TEMPS DU VECTEUR VITESSE+PRESSION dans ./incl/ctemps.inc
C NUTYEL : 13 TRIANGLE  BREZZI FORTIN 2D
C          15 TRIANGLE  TAYLOR-HOOD   2D
C          19 TETRAEDRE BREZZI FORTIN 3D
C          20 TETRAEDRE TAYLOR-HOOD   3D
C NDIM   : DIMENSION 2 OU 3 DE L'ESPACE DE L'OBJET
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NODDL  : TABLEAU DU NUMERO DU DERNIER D.L. DE CHAQUE NOEUD
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES PRESSIONS
C NCAS   : NUMERO DE LA CARTE DES VITESSES_PRESSIONS A AFFICHER
C VXYZPN : LES CARTES DES VITESSES-PRESSIONS AUX NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C NBVPFILE:NOMBRE D'APPELS ficvipr a la FIN c-a-d
C          NOMBRE DE FICHIERS VECTEURS VITESSES-PRESSIONS ECRITS
C MISE SUR LE FICHIER KNOMFIC DES VITESSES PRESSIONS AU TEMPS ACTUEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Veulettes sur mer               Juillet 2020
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

      CHARACTER*(*)     KNOMOB
      CHARACTER*128     KNOMFIC
      LOGICAL           LEXIST, LOPEN
      INTEGER           NODDL(0:NBNOVI)
      DOUBLE PRECISION  VXYZPN(NTDLVP,NCAS)
      DOUBLE PRECISION  VITNOR, VITMIN, VITMAX, VITMOY,
     %                  PREMIN, PREMAX, PREMOY, PRECALC,
     %                  XP, YP, ZP
      INTRINSIC         SQRT

C     RECHERCHE D'UNE UNITE OPEN DE FICHIER
      CALL TRUNIT( NFUNIT )

C     LE NOM DU FICHIER du VECTEUR VITESSE-PRESSION
      NBK = NUDCNB( KNOMOB )
      KNOMFIC = KNOMOB(1:NBK) // '_VitPres_                            '
      NBK = NUDCNB( KNOMFIC )
      KNOMFIC = KNOMFIC(1:NBK) // '_T_                    '
      WRITE(KNOMFIC(NBK+4:NBK+18),'(G15.7)') TEMPS
      CALL SANSBL( KNOMFIC, NBK )
      KNOMFIC = KNOMFIC(1:NBK) // '_                    '           
      WRITE(KNOMFIC(NBK+2:NBK+10),'(I9)') NBVPFILE
      CALL SANSBL( KNOMFIC, NBK )

      PRINT*
      PRINT*,'ficvipr: au TEMPS=',TEMPS,
     %       ' Ouverture du fichier ',KNOMFIC(1:NBK)

C     SI LE FICHIER KNOMFIC EXISTE ALORS IL EST DETRUIT PUIS RECONSTRUIT
      INQUIRE( FILE=KNOMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( LEXIST ) THEN
C        LE FICHIER KNOMFIC EXISTE
         IF( .NOT. LOPEN ) THEN
C           OUVERTURE DU FICHIER KNOMFIC
            CALL TRUNIT( NFUNIT )
            OPEN( FILE=KNOMFIC, UNIT=NFUNIT, STATUS='OLD' )
         ENDIF
C        DESTRUCTION DU FICHIER
         CLOSE( NFUNIT, STATUS='DELETE' )
      ENDIF

C     CREATION DU FICHIER KNOMFIC
      CALL TRUNIT( NFUNIT )
      OPEN( UNIT=NFUNIT, ERR=9998, STATUS='NEW',
     %      FILE=KNOMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )

C     ======================================================================
C     LES COMPOSANTES DE LA VITESSE ET DE LA PRESSION AUX NOEUDS DU MAILLAGE
C     ======================================================================
      NOEUD1 = 1
      NOEUD2 = NBNOVI

      IF( LANGAG .EQ. 0 ) THEN

C        ECRITURE EN FRANCAIS
         WRITE(NFUNIT,10019) TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
10019    FORMAT(/'Au Temps ',G15.7,': les VITESSES et PRESSIONS des ',
     %    'NOEUDS',I8,' a ',I8,'/',I9,' NOEUDS (Au total',I9,' DL):')
C
         DO 20 I = NOEUD1, NOEUD2

C           L'ADRESSE MCN DES COORDONNEES DES NOEUDS VITESSE
C           ICI LE TMS XYZNOEUD A ETE ENRICHI DES BARYCENTRES DES EF
C           DANS LES CAS DES EF DE BREZZI-FORTIN
            MN  = MNXYZN + WYZNOE - 3 + 3*I

C           LES DEGRES DE LIBERTE AU NOEUD I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL
            IF( NUTYEL .EQ. 19 .OR.  NUTYEL .EQ. 20 ) THEN
C
C              3D BF ou TH : VITESSE"PRESSION
               IF( ND .GT. NDIM ) THEN
                  WRITE(NFUNIT,10023) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K,NCAS),K=1,ND)
               ELSE
                  WRITE(NFUNIT,10033) I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K,NCAS),K=1,ND)
               ENDIF
c
            ELSE
C
C              2D: NOMBRE DE DL EN LE NOEUD I
               IF( ND .GT. NDIM ) THEN
                  WRITE(NFUNIT,10022) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K,NCAS),K=1,ND)
               ELSE
                  WRITE(NFUNIT,10032) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K,NCAS),K=1,ND)
               ENDIF
            ENDIF
 20      ENDDO
C
10022 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' PRESSION=',G15.7)
C
10032 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7)
C
10023 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,
     %'  PRESSION=',G15.7)
C
10033 FORMAT('Noeud',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VITESSE VX=',G15.7,' VY=',G15.7, ' VZ=',G15.7 )

      ELSE

C        ENGLISH PRINTING
         WRITE(NFUNIT,20019) TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
20019 FORMAT(/'At Time ',G15.7,': VELOCITIES and PRESSURES of nodes'
     %,I8,' to ',I8,'/',I8,' NODES (Total:',I9,' DoF):')

         DO 30 I = NOEUD1, NOEUD2
            MN  = MNXYZN + WYZNOE - 3 + 3*I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL

C           ECRITURE DES COMPOSANTES DE LA VITESSE ET
C           EVENTUELLEMENT DE LA PRESSION
            IF( NDIM .EQ. 3 ) THEN

C              ELEMENT FINI 3D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  WRITE(NFUNIT,20023) I, (RMCN(MN+K),K=0,2),
     %                          (VXYZPN(NDL+K,NCAS),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  WRITE(NFUNIT,20033) I, (RMCN(MN+K),K=0,2),
     %                          (VXYZPN(NDL+K,NCAS),K=1,ND)
               ENDIF

            ELSE

C              ELEMENT FINI 2D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  WRITE(NFUNIT,20022) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K,NCAS),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  WRITE(NFUNIT,20032) I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K,NCAS),K=1,ND)
               ENDIF
            ENDIF
 30      ENDDO

      ENDIF

20022 FORMAT('NODE',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,'  PRESSURE=',G15.7)

20032 FORMAT('NODE',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7)

20023 FORMAT('NODE',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,' VZ=',G15.7,
     %'  PRESSURE=',G15.7)

20033 FORMAT('NODE',I8,': X=',G13.5,' Y=',G13.5,' Z=',G13.5,
     %'  VELOCITY VX=',G15.7,' VY=',G15.7,' VZ=',G15.7)

C     ========================================================================
C     MIN MAX MOYENNE des |Vitesse| et PRESSIONS CALCULEES:
C     VITMIN : MINIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C     VITMAX : MAXIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C     VITMOY : MOYENNE DE LA NORME DE LA VITESSE CALCULEE EN TOUS LES NOEUDS
C     PREMIN : MINIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C     PREMAX : MAXIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C     PREMOY : MOYENNE DE LA PRESSION CALCULEE EN TOUS LES SOMMETS DU MAILLAGE
C     ========================================================================
      VITMOY =  0D0
      VITMIN =  1D100
      VITMAX = -1D100
      PREMIN =  1D100
      PREMAX = -1D100
      PREMOY = 0D0
      NBST   = 0

      DO 90 I=1,NBNOVI
         NDL = NODDL(I-1)

C        MODULE DE LA VITESSE CALCULEE AU NOEUD I
         IF( NDIM .EQ. 2 ) THEN
            VITNOR = SQRT( VXYZPN(NDL+1,NCAS)**2
     %                   + VXYZPN(NDL+2,NCAS)**2 )
         ELSE
            VITNOR = SQRT( VXYZPN(NDL+1,NCAS)**2
     %                   + VXYZPN(NDL+2,NCAS)**2
     %                   + VXYZPN(NDL+3,NCAS)**2 )
         ENDIF

C        NORME MOYENNE DU MODULE DE LA VITESSE CALCULEE
         VITMOY = VITMOY + VITNOR

C        NORME MIN et MAX DU MODULE DE LA VITESSE CALCULEE
         IF( VITNOR .LT. VITMIN ) VITMIN = VITNOR
         IF( VITNOR .GT. VITMAX ) THEN
            VITMAX = VITNOR
            NOVMAX = I
         ENDIF

C        PRESSION CALCULEE AU NOEUD I
         IF( NODDL(I)-NDL .GT. NDIM ) THEN
            NBST    = NBST + 1
            PRECALC = VXYZPN(NODDL(I),NCAS)
            PREMOY  = PREMOY + PRECALC
C           PRECALC CALCULEE MIN et MAX EN UN NOEUD
            IF( PRECALC .LT. PREMIN ) PREMIN=PRECALC
            IF( PRECALC .GT. PREMAX ) PREMAX=PRECALC
         ENDIF
 90   ENDDO

C     MODULE MOYEN AUX NOEUDS DE LA VITESSE
      VITMOY = VITMOY / NBNOVI

C     PRESSION MOYENNE AUX SOMMETS
      PREMOY = PREMOY / NBST

C     COORDONNEES DU NOEUD DE VITESSE MAXIMALE
      MN = MNXYZN + WYZNOE - 3 + 3*NOVMAX
      XP = RMCN(MN)
      YP = RMCN(MN+1)
      ZP = RMCN(MN+2)

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NFUNIT,10090)
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMIN, PREMAX, PREMAX-PREMIN
      ELSE
         WRITE(NFUNIT,20090)
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMIN, PREMAX, PREMAX-PREMIN
      ENDIF

10090 FORMAT(/'Au Temps',G13.5,
     %' |VITESSE|Moyenne=',G14.6,
     %' |VITESSE|Max=',G14.6,' au noeud',I9,' XYZ=',3G14.6/
     %        'Au Temps',G13.5,
     %' PRESSION Moyenne=',G14.6,' PRESSION Min=',G14.6,
     %' PRESSION Max=',G14.6,' PRESSION Max-Min=',G14.6 )

20090 FORMAT(/'At Time',G13.5,
     %' |VELOCITY|Mean=',G14.6,
     %' |VELOCITY|Max=',G14.6,' at node',I9,' XYZ=',3G14.6/
     %        'At Time',G13.5,
     %' PRESSURE  Mean=',G14.6,' PRESSURE  Min=',G14.6,
     %' PRESSURE  Max=',G14.6,' PRESSURE Max-Min=',G14.6 )

C     FERMETURE DU FICHIER KNOMFIC
      CLOSE( NFUNIT )
      PRINT*,'ficvipr: TEMPS=',TEMPS,' Fermeture du fichier ',
     %        KNOMFIC(1:NBK)
      NBVPFILE = NBVPFILE + 1
      GOTO 9999

C     TRAITEMENT DES ERREURS
 9998 PRINT*,'ficvipr: Le FICHIER ',NOMFIC,' UNITE',NFUNIT,
     %          ' NE PEUT PAS ETRE DECLARE'
      IERR = 1

 9999 RETURN
      END
