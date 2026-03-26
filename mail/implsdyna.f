      SUBROUTINE IMPLSDYNA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPORTER 2 FICHIERS de CARACTERES ASCII 
C ----- NomSurface.xyzlsd et NomSurface.nselsd ISSUS du LOGICIEL LSDYNA
C       contenant les XYZ des SOMMETS  (->'XYZSOMMET')
C       les NUMEROS DES SOMMETS DU MAILLAGE (->'NSEF')
C       d'une Surface POUR CREER la SURFACE et son MAILLAGE
C       dans le PROJET MEFISTO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY           Janvier 2021
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      CHARACTER*24   KNOMSURF
      CHARACTER*1    KSUFIX(4)
      DATA  KSUFIX / 'p', 'l', 's', 'v' /

C     LECTURE DU NOM De la SURFACE a CREER
 3    CALL INVITE( 42 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMSURF )
      IF( NCVALS .EQ. -1 ) GOTO 9999

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'Construction de la SURFACE: ',KNOMSURF
      ELSE
         PRINT *,'Making of SURFACE: ',KNOMSURF
      ENDIF

C     RECHERCHE DU NOM DE LA SURFACE DANS LE LEXIQUE DES SURFACES
C     -----------------------------------------------------------
      CALL MAJUSC( KNOMSURF )
      CALL LXLXOU( NTSURF, KNOMSURF, NTLXSU, MNLXSU )
C     S'IL EXISTE ALORS IL EST DETRUIT
      IF( NTLXSU .GT. 0 ) THEN
         CALL LXLXDS( NTSURF, KNOMSURF )
      ENDIF

C     CREATION DE LA SURFACE DANS LE REPERTOIRE DES SURFACES
C     ------------------------------------------------------
C     LA SURFACE EST CREEE ET OUVERTE
      CALL LXLXDC( NTSURF, KNOMSURF, 24, 8 )
      CALL LXLXOU( NTSURF, KNOMSURF, NTLXSU, MNLXSU )
      IF( NTLXSU .LE. 0 ) GOTO 3
      CALL NUOBNM( 'SURFACE', KNOMSURF, NUMSURF )

C     LIRE SUR 2 FICHIERS NomSurface.xyzlsd et NomSurface.nselsd
C     UN MAILLAGE D'UNE SURFACE ISSU DU MAILLEUR DU LOGICIEL LS-DYNA
C     LES 3 COORDONNEES DES SOMMETS
C     LES NUMEROS DES SOMMETS DES EF D'UNE SURFACE
C     ET EN CONSTITUER UNE SURFACE MEFISTO
C     --------------------------------------------------------------
      CALL SUEXLSDY( NUMSURF, NTLXSU,
     %               NTNSEF,  MNNSEF, NTXYZS, MNXYZS, IERR )
      IF( IERR .NE. 0 ) THEN
         CALL LXLXDS( NTSURF, KNOMSURF )
         GOTO 9999
      ENDIF

C     LE TMS DEFINITION DU P L S V A UN TMS DEFINITION
C     A PARTIR DES TMS XYZSOMMET et NSEF
C     10: 'DONNEE DIRECTE DES TMS XYZSOMMET ET NSEF'
C     ------------------------------------------------
      CALL LXTNDC( NTLXSU, 'DEFINITION', 'MOTS', WUTSSS+1 )
      CALL LXTSOU( NTLXSU, 'DEFINITION', NTDEFI, MNDEFI )
C     TRANSFORMATION (I pour IDENTITE)
      MCN( MNDEFI + WTYTRS ) = 1

C     TYPE DE DEFINITION DU MAILLAGE DE LA SURFACE
C     10: 'DONNEE DIRECTE DES TMS XYZSOMMET ET NSEF'
      MCN( MNDEFI + WUTYSU ) = 10

C     NUTSOV 'tableau MS xyzsommet' tms ~>>>XYZSOMMET
      MCN( MNDEFI + WUTSOS ) = NTXYZS

C     NUTSSV 'tableau MS no des sommets des EF' tms ~>>>NSEF
      MCN( MNDEFI + WUTSSS ) = NTNSEF

C     AJOUT DE LA DATE du TMS 'DEFINITION'
      CALL ECDATE( MCN(MNDEFI) )

C     AJOUT DE LA DATE du TMS 'XYZSOMMET' POSTERIEURE a DEFINITION
      CALL ECDATE( MCN(MNXYZS) )

C     MISE A JOUR DE LA DIMENSION DES COORDONNEES
      CALL DIMCOO( MCN(MNXYZS+WNBSOM), MCN(MNXYZS+WYZSOM), NDIMLI )

C     MISE A JOUR DU CADRE MINMAX DES NBCOOR COORDONNEES DES SOMMETS
      CALL MAJEXT( MNXYZS )

      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'Creation ','SURFACE ',KNOMSURF,' CORRECTE'
      ELSE
         PRINT *,'CORRECT CREATION of ','SURFACE ',KNOMSURF
      ENDIF
      PRINT*

C     TRACE du MAILLAGE du P L S V
      IF( INTERA .GE. 1 .AND. MNXYZS .GT. 0 ) THEN
         INIEXT = 0
         NOTYVI = 0
         CALL MAJEXT( MNXYZS )
         CALL T1SOBJ( 3,      KNOMSURF, NUMSURF,
     %                NTNSEF, MNNSEF,   NTXYZS, MNXYZS )
      ENDIF

 9999 RETURN
      END
