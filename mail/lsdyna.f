      SUBROUTINE LSDYNA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SORTIR SUR UN FICHIER LES COORDONNEES DES SOMMETS
C ----- ET LES NUMEROS DES SOMMETS DU MAILLAGE D'UNE SURFACE + UN VOLUME
C       SELON LE FORMAT DE LECTURE DU MAILLAGE DU LOGICIEL LSDYNA

C SORTIE :
C --------
C LE FICHIER lsdyna_mesh DANS LE REPERTOIRE DU PROJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Veulettes sur mer    Octobre 2015
C MODIFS : ALAIN PERRONNET LJLL UPMC & Veulettes sur mer    Mars    2017
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      COMMON / UNITES /  LECTEU, IMPRIM, NUNITE(30)
      CHARACTER          KDATE*8,KCREAT*24,KINFO*80,KHEURE*8

      CHARACTER*80       KAUX
      CHARACTER*24       NMSURF, NMVOLU
      CHARACTER*34       KNMFIC
      LOGICAL            LEXIST, LOPEN
      INTEGER            NOSOEF(8)

C     NOMBRE DE SURFACE A TRAITER
      NBSF = 0

C     NOMBRE DE VOLUME A TRAITER
      NBVL = 0

C     NUMERO DU MATERIAU
      NUMATE = 0

C     NOMBRE DE SOMMETS
      NBSOMM = 0

C     NOMBRE D'EF
      NBEF  = 0
      NBEFS = 0
      NBEFV = 0

C     NOM DU TYPE DU SV A TRAITER
C     ---------------------------
 10   CALL LIMTCL( 'typ_sv', NUTYSV )
      IF( NUTYSV .LT. 0 ) RETURN
C
C     NOM DU PLSV A TRAITER
      GOTO( 11, 12, 13, 20, 100 ),NUTYSV

 11   PRINT *,'lsdyna: TYPE POINT NON TRAITE'
      NBLGRC(NRERR) = 1
      KERR(1) ='ERREUR: TYPE POINT NON TRAITE'
      CALL LEREUR
      GOTO 10

 12   PRINT *,'lsdyna: TYPE LIGNE NON TRAITE'
      NBLGRC(NRERR) = 1
      KERR(1) ='ERREUR: TYPE LIGNE NON TRAITE'
      CALL LEREUR
      GOTO 10

C     LECTURE DU NOM DE LA SURFACE MATERIAU 1 POUR LS-DYNA
C     ----------------------------------------------------
 13   CALL INVITE( 42 )
      IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, NMSURF )
      IF( NCVALS .EQ. -1 ) RETURN

C     RECHERCHE DU NOM DE LA SURFACE DANS LE LEXIQUE DES SURFACES
      CALL LXLXOU( NTSURF, NMSURF, NTLXSU , MNLXSU )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMSURF
         KERR(2) ='ERREUR: SURFACE INCONNUE'
         CALL LEREUR
         CALL LXIM( NTSURF )
         GOTO 10
      ENDIF
C
C     RECHERCHE DU TABLEAU XYZSOMMET de la SURFACE
      CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTXYZSU, MNXYZSU )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE LA SURFACE
      IF( NTXYZSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMSURF
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'PAS DE TMS XYZSOMMET de la SURFACE'
         ELSE
            KERR(2) = 'TMS XYZSOMMET UNKNOWN of the SURFACE'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     LE NOMBRE DE SOMMETS DU MAILLAGE de la SURFACE
      NBSOSU = MCN( MNXYZSU + WNBSOM )
C
C     RECHERCHE DU TABLEAU NSEF de la SURFACE
      CALL LXTSOU( NTLXSU, 'NSEF', NTNSEFSU, MNNSEFSU )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE LA SURFACE
      IF( NTNSEFSU .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMSURF
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'PAS DE TMS NSEF de la SURFACE'
         ELSE
            KERR(2) = 'TMS NSEF UNKNOWN of the SURFACE'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     LE NOMBRE D'EF du MAILLAGE de la SURFACE
      NBEFSU = MCN( MNNSEFSU + WBEFOB )

C     NOMBRE DE SURFACES
      NBSF = NBSF + 1
      GOTO 30

C     LECTURE DU NOM DU VOLUME MATERIAU 2 POUR LS-DYNA
C     ------------------------------------------------
 20   CALL INVITE( 60 )
      IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, NMVOLU )
      IF( NCVALS .EQ. -1 ) RETURN

C     RECHERCHE DU NOM DU VOLUME DANS LE LEXIQUE DES VOLUMES
      CALL LXLXOU( NTVOLU, NMVOLU, NTLXVO , MNLXVO )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         KERR(2) ='ERREUR: VOLUME INCONNU'
         CALL LEREUR
         CALL LXIM( NTVOLU )
         GOTO 20
      ENDIF
C
C     RECHERCHE DU TABLEAU XYZSOMMET DU VOLUME
      CALL LXTSOU( NTLXVO, 'XYZSOMMET', NTXYZVO, MNXYZVO )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DU VOLUME
      IF( NTXYZVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'PAS DE TMS XYZSOMMET du VOLUME'
         ELSE
            KERR(2) = 'TMS XYZSOMMET UNKNOWN of the VOLUME'
         ENDIF
         CALL LEREUR
         GOTO 20
      ENDIF

C     LE NOMBRE DE SOMMETS DU MAILLAGE DU VOLUME
      NBSOVO = MCN( MNXYZVO + WNBSOM )
C
C     RECHERCHE DU TABLEAU NSEF du VOLUME
      CALL LXTSOU( NTLXVO, 'NSEF', NTNSEFVO, MNNSEFVO )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DU VOLUME
      IF( NTNSEFVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'PAS DE TMS NSEF du VOLUME'
         ELSE
            KERR(2) = 'TMS NSEF UNKNOWN of the VOLUME'
         ENDIF
         CALL LEREUR
         GOTO 20
      ENDIF

C     LE NOMBRE D'EF du MAILLAGE du VOLUME
      NBEFVO = MCN( MNNSEFVO + WBEFOB )

C     NOMBRE DE VOLUMES
      NBVL = NBVL + 1

C     TRAITEMENT DU MATERIAU NOMATE
C     -----------------------------
 30   NUMATE = NUMATE + 1

      IF( NUMATE .EQ. 1 ) THEN

C        DECLARATION OUVERTURE DU FICHIER lsdyna_mesh
C        --------------------------------------------
         KNMFIC = 'lsdyna_mesh'

C        SI LE FICHIER KNMFIC EXISTE ALORS IL EST DETRUIT PUIS RECONSTRUIT
         INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
         IF( LEXIST ) THEN
C           LE FICHIER KNMFIC EXISTE
            IF( .NOT. LOPEN ) THEN
C           OUVERTURE DU FICHIER KNMFIC
               CALL TRUNIT( NF )
               OPEN( FILE=KNMFIC, UNIT=NF, STATUS='OLD' )
            ENDIF
C           DESTRUCTION DU FICHIER
            CLOSE( NF, STATUS='DELETE' )
         ENDIF

C        CREATION DU FICHIER KNMFIC
         CALL TRUNIT( NF )
         OPEN( UNIT=NF, ERR=9999, STATUS='NEW',
     %      FILE=KNMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )

C        ECRITURE DE LA PREMIERE LIGNE DU FICHIER lsdyna
C        DATE
         KAUX = KINFO( 'DATE' )
         KDATE(1:8) = KAUX(1:8)

C        HEURE
         KAUX = KINFO( 'HEURE' )
         KHEURE(1:8) = KAUX(1:8)

C        NOM DE L'UTILISATEUR
         KAUX = KINFO( 'UTILISATEUR' )
         KCREAT(1:24) = KAUX(1:24)

10001    FORMAT('$Date ',A8,'   ',A2,'h ',A2,'m ',A2,'s  ',A)
         WRITE(NF,10001) KDATE(1:8), KHEURE(1:2), KHEURE(3:4),
     %                   KHEURE(5:6), KCREAT

         WRITE(IMPRIM,*)
         L = NUDCNB( KNMFIC )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'En tete du fichier ',KNMFIC(1:L)
         ELSE
            WRITE(IMPRIM,*)'Header of file ',KNMFIC(1:L)
         ENDIF
         WRITE(IMPRIM,10001) KDATE(1:8), KHEURE(1:2), KHEURE(3:4),
     %                       KHEURE(5:6), KCREAT

C        SANS les TANGENTES dans le FICHIER lsdyna_mesh
         NOTG = 0

      ENDIF

      WRITE(IMPRIM,*)

      IF( NUTYSV .EQ. 3 ) THEN

C        SURFACE
C     ------------------------------------------------------------------
      WRITE( NF,   '(A5)' ) '*NODE'
      WRITE( NF,   '(A7,I8,A21,A24)')'$XYZ of',NBSOSU,
     %     ' vertices of SURFACE ',NMSURF
      WRITE(IMPRIM,'(A7,I8,A21,A24)')'$XYZ of',NBSOSU,
     %     ' vertices of SURFACE ',NMSURF

C     LES XYZ DES SOMMETS DE LA SURFACE
      MN = MNXYZSU + WYZSOM -1
      DO N = 1, NBSOSU
         WRITE( NF, '(I8,3(G16.8),I8,I8)' )N+NBSOMM, (RMCN(MN+L),L=1,3),
     %                                     0, 0
         MN = MN + 3
      ENDDO

C     ECRITURE DE LA PREMIERE LIGNE DES NUMEROS DES SOMMETS DES EF DE LA SURFACE
C     LEUR NUMERO EST TRANSLATE DU NOMBRE DE SOMMETS des PRECEDENTS
      WRITE( NF,   '(A14)' ) '*ELEMENT_SHELL'
      WRITE( NF,   '(A9,I8,A15,A24)' ) '$Nodes of',NBEFSU,
     %             ' FE of SURFACE ',NMSURF
      WRITE(IMPRIM,'(A9,I8,A15,A24)' ) '$Nodes of',NBEFSU,
     %             ' FE of SURFACE ',NMSURF

C     LES NUMEROS DES SOMMETS DE LA SURFACE DE NUMERO MATERIAU 1 POUR LSDYNA
      MN = MNNSEFSU + WUSOEF -1
      DO N = 1, NBEFSU
         DO L=1,3
            MCN(MN+L) = MCN(MN+L) + NBSOMM
         ENDDO
         NS4 = MCN(MN+4)
         IF( NS4 .EQ. 0 ) THEN
C           LE 4-EME SOMMET EST IDENTIQUE AU 3-EME POUR LSDYNA
            NS4 = MCN(MN+3)
         ELSE
            NS4 = MCN(MN+4)+NBSOMM
         ENDIF
         WRITE( NF, '(10I8)' ) N+NBEF, NUMATE, (MCN(MN+L),L=1,3), NS4,
     %                         0, 0, 0, 0
         MN = MN + 4
      ENDDO

      PRINT *
      L = NUDCNB( KNMFIC )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'Sur le fichier ',KNMFIC(1:L),' sont AJOUTES'
         PRINT *,  NBSOSU,' SOMMETS DE LA SURFACE ',NMSURF
         PRINT *,  NBEFSU,' QUADRANGLES'
      ELSE
         PRINT *, 'On the file ',KNMFIC(1:L),' are ADDED'
         PRINT *,  NBSOSU,' VERTICES of the SURFACE ',NMSURF
         PRINT *,  NBEFSU,' QUADRANGLES'
      ENDIF

      NBSOMM = NBSOMM + NBSOSU
      NBEFS  = NBEFS  + NBEFSU

      ELSE

C        VOLUME
C     ----------------------------------------------------------------
C     LES SOMMETS DU VOLUME SONT AJOUTES APRES CEUX 
C     ET LEUR NUMERO EST TRANSLATE DU NOMBRE DE SOMMETS des PRECEDENTS

C     LES XYZ DES SOMMETS DU VOLUME
      WRITE( NF,   '(A5)' ) '*NODE'
      WRITE( NF,   '(A7,I8,A20,A24)') '$XYZ of',NBSOVO,
     %             ' vertices of VOLUME ',NMVOLU
      WRITE(IMPRIM,'(A7,I8,A20,A24)') '$XYZ of',NBSOVO,
     %             ' vertices of VOLUME ',NMVOLU

      MN = MNXYZVO + WYZSOM -1
      DO N = 1, NBSOVO
         WRITE( NF, '(I8,3(G16.8),I8,I8)' ) N+NBSOMM,(RMCN(MN+L),L=1,3),
     %                                      0, 0
         MN = MN + 3
      ENDDO

C     ECRITURE DE LA PREMIERE LIGNE DES NUMEROS DES SOMMETS DES EF DU VOLUME
C     LEUR NUMERO EST TRANSLATE DU NOMBRE DE SOMMETS de la SURFACE
      WRITE( NF,   '(A14)' ) '*ELEMENT_SOLID'
      WRITE( NF,   '(A9,I8,A14,A24)' ) '$Nodes of',NBEFVO,
     %             ' FE of VOLUME ',NMVOLU
      WRITE(IMPRIM,'(A9,I8,A14,A24)' ) '$Nodes of',NBEFVO,
     %             ' FE of VOLUME ',NMVOLU

C     LES NUMEROS DES SOMMETS DE LA SURFACE DE NUMERO MATERIAU 2 POUR LSDYNA
      MN = MNNSEFVO + WUSOEF -1
      DO N = 1, NBEFVO

C        LS-DYNA:  LS-DYNA_manual_Vol_I_R7.1.pdf   page 1446
C                  NUMEROTATION DES SOMMETS DANS LES EF
         IF( MCN(MN+5) .EQ. 0 .AND. MCN(MN+6) .EQ. 0 .AND.
     %       MCN(MN+7) .EQ. 0 .AND. MCN(MN+8) .EQ. 0 ) THEN

C           TETRAEDRE Mefisto: NO du 5-EME au 8-EME SOMMET =0 SI TETRAEDRE
C           4-noded tetrahedron N1, N2, N3, N4, N4, N4, N4, N4, 0, 0
C           MEME NUMEROTATION DES 4 SOMMETS DANS LS-DYNA et Mefisto
            DO L=1,4
               NOSOEF(L)   = MCN(MN+L)
               NOSOEF(L+4) = MCN(MN+4)
            ENDDO

         ELSE IF( MCN(MN+7) .EQ. 0 .AND. MCN(MN+8) .EQ. 0 ) THEN

C          PENTAEDRE Mefisto: 7-EME et 8-EME SOMMET =0 SI PENTAEDRE
C          LS-DYNA: 6-noded pentahedron N1, N2, N3, N4, N5, N5, N6, N6
C          5_________6
C         /|        /|
C        1_________4 |
C         \|       \ |
C          \2_______\3

C          Mefisto: LE NUMERO du 7-EME et 8-EME SOMMET =0 SI PENTAEDRE, >0 SI HEXAEDRE
C          3_________6
C         /|        /|
C        1_________4 |
C         \|       \ |
C          \2_______\5

            NOSOEF(1) = MCN(MN+1)
            NOSOEF(2) = MCN(MN+2)
            NOSOEF(3) = MCN(MN+5)
            NOSOEF(4) = MCN(MN+4)
            NOSOEF(5) = MCN(MN+3)
            NOSOEF(6) = MCN(MN+3)
            NOSOEF(7) = MCN(MN+6)
            NOSOEF(8) = MCN(MN+6)

         ELSE

C           HEXAEDRE: NUMEROTATION IDENTIQUE DANS LS-DYNA et Mefisto
            DO L=1,8
               NOSOEF(L) = MCN(MN+L)
            ENDDO

         ENDIF

C        ECRITURE DES 8 NUMEROS DES SOMMETS DE l'HEXAEDRE
         DO L=1,8
            NOSOEF(L) = NOSOEF(L) + NBSOMM
         ENDDO
         WRITE( NF,'(10I8)' ) N+NBEF, NUMATE, (NOSOEF(L),L=1,8)
         MN = MN + 8

      ENDDO

      PRINT *
      L = NUDCNB( KNMFIC )
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'Sur le fichier ',KNMFIC(1:L),' sont AJOUTES'
         PRINT *,  NBSOVO,' SOMMETS DU VOLUME ', NMVOLU
         PRINT *,  NBEFVO,' HEXAEDRES'
      ELSE
         PRINT *, 'On the file ',KNMFIC(1:L),' are ADDED'
         PRINT *,  NBSOVO,' VERTICES of the VOLUME ', NMVOLU
         PRINT *,  NBEFVO,' HEXAEDRA'
      ENDIF

      NBSOMM = NBSOMM + NBSOVO
      NBEFV  = NBEFV  + NBEFVO
      NBEF   = NBEFS  + NBEFV

      ENDIF

      GOTO 10

C     LA FIN DE FICHIER POUR LS-DYNA
C     ------------------------------
 100  WRITE( NF, '(A4)' ) '*END'

C     SORTIE CORRECTE
      CLOSE( NF )

C     AU TOTAL
      L = NUDCNB( KNMFIC )
      PRINT *
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'Sur le fichier', KNMFIC(1:L),' ont ete STOCKES'
         PRINT *,  NBSF,   ' SURFACES'
         PRINT *,  NBVL,   ' VOLUMES'
         PRINT *,  NBSOMM, ' SOMMETS'
         PRINT *,  NBEFS,  ' QUADRANGLES'
         PRINT *,  NBEFV,  ' HEXAEDRES'
      ELSE
         PRINT *, 'On the file ',KNMFIC(1:L),' are STORED'
         PRINT *,  NBSF,   ' SURFACES'
         PRINT *,  NBVL,   ' VOLUMES'
         PRINT *,  NBSOMM, ' VERTICES'
         PRINT *,  NBEFS,  ' QUADRANGLES'
         PRINT *,  NBEFV,  ' HEXAEDRA'
      ENDIF

      RETURN
C
C     PROBLEME A L'OUVERTURE DU FICHIER KNMFIC
 9999 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'PB A L''OUVERTURE DU FICHIER ' // KNMFIC
      ELSE
         KERR(1) = KNMFIC // 'FILE CAN NOT BE OPENED'
      ENDIF
      CALL LEREUR
      GOTO 10

      END
