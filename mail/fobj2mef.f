      SUBROUTINE FOBJ2MEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LE NOM D'UN FICHIER.obj CONTENANT LE MAILLAGE d'une SURFACE
C ----- LIRE LE MAILLAGE SELON LE FORMAT D'UN FICHIER.obj
C       CONSTRUIRE LES TMS XYZSOMMET et NSEF de la SURFACE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray             Avril 2020
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"

      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*96      KNMFIC, KTEXT
      CHARACTER*24      KNOMSF
      CHARACTER*1       KAR
      LOGICAL           LOPEN, LEXIST
      INTEGER           NUSTFA(3,4), NV(24), NBARXF(6)
      DATA              NUSTFA / 1,2,3, 2,3,4, 3,4,1, 4,1,2 /

      IERR = 0

C     LECTURE DU NOM_DU_FICHIER.obj  LE CARACTERE @ POUR FINIR
 5    CALL INVITE( 49 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNMFIC )
      IF( NCVALS .EQ. -1 ) GOTO 9999

C     PASSAGE EN MINUSCULES
      CALL MINUSC( KNMFIC )
      N = INDEX( KNMFIC, '.obj' )
      IF( N .LE. 0 ) THEN 
         NBLGRC(NRERR) = 3
         KERR(1) = KNMFIC
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'Ce FICHIER est INCONNU. Est il en minuscules'
      KERR(3) = 'Redonner le NOM en minuscules du FICHIER.obj'
         ELSE
            KERR(2) = 'This FILE is UNKNOWN. Is it a lower case NAME?'
            KERR(3) = 'GIVE AGAIN the lower case FILE.obj NAME'
          ENDIF
         CALL LEREUR
         GOTO 5
      ENDIF

C     IMPRESSION du NOM DU FICHIER.obj
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*, 'LECTURE du FICHIER.obj de NOM: ', KNMFIC
      ELSE
         PRINT*, 'READING a FILE.obj NAME: ', KNMFIC
      ENDIF

C     SI LE FICHIER N'EXISTE PAS, ALORS ERREUR
      INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( .NOT. LEXIST ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMFIC
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'Le FICHIER N''EXISTE PAS'
         ELSE
            KERR(2) = 'The FILE DOES NOT EXIST'
         ENDIF
         CALL LEREUR
         GOTO 5
      ENDIF
C
C     LE FICHIER KNMFIC EXISTE
      IF( .NOT. LOPEN ) THEN
C        OUVERTURE DU FICHIER
         CALL TRUNIT( NOFILE )
         OPEN( FILE=KNMFIC, UNIT=NOFILE, STATUS='OLD' )
      ENDIF

C     LE NOM de la SURFACE EST REDUIT A AU PLUS 24 CARACTERES
      N = NUDCNB( KNMFIC )
      N = INDEX( KNMFIC(1:N), '.obj' )
      N = N - 1
C     N DERNIER CARACTERE AVANT LE . DE .obj
      IF( N .GT. 24 ) THEN
         N = 24
      ENDIF
C     PASSAGE EN MAJUSCULES
      KNOMSF = KNMFIC(1:N)
      CALL MAJUSC( KNOMSF )

C     RECHERCHE DU NOM DE LA SURFACE DANS LE LEXIQUE DES SURFACES
      IF( NTSURF .LE. 0 ) RETURN
      CALL LXLXOU( NTSURF, KNOMSF, NTLXSF, MNLXSF )

C     S'IL N'EXISTE PAS IL EST CREE ET OUVERT
 8    IF( NTLXSF .LE. 0 ) THEN
         CALL LXLXDC( NTSURF, KNOMSF, 24, 8 )
         CALL LXLXOU( NTSURF, KNOMSF, NTLXSF, MNLXSF )
         ILEXIS = 0
      ELSE
C        S'IL EXISTE CHOIX DE DESTRUCTION OU NON
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(3) = 'L''SURFACE ' // KNOMSF
            L = NUDCNB(KERR(3))
            KERR(1) = KERR(3)(1:L)  // ' EXISTE DEJA'
            KERR(2) = 'VOULEZ VOUS LE DETRUIRE?'
         ELSE
            KERR(3) = 'The OBJECT ' // KNOMSF
            L = NUDCNB(KERR(3))
            KERR(1) = KERR(3)(1:L)  // ' ALREADY EXISTS'
            KERR(2) = 'DO YOU WANT DESTROY IT?'
         ENDIF
         CALL LERESU
         CALL LIMTCL( 'non_oui', N )
         IF( N .EQ. 0 ) GOTO 5
C        DESTRUCTION DE LA SURFACE
         CALL LXLXDS( NTSURF, KNOMSF )
         NTLXSF = 0
         GOTO 8
      ENDIF

C     TRACE DU NOM DE LA SURFACE
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'SURFACE: ' // KNOMSF
      PRINT*, 'CREATION ' // KHIST(1)
      CALL LHISTO
      ILEXIS = 1


C     LE TABLEAU XYZ DES SOMMETS DE LA SURFACE
C     ----------------------------------------
      MXSOMM = 4096*16
      MNXYZS = 0
      CALL TNMCDC( 'REEL', 3*MXSOMM, MNXYZS )
      MNST = MNXYZS

C     LE TABLEAU DES NBSOEF=4 NO DE SOMMETS DES MXEF QTANGLES
C     -------------------------------------------------------
C     NOMBRE DE NO STOCKES POUR UN QTANGLE
      NBSOEF = 4

      MXEF   = 2 * MXSOMM
      MNNSEF = 0
      CALL TNMCDC( 'ENTIER', NBSOEF*MXEF, MNNSEF )
      MNEF = MNNSEF

      NBV = 0
      NBF = 0
      NBTRIA = 0
      NBQUAD = 0
      NBERR  = 0

C     LECTURE SUR LE FICHIER.obj LIGNE par LIGNE
C     ==========================================
 10   READ( NOFILE, '(A)', END=100, ERR=9000 ) KTEXT
ccc      PRINT*, 'fobj2mef:',KTEXT

      NCDER = NUDCNB( KTEXT )
      IF( NCDER .GE. 96 ) THEN
         NCDER = 96
      ELSE
         NCDER = NCDER + 1
      ENDIF
C     AJOUT D'UN BLANC
      KTEXT( NCDER:NCDER) = ' '

C     RECHERCHE d'un vertex
      M = INDEX( KTEXT, 'v ' )
      IF( M .GT. 0 ) THEN

C        UN VERTEX NBV de COORDONNEES X Y Z A RECOLTER
         NBV = NBV + 1
         READ( KTEXT(M+2:NCDER), * ) X, Y, Z
ccc      PRINT*, 'fobj2mef: Vertex=',NBV,': X=',X, ' Y=',Y, ' Z=',Z

         IF( NBV .GT. MXSOMM ) THEN
            CALL TNMCAU( 'REEL', 3*MXSOMM, 2*3*MXSOMM, 3*MXSOMM, MNXYZS)
            MXSOMM = 2*MXSOMM
            MNST   = MNXYZS + 3 * NBV - 3
         ENDIF

         RMCN( MNST     ) = X
         RMCN( MNST + 1 ) = Y
         RMCN( MNST + 2 ) = Z
         MNST = MNST + 3

         GOTO 10

      ENDIF

C     RECHERCHE d'une face ou QTANGLE
      NCD = INDEX( KTEXT, 'f ' )
      IF( NCD .GT. 0 ) THEN

C        FACE v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...
         NBF = NBF + 1

C        RECHERCHE des TRIPLETS  vi/vti/vni
C        NBSF = NOMBRE DE SOMMETS de la FACE
         NBSF = 0
         NCD  = NCD + 2

 18      NCFB = INDEX( KTEXT(NCD:NCDER), ' ' )
         IF( NCFB .GT. 0 ) THEN

C           UN SOMMET DE PLUS DE LA FACE
            NBSF = NBSF + 1
            NCS  = INDEX( KTEXT(NCD:NCDER), '/' )
            NCF  = NCD+NCS-2
            READ( KTEXT(NCD:NCF), * ) NV(NBSF)
            NCD = NCD + NCFB
            IF( NCD .LT. NCDER ) GOTO 18

         ENDIF
ccc         PRINT*, 'fobj2mef: Face',NBF,' Vertices:',(NV(K),K=1,NBSF)

C        STOCKAGE DU NO DES NBSF SOMMETS DU QTANGLE
         IF( NBSF .GT. 4 ) THEN
            PRINT*
            PRINT*,'fobj2mef: PB Face',NBF,' AVEC PLUS de 4 SOMMETS'
            PRINT*,'fobj2mef: Face',NBF,' Vertices:',(NV(K),K=1,NBSF)
            PRINT*,'fobj2mef: Sous Programme A AMELIORER'
            NBERR = NBERR + 1
C           PERTE DU NO DES SOMMETS AU DELA DU 4-EME
            NBSF = 4
         ENDIF

         IF( NBF .GT. MXEF ) THEN
            CALL TNMCAU( 'ENTIER', 4*MXEF, 2*4*MXEF, 4*MXEF, MNNSEF )
            MXEF = 2*MXEF
            MNEF = MNNSEF + NBSOEF * NBF - NBSOEF
         ENDIF

         GOTO( 10, 10, 13, 14 ), NBSF
 13      NBTRIA = NBTRIA + 1
         GOTO 15
 14      NBQUAD = NBQUAD + 1

 15      DO K=1,NBSF
            MCN( MNEF-1+K ) = NV( K )
         ENDDO
         DO K=NBSF+1,NBSOEF
            MCN( MNEF-1+K ) = 0
         ENDDO
         MNEF = MNEF + NBSOEF

         GOTO 10

      ENDIF

      M = INDEX( KTEXT, 'vt' )
      IF( M .EQ. 1 ) THEN
C        VERTEX TEXTURE NON EXPLOITE
         GOTO 10
      ENDIF

      M = INDEX( KTEXT, 'vn' )
      IF( M .EQ. 1 ) THEN
C        VERTEX NORMAL NON EXPLOITE
         GOTO 10
      ENDIF

      M = INDEX( KTEXT, 'vp' )
      IF( M .EQ. 1 ) THEN
C        VERTEX PARAMETER (u,v,...)  NON EXPLOITE
         GOTO 10
      ENDIF

      M = INDEX( KTEXT, 'g ' )
      IF( M .EQ. 1 ) THEN
C        group  NON EXPLOITE
         GOTO 10
      ENDIF

      M = INDEX( KTEXT, 's ' )
      IF( M .EQ. 1 ) THEN
C        s  NON EXPLOITE
         GOTO 10
      ENDIF

      KAR = KTEXT(1:1)
 17   IF( KAR .EQ. '#' .OR. KAR .EQ. ' ' ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*, 'LIGNE NON EXPLOITEE: ',KTEXT
         ELSE
            PRINT*, 'NON EXPLOITED LINE: ',KTEXT
         ENDIF
         GOTO 10
      ENDIF

C     KTEXT INEXPLOITABLE
 19   KAR = ' '
      GOTO 17


C     ERREUR EN LECTURE DE KTEXT
 9000 PRINT*,'fobj2mef: ERREUR en LECTURE sur le FICHIER.obj',KTEXT
      GOTO 19


C     FERMETURE DU FICHIER.obj
C     ------------------------
 100  CLOSE( NOFILE )

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10110) KNOMSF, NBV, NBTRIA, NBQUAD, NBF, NBERR,
     %                       NBSOEF
      ELSE
         WRITE(IMPRIM,20110) KNOMSF, NBV, NBTRIA, NBQUAD, NBF, NBERR,
     %                       NBSOEF
      ENDIF
10110 FORMAT(' LECTURE de la SURFACE du FICHIER.obj: ', A /
     %' NOMBRE de SOMMETS       =',I12/
     %' NOMBRE de TRIANGLES     =',I12/
     %' NOMBRE de QUADRANGLES   =',I12/
     %' NOMBRE de QUA-TRIANGLES =',I12/
     %' NOMBRE de EF avec ERREUR=',I12/
     %' NOMBRE de SOMMETS par QT=',I12/)
20110 FORMAT(' READING of SURFACE.obj: ', A /
     %' NUMBER of VERTICES      =',I12/
     %' NUMBER of TRIANGLES     =',I12/
     %' NUMBER of QUADRANGLES   =',I12/
     %' NUMBER of QUA-TRIANGLES =',I12/
     %' NUMBER of FE with ERROR =',I12/
     %' NUMBER of VERTICES by QT=',I12/)

C     LES XYZ EXTREMES POUR LES TRACES
      COOEXT(1,1) = RMCN(MNXYZS)
      COOEXT(1,2) = RMCN(MNXYZS)
      COOEXT(2,1) = RMCN(MNXYZS+1)
      COOEXT(2,2) = RMCN(MNXYZS+1)
      COOEXT(3,1) = RMCN(MNXYZS+2)
      COOEXT(3,2) = RMCN(MNXYZS+2)
      MNS = MNXYZS + 3
      DO N=2, NBV
         COOEXT(1,1) = MIN( COOEXT(1,1), RMCN(MNS) )
         COOEXT(1,2) = MAX( COOEXT(1,2), RMCN(MNS) )
         COOEXT(2,1) = MIN( COOEXT(2,1), RMCN(MNS+1) )
         COOEXT(2,2) = MAX( COOEXT(2,2), RMCN(MNS+1) )
         COOEXT(3,1) = MIN( COOEXT(3,1), RMCN(MNS+2) )
         COOEXT(3,2) = MAX( COOEXT(3,2), RMCN(MNS+2) )
         MNS = MNS + 3
      ENDDO

      PRINT*,'fobj2mef: XMIN=', COOEXT(1,1),' XMAX=', COOEXT(1,2)
      PRINT*,'fobj2mef: YMIN=', COOEXT(2,1),' YMAX=', COOEXT(2,2)
      PRINT*,'fobj2mef: ZMIN=', COOEXT(3,1),' ZMAX=', COOEXT(3,2)

C     CONSTRUCTION DE LA SURFACE SELON MEFISTO
C     ACTUELLEMENT SEULE UNE SURFACE EST PRISE EN COMPTE...
C     -----------------------------------------------------
C     LE TMS DEFINITION DE LA SURFACE
      CALL LXTNDC( NTLXSF, 'DEFINITION', 'MOTS', WUTSSS+1 )
      CALL LXTSOU( NTLXSF, 'DEFINITION', NTDFSF, MNDFSF )

C     LE TMS XYZSOMMET DE LA SURFACE
      CALL LXTNDC( NTLXSF, 'XYZSOMMET', 'MOTS', WYZSOM+NBV*3 )
      CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTXYZSF, MNXYZSF )

C     LE TMS NSEF DE LA SURFACE
      CALL LXTNDC( NTLXSF,'NSEF','MOTS', WUSOEF+NBSOEF*NBF )
      CALL LXTSOU( NTLXSF,'NSEF', NTNSESF, MNNSESF )

C     SURFACE DEFINIE PAR SES TMS XYZSOMMET NSEF
      MCN( MNDFSF + WTYTRS ) = 1
      MCN( MNDFSF + WUTYSU ) = 10
      MCN( MNDFSF + WUTSOS ) = NTXYZSF
      MCN( MNDFSF + WUTSSS ) = NTNSESF
C
C     LA DATE DU TMS DEFINITION DE LA SURFACE
      CALL ECDATE( MCN(MNDFSF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNDFSF + MOTVAR(6) ) = NONMTD( '~>SURFACE>>DEFINITION' )


C     RECHERCHE DES FACES ISOLEES ET SUPPRESSION
C     ------------------------------------------
C     CONSTRUCTION DU TABLEAU DES ARETES DES FACES DE LA SURFACE
      MXFAAR = 6
      MNARFA = 0
      CALL GEARFA( MCN(MNXYZS), NBSOEF, NBF, MCN(MNNSEF), MXFAAR,
     %             L1ARFA, L2ARFA, MNARFA, NBARXF, IERR )
      IF( IERR .LT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'fobj2mef: TABLEAU des ARETES NON CALCULABLE'
            PRINT*,'fobj2mef: LA SURFACE Mefisto',KNOMSF,
     %             ' N''EST PAS CREEE'
         ELSE
            PRINT*,'fobj2mef: the EDGE ARRAY IS NOT COMPUTABLE'
            PRINT*,'fobj2mef: the Mefisto ',KNOMSF,
     %             ' IS NOT CREATED'
         ENDIF
C        DESTRUCTION DU LEXIQUE DE LA SURFACE KNOMSF
         CALL LXLXDS( NTSURF, KNOMSF )
         GOTO 9999
      ENDIF
      CALL TUERTA1F( NBF, MCN(MNNSEF), L1ARFA, L2ARFA, MCN(MNARFA),
     %               NBFSUP )
      IF( MNARFA.GT.0 ) CALL TNMCDS('ENTIER',L1ARFA*L2ARFA,MNARFA)


C     NETTOYAGE DU MAILLAGE LU Sur le fichier.obj
C     -------------------------------------------
C     MISE A JOUR DU TABLEAU XYZSOM ET NSEF EN RENUMEROTANT
C     LES SOMMETS ET EN ELIMINANT LES EF DESACTIVES ou
C     AYANT 2 SOMMETS DE MEME NUMERO
C     NEF( NBSOEF:L1SOEF, NoEF ) SONT MIS A ZERO
C     NEWS EST UN TABLEAU AUXILLIAIRE
      MXS    = 1+NBV
      MNNEWS = 0
      CALL TNMCDC( 'ENTIER', MXS, MNNEWS )
      CALL MAJXYZNSE( NBV,    MCN(MNXYZS), MCN(MNNEWS),
     %                NBSOEF, NBF,         MCN(MNNSEF) )
      CALL TNMCDS( 'ENTIER', MXS, MNNEWS )

C     TMS xyzsommet de la SURFACE
      MCN( MNXYZSF + WNBSOM ) = NBV
      MCN( MNXYZSF + WNBTGS ) = 0
      MCN( MNXYZSF + WBCOOR ) = 3
C     LES XYZ DES NBV SOMMETS
      CALL TRTATA( RMCN(MNXYZS), RMCN(MNXYZSF+WYZSOM), 3*NBV )
C     LA DATE DU TMS XYZSOMMET DE LA SURFACE NSF
      CALL ECDATE( MCN(MNXYZSF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZSF + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C     TMS nsef de la SURFACE
C     TYPE SURFACE
      MCN( MNNSESF + WUTYOB ) = 3
C     FERMETURE INCONNUE du MAILLAGE
      MCN( MNNSESF + WUTFMA ) = -1
C     Nombre de sommets par EF
      MCN( MNNSESF + WBSOEF ) = NBSOEF
C     Nombre de tangentes par EF
      MCN( MNNSESF + WBTGEF ) = 0
C     Nombre des QTANGLES de la SURFACE
      MCN( MNNSESF + WBEFOB ) = NBF
C     Nombre des EF avec TG
      MCN( MNNSESF + WBEFTG ) = 0
C     Nombre des EF avec POINTEUR sur EF
      MCN( MNNSESF + WBEFAP ) = 0
C     Numero de type du maillage NON STRUCTURE
      MCN( MNNSESF + WUTYMA ) = 0
C     LES NBSOEF NO DES SOMMETS DES NBF QTANGLES
      CALL TRTATA( RMCN(MNNSEF), RMCN(MNNSESF+WUSOEF), NBSOEF*NBF )
C     LA DATE DU TMS NSEF DE LA SURFACE NSF
      CALL ECDATE( MCN(MNNSESF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSESF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C     DESTRUCTIONS DES TABLEAUX MCN AUXILIAIRES
      CALL TNMCDS( 'ENTIER', NBSOEF*MXEF, MNNSEF )
      CALL TNMCDS( 'REEL',   3*MXSOMM,    MNXYZS )

C     AFFICHAGE DE LA QUALITE du MAILLAGE FINAL
      CALL IMPQUA( 3, KNOMSF, MNNSESF, MNXYZSF, NBEFMQ, QUAMIN,
     %             SURVOLEF )

C     SORTIE
 9999 RETURN
      END
