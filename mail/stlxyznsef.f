      SUBROUTINE STLXYZNSEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   EXPORTER SUR UN FICHIER .stl en CARACTERES ASCII
C -----   la TRIANGULATION d'une SURFACE 3D FERMEE
C         c-a-d Le nom de la surface fermee
C               le nombre de triangles
C               Pour chaque triangle
C                   XYZN du VECTEUR NORMAL dirige vers l'exterieur
C                   XYZ1 du SOMMET 1 du triangle
C                   XYZ1 du SOMMET 2 du triangle
C                   XYZ1 du SOMMET 3 du triangle dans le sens direct
C
C         obtenus a partir des TMS XYZSOMMET et NSEF Mefisto
C         de la triangulation de la SURFACE 3D FERMEE
C
C SORTIE :
C --------
C         LE FICHIER ASCII  KNMSFFR.stl  est ECRIT dans le
C         REPERTOIRE du PROJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY          Decembre 2024
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
      CHARACTER*10       NMTYOB
      CHARACTER*24       KNMSFFR
      CHARACTER*1        KSUFIX(4)
      DATA               KSUFIX / 'p', 'l', 's', 'v' /

C     NOM DE LA SURFACE FERMEE A TRAITER
 10   NUTYOB = 3
      CALL INVITE( 42 )
      IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNMSFFR )
      IF( NCVALS .EQ. -1 ) GOTO 9999

C     RECHERCHE DU NOM DE LA SURFACE FERMEE DANS LE LEXIQUE DES SURFACES
      CALL LXLXOU( NTMN(NUTYOB), KNMSFFR, NTLXSFFR, MNLXSFFR )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXSFFR .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMSFFR
         KERR(2) ='ERREUR:' // NMTYOB(NUTYOB) // 'INCONNUE'
         CALL LEREUR
         CALL LXIM( NTMN(NUTYOB) )
         GOTO 10
      ENDIF
C     RECHERCHE DU NUMERO DE LA SURFACE DANS LE LEXIQUE des SURFACES
      CALL LXNMNU( NTMN(NUTYOB), KNMSFFR, N, NUSFFR, MNLXSFFR )
      PRINT*,'stlxyznsef: EXPORT SURFACE ', KNMSFFR, 'No=',NUSFFR

C     RECHERCHE DU TABLEAU XYZSOMMET
      CALL LXTSOU( NTLXSFFR, 'XYZSOMMET', NTXYZS, MNXYZS )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE LA SURFACE
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMSFFR
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'PAS de TMS XYZSOMMET'
         ELSE
            KERR(2) = 'TMS XYZSOMMET UNKNOWN'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     LA SURFACE EST ELLE 3D?
C     DIMENSION NDIM DE L'ESPACE DES COORDONNEES DE LA SURFACE
      CALL DIMCOO( NBSOM, MCN(MNXYZS+WYZSOM), NDIM )
      IF( NDIM .NE. 3 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = KNMSFFR
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE NON 3D'
         ELSE
            KERR(2) = 'N0 3D SURFACE'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

      NBSOM = MCN(MNXYZS+WNBSOM)
      IF( NBSOM .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'INCORRECT Nombre Sommets=',NBSOM,' de ',KNMSFFR
         ELSE
            PRINT*,'INCORRECT Vertice Number=',NBSOM,' of ',KNMSFFR
         ENDIF
         GOTO 10
      ENDIF

C     RECHERCHE DU TABLEAU NSEF DE LA SURFACE FERMEE
      CALL LXTSOU( NTLXSFFR, 'NSEF', NTNSEF, MNNSEF )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE LA SURFACE
      IF( NTNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMSFFR
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'ABSENCE de TMS NSEF de la SURFACE'
         ELSE
            KERR(2) = 'NSEF SURFACE TMS is UNKNOWN'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     LA SURFACE EST ELLE FERMEE ?
C     SON TYPE DE FERMETURE EST FORCE INCONNU => MISE A JOUR
      MCN( MNNSEF + WUTFMA ) = -1
cccC     SANS AFFICHAGE DES ARETES SIMPLES (DANS 1 SEULE FACE)
ccc      CALL OBJFER( NUTYOB, NUSFFR, 0, NOFERM )
C     AVEC AFFICHAGE DES ARETES SIMPLES (DANS 1 SEULE FACE)
      CALL OBJFER( NUTYOB, NUSFFR, 1, NOFERM )
C     NOFERM : =1 SI L'OBJET EST FERME
C              =0 SI L'OBJET N'EST PAS FERME
C              <0 SI SATURATION DU TABLEAU SERVANT AU HACHAGE
      IF( NOFERM .NE. 1 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMSFFR
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'SURFACE NON FERMEE'
         ELSE
            KERR(2) = 'UNCLOSED SURFACE'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     DECLARATION OUVERTURE ECRITURE du FICHIER NOM_SFFR.stl
C     ======================================================
      CALL STLXYZNSEF1( KNMSFFR, MNXYZS, MNNSEF, IERR )
      GOTO 10

 9999 RETURN
      END
