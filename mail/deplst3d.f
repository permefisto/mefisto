      SUBROUTINE DEPLST3D( NBSOM,  XYZSOM, NBFACE, NOSOEF,
     %                     NUMST0, XYZPT0, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   DEPLACER UN SOMMET 3D DANS UN MAILLAGE XYZSOMMET D'UNE SURFACE
C -----
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE DE LA SURFACE
C XYZSOM : XYZ 3 COORDONNEES CARTESIENNES DES SOMMETS DU MAILLAGE
C NBFACE : NOMBRE DE FACES TRIANGLES OU QUADRANGLES DU MAILLAGE
C NOSOEF : NO DES 4 SOMMETS DES TRIANGLES ET/OU QUADRANGLES DU MAILLAGE

C SORTIES:
C --------
C NUMST0 : >0 NUMERO DU SOMMET VISIBLE CLIQUE A DEPLACER et DEPLACE
C          =0 SI AUCUN TRIANGLE OU QUADRANGLE NE CONTIENT LE POINT CLIQUE
C             ou ABANDON DEMANDE
C XYZPT0 : XYZ DU SOMMET CLIQUE NUMST0 AVANT DEPLACEMENT
C IERR   : =0 SI PAS D'ERREUR
C          =1 SI ERREUR ou ABANDON DEMANDE NUMST0=0 XYZPT0 NON INITIALISE

C REMARQUE: XYZSOM(NUMST0) en SORTIE SONT les XYZ du POINT CLIQUE
C           dans la FACE si NUMST0>0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C MODIFS : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      REAL          XYZSOM(3,NBSOM), XYZPT0(3), XYZPT1(3)
      INTEGER       NOSOEF(4,NBFACE), NOST(4)
      CHARACTER*72  KTITRE

      IERR   = 0
      NUMST0 = 0

C     SELECTIONNER A L'AIDE D'UN CLIC SOURIS UN SOMMET DU MAILLAGE
C     VISIBLE DANS LA FENETRE GRAPHIQUE SUR L'ECRAN C-A-D
C     RETOURNER LE NUMERO NUMST0 DU SOMMET VISIBLE LE PLUS PROCHE
C     DU POINT CLIQUE EN NX, NY
C     ------------------------------------------------------------
      CALL SESTCLIC( XYZSOM,  ITST0, NUMST0 )
      IF( NUMST0 .LE. 0 ) GOTO 9990

C     PROTECTION DES ANCIENNES COORDONNEES DU SOMMET NUMST0 A DEPLACER
C     ----------------------------------------------------------------
      DO K = 1, 3
         XYZPT0( K ) = XYZSOM( K, NUMST0 )
      ENDDO

C     TRACE DU MAILLAGE + LA FACE CLIQUEE + LE SOMMET CLIQUE a DEPLACER
C     -----------------------------------------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         KTITRE='CLIQUER la NOUVELLE POSITION du SOMMET              dan
     %s la FACE CLIQUEE'
         WRITE( KTITRE(40:48), '(I9)' ) NUMST0
      ELSE
         KTITRE='CLICK the NEW POSITION of the VERTEX            inside 
     %the CLICKED FACE'
         WRITE( KTITRE(38:46), '(I9)' ) NUMST0
      ENDIF
      CALL SANSDBL( KTITRE, NBC )
      CALL TRFINS( KTITRE(1:NBC) )


C     SELECTIONNER A L'AIDE D'UN CLIC SOURIS UNE FACE DU MAILLAGE
C     PRESENTE DANS LA FENETRE GRAPHIQUE C-A-D
C     RETOURNER LE NUMERO NOFACL DE LA FACE VISIBLE DONT LE CLIC
C     EST INTERNE POUR UN MAILLAGE d'une SURFACE 2D ou 3D
C     -----------------------------------------------------------
      NBNOST = 1
      NOST(NBNOST) = NUMST0
      CALL SEFACLIC( XYZSOM, NCJAUN, NBNOST, NOST, NOSOEF, NOFACL )
      IF( NOFACL .LE. 0 ) GOTO 9990


C     CLIC D'UNE POSITION NOUVELLE DU SOMMET NUMST0 DANS LA FACE NOFACL
C     -----------------------------------------------------------------
      CALL PTFACLIC( 3, XYZSOM, NOFACL, NOSOEF,  XYZPT1, IERR )
      IF( IERR .NE. 0 ) GOTO 9990


C      MODIFICATION DES COORDONNEES OBJET DU SOMMET NUMST0 CLIQUE
C      PAR LES 3 COORDONNEES DU POINT CLIQUE
C      ----------------------------------------------------------
       IF( LANGAG .EQ. 0 ) THEN
          PRINT*,'deplst3d: DEPLACEMENT du SOMMET',NUMST0,' de',
     %           (XYZSOM(K,NUMST0),K=1,3),' a',(XYZPT1(K),K=1,3)
       ELSE
          PRINT*,'deplst3d: DISPLACEMENT of VERTEX',NUMST0,' from',
     %           (XYZSOM(K,NUMST0),K=1,3),' to',(XYZPT1(K),K=1,3)
       ENDIF

C      MODIFICATION DES XYZ DU SOMMET NUMST0 EN LES XYZ DU POINT CLIQUE
C      ----------------------------------------------------------------
       DO K=1,3
          XYZSOM(K,NUMST0) = XYZPT1(K)
       ENDDO
       GOTO 9999


C      PAS DE POINT CLIQUE
 9990  NUMST0 = 0
       IERR   = 1


 9999  RETURN
       END
