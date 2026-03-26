      SUBROUTINE SEFACLIC( XYZSOM, NCNOST, NBNOST, NOST, NOSOEF, NUFACL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DES ITEMS FACES (TRIANGLES ou QUADRANGLES) 3D
C -----    VISIBLES SUR L'ECRAN et DES XYZ DE LEUR BARYCENTRE
C          SELECTIONNER UNE FACE A L'AIDE D'UN CLIC SOURIS INTERNE
C          RETOURNER LE NUMERO NUFACL DE CETTE FACE VISIBLE CLIQUEE

C ENTREES:
C --------
C XYZSOM : XYZ 3 COORDONNEES CARTESIENNES DES SOMMETS DU MAILLAGE
C NCNOST : >=0 NUMERO DE LA COULEUR DE TRACE DU NUMERO DES SOMMETS NOST
C          <0  PAS DE TRACE DE '+NOST' (GENERALEMENT LE DERNIER CLIQUE)
C NBNOST : NOMBRE DE SOMMETS VISIBLES DE NUMERO A TRACER
C NOST   : NUMERO XYZSOM DU SOMMET A TRACER AU DESSUS DE LA FACE CLIQUEE
C NOSOEF : NO DES 4 SOMMETS DES TRIANGLES ET/OU QUADRANGLES DU MAILLAGE

C SORTIES:
C --------
C NUFACL : >0 NUMERO DE LA FACE VISIBLE LA PLUS PROCHE DU POINT CLIQUE
C             DANS LE PLAN AXONOMETRIQUE => SELECTIONNE
C          =0 SI LE POINT CLIQUE EST HORS MAILLAGE ou ABANDON DEMANDE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      REAL          RMCN(MOTMCN)
      COMMON         MCN(MOTMCN)
      EQUIVALENCE  (RMCN(1),MCN(1))
      COMMON /LECLXQ/LECLEX
      INTEGER        NOSOEF(4,*),NOST(NBNOST)
      REAL           XYZSOM(3,*)
      CHARACTER*40   KTITRE

      NUITFA = 0
      NUFACL0 = 0

C     LES ITEMS FACES VISIBLES SONT DEJA TRACES
C     -----------------------------------------
      LECLEX = 8

      IF( LANGAG .EQ. 0 ) THEN
         KTITRE='CLIQUER une FACE Quadrangle ou Triangle'
      ELSE
         KTITRE='CLICK a Quadrangle or Triangle FACE'
      ENDIF
      CALL SANSBL( KTITRE, NBC )
      CALL TRFINS( KTITRE(1:NBC) )


C     ACTIVATION DE LA SAISIE A LA SOURIS ou FRAPPE AU CLAVIER
C     SAISIE D'UN CLIC (NX,NY) DE LA SOURIS 
C     ou UN DEPLACEMENT du POINTEUR de la SOURIS
C     ou d'UN CARACTERE AU CLAVIER
C     --------------------------------------------------------
 10   CALL XVSOURIS( NOTYEV, NOBOCA, NX, NY )
C     RETOURNE NOTYEV QUI VAUT
C     = 0 Si ABANDON demande par clic du bouton 2 de la souris
C         ou par frappe de la touche Echappement ou @
C     = 1 Si CLIC ENFONCE et RELACHE D'UN BOUTON DE LA SOURIS
C         => NX, NY et NOBOCA=NUMERO DU BOUTON ENFONCE RELACHE
C     =-1 Si CLIC SEULEMENT ENFONCE  D'UN BOUTON DE LA SOURIS
C         => NX, NY et NOBOCA=NUMERO DU BOUTON ENFONCE
C     =-2 Si le POINTEUR de la SOURIS A BOUGE SANS BOUTON ENFONCE
C         => NX, NY
C     = 2 Si FRAPPE D'UN CARACTERE AU CLAVIER
C         => NOBOCA=NUMERO DU CARACTERE DANS LA TABLE ASCII

      IF( NOTYEV .EQ. 2 ) THEN
C        CARACTERE de NUMERO ASCII NOBOCA TAPE AU CLAVIER
         IF( NOBOCA .EQ. 27 ) THEN
cccC        LE CARACTERE 'ECHAPPEMENT' DEVIENT LE CARACTERE '@'
ccc         NOBOCA = 64
            NOTYEV = 0
            GOTO 12
         ELSE
            GOTO 10
         ENDIF
      ENDIF

 12   IF( NOTYEV .EQ. 0 ) THEN
C        ABANDON DE LA RECHERCHE d'une FACE CLIQUEE
         NUFACL = 0
         GOTO 9999
      ENDIF

C     Ici: NOTYEV=1 ou -1 ou -2
C     CLIC DANS L'ECRAN GRAPHIQUE EN (NX,NY) PIXELS
C     RECHERCHE DE LA PLUS PROCHE FACE DU CLIC PARMI LES ITEMFA
C     ---------------------------------------------------------
      CALL SEITCLIC( NX, NY, 8, MNITFA, NUITFA, NUFACL )
      IF( NUITFA .EQ. 0 ) GOTO 10

C     PAS DE TRACE DU NUMERO DES SOMMETS DE LA FACE CLIQUEE
      NCNS = -1

      IF( NUFACL0 .GT. 0 .AND. NUFACL0 .NE. NUFACL ) THEN

C        RESTAURATION DES PIXELS DE LA FENETRE PIXELS MEFISTO
C        SAUVEGARDEE JUSTE AVANT LES MODIFICATIONS
         CALL RESTAUREMEMPX
         CALL MEMPXFENETRE
         CALL XVVOIR

      ENDIF

C     TRACE EN BLANC DE L'ITEM FACE '+NUFACL' CLIQUE
C     ----------------------------------------------
      CALL TIFACE3D( NCBLAN, NCMAGE, NCGRIS, NCNS,
     %               NUFACL, NOSOEF, XYZSOM )

C     TRACE DES NBNOST SOMMETS CLIQUES NOST
C     -------------------------------------
      IF( NCNOST .GE. 0 ) THEN
         DO K = 1, NBNOST
C           TRACE DU NO DE SOMMET '+NOST'
            CALL TRST3D( NCNOST, NOST(K), XYZSOM )
         ENDDO
      ENDIF

      CALL TRFINS( KTITRE(1:NBC) )
      NUFACL0 = NUFACL

      IF( NOTYEV .LT. 0 ) THEN
         GOTO 10
      ENDIF

C     Ici: NOTYEV = 1  NOBOCA EST LE BOUTON CLIQUE et RELACHE
C                 =>   la FACE NUFACL CLIQUEE EST SELECTIONNEE
C     ========================================================
C     L'ADRESSE MCN DU DEBUT DE L'ITEMFA NUITFA
      MNIT = MNITFA + MCN(MNITFA) * NUITFA

C     NOMBRE DE SOMMETS DE LA FACE (A PRIORI 3:TRIANGLE ou 4:QUADRANGLE)
      IF( NOSOEF(4,NUFACL) .EQ. 0 ) THEN
         NBS = 3
      ELSE
         NBS = 4
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'sefaclic: FACE CLIQUEE: No FACE=',NUFACL,
     %          ' Barycentre:',(RMCN(MNIT+K),K=0,2),
     %          ' Sommets:',(NOSOEF(K,NUFACL),K=1,NBS)
      ELSE
         PRINT*,'sefaclic: CLICKED FACE: No FACE=',NUFACL,
     %          ' Barycenter:',(RMCN(MNIT+K),K=0,2),
     %          ' Vertices:',(NOSOEF(K,NUFACL),K=1,NBS)
      ENDIF

 9999 RETURN
      END
