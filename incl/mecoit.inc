C     TOUTES LES OPTIONS DES TRACES GRAPHIQUES
C     fichier mecoit.inc
      PARAMETER        (MLATTR=54)

C     LA GESTION des ITEMS A TRACER   -----------------------------------
      COMMON / MECOI1 / MNITEP, MNITEL, MNITES, MNITEV, MNITEO,
     %                  MNITST, MNITAR, MNITFA, MNITEF
      INTEGER           MNITEM(9)
      EQUIVALENCE      (MNITEM(1),MNITEP)
C     MNITEP,...,MNITEF : ADRESSE MCN DU TABLEAU DES ITEMS
C                         1:POINTS,  2:LIGNES, 3:SURFACES, 4:VOLUMES, 5:OBJETS,
C                         6:SOMMETS, 7:ARETES, 8:FACES,    9:EF3D

C     MCN(MNITE.  ) = NOMBRE DE MOTS PAR ITEM (=4 ou PLUS)
C     MCN(MNITE.+1) = NOMBRE MAXIMAL D'ITEMS STOCKABLES DANS CE TABLEAU MCN
C     MCN(MNITE.+2) = NOMBRE ACTUEL  D'ITEMS

C     RMCN(MNITE.+3*I  ) = X ABSCISSE DU BARYCENTRE DU POINT ou ... ou EF3D I
C     RMCN(MNITE.+3*I+1) = Y ORDONNEE DU BARYCENTRE DU POINT ou ... ou EF3D I
C     RMCN(MNITE.+3*I+2) = Z COTE     DU BARYCENTRE DU POINT ou ... ou EF3D I
C                        = Cte SI NDIM=2. EN FAIT LA COTE Z EST INACTIVE
C      MCN(MNITE.+3*I+3) = NUMERO DU POINT ou ... ou EF3D I
C                   I VARIANT DE 1 A MCN(MNIT..+2)

C     LA GESTION DES COULEURS A TRACER -----------------------
      COMMON / MECOI2 / LCRITR, LATTRI(MLATTR)
      REAL              RATTRI(1:MLATTR) 
      EQUIVALENCE      (RATTRI( 1),LATTRI(1),NBRCOU),
     %                 (LATTRI( 2),NCOFON),(LATTRI( 3),IAVARE),
     %                 (LATTRI( 4),NPREDA),(LATTRI( 5),NEPARL),
     %                 (LATTRI( 6),NCOUAL),(LATTRI( 7),IAVFAC),
     %                 (LATTRI( 8),NPREDF),(LATTRI( 9),NEPARF),
     %                 (LATTRI(10),NCOUFA),(LATTRI(11),NCOUAF),
     %                 (LATTRI(12),IAVELO),(LATTRI(13),IAVNSO),
     %                 (LATTRI(14),NCONSO),(LATTRI(15),IAVNEF), 
     %                 (LATTRI(16),NCONEF),(LATTRI(17),IAVTIT),
     %                 (LATTRI(18),LPPOIN),(LATTRI(19),LPLIGN),
     %                 (LATTRI(20),LPSURF),(LATTRI(21),LPVOLU),
     %                 (LATTRI(22),LPOBJE),(LATTRI(23),NBSUAR),
     %                 (LATTRI(24),NBSUAF),(LATTRI(25),IAVTGA),
     %                 (LATTRI(26),NCOTGA),(LATTRI(27),NTRTGA),
     %                 (LATTRI(28),IAVTGF),(LATTRI(29),NCOTGF),
     %                 (LATTRI(30),NTRTGF),(LATTRI(31),NCOPOI),
     %                 (LATTRI(32),NCOLIG),(LATTRI(33),NCOSUR),
     %                 (LATTRI(34),NCOVOL),(LATTRI(35),NCOOBJ),
     %                 (LATTRI(36),NTY3PL),(LATTRI(37),NEP3PL),
     %                 (LATTRI(38),NC13PL),(LATTRI(39),NC23PL),
     %                 (LATTRI(40),NC33PL),(LATTRI(41),MAR3PL),
     %                 (LATTRI(42),NCOUFL),(LATTRI(43),NCOUAD),
     %                 (LATTRI(44),NCOAPL),(LATTRI(45),NTLAPL),
     %                 (LATTRI(46),NCOAFR),(LATTRI(47),NTLAFR),
     %                 (LATTRI(48),NTRAXE),(LATTRI(49),NEPFLE),
     %                 (LATTRI(50),NTRAXZ),(LATTRI(51),IAVNRF),
     %                 (LATTRI(52),NCONRF),(LATTRI(53),NTRNRF),
     %                 (LATTRI(54),NPAFLE)

C     LCRITR : CHOIX DU CRITERE A TRACER
C              0 : TRACE SANS la QUALITE DES ELEMENTS FINIS
C              1 : TRACE   de la QUALITE DES ELEMENTS FINIS
C              2... : TEMPERATURE, CONTRAINTES, ... )
C             -2 : PRISE EN COMPTE DES COULEURS AUX SOMMETS avec les
C                  valeurs ZZZMIN ZZZMAX: VALEUR MINIMALE et MAXIMALE
C                  de Z utilisee pour le CALCUL de la COULEUR aux SOMMETS

C     LATTRI : LES ATTRIBUTS DU TRACE DU MAILLAGE DE L'OBJET
C              OBTENUS SOIT PAR DEFAUT SOIT PAR LE TABLEAU ~/TD/d/a___trace 
C     NBRCOU : >0  NOMBRE DE COULEURS DES PALETTES DE COULEURS
C              =0  SI ECRAN EN NOIR ET BLANC
C     NCOFON : COULEUR du FOND de l'ECRAN

C     IAVARE : avec 0 ou non 1  trace des aretes des lignes
C     NPREDA : pourcentage de reduction des ARETES des lignes OBSOLETE (INUTILISE)
C     NEPARL : nombre d'epaisseurs des ARETES des lignes
C     NCOUAL : COULEUR des ARETES des LIGNES

C     IAVFAC : avec ou non trace des faces
C     NPREDF : pourcentage de reduction des faces OBSOLETE (INUTILISE)
C     NEPARF : nombre d'epaisseurs des ARETES des FACES 
C     NCOUFA : COULEUR des FACES  
C     NCOUAF : COULEUR des ARETES des FACES

C     IAVELO : avec ou sans eloignement
C     IAVNSO : avec ou non trace du numero des SOMMETS 
C     NCONSO : COULEUR des NUMEROS des SOMMETS des EF

C     IAVNEF : avec ou non trace du numero des EF
C     NCONEF : COULEUR des NUMEROS des EF

C     IAVTIT : avec ou non trace du TITRE 
C     LPPOIN : LA POIGNEE ou NON DES POINTS
C     LPLIGN : LA POIGNEE ou NON DES LIGNES 
C     LPSURF : LA POIGNEE ou NON DES SURFACES
C     LPVOLU : LA POIGNEE ou NON DES VOLUMES 
C     LPOBJE : LA POIGNEE ou NON DES OBJETS

C     NBSUAR : NOMBRE DE SUBDIVISIONS DES ARETES POLYNOMES DE DEGRE 3
C     NBSUAF : NOMBRE DE SUBDIVISIONS DES ARETES DES FACES POLYNOMES DE DEGRE 3

C     IAVTGA : avec ou non trace des 2 tangentes des aretes degre 3
C     NCOTGA : COULEUR des TANGENTES des ARETES  de DEGRE 3
C     NTRTGA : Numero du TYPE de TRAIT des TANGENTES des ARETES  de DEGRE 3
C     IAVTGF : avec ou non trace des 8 ou 6 tangentes des aretes des faces degre 3
C     NCOTGF : COULEUR des TANGENTES des ARETES des FACES de DEGRE 3
C     NTRTGF : Numero du TYPE de TRAIT des TANGENTES des FACES de DEGRE 3

C     IAVNRF : avec ou non trace de la NORMALE au barycentre des faces
C     NCONRF : COULEUR de la NORMALE au barycentre des faces
C     NTRNRF : Numero du TYPE de TRAIT de la NORMALE au barycentre des faces

C     NCOPOI : COULEUR de trace du NOM des POINTS
C     NCOLIG : COULEUR de trace du NOM des LIGNES
C     NCOSUR : COULEUR de trace du NOM des SURFACES
C     NCOVOL : COULEUR de trace du NOM des VOLUMES
C     NCOOBJ : COULEUR de trace du NOM des OBJETS

C     NTY3PL : TYPE du trace des 3PLANS de FOND des PLSVO
C              =<0 PAS DE TRACE
C              =1  TRACE SIMPLE DE LA GRILLE AVEC LA COULEUR NC13PL
C              =2  TRACE D'UN DAMIER AVEC LES 2 COULEURS NC23PL ET NC33PL
C     MAR3PL : MARGE EN CARREAUX DU DAMIER au dela de l'hexaedre englobant
C     NEP3PL : NOMBRE D'EPAISSEURS de trace de la GRILLE des 3 PLANS
C     NC13PL : COULEUR de trace de la grille des 3 PLANS
C     NC23PL : COULEUR 1 de trace du damier  des 3 PLANS
C     NC33PL : COULEUR 2 de trace du damier  des 3 PLANS

C     NCOUFL : COULEUR des FLECHES a TRACER
C     NEPFLE : NOMBRE D'EPAISSEURS DE TRAIT D'UNE FLECHE
C     NPAFLE : NOMBRE DE PAS DE TRACE DES FLECHES VITESSE (1 sur NPAFLE)
C     NCOUAD : COULEUR des ARETES DEFORMEES

C     NCOAPL : COULEUR des ARETES DANS UN PLAN DE SECTION ou
C              COULEUR des ARETES des SURFACES ISOTHERMES ou ISOVALEURS
C     NTLAPL : NUMERO DU TYPE DE TRACE DES ARETES DE LA FRONTIERE
C              des ARETES des SURFACES ISOTHERMES ou ISOVALEURS
C              0:CONTINU, 1:TIRETE, 2:TIRETE DOUBLE
C
C     NCOAFR : COULEUR des ARETES DU MAILLAGE DE LA FRONTIERE
C     NTLAFR : NUMERO DU TYPE DE TRACE DES ARETES DE LA FRONTIERE
C              0:CONTINU, 1:TIRETE, 2:TIRETE DOUBLE

C     LA GESTION DES AXES A TRACER ----------------------------
C     NTRAXE : 0: PAS DE TRACE DES AXES
C              >0: CODE DE TRACE DES AXES
C     NTRAXZ : 0 PAS DE TRACE SPECIAL DES GRADUATIONS EN Z
C              >0: CODE DE TRACESPECIAL  DES GRADUATIONS EN Z
C     ZMIAXZ : REEL MIN DE LA GRADUATION SPECIALE DE L'AXE Z
C     ZMXAXZ : REEL MAX DE LA GRADUATION SPECIALE DE L'AXE Z
      REAL              ZMIAXZ, ZMXAXZ
      COMMON / MECOI3 / ZMIAXZ, ZMXAXZ

C     LES REDUCTIONS DE TRACE DES FACES ARETES
C     PREDUA : pourcentage de reduction des ARETES des lignes
C     PREDUF : pourcentage de reduction des faces  
      REAL              PREDUA, PREDUF
      COMMON / MECOI4 / PREDUA, PREDUF

C     LA GESTION des COULEURS dans t31fco.f et tritemv012.f
C     VALEURS UTILISEES SI LCRITR=-2
C     ZZZMIN : VALEUR MINIMALE DE Z utilisee pour le CALCUL 
C               de la COULEUR a TRACER
C     ZZZMAX : VALEUR MAXIMALE DE Z utilisee pour le CALCUL 
C               de la COULEUR a TRACER
      REAL              ZZZMIN, ZZZMAX
      COMMON / MECOI5 / ZZZMIN, ZZZMAX
