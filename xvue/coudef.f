      SUBROUTINE COUDEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISER LES COULEURS PAR DEFAUT DES ARETES , FACES , ...
C -----                 VERSION xvue
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / T3PLAN / XYZTRI(3,2), H(3), NDIV(3), NOFA3P(3), NOET3P
C
C     A PRIORI PAS DE VIDEO DEMANDEE
      LVIDEO = 0
C
C     OPTIONS PAR DEFAUT POUR LES TRACES
C     ----------------------------------
C     INCLINAISON DES TEXTES DE 0.0 DEGRES
      DRGCAR = 0.0
C     IAVARE 'avec ou non trace des aretes'
      IAVARE = 1
C     IAVFAC 'avec ou non trace des faces'
      IAVFAC = 1
C     PREDUF 'pourcentage de reduction des aretes'
      PREDUA = 0.0
C     NEPARL 'nombre d'epaisseurs des aretes des lignes'
      NEPARL = 0
C     NEPARF 'nombre d'epaisseurs des aretes des faces'
      NEPARF = 0
C     PREDUF 'pourcentage de reduction des faces'
      PREDUF = 0.0
C     IAVELO 'avec ou sans eloignement'
      IAVELO = 1
C     IAVNSO 'avec ou non trace du numero des SOMMETS'
      IAVNSO = 0
C     IAVNEF 'avec ou non trace du numero des EF'
      IAVNEF = 0
C
C     NCOFON  'COULEUR du FOND'
      NCOFON = NCFOND()
C
C     NCOUAL  'COULEUR des ARETES DES LIGNES'
CCC      NCOUAL = NCGRIC
CCC      NCOUAL = NCBLAN
      NCOUAL = NCNOIR

C     NCOUAF  'COULEUR des ARETES DES FACES'
C     ATTENTION POUR LE POSTSCRIPT LE BLANC DEVIENT NOIR ET NOIR DEVIENT BLANC
CCC      NCOUAF = NCBLAN
CCC      NCOUAF = NCNOIR
      NCOUAF = NCGRIM

C     NCOUFA  'COULEUR des FACES'
CCC      NCOUFA = NCBLAN
      NCOUFA = NCGRIC
C
C     IAVTIT 'avec ou non trace du titre'
      IAVTIT = 1
C
C     PALETTE DES COULEURS DU TRACE DU MAILLAGE
      NOPACL = 0
C
C     TRACE DU MAILLAGE SELON SA QUALITE
      LCRITR = 1
C     PALETTE DES COULEURS DE LA QUALITE DU MAILLAGE
      CALL PALCDE( 12 )
C
C     TRACE DEMANDE DE LA POIGNEE DES POINTS, LIGNES, ...
      LPPOIN = 1
      LPLIGN = 1
      LPSURF = 1
      LPVOLU = 1
      LPOBJE = 1
C
C     IAVTGA : avec ou non trace des 2 tangentes des aretes degre 3
      IAVTGA = 0
C     NCOTGA : COULEUR des TANGENTES des ARETES de DEGRE 3
      NCOTGA = NCMAGE
C     IAVTGF : avec ou non trace 8 ou 6 tangentes des aretes des faces degre 3
      IAVTGF = 0
C     NCOTGF : COULEUR des TANGENTES des ARETES des FACES de DEGRE 3
      NCOTGF = NCBLAN
C
C     IAVNRF : avec ou non trace du VECTEUR NORMAL au barycentre des faces
      IAVNRF = 0
C     NCONRF : COULEUR des TANGENTES des ARETES des FACES de DEGRE 3
      NCONRF = NCBLAN
C
C     NBSUAR : MAXimum DE SUBDIVISIONS DES ARETES POLYNOMES DE DEGRE 3
C     NBSUAF : MAXimum DE SUBDIVISIONS DES ARETES DES FACES POLYNOMES DE DEGRE 3
      NBSUAR = 16
      NBSUAF = 16
C
C     NTRTGA : Numero du TYPE de TRAIT des TANGENTES des ARETES de DEGRE 3
C     NTRTGF : Numero du TYPE de TRAIT des TANGENTES des FACES  de DEGRE 3
      NTRTGA = 0
      NTRTGF = 0
C     NTRTGF : Numero du TYPE de TRAIT du VECTEUR NORMAL aux FACES
      NTRNRF = 0
C
C     COULEUR DES NUMEROS DES EF ET DES NUMEROS DES SOMMETS
      NCONEF = NCORAN
      NCONSO = NCVERT
C
C     NOMBRE DE COULEURS DES PALETTES DANS ~/incl/mecoit.inc
      NBRCOU = NDCOUL - N1COUL
C
C     NCOPOI : COULEUR de trace du NOM des POINTS
C     NCOLIG : COULEUR de trace du NOM des LIGNES
C     NCOSUR : COULEUR de trace du NOM des SURFACES
C     NCOVOL : COULEUR de trace du NOM des VOLUMES
C     NCOOBJ : COULEUR de trace du NOM des OBJETS
      NCOPOI = NCNOIR
      NCOLIG = NCORAN
      NCOSUR = NCSAUM
      NCOVOL = NCMAGE
      NCOOBJ = NCNOIR
C
C     NTY3PL : TYPE du trace des 3PLANS de FOND des PLSVO
C              =<0 PAS DE TRACE
C              =1  TRACE SIMPLE DE LA GRILLE AVEC LA COULEUR NC13PL
C              =2  TRACE D'UN DAMIER AVEC LES 2 COULEURS NC23PL ET NC33PL
C     MAR3PL : MARGE en CARREAUX DU DAMIER 0 a 10
C     NEP3PL : NOMBRE D'EPAISSEURS de trace de la GRILLE des 3 PLANS
C     NC13PL : COULEUR de trace de la grille des 3 PLANS
C     NC23PL : COULEUR 1 de trace du damier  des 3 PLANS
C     NC33PL : COULEUR 2 de trace du damier  des 3 PLANS
      NTY3PL = 0
      MAR3PL = 0
      NEP3PL = 2
      NC13PL = NCROUG
      NC23PL = NCBLAN
      NC33PL = NCNOIR
C     LES TABLEAUX DU COMMON/T3PLAN/ NE SONT PAS INITIALISES
      NOET3P = 0
      NDIV(1) = 0
C
C     PAS DE LAMPE PAR DEFAUT
      NBLAMP = 0
C
C     COULEUR DES FLECHES
      NCOUFL = NCORAN
C     NOMBRE D'EPAISSEURS DE TRAIT D'UNE FLECHE
      NEPFLE = 1
C     PAS DE TRACE DES FLECHES VITESSE ( 1 sur NPAFLE, 1=>TOUTES)
      NPAFLE = 1
C
C     COULEUR DES ARETES DEFORMEES
      NCOUAD = NCGRIS
C
C     COULEUR des ARETES DANS UN PLAN DE SECTION: INVISIBLE ou GRIS CLAIR
ccc      NCOAPL = NCGRIC
      NCOAPL = -2
C
C     COULEUR des ARETES DU MAILLAGE DE LA FRONTIERE
      NCOAFR = NCGRIC
C     NTLAFR NUMERO DU TYPE DE TRACE DES ARETES DE LA FRONTIERE
C            0: CONTINU, 1:TIRETE, 2: TIRETE DOUBLE
      NTLAFR = 0
C
C     TRACE DES AXES en 3D: X Y //X //Y Z, en 2D: X Y //X //Y
      NTRAXE = 3
C     NTRAXZ : 0 PAS DE TRACE SPECIAL DES GRADUATIONS EN Z
      NTRAXZ = 0
      ZMIAXZ = 0
      ZMXAXZ = 0
C
C     VALEURS UTILISEES SI LCRITR=-2 et DANS le sp mail/t31fco.f
      ZZZMIN = 0
      ZZZMAX = 0
C
      RETURN
      END
