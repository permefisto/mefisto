      SUBROUTINE TRFINS( KTITRE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LE TITRE D'UN DESSIN
C -----
C ENTREE :
C --------
C KTITRE : LA CHAINE DE CARACTERES A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1994
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/traaxe.inc"
      include"./incl/xvfontes.inc"
      include"./incl/xyzext.inc"
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     KTITRE
      CHARACTER*256     KBLANC

      IF( INTERA .LE. 0 ) RETURN
C
C     SAUVEGARDE DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      NOFONT0 = NOFONT
C
C     TRACE DES 3 ARETES VUES DU TRIEDRE
      CALL T3PLAP
C
C     RETOUR DIRECT SI PAS DE TITRE A TRACER
      IF( IAVTIT .EQ. 0 ) GOTO 1000
C
C     TRACE DE LA 1-ERE LIGNE EN HAUT A GAUCHE DE LA FENETRE DE TRACE
C     AVEC 'MEFISTO', NomUtilisateur, DATE, ...
      CALL TIT1LG
C
C     CHANGEMENT DE POLICE DE CARACTERES POUR UNE DE 20 PIXELS DE HAUT
      LHPXCA = 20
      CALL CHOIXFONTE( LHPXCA )
C
C     UNE LIGNE ENTRE LES 2 EST DISPONIBLE POUR AFFICHER LE TITRE du TRACE
ccc      CALL XVCOULEUR( NCBLAN )
      L = NUDCNB( KTITRE )

C     EFFACEMENT DE LA LIGNE - Phase 9.1 fix 2026-05-06: l'ancien code
C     dessinait L espaces en NCROUG par-dessus l'ancien titre. XVTEXTE
C     n'efface pas les pixels du glyphe precedent (les espaces ne
C     dessinent rien). Resultat: TRFINS appele plusieurs fois avec des
C     titres differents -> tous les titres rouges visibles superposes.
C     Solution: rectangle plein de la couleur de fond avant le trace.
      CALL XVCOULEUR( NCNOIR )
      CALL XVRECTANGLE( 0, 3*LHPXCA-NPHACA-2, LAPXFE, NPHACA+8 )

C     TRACE DU TITRE
      CALL XVCOULEUR( NCROUG )
      CALL XVTEXTE( KTITRE(1:L), L, 20, 3*LHPXCA+9 )
C
C     COPIE DE MEMPX DANS FENETRE
 1000 CALL MEMPXFENETRE
C
C     VIDE LE BUFFER DE X11
      CALL XVVOIR
C
C     RETOUR AU TRACE CONTINU DES TRAITS
      CALL XVTYPETRAIT( 0 )

C     PAS D'EPAISSEUR SUPPLEMENTAIRE DES TRAITS
      CALL XVEPAISSEUR( 0 )
C
C     TRACE DES AXES EFFECTUE
      NETAXE = 0
C
C     RESTAURATION DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      CALL CHARGEFONTE( NOFONT0 )

      RETURN
      END
