      SUBROUTINE MXFLSD( NTLXOB, NBTYEL, MNELEM, FLUXMX, SDCOIN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE FLUX DE CHALEUR MAXIMAL ET
C -----    L'HEXAEDRE ENGLOBANT TOUS LES POINTS (SOUS-DOMAINES)
C
C ENTREES :
C ---------
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET A TRAITER
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNELEM : ADRESSE MCN DU TABLEAU DES ADRESSES MCN DES TABLEAUX ELEMENTS
C
C RESULTATS :
C -----------
C SDCOIN : MIN ET MAX DES COORDONNEES DES POINTS OU SONT CALCULEES
C          LES GRADIENTS
C FLUXMX : VALEUR ABSOLUE DU FLUX MAXIMAL
C IERR   : 0 SI PAS D'ERREUR >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              FLUXMX, SDCOIN(6,2)
      CHARACTER*160     KNOM
      CHARACTER*4       NOMELE(2)
C
C     INITIALISATIONS
      RMAX   = RINFO( 'GRAND' )
      DO 10 I=1,3
C        LE MINIMUM
         SDCOIN(I,1) =  RMAX
C        LE MAXIMUM
         SDCOIN(I,2) = -RMAX
 10   CONTINUE
      FLUXMX = 0.
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS DU MAILLAGE
      DO 20 I = 0 , NBTYEL-1
C        L'ADRESSE MCN DU DEBUT DU TABLEAU ELEMENTS
         MNELE = MCN( MNELEM + I )
C        LE NOM DU TABLEAU FLUX ASSOCIE
         NUTYEL = MCN( MNELE + WUTYEL )
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL , NOMELE )
         KNOM = 'FLUXPT"' // NOMELE(2)
C        OUVERTURE DU TABLEAU
         CALL LXTSOU( NTLXOB , KNOM , NTFLUX , MNFLUX )
         IF( NTFLUX .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'OBJET SANS TABLEAU FLUXPT'
            CALL LEREUR
            IERR = 1
            GOTO 20
         ENDIF
C
C        RECHERCHE DU MIN MAX DES COORDONNEES ET DES FLUX NORMAUX
         NBPNFX = MCN( MNFLUX + WBPNFX )
         NBELFX = MCN( MNFLUX + WBELFX )
         NBCAFX = MCN( MNFLUX + WBCAFX )
         MNFLNP = MNFLUX + WLUXNP
         MNCOPN = MNFLNP + MOTVAR(6) * NBPNFX * NBELFX * NBCAFX
         CALL MXFLU1( NBPNFX, NBELFX, NBCAFX,
     %                MCN(MNCOPN), MCN(MNFLNP),
     %                SDCOIN, FLUXMX )
 20   CONTINUE
C
      IF( SDCOIN(3,1) .EQ. RMAX ) THEN
C        DIMENSION 2
         SDCOIN(3,1) = 0
         SDCOIN(3,2) = 0
         NDIMLI = 2
      ELSE
         NDIMLI = 3
      ENDIF
C
      END
