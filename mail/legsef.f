      SUBROUTINE LEGSEF( NUTYOB, NOMOBJ, NBST, NBEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRACER LE NOM DU LSV, LE NOMBRE DE SOMMETS ET D'EF DU MAILLAGE
C -----
C
C ENTREES :
C ---------
C NUTYOB  : NUMERO DU TYPE DE LSV ( 2:LIGNE, 3:SURFACE, 4:VOLUME )
C NOMOBJ  : NOM DU LSV
C NBST    : NOMBRE DE SOMMETS DU MAILLAGE
C NBEF    : NOMBRE DES ELEMENTS FINIS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1994
C....................................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      CHARACTER*(*)  NOMOBJ
      CHARACTER*34   KNOM
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
C
      IF( INTERA .GE. 1 ) THEN
C
C        LE NOM DU LSV
C        =============
C        LE TYPE DE LSV
         IF( NUTYOB .EQ. 1 ) THEN
             L = 7
             KNOM(1:L) = 'POINT :'
         ELSE IF( NUTYOB .EQ. 2 ) THEN
             L = 7
             KNOM(1:L) = 'LIGNE :'
         ELSE IF( NUTYOB .EQ. 3 ) THEN
             L = 9
             KNOM(1:L) = 'SURFACE :'
         ELSE IF( NUTYOB .EQ. 4 ) THEN
             L = 8
             KNOM(1:L) = 'VOLUME :'
         ELSE
             L = 0
         ENDIF
C        AJOUT DU NOM DU LSV
         LN   = NUDCNB( NOMOBJ )
         KNOM = KNOM(1:L) // ' ' // NOMOBJ(1:LN)
         L    = L + 1 + LN
C        NOMBRE DE PIXELS EN LARGEUR ET HAUTEUR
         CALL XVNBPIXELTEXTE( KNOM(1:L), L, NBPXLA, NBPXHA )
C
C        POSITION INITIALE
         NX = MIN( LAPXFE - 170, LAPXFE - NBPXLA - 10 )
         NY = 750
C
C        EFFACEMENT DE LA LIGNE DU NOM DU LSV
         CALL XVCOULEUR( NCOFON )
         LXBORD = LAPXFE - 10 - NX + 4 * NPLACA
         NX0    = NX - 4 * NPLACA
         CALL XVRECTANGLE( NX0, NY-NPHACA-2, LXBORD, NPHACA+4 )
C
C        TRACE DU NOM DU LSV
         CALL XVCOULEUR( NCBLAN )
         CALL XVTEXTE( KNOM(1:L), L, NX, NY )
         NY = NY + 5
         CALL XVTRAIT( NX, NY, NX+NBPXLA, NY )
C
C        LE NOMBRE DE SOMMETS
         NY = NINT( NY + 2.2 * NPHACA )
         CALL XVCOULEUR( NCOFON )
         CALL XVRECTANGLE( NX0, NY-NPHACA-2, LXBORD, NPHACA+4 )
         CALL XVCOULEUR( NCBLAN )
         CALL XVTEXTE( 'ST:', 3, NX, NY )
         CALL ENTIERPX( NCBLAN, NX+50, NY, NBST )
C
C        LE NOMBRE D'EF
         NY = NY + 2 * NPHACA
C        EFFACEMENT DE LA LIGNE
         CALL XVCOULEUR( NCOFON )
         CALL XVRECTANGLE( NX0, NY-NPHACA-2, LXBORD, NPHACA+4 )
         CALL XVCOULEUR( NCBLAN )
         CALL XVTEXTE( 'EF:', 3, NX, NY )
         CALL ENTIERPX( NCBLAN, NX+50, NY, NBEF )
C
      ENDIF
      END
