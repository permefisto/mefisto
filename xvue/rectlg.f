      SUBROUTINE RECTLG( NUR , NETEXT , NATEXT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DES LIGNES ENTRE 2 LIGNES DE TEXTE
C -----    VERSION xvue
C
C ENTREE :
C --------
C NUR    : LE NUMERO DU RECTANGLE A TRACER
C NETEXT : 0 PAS DE TRACE DES INTERLIGNES
C          1 TRACE DES INTERLIGNES SANS PREOCCUPATION DE NATEXT
C          2 TRACE DES INTERLIGNES SI NA ENTRE CES 2 LIGNES
C                                        SONT DIFFERENTS
C NATEXT : LE NUMERO ASSOCIE A CHAQUE LIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       MARS 1990
C23456---------------------------------------------------------------012
      include"./incl/epombr.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      INTEGER        NATEXT(1:*)
C
      IF( INTERA .LT. 3 ) RETURN
      IF( NETEXT .LE. 0 ) RETURN
C
C     TRACE DES LIGNES EN CONTINU
      CALL XVTYPETRAIT( 0 )
C
      NX1 = XRECT( NUR )
      NX2 = NX1 + DXRECT( NUR )
      NY1 = YRECT( NUR )
C
      IF( NETEXT .EQ. 1 ) THEN
C
         DO 10 I = 1,NBLGRC( NUR ) - 1
            NY1 = NY1 + DYLGRC( NUR )
C           TRACE SI LES OPTIONS ENTRE 2 LIGNES DIFFERENT
C           LES EPAISSEURS DE LA LIGNE
            CALL XVEPAISSEUR( NBEPO2 )
C           LA COULEUR DE LA LIGNE BORD DU RECTANGLE
            CALL XVCOULEUR( NBRECT( NUR ) )
            CALL XVFTRAIT( NX1, NY1, NX2, NY1 )
C
C           AU DESSUS UN TRAIT BLANC
            CALL XVEPAISSEUR( 1 )
            CALL XVCOULEUR( NCBLAN )
            CALL XVFTRAIT( NX1+NBEPO4, NY1-NBEPO4,
     %                     NX2-NBEPO4, NY1-NBEPO4 )
C           AU DESSOUS UN TRAIT NOIR
            CALL XVCOULEUR( NCNOIR )
            CALL XVFTRAIT( NX1+NBEPO4, NY1+NBEPO4,
     %                     NX2-NBEPO4, NY1+NBEPO4 )
 10      CONTINUE
C
      ELSE
C
C        LE TRACE DES INTERLIGNES LIGNES SOUS CONDITION
         DO 100 I = 1,NBLGRC( NUR ) - 1
            NY1 = NY1 + DYLGRC( NUR )
C           TRACE SI LES OPTIONS ENTRE 2 LIGNES DIFFERENT
            IF( NATEXT(I) .NE. NATEXT(I+1) ) THEN
C              LES EPAISSEURS DE LA LIGNE
               CALL XVEPAISSEUR( NBEPO2 )
C              LA COULEUR DE LA LIGNE BORD DU RECTANGLE
               CALL XVCOULEUR( NBRECT( NUR ) )
               CALL XVFTRAIT( NX1, NY1, NX2, NY1 )
C
C              AU DESSUS UN TRAIT BLANC
               CALL XVEPAISSEUR( 1 )
               CALL XVCOULEUR( NCBLAN )
               CALL XVFTRAIT( NX1+NBEPO4, NY1-NBEPO4,
     %                        NX2-NBEPO4, NY1-NBEPO4 )
C              AU DESSOUS UN TRAIT NOIR
               CALL XVCOULEUR( NCNOIR )
               CALL XVFTRAIT( NX1+NBEPO4, NY1+NBEPO4,
     %                        NX2-NBEPO4, NY1+NBEPO4 )
            ENDIF
 100     CONTINUE
      ENDIF
      END
