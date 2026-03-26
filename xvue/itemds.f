      SUBROUTINE ITEMDS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DESTRUCTION DES TABLEAUX DE GESTION DES ITEMS
C -----
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain Saint PIERRE du PERRAY             Avril 2020
C ......................................................................
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)

      DO I = 1, 9

C        NOMBRE DE MOTS PAR ITEM =3
         NBMOIT = MCN( MNITEM(I) )

C        NOMBRE MAXIMAL D'ITEMS DECLARES
         MXITEM = MCN( MNITEM(I) + 1 )

C        DESTRUCTION DU TABLEAU DE CE TYPE D'ITEM
         CALL TNMCDS( 'ENTIER', NBMOIT*MXITEM+NBMOIT, MNITEM(I) )

      ENDDO

      RETURN
      END
