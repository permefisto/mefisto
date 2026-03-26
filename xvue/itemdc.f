      SUBROUTINE ITEMDC
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DECLARE LES 9 TABLEAUX DE GESTION DES ITEMS VISIBLES
C -----    sur la FENETRE GRAPHIQUE lors de CERTAINS TRACES
C          1 POINT  2:LIGNE  3:SURFACE  4:VOLUME  5:OBJET
C          6:SOMMET 7:ARETE  8:FACE EF3D 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain ANALYSE NUMERIQUE UPMC PARIS        JUIN 1994
C MODIFS : PERRONNET Alain Saint PIERRE du PERRAY             Avril 2020
C MODIFS : PERRONNET Alain Saint PIERRE du PERRAY           Janvier 2021
C ......................................................................
      PARAMETER     (MXITEM=1024-1)
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)

      DO I = 1, 9

C        DECLARATION DU TABLEAU DU TYPE I D'ITEM
         GOTO( 10, 10, 10, 10, 10, 60, 70, 80, 10 ), I

C        ITEMS VISIBLES POINT LIGNE SURFACE VOLUME OBJET et EF3D
C        STOCKAGE de 1:ADRESSE MCN, 2:Nb MOTS/ITEM, 3:MXITEMS
 10      NbMOTIT = 4
         GOTO 100

C        ITEMS SOMMET VISIBLE: ADRESSE MCN, Nb MOTS/ITEM, MXITEMS,...
C                              XYZBAR, No Sommet dans XYSOMMET
 60      NbMOTIT = 4
         GOTO 100

C        ITEMS ARETE VISIBLE: ADRESSE MCN, Nb MOTS/ITEM, MXITEMS,...
C                             XYZBAR, No Arete dans LARETE
 70      NbMOTIT = 4
         GOTO 100

C        ITEMS FACE VISIBLE: ADRESSE MCN, Nb MOTS/ITEM, MXITEMS,...
C                            XYZBAR, No FACE dans NSEF
 80      NbMOTIT = 4


C        DECLARATION MCN DU TABLEAU DES ITEMS I
 100     CALL TNMCDC( 'ENTIER', NbMOTIT*(1+MXITEM), MNITEM(I) )

C        NOMBRE DE MOTS PAR ITEM
         MCN( MNITEM(I)     ) = NbMOTIT

C        NOMBRE MAXIMAL D'ITEMS DECLARABLES
         MCN( MNITEM(I) + 1 ) = MXITEM

C        AUCUN ITEM N'EST VISIBLE POUR L'INSTANT
         MCN( MNITEM(I) + 2 ) = 0

      ENDDO

      RETURN
      END
