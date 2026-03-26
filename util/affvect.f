         SUBROUTINE AFFVECT( TITRE, NB, V )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AFFICHE UN TITRE ET LES NB REEL DOUBLE PRECISION DU VECTEUR V
C -----

C ENTREES:
C --------
C TITRE : LE TITRE A AFFICHER
C NB    : NOMBRE DE COMPOSANTES DU VECTEUR V
C V     : VECTEUR REEL DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2008
C23456---------------------------------------------------------------012
         CHARACTER*(*)     TITRE
         DOUBLE PRECISION  V(NB)

         PRINT *, TITRE
         PRINT 10000,(I,V(I),I=1,NB)
10000    FORMAT(5('  V(',I4,')=',G14.7))

ccc         CALL NUDCNB( TITRE, N )
ccc         PRINT 10010,TITRE(1:N),(I,V(I),I=1,NB)
ccc10010    FORMAT(A,5('  V(',I4,')=',G14.7))

         RETURN
         END
