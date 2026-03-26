      SUBROUTINE NMGRAF( KNOM , KTANOM , NO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFECTER LE NOM KNOM EN POSITION NO DU TABLEAU KTANOM
C -----
C ENTREES :
C ---------
C KNOM   : NOM A AFFECTER
C NO     : NO DANS LE TABLEAU KTANOM ( 0=< NO =< MXNOMS )
C ENTREE ET SORTIE :
C ------------------
C KTANOM : TABLEAUX DE NOMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE PARIS  OCTOBRE 1984
C.......................................................................
      CHARACTER*(*) KNOM
      CHARACTER*24  KTANOM(0:*)
C
      KTANOM( NO ) = KNOM
      END
