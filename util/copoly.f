      SUBROUTINE COPOLY(NOM,IA,L,NDIM,NDEGRE,NPOUQ,NBPOLY,POLYNO)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : COPIER UN TABLEAU DES COEFFICIENTS DES POLYNOMES DANS UN TABLEAU
C ----- GERE DANS MCN ET LUI AJOUTER EN TETE 4 INFORMATIONS

C PARAMETRES D ENTREE :
C ---------------------
C NOM    : NOM DU TABLEAU DANS LE SUPER TABLEAU M
C NDIM   : NOMBRE DE VARIABLES DES POLYNOMES(1 OU 2 OU 3)
C NDEGRE : DEGRE DES POLYNOMES
C NPOUQ  : 0 => POLYNOME EN (X,Y,Z)
C          1 => POLYNOME EN (X) , (Y) , (Z)
C          2 => POLYNOME EN (X,Y) , (Z)
C NBPOLY : NOMBRE DE POLYNOMES
C POLYNO : LES COEFFICIENTS(NDEGRE+1,NDEGRE+1,NBPOLY) DES POLYNOMES

C PARAMETRE RESULTAT :
C --------------------
C IA     : ADRESSE DANS M DU 1-ER MOT DU TABLEAU COPIE
C L      : NOMBRE DE MOTS DU TABLEAU APRES COPIE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN 189 PARIS ET INRIA  SEPTEMBRE 1981
C ......................................................................
      DOUBLE PRECISION POLYNO(*)
       include"./incl/pppoba.inc"
      COMMON     MCN( MOTMCN )

      PRINT*,'Execution de copoly sur le tableau ',NOM

C     DECLARATION DU TABLEAU DANS LE SUPER-TABLEAU
C     --------------------------------------------
      IA = 0
      L  = 2 * (NDEGRE+1) ** NDIM * NBPOLY + 4
      CALL TNMCDC( 'REEL2' , L , IA )

C     ADJONCTION DES RENSEIGNEMENTS SUPPLEMENTAIRES
C     ---------------------------------------------
      MCN( IA ) = NDIM
      MCN(IA+1) = NDEGRE + 1
      MCN(IA+2) = NPOUQ
      MCN(IA+3) = NBPOLY

C     COPIE DU TABLEAU POLYNO DANS MCN(IA +....)
C     ----------------------------------------
      CALL TRTATA( POLYNO , MCN(IA + 4) , L - 4 )

      RETURN
      END
