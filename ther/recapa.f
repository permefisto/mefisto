      SUBROUTINE RECAPA( NYOBJT, NUOBJT, NBCOOR, XYZPI,
     %                   MNMASS, MNCHMA, CAPACT  )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT : CALCULER LA CAPACITE CALORIFIQUE = CHALEUR MASSIQUE Cp x MASSE Rho
C ----- ( ou encore NOMMEE CAPACITE THERMIQUE ) EN UN POINT A L'INSTANT TEMPS

C REMARQUE: CHALEUR MASSIQUE ou CHALEUR SPECIFIQUE notee Cp dans les FORMULES
C           est maintenant NOMMEE  CAPACITE THERMIQUE MASSIQUE
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET 1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE

C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C          EN 1D FOURNIR NBCOOR=3 ET LES VALEURS XYZPI(2:3)=0
C          EN 2D FOURNIR NBCOOR=3 ET LA VALEUR XYZPI(3)=0
C XYZPI  : LES NBCOOR COORDONNEES DU POINT DE CALCUL DE LA CHALEUR MASSIQUE
C MNMASS : ADRESSE MCN DU TABLEAU 'MASSE'
C MNCHMA : ADRESSE MCN DU TABLEAU 'CHALEUR MASSIQUE'

C SORTIE :
C --------
C CAPACT : CAPACITE = Rho DENSITE MASSE(Kg/m3)  *
C                     Cp CAPACITE THERMIQUE MASSIQUE(J/Kelvin/kg)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1999
C....6...............................................................012
      include"./incl/a___chaleurmassique.inc"
      include"./incl/a___masse.inc"
      include"./incl/donthe.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      DOUBLE PRECISION XYZPI(NBCOOR), Cp, CAPACT, DMASSE
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))

      IF( MNCHMA .LE. 0 ) THEN
C        ATTENTION: VALEUR FIXEE A 1 POUR L'EQUATION DES ONDES
C        POUR LAQUELLE ON A LE TERME:  ROT d2 THETA / dT2 - ...
         CAPACT = 1D0
         GOTO 9999
      ENDIF

C     RECHERCHE DE LA CHALEUR MASSIQUE,  Notee Cp dans les FORMULES
C     (DITE maintenant CAPACITE THERMIQUE MASSIQUE en J/Kelvin/kg)
      CALL RECATHMA( NYOBJT, NUOBJT, NBCOOR, XYZPI, MNCHMA, Cp )

C     RECHERCHE DE LA DENSITE DE MASSE EN CE POINT
      CALL REMASS( NYOBJT, NUOBJT, NBCOOR, XYZPI, MNMASS,  DMASSE )

C     LA CAPACITE CALORIFIQUE  RHO * Cp  EN CE POINT
      CAPACT = DMASSE * Cp

 9999 RETURN
      END
