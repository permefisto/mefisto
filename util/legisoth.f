      SUBROUTINE LEGISOTH( KNOMOB, NCAS, TMIN, TMAX, NBISO, MNVISO,
     %                     MODECO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LA LEGENDE DU TRACE DES ISOTHERMES
C -----
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NCAS   : NUMERO DU CAS A TRACER PARMI LES CAS
C TMIN   : VALEUR MINIMALE DE LA SOLUTION NCAS
C TMAX   : VALEUR MAXIMALE DE LA SOLUTION NCAS
C NBISO  : NOMBRE D'ISOTHERMES A TRACER
C MNVISO : ADRESSE MCN DES VALEURS DES ISO-TEMPERATURES
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT DES TEMPERATURES
C          =2 CE SONT DES MODES PROPRES
C          =3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT en 2D P2 SOIT P1+BULLE P3 OU en 3D P1 OU en 3D P2
C          =4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C          =5 CAS PARTICULIER DU 1D TRACE DES VECTEURS TEMPERATURES
C          =6 CE SONT LES NORMES D'UNE VITESSE
C          =7 FONCTION COURANT D'UN FLUIDE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/ctemps.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*120     KNOM
C
C     EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
        IF ( LASOPS .EQ. 1 ) THEN
          LASOPS = -11
        ELSE
          IF ( LASOPS .EQ. 2 ) THEN
            LASOPS = -12
          ELSE
            LASOPS = 0
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LEGISOTH : MAUVAISE VALEUR de LASOPS'
               KERR(2) = '           ARRET du TRACE POSTSCRIPT'
            ELSE
               KERR(1) = 'LEGISOTH : BAD VALUE of LASOPS'
               KERR(2) = '           STOP of the POSTSCRIPT DRAWING'
            ENDIF
            CALL LEREUR
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT(LASOPS)
        LASOPS = - LASOPS
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
C     AFFICHAGE DE LA LEGENDE D'AU MAXIMUM 20 ISOS
      CALL LEGCOULIS( NBISO,  RMCN(MNVISO) )
C
C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
      IF ( LASOPS.NE.0 ) THEN
        LASOPS = LASOPS -10
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
C     TRACE DE LA 2-EME LIGNE DU TITRE DU TRACE
C     -----------------------------------------
      CALL TIT2LG( KNOMOB, MODECO )
C
C     DEFINITION DE LA 3-EME LIGNE DU TITRE
C     -------------------------------------
      CALL TIT3LG( MODECO, NCAS, TEMPS, TMIN, TMAX, KNOM )
C
C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
C
C     TRACE DU TITRE
C     --------------
      CALL TRFINS( KNOM )
C
      RETURN
      END
