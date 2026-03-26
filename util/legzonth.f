      SUBROUTINE LEGZONTH( KNOMOB, NCAS, NOPROJ,  MODECO,
     %                     VMIN,   VMAX, VMINCAS, VMAXCAS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE EFFECTIF avec LA LEGENDE DU TRACE DES COULEURS
C -----                        SELON LE TYPE DE LA SOLUTION
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NCAS   : NUMERO DU CAS A TRACER PARMI LES NDSM CAS
C NOPROJ : TYPE DE PROJECTION 0 CI-DESSOUS FIXE LA COORDONNEE A ZERO
C         <=0 : PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1  CE SONT DES TEMPERATURES
C          =2  CE SONT DES VECTEURS PROPRES
C          =3  CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C              SOIT P2 SOIT P1+BULLE P3
C          =4  CE SONT DES ERREURS PONCTUELLES SOLEX-SOLCAL
C              => TRACE DE L'ERREUR ENTRE TEMPERATURE_EXACTE
C                 ET LA TEMPERATURE CALCULEE
C          =5  C'EST LE PARCOURS DE PARTICULES SUIVANT L'ECOULEMENT DU FLUIDE
C          =6  CE SONT DES NORMES DE VITESSE AUX NOEUDS D'UN MAILLAGE
C          =7  FONCTION COURANT D'UN FLUIDE
C          =8  CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C              DU MODULE D'UNE ONDE COMPLEXE
C             (CALCUL DU MODULE DE L'ERREUR COMPLEXE A FAIRE ICI)
C          =9  PARTIE REELLE     D'UNE ONDE COMPLEXE NLSE
C          =10 PARTIE IMAGINAIRE D'UNE ONDE COMPLEXE NLSE
C          =11 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C              DU MODULE DE LA VITESSE D'UN FLUIDE
C             (CALCUL DE L'ERREUR(VITESSE EXACTE - CALCULEE) DEJA CALCULE)
C          =12 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C              DE LA PRESSION P1 DANS UN FLUIDE
C             (CALCUL DE L'ERREUR(PRESSION EXACTE-CALCULEE) DEJA CALCULE)
C VMIN   : SOLUTION MINIMALE DES NDSM CAS
C VMAX   : SOLUTION MAXIMALE DES NDSM CAS
C VMINCAS: SOLUTION MINIMALE DU  CAS NCAS
C VMAXCAS: SOLUTION MAXIMALE DU  CAS NCAS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
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

C     EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
C     ---------------------------------------
      IF( LASOPS .NE. 0 ) THEN
        IF( LASOPS .EQ. 1 ) THEN
          LASOPS = -11
        ELSE
          IF( LASOPS .EQ. 2 ) THEN
            LASOPS = -12
          ELSE
            LASOPS = 0
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'LEGZONTH : MAUVAISE VALEUR de LASOPS'
               KERR(2) = '           ARRET du TRACE POSTSCRIPT'
            ELSE
               KERR(1) = 'LEGZONTH : BAD VALUE of LASOPS'
               KERR(2) = '           STOP of the POSTSCRIPT DRAWING'
            ENDIF
            CALL LEREUR
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT(LASOPS)
        LASOPS = - LASOPS
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF

C     LE TRACE DE LA LEGENDE: COULEURS => VALEURS
C     -------------------------------------------
      CALL LEGCOULSO( VMIN, VMAX )

C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
C     --------------------------------------
      IF ( LASOPS.NE.0 ) THEN
        LASOPS = LASOPS - 10
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF

C     TRACE DE LA 2-EME LIGNE DU TITRE DU TRACE
C     -----------------------------------------
      CALL TIT2LG( KNOMOB, MODECO )

C     DEFINITION DE LA 3-EME LIGNE DU TITRE
C     -------------------------------------
      CALL TIT3LG( MODECO, NCAS, TEMPS, VMINCAS, VMAXCAS, KNOM )

C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )

      CALL PROJ6C( NOPROJ, KNOM )

C     TRACE DU TITRE
C     --------------
      CALL TRFINS( KNOM )

      RETURN
      END
