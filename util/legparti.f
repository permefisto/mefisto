      SUBROUTINE LEGPARTI( KNOMOB, NBPART, NCAS, VMIN, VMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE EFFECTIF avec LA LEGENDE DU TRACE DU PARCOURS DES
C -----                        PARTICULES
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NBPART : NOMBRE DE PARTICULES A TRACER LE PARCOURS
C NCAS   : NUMERO DU 1-ER CAS DU PARCOURS
C TEMPS  : TEMPS de DEPART des PARCOURS (Cf incl/ctemps.inc)
C VMIN   : VITESSE MINIMALE D'UNE PARTICULE DURANT LES PARCOURS
C VMAX   : VITESSE MAXIMALE D'UNE PARTICULE DURANT LES PARCOURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY         decembre 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
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
               KERR(1) = 'legparti: MAUVAISE VALEUR de LASOPS'
               KERR(2) = '          ARRET du TRACE POSTSCRIPT'
            ELSE
               KERR(1) = 'legparti: BAD VALUE of LASOPS'
               KERR(2) = '          STOP of the POSTSCRIPT DRAWING'
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

C     SAUVEGARDE DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      NOFONT0 = NOFONT

C     CHANGEMENT DE LA POLICE DE CARACTERES
      LHPXCA = 20
      CALL CHOIXFONTE( LHPXCA )

C     LA 2-EME LIGNE DU TITRE DU TRACE
C     --------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'PARCOURS de         PARTICULES dans le FLUIDE: '
     %           // KNOMOB
         WRITE( KNOM(13:18), '(I6)' ) NBPART
      ELSE
         KNOM = '         PARTICLE TRIPS in FLUID: ' // KNOMOB
         WRITE( KNOM(1:6), '(I6)' ) NBPART
      ENDIF
      CALL XVCOULEUR( NCBLEU )

C     SUPPRESSION DES DOUBLE BLANCS
      CALL SANSDBL( KNOM, L )

C     AFFICHAGE DU TEXTE
      CALL XVTEXTE( KNOM(1:L), L, 20, 2*LHPXCA+6 )

C     LA 3-EME LIGNE DU TITRE
C     -----------------------
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'PARCOURS a partir du CAS        au TEMPS              '
         WRITE( KNOM(26:30), '(I5)'    ) NCAS
         WRITE( KNOM(42:55), '(G14.6)' ) TEMPS
      ELSE
         KNOM = 'TRIP from CASE        at TIME                         '
         WRITE( KNOM(16:20), '(I5)'    ) NCAS
         WRITE( KNOM(31:44), '(G14.6)' ) TEMPS
      ENDIF

C     FIN de la3-EME LIGNE
      I = NUDCNB( KNOM )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM(I+1:I+14) = ' Vitesse MAX= '
         WRITE( KNOM(I+15:I+28), '(G14.6)' ) VMAX
      ELSE
         KNOM(I+1:I+15) = ' Velocity MAX= '
         WRITE( KNOM(I+16:I+29), '(G14.6)' ) VMAX
      ENDIF

C     SUPPRESSION DES DOUBLE BLANCS
      CALL SANSDBL( KNOM, L )

C     TRACE DU TITRE EN 3-EME LIGNE EN HAUT A GAUCHE
C     ----------------------------------------------
      CALL TRFINS( KNOM )

C     RESTAURATION DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      CALL CHARGEFONTE( NOFONT0 )

      RETURN
      END
