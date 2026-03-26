      SUBROUTINE LEGVIT6( KNOMOB, NCAS, VITMOY, VITMIN, VITMAX, CMVITE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA LEGENDE DE LA NORME DE LA VITESSE pour BF6 ou TH6
C------
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NCAS   : NUMERO DU CAS A TRACER PARMI LES NDSM CAS
C VITMOY : SOLUTION MOYENNE  DES NDSM CAS
C VITMIN : SOLUTION MINIMALE DES NDSM CAS
C VITMAX : SOLUTION MAXIMALE DES NDSM CAS
C CMVITE : CM par UNITE de MODULE de la VITESSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint Pierre du Perray          Fevrier 2021
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
      include"./incl/ctemps.inc"
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION  VITMIN, VITMAX, VITMOY
      REAL              CMVITE
      CHARACTER*(*)     KNOMOB
      CHARACTER*128     KNOM
      CHARACTER*16      TEXT
      INTRINSIC         REAL
C
C     TRACE DE LA LEGENDE DES QUALITES
C     --------------------------------
      IF( INTERA .LE. 0 ) RETURN
C
C     EFFACEMENT DE LA LEGENDE DE LA QUALITE SUR POSTSCRIPT
      IF ( LASOPS.NE.0 ) THEN
         IF ( LASOPS.EQ.1 ) THEN
            LASOPS = -11
         ELSE
            IF ( LASOPS.EQ.2 ) THEN
               LASOPS = -12
            ELSE
               LASOPS = 0
               NBLGRC(NRERR) = 2
               KERR(1) = 'LEGVIT: MAUVAISE VALEUR DE LASOPS'
               KERR(2) = '        ARRET DU POSTSCRIPT'
               CALL LEREUR
            ENDIF
         ENDIF
         CALL XVPOSTSCRIPT(LASOPS)
         LASOPS = - LASOPS
         CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
C     TRACE DE LA LEGENDE DE LA NORME DE LA VITESSE
C     ---------------------------------------------
      CALL XVCOULEUR( NCBLEU )
      NOFONT0 = NOFONT
      LHPXCA = 20
      CALL CHOIXFONTE( LHPXCA )

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'Vecteurs VITESSES dans ' // KNOMOB
      ELSE
         KNOM = 'VELOCITY VECTORS in ' // KNOMOB
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 20, 45 )
C
C     LE TRACE DE LA LEGENDE : COULEURS => VALEURS
      VMIN = REAL( VITMIN )
      VMAX = REAL( VITMAX )
      CALL LEGCOULSO( VMIN, VMAX )
C
C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
         LASOPS = LASOPS - 10
         CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     FIN DU TRACE
      IAVTIT = 1
C     FIN DU TRACE
      WRITE( TEXT, '(I4)' ) NCAS
      L = NUDCNB( TEXT )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'Cas ' // TEXT(1:L)
      ELSE
         KNOM = 'Case '// TEXT(1:L)
      ENDIF
      M = NUDCNB( KNOM )
C
C     LE TEMPS
      WRITE( TEXT, '(G14.6)' ) TEMPS
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = KNOM(1:M) // '  VITESSES m/s au TEMPS=' // TEXT(1:L)
      ELSE
         KNOM = KNOM(1:M) // '  VELOCITIES m/s at TIME=' // TEXT(1:L)
      ENDIF
      M = NUDCNB( KNOM )

C     LA VITESSE MOYENNE
      WRITE( TEXT, '(G14.6)' ) VITMOY
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = KNOM(1:M) // '  VITESSE MOYENNE=' // TEXT(1:L)
      ELSE
         KNOM = KNOM(1:M) // '  MEAN VELOCITY=' // TEXT(1:L)
      ENDIF
      M = NUDCNB( KNOM )

C     LA VITESSE MAX
      WRITE( TEXT, '(G14.6)' ) VITMAX
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = KNOM(1:M) // '  VITESSE MAX=' // TEXT(1:L)
      ELSE
         KNOM = KNOM(1:M) // '  MAX VELOCITY=' // TEXT(1:L)
      ENDIF
      M = NUDCNB( KNOM )

C     CMVITE : NOMBRE DE CM PAR UNITE DE VITESSE
      WRITE( TEXT, '(G14.6)' ) CMVITE
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = KNOM(1:M) //' (fleche 1 m/s='//TEXT(1:L)
     %                    //' CM)'
      ELSE
         KNOM = KNOM(1:M) //'  (1 m/s arrow='//TEXT(1:L)
     %                    //' CM)'
      ENDIF
      M = NUDCNB( KNOM )
C
C     TRACE DANS LA FENETRE
      CALL TRFINS( KNOM(1:M) )

C     REMISE EN ETAT DE LA FONTE COURANTE
      CALL CHARGEFONTE( NOFONT0 )
C
      RETURN
      END
