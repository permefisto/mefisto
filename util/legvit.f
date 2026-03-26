      SUBROUTINE LEGVIT( KNOMOB, NCAS, VITMOY, VITMIN, VITMAX, CMVITE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA LEGENDE DE LA NORME DE LA VITESSE
C------
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C NCAS   : NUMERO DU CAS A TRACER PARMI LES NDSM CAS
C VITMOY : SOLUTION MOYENNE  DES NDSM CAS
C VITMIN : SOLUTION MINIMALE DES NDSM CAS
C VITMAX : SOLUTION MAXIMALE DES NDSM CAS
C CMVITE : CM par UNITE de MODULE de la VITESSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Juillet 2010
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

C     TRACE DE LA LEGENDE DES QUALITES
C     --------------------------------
      IF( INTERA .LE. 0 ) RETURN

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

C     TRACE DE LA LEGENDE DE LA NORME DE LA VITESSE
C     ---------------------------------------------
      CALL XVCOULEUR( NCBLEU )
      NOFONT0 = NOFONT
      LHPXCA = 20
      CALL CHOIXFONTE( LHPXCA )

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'MODULE de la VITESSE de ' // KNOMOB
      ELSE
         KNOM = 'VELOCITY MAGNITUDE of ' // KNOMOB
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVTEXTE( KNOM(1:I), I, 20, 45 )

C     LE TRACE DE LA LEGENDE : COULEURS => VALEURS
      VMIN = REAL( VITMIN )
      VMAX = REAL( VITMAX )
      CALL LEGCOULSO( VMIN, VMAX )

C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
         LASOPS = LASOPS - 10
         CALL XVPOSTSCRIPT( LASOPS )
      ENDIF

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

C     LE TEMPS
      WRITE( TEXT, '(G14.6)' ) TEMPS
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = KNOM(1:M) // '  VITESSES au TEMPS=' // TEXT(1:L)
      ELSE
         KNOM = KNOM(1:M) // '  VELOCITIES at TIME=' // TEXT(1:L)
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

C     CMVITE : NOMBRE DE CM PAR UNITE DE VITESSE
      WRITE( TEXT, '(G14.6)' ) CMVITE
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )

      IF( LANGAG .EQ. 0 ) THEN
         KNOM = KNOM(1:M) //'  '//TEXT(1:L)//' CM par UNITE de VITESSE'
      ELSE
         KNOM = KNOM(1:M) //'  '//TEXT(1:L)//' CM per ONE VELOCITY UNIT'
      ENDIF
      M = NUDCNB( KNOM )

C     TRACE DANS LA FENETRE
      CALL TRFINS( KNOM(1:M) )

C     REMISE EN ETAT DE LA FONTE COURANTE
      CALL CHARGEFONTE( NOFONT0 )

      RETURN
      END
