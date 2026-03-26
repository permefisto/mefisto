      SUBROUTINE TIT2LG( KNOMOB, MODECO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LA 2-EME LIGNE DU TITRE  DU TRACE SELON MODECO
C -----    CHOISIR LES COULEURS ET LA FONTE DU TRACE
C          TRACER CETTE 2-EME LIGNE DU TITRE
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET EN COURS DE TRAITEMENT
C MODECO = 1  CE SONT DES TEMPERATURES
C        = 2  CE SONT DES VECTEURS PROPRES
C        = 3  CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT P2 SOIT P1+BULLE P3
C        = 4  CE SONT DES ERREURS PONCTUELLES SOLEX-SOLCAL
C             => TRACE DE L'ERREUR ENTRE TEMPERATURE_EXACTE
C                ET LA TEMPERATURE CALCULEE
C        = 5  C'EST LE SUIVI DE PARTICULES SUIVANT L'ECOULEMENT DU FLUIDE
C        = 6  CE SONT DES MODULES DE VITESSE AUX NOEUDS D'UN MAILLAGE
C        = 7  FONCTION COURANT D'UN FLUIDE
C        = 8  CE SONT LES MODULES D'UNE ONDE COMPLEXE AUX NOEUDS DU MAILLAGE
C        = 9  CE SONT DES PARTIES REELLES     D'UNE ONDE COMPLEXE
C        =10  CE SONT DES PARTIES IMAGINAIRES D'UNE ONDE COMPLEXE
C        =11  CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE DE LA VITESSE D'UN FLUIDE
C        =12  CE SONT DES ERREURS DE LA PRESSION AUX NOEUDS D'UN MAILLAGE
C
C SORTIE :
C --------
C KTXT2L : LE TEXTE DE LA 2-EME LIGNE DU TITRE A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      INTEGER         LHPXCA
      PARAMETER      (LHPXCA=20)
C
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
      CHARACTER*(*)  KNOMOB
      CHARACTER*120  KTXT2L
C
C     LE NOM DE L'OBJET
      IF( MODECO .EQ. 12 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L='ERREUR(PRESSION(t,X)) de l''OBJET: '// KNOMOB
         ELSE
            KTXT2L='ERROR(PRESSURE(t,X)) on the OBJECT: '//KNOMOB
         ENDIF
         CALL XVCOULEUR( NCROUG )
C
      ELSE IF( MODECO .EQ. 11 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L='ERREUR(|VITESSE(t,X)|) de l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L='ERROR(|VELOCITY(t,X)|) on the OBJECT: '//KNOMOB
         ENDIF
         CALL XVCOULEUR( NCROUG )
C
      ELSE IF( MODECO .EQ. 10 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'Partie Imaginaire de l''ONDE dans l''OBJET: '
     %               // KNOMOB
         ELSE
            KTXT2L = 'Wave Imaginary Part in the OBJECT: '// KNOMOB
         ENDIF
         CALL XVCOULEUR( NCROUG )
C
      ELSE IF( MODECO .EQ. 9 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'Partie Reelle de l''ONDE dans l''OBJET: ' //KNOMOB
         ELSE
            KTXT2L = 'Wave Real Part in the OBJECT: '// KNOMOB
         ENDIF
         CALL XVCOULEUR( NCVERT )
C
      ELSE IF( MODECO .EQ. 8 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = '|ONDE(t,xyz)| de l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = '|WAVE(t,xyz)| on the OBJECT: '// KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE IF( MODECO .EQ. 7 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'FONCTION COURANT dans l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = 'STREAM FUNCTION in the OBJECT: ' // KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE IF( MODECO .EQ. 6 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = '|VITESSE(t,X)| de l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = '|VELOCITY(t,X)| on the OBJECT: ' // KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE IF( MODECO .EQ. 5 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'SUIVI de PARTICULES dans l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = 'PARTICLES RUN in the OBJECT: ' // KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE IF( MODECO .EQ. 4 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'ERREUR sur la SOLUTION de l''OBJET: '// KNOMOB
         ELSE
            KTXT2L = 'SOLUTION ERROR of the OBJECT: '// KNOMOB
         ENDIF
         CALL XVCOULEUR( NCROUG )
C
      ELSE IF( MODECO .EQ. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'PRESSION(t,X) de l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = 'PRESSURE(t,X) on the OBJECT: ' // KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE IF( MODECO .EQ. 2 ) THEN
         IF( LANGAG .EQ. 2 ) THEN
            KTXT2L = 'VECTEURPROPRE(X) de l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = 'EIGENVECTOR(X) on the OBJECT: ' // KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE IF( MODECO .EQ. 1 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT2L = 'TEMPERATURE(t,X) de l''OBJET: ' // KNOMOB
         ELSE
            KTXT2L = 'TEMPERATURE(t,X) on the OBJECT: ' // KNOMOB
         ENDIF
         CALL XVCOULEUR( NCBLEU )
C
      ELSE
         KTXT2L = KNOMOB
      ENDIF
C
C     SAUVEGARDE DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      NOFONT0 = NOFONT
C
C     CHANGEMENT DE POLICE DE CARACTERES
      CALL CHOIXFONTE( LHPXCA )
C
C     AFFICHAGE DU TEXTE
      L = NUDCNB( KTXT2L )
      CALL XVTEXTE( KTXT2L(1:L), L, 20, 2*LHPXCA+6 )
C
C     RESTAURATION DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      CALL CHARGEFONTE( NOFONT0 )
C
      RETURN
      END
