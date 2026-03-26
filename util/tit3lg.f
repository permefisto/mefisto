      SUBROUTINE TIT3LG( MODECO, NCAS, TEMPS, VALMIN, VALMAX,
     %                   KTXT3L )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LA 3-EME LIGNE DU TITRE  DU TRACE SELON MODECO
C -----    RESTE A FAIRE UN CALL TRFINS( KTXT3L ) POUR LE TRACE EFFECTIF
C
C ENTREES:
C --------
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
C        = 8  CE SONT LES MODULES D'UNE ONDE COMPLEXE
C             AUX NOEUDS D'UN MAILLAGE
C        = 9  CE SONT DES PARTIES REELLES     D'UNE ONDE COMPLEXE
C        =10  CE SONT DES PARTIES IMAGINAIRES D'UNE ONDE COMPLEXE
C        =11  CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE DE LA VITESSE D'UN FLUIDE
C        =12  CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C        =13  CE SONT DES FLECHES REPRESENTANT DES VECTEURS VITESSE
C        =14  CE SONT DES ROTATIONNEL DE Vitesses ou TOURBILLONS ou VORTICITES
C
C NCAS   : NUMERO DU CAS TRACE
C TEMPS  : TEMPS DU CALCUL OU CETTE SOLUTION A ETE CALCULEE
C VALMIN : MINIMUM DE LA VALEUR A TRACER
C VALMAX : MAXIMUM DE LA VALEUR A TRACER
C
C SORTIE :
C --------
C KTXT3L : LE TEXTE DE LA 3-EME LIGNE DU TITRE A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
      CHARACTER*(*)  KTXT3L
      REAL           TEMPS, VALMIN, VALMAX
      INTEGER        NCAS
C
C     3-LIGNE DU TITRE
      KTXT3L = '         '
      IF( MODECO .EQ. 1 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      au TEMPS                : '
         ELSE
            KTXT3L = 'Case     at TIME                 : '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(19:33), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'VALEUR PROPRE               : '
         ELSE
            KTXT3L = 'EIGENVALUE                  : '
         ENDIF
         WRITE( KTXT3L(15:18),   '(I4)'  ) NCAS
         WRITE( KTXT3L(24:38), '(G15.7)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Les PRESSIONS au TEMPS  '
         ELSE
            KTXT3L = 'Case     The PRESSURES at TIME   '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 4 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Les ERREURS au TEMPS  '
         ELSE
            KTXT3L = 'Case     The ERRORS at TIME   '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 5 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Parcours de PARTICULES au TEMPS  '
            WRITE( KTXT3L(42:55), '(G14.6)' ) TEMPS
         ELSE
            KTXT3L = 'Case     Particles RUN at TIME        '
            WRITE( KTXT3L(32:45), '(G14.6)' ) TEMPS
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
C
      ELSE IF( MODECO .EQ. 6 .OR. MODECO .EQ. 13 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas       Les VITESSES  au TEMPS  '
         ELSE
            KTXT3L = 'Case      The VELOCITIES at TIME  '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 7 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      FONCTION COURANT au TEMPS  '
         ELSE
            KTXT3L = 'Case     STREAM FUNCTION at TIME  '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 8 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      MODULES au TEMPS  '
         ELSE
            KTXT3L = 'Case     MAGNITUDES at TIME  '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 9 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Partie Reelle de l''Onde au TEMPS  '
         ELSE
            KTXT3L = 'Case     Wave Real Part at TIME   '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 10 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Partie Imaginaire de l''Onde au TEMPS  '
         ELSE
            KTXT3L = 'Case     Wave Imaginary Part at TIME            '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(35:48), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 11 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Les ERREURS |Vitesse| au TEMPS  '
         ELSE
            KTXT3L = 'Case     The |Velocity| ERRORS at TIME   '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(42:55), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 12 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      ERREURS sur la PRESSION au TEMPS  '
         ELSE
            KTXT3L = 'Case     The PRESSURE ERRORS at TIME       '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(41:54), '(G14.6)' ) TEMPS
C
      ELSE IF( MODECO .EQ. 14 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KTXT3L = 'Cas      Rot des VITESSES au TEMPS  '
         ELSE
            KTXT3L = 'Case     Rot of VELOCITIES at TIME  '
         ENDIF
         WRITE( KTXT3L(5:8),   '(I4)'    ) NCAS
         WRITE( KTXT3L(38:51), '(G14.6)' ) TEMPS
      ENDIF
C
C     3-EME LIGNE FIN
      I = NUDCNB( KTXT3L )
      KTXT3L(I+1:I+6) = ' MIN= '
      WRITE( KTXT3L(I+7:I+20), '(G14.6)' ) VALMIN
      I = I + 21
      KTXT3L(I:I+6) = ' MAX= '
      WRITE( KTXT3L(I+7:I+20), '(G14.6)' ) VALMAX
C
CCCC     IL RESTE A FAIRE LE TRACE FINAL DE LA LIGNE 3 DU TITRE PAR
CCC      CALL TRFINS( KTXT3L )
C
      RETURN
      END
