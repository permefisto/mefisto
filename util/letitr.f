      SUBROUTINE LETITR( NOPROJ, MODECO, NCAS, TEMPS, KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     DEFINIR LE TITRE LEGENDE D'UN TRACE
C -----
C
C ENTREES:
C --------
C NOPROJ  : TYPE DE PROJECTION 0 CI-DESSOUS FIXE LA COORDONNEE A ZERO
C          <=0 PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          = 1 CE SONT DES ISO-SOLUTIONS
C          = 2 CE SONT DES MODES PROPRES
C          = 3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C              SOIT en 2D P2 SOIT P1+BULLE P3 OU en 3D P1 OU en 3D P2
C          = 4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C          = 5 CE SONT DES SOLUTIONS         AUX NOEUDS D'UN MAILLAGE
C          = 6 CE SONT DES NORMES DE VITESSE AUX NOEUDS D'UN MAILLAGE
C          = 7 FONCTION COURANT  DANS UN FLUIDE
C          = 8 MODULE            D'UNE ONDE COMPLEXE NLSE
C          = 9 PARTIE REELLE     D'UNE ONDE COMPLEXE NLSE
C          =10 PARTIE IMAGINAIRE D'UNE ONDE COMPLEXE NLSE
C          =11 ERREURS DU MODULE DE LA VITESSE     D'UN FLUIDE
C          =12 ERREURS DU MODULE DE LA PRESSION P1 D'UN FLUIDE
C          =13 FLECHES REPRESENTANT DES VECTEURS VITESSE
C          =14 ROTATIONNEL DE VITESSES ou TOURBILLONS ou VORTICITES
C
C NCAS   : NUMERO DU VECTEUR A TRAITER PARMI LES VECTEURS
C TEMPS  : TEMPS ACTUEL
C
C SORTIE :
C --------
C KNOM   : CONTENU DU TITRE A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Fevrier 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      CHARACTER*(*)  KNOM
      CHARACTER*16   TEXT
C
      WRITE( TEXT, '(G15.7)' ) TEMPS
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, LTEXT )
C
      IF( MODECO .EQ. 1 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      Les ISOVALEURS au TEMPS  '
           ELSE
              KNOM = 'CASE     The ISOVALUES at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(33:32+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'VALEUR PROPRE                 : Le VECTEUR PROPRE'
         ELSE
         KNOM = 'EIGENVALUE                    : The EIGENVECTOR'
         ENDIF
         KNOM(15:14+LTEXT) = TEXT(1:LTEXT)
         I = NUDCNB( KNOM )
         WRITE( KNOM(I+1:I+5), '(I5)' ) NCAS
C
      ELSE IF( MODECO .EQ. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = 'CAS      Les PRESSIONS au TEMPS  '
         ELSE
            KNOM = 'CASE     The PRESSURES at TIME   '
         ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(35:34+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 4 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      Les ERREURS au TEMPS  '
           ELSE
              KNOM = 'CASE     The ERRORS at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(33:32+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 5 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      Les SOLUTIONS au TEMPS  '
           ELSE
              KNOM = 'CASE     The SOLUTIONS at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(33:32+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 6 .OR. MODECO .EQ. 13 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      Les VITESSES au TEMPS     '
           ELSE
              KNOM = 'CASE     The VELOCITIES at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(34:33+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 7 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      La FONCTION COURANT au TEMPS     '
           ELSE
              KNOM = 'CASE     The STREAM FUNCTION at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(38:37+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 8 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      MODULE de l''Onde NLSE au TEMPS  '
           ELSE
              KNOM = 'CASE     MODULE of NLSE WAVE at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(38:37+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 9 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      PARTIE REELLE Onde NLSE au TEMPS  '
           ELSE
              KNOM = 'CASE     REAL PART of NLSE WAVE at TIME    '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(38:37+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 10 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      PARTIE IMAGINAIRE Onde NLSE au TEMPS '
           ELSE
              KNOM = 'CASE     IMAGINARY PART of NLSE WAVE at TIME  '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(38:37+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 11 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
              KNOM = 'CAS      ERREUR du MODULE de la VITESSE au TEMPS '
           ELSE
              KNOM = 'CASE     VELOCITY MAGNITUDE ERROR at TIME  '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(38:37+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 12 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
             KNOM = 'CAS      ERREUR de la PRESSION du FLUIDE au TEMPS '
           ELSE
             KNOM = 'CASE     FLUID PRESSURE ERROR at TIME '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(38:37+LTEXT) = TEXT(1:LTEXT)
C
      ELSE IF( MODECO .EQ. 14 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
             KNOM = 'CAS      Rot Vitesses au TEMPS  '
           ELSE
             KNOM = 'CASE     Rot Velocities at TIME '
           ENDIF
         WRITE( KNOM(5:8),   '(I4)'    ) NCAS
         KNOM(34:33+LTEXT) = TEXT(1:LTEXT)
C
      ENDIF
C
      CALL PROJ6C( NOPROJ, KNOM )
      CALL TRFINS( KNOM )
C
      RETURN
      END
