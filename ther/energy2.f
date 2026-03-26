      SUBROUTINE ENERGY2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VERSION 2D avec TRIANGULATIONS SUBDIVISES 1TRIANGLE=>4SOUS-TRIANGLES
C ----- COMPARAISONS DES ENERGIES CINETIQUE ET POTENTIELLES
C       TRACER LA COURBE Log( |-2KE + 4PE_1 + 2 PE_2 | )
C       CALCULES PAR SSPACE en FONCTION de -Log(h)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Paris & Veulettes 76      Aout 2008
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*120 TITRE
C
C     LES DONNEES CALCULEES PAR SSPACE
      PARAMETER (NBI2=9,NBH=3)
      INTEGER    NBTRIA(NBH), NBNODE(NBH)
      REAL       H(NBH)
      REAL       VALEUR(4,NBI2,NBH)
C     VALEUR(1,noEV,h) EIGENVALUE
C     VALEUR(2,noEV,h) KINETIC ENERGY
C     VALEUR(3,noEV,h) POTENTIAL1 ENERGY
C     VALEUR(4,noEV,h) POTENTIAL2 ENERGY
      REAL ENER(15,5)
C     ENER(noEV,h) = -2KE + 4PE_1 + 2 PE_2
C
C     TABLEAUX AUXILIAIRES
      REAL       PT(2,NBH), LNH0, LNH1, V0, V1
      DOUBLE PRECISION  DROITE(3)
C
C     INITIALISATIONS DES DONNEES
      DATA H      /  2.0, 1.0,   0.5  /
      DATA NBNODE / 1245, 4237, 19665 /
      DATA NBTRIA /  612, 2448,  9792 /
C
      DATA VALEUR /
     %  0.650319159,  0.458237639,  0.548981332,   -0.0824091556,
     %  1.38342834,   1.06795683,   1.56991814,    -0.469487572,
     %  1.6557976,    0.958412871,  0.546582008,    0.424093674,
     %  2.39979124,   1.57222976,   1.56117664,     0.0469731237,
     %  2.65909529,   1.97795513,   2.61212871,    -0.624924186,
     %  2.67525339,   1.45106025,   0.551005037,    0.948690539,
     %  3.43010044,   2.05845834,   1.5930242,      0.575129938,
     %  3.69279265,   1.96953592,   0.584140899,    1.43118631,
     %  3.71481848,   2.45405067,   2.69875039,    -0.0886073198,
     %  0.647447705,  0.459343495,  0.543184708,  -0.0834881699,
     %  1.37321007,   1.06994379,   1.53664464,   -0.465056054,
     %  1.64789271,   0.959216824,  0.543316866,   0.417017496,
     %  2.37433624,   1.56993058,   1.53720899,    0.0358011914,
     %  2.63166428,   1.96129201,   2.59268797,   -0.625971805,
     %  2.64922643,   1.45922951,   0.54353956,    0.918227189,
     %  3.37663364,   2.0713881,    1.53844844,    0.536021255,
     %  3.63442063,   2.46166638,   2.59526178,   -0.124876632,
     %  3.65147376,   1.96207754,   0.543582104,   1.41760518,
     %  0.647249043,  0.459374279,  0.543045227,  -0.0836478766,
     %  1.37232399,   1.06988309,   1.53510546,   -0.465111798,
     %  1.64727843,   0.959365316,  0.543053875,   0.416386204,
     %  2.37239838,   1.57000085,   1.53503136,    0.0348819569,
     %  2.62822628,   1.96068664,   2.58629343,   -0.625607028,
     %  2.64736724,   1.45977624,   0.543138773,   0.916021597,
     %  3.37257433,   2.07055125,   1.5350959,     0.534475185,
     %  3.62858033,   2.46202243,   2.58626376,   -0.126573967,
     %  3.64790845,   1.96332662,   0.544548191,   1.41230775
     % /
C     LA FENETRE DE TRACE
      print *
      HMIN =  1E20
      HMAX = -1E20
      RMIN =  1E20
      RMAX = -1E20
      DO 3 K = 1, NBH
         HH = -Log( H(K) )
         HMIN = MIN( HMIN, HH )
         HMAX = MAX( HMAX, HH )
         DO 2 J = 1, NBI2
            ENER(J,K) = -2*VALEUR(2,J,K)+4*VALEUR(3,J,K)+2*VALEUR(4,J,K)
            print *,H(K),'h: eigv=',VALEUR(1,J,K),
     %      ' KE =',VALEUR(2,J,K),
     %      ' PE1=',VALEUR(3,J,K),
     %      ' PE2=',VALEUR(4,J,K),
     %      '  -2KE+ 4PE_1+ 2 PE_2=',ENER(J,K)
            V1 = Log( ABS( ENER(J,K) ) )
            RMIN = MIN( RMIN, V1 )
            RMAX = MAX( RMAX, V1 )
 2       CONTINUE
 3    CONTINUE
      print *,'HMIN=',HMIN,' HMAX=',HMAX,' RMIN=',RMIN,'  RMAX=',RMAX
      HH = ( HMAX - HMIN ) / 15
      RH = ( RMAX - RMIN ) / 15
C
C     MISE SUR FICHIER.eps du TRACE
      CALL xvinitierps( 1 )
      CALL EFFACE
      CALL FENETRE( HMIN-HH, HMAX+HH, RMIN-RH, RMAX+2*RH )
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
      CALL TEXTE2D( NCNOIR, HMIN, RMAX+RH/3,
     %             'Log( |-2 KE + 4 PE_1 + 2 PE_2| )')
      CALL TEXTE2D(NCNOIR, HMAX, RMIN, '-Log(h)' )
C
C     LE TRACE DES AXES 2D
      CALL TRAXE2
C
C     LE TRACE DES ( -Log(hi), Log( |-2 KE + 4 PE_1 + 2 PE_2| ) ) (-Log(hi))
      DO 20 J = 1, NBI2
C
C        CHOIX DE LA COULEUR DE TRACE
         NC = J
         IF( NC .GT. 9 ) NC = NC + 2
         print *
C
C        LE POINT A TRACER ( -Log(H(1)), Log( ABS( |-2 KE + 4 PE_1 + 2 PE_2| ) )
         LNH0 = -Log(H(1))
         V0  =  Log( ABS( ENER(J,1) ) )
C        LE NO DE LA VALEUR PROPRE A GAUCHE
         CALL ENTIER2D(  NC, HMAX+HH/7, V0, J )
         CALL SYMBOLE2D( NC, LNH0, V0, '*' )
C
C        DROITE PAR MOINDRES CARRES
         PT(1,1) = LNH0
         PT(2,1) = V0
C
C        BOUCLE SUR LA TAILLE DES ARETES H
         DO 10 K = 2, NBH
C
C           LE POINT A TRACER
            LNH1 = -Log(H(K))
            V1   =  Log( ABS( ENER(J,K) ) )
            CALL XVEPAISSEUR( 0 )
            CALL TRAIT2D(   NC, LNH0, V0,  LNH1, V1 )
            CALL SYMBOLE2D( NC, LNH1, V1, '*' )
C
C           DROITE PAR MOINDRES CARRES
            PT(1,K) = LNH1
            PT(2,K) = V1
C
            print *,'SLOPE j=',j,' k=',K,'=',
     %              (V1-V0)/(-Log(H(K))+Log(H(K-1)))
C
            LNH0 = LNH1
            V0  = V1
 10      CONTINUE
C
C        LE NO DE LA VALEUR PROPRE A DROITE
         CALL ENTIER2D( NC, HMIN-HH/2, V1, J )
C
C        DROITE PAR MOINDRES CARRES
         CALL DRDIMI( NBH, PT, DROITE, IERR )
C
         print *
         print *,'EIGV',J,': Least square LINE(',J,')=',DROITE
         print *,'Y=',-DROITE(1)/DROITE(2),' * X +',-DROITE(3)/DROITE(2)
         print *,'Least square LINE SLOPE',J,'=',-DROITE(1)/DROITE(2)
C
C        TRACE DE LA DROITE
         LNH0 = -Log(H(1))
         V0  = REAL( - ( DROITE(1) * LNH0 + DROITE(3) ) / DROITE(2) )
         LNH1 = -Log(H(NBH))
         V1  = REAL( - ( DROITE(1) * LNH1 + DROITE(3) ) / DROITE(2) )
         CALL XVEPAISSEUR( 2 )
         CALL TRAIT2D( NC, LNH0, V0,  LNH1, V1 )
C
C        LE NO DE LA DROITE A GAUCHE
         CALL ENTIER2D(  NC, LNH0, V0, J )
C        LE NO DE LA DROITE A DROITE
         CALL ENTIER2D( NC, LNH1, V1, J )
C
 20   CONTINUE
C
C     LE TITRE DU GRAPHIQUE
      CALL CHOIXFONTE( 25 )
      TITRE = '[(-1/2)Laplacian +(x^4-2x^2+2y^2)/4] psi(x,y)= E*psi(x,y)
     %: Log( |-2 KE + 4 PE_1 + 2 PE_2|) (-Log(h))'
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+HH, RMAX+RH, TITRE(1:L) )
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
C     MISE SUR FICHIER energyh.eps du TRACE
C     ATTENTION PASSAGE PAR VARIABLE OBLIGATOIRE
      NBC = 8
      CALL xvsauverps( 'energy2d', NBC )
      print *, 'NBC=',NBC
C
C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO
C
      CALL XVEPAISSEUR( 1 )
      RETURN
      END
