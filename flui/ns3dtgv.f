      SUBROUTINE NS3DTGV
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES RESULTATS DES TESTS SUR LE TOURBILLON DE TAYLOR-GREEN
C ----- OU LA SOLUTION EXACTE EST CONNUE ET POUR LE CAS WATER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  TEXAS A & M UNIVERSITY at QATAR Fevrier 2012
C23456---------------------------------------------------------------012
      PARAMETER (NBMETH=9, NBCAS=29, NBK=19)
C
      include"./incl/langue.inc"
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      include"./incl/traaxe.inc"
      include"./incl/mecoit.inc"
C
C     SYMBOLE <=> METHODE UTILISEE
      CHARACTER*48  METHOD(NBMETH)
      DATA          METHOD /
     %     '# Mefisto Scheme 2      dt=0.01',
     %     '§ Mefisto Scheme 2+Lagr dt=0.01',
     %     '+ Mefisto Scheme 3 PISO dt=0.01',
     %     '@ OpenFOAM RAS Model    dt=0.001',
     %     '& OpenFOAM Dirichlet    dt=0.2',
     %     '§ OpenFOAM Dirichlet    dt=0.01',
     %     '£ OpenFOAM Dirichlet    dt=0.001',
     %     '$ IcoFOAM  Laminar Flow dt=0.005',
     %     'O IcoFOAM  Dirichlet    dt=0.01'  /
C
C     COULEUR = NUMERO DE LA METHODE
      INTEGER COULMETH(NBCAS)
      DATA    COULMETH /
     %        1, 1, 3, 1, 2, 3,
     %        4, 4, 4,
     %        5, 5, 5,
     %        6, 7,
     %        8, 8,
     %        9, 9, 9,
     %        1, 1, 1,
     %        4, 5,
     %        8, 8,
     %        9, 9, 9  /
C
C     NOMBRE D'ELEMENTS EQUIVALENTS CUBES DU MAILLAGE
      REAL  NBELT(NBCAS)
      DATA  NBELT /  1000,   8000,    8000, 64000, 64000, 64000,
     %               1000,   8000,   64000,
     %               1000,   8000,   64000, 64000, 64000,
     %              15625, 125000,
     %               1000,   8000,   64000,
     %               1000,   8000,   64000,
     %               8000,   8000,
     %              15625, 125000,
     %               1000,   8000,   64000 /
C
C     MeanVelocityError=  %
      REAL  WMVE(NBCAS)
      DATA  WMVE /
     %        2.84433,
     %        1.04473,
     %        1.14438,
     %        0.412965,
     %        0.412966,
     %        0.501447,
     %        0.623138,
     %        0.280513,
     %        0.127972,
     %        1.92489,
     %        0.3834,
     %        0.0849529,
     %        0.155654,
     %        0.536506,
     %        0.960991,
     %        0.418052,
     %        3.29929,
     %        1.32149,
     %        0.380152,
     %        3.50866,
     %        1.03251,
     %        0.397715,
     %        0.280378,
     %        0.38324,
     %        0.960751,
     %        0.417437,
     %        3.29883,
     %        1.32018,
     %        0.378208 /
C
C     MeanPressureError=  %
      REAL  WMPE(NBCAS)
      DATA  WMPE /
     %      1.58405,
     %      0.688216,
     %      0.334859,
     %      0.178636,
     %      0.178632,
     %      0.342025,
     %      1.2521,
     %      0.421843,
     %      0.137297,
     %      0.845348,
     %      0.165799,
     %      0.0469224,
     %      0.0626353,
     %      0.267623,
     %      1.39141,
     %      0.534647,
     %      1.44135,
     %      0.55885,
     %      0.131764,
     %      1.58600,
     %      0.684440,
     %      0.169819,
     %      0.477805,
     %      0.165827,
     %      1.38899,
     %      0.532008,
     %      1.44075,
     %      0.557838,
     %      0.130068 /
C
C     CPU on Sony=   seconds  (Uses an equivalent ratio coefficient)
      REAL   WCPU(NBCAS)
      DATA   WCPU /
     %   68,
     %   499,
     %   715,
     %   4472,
     %   4394,
     %   6025,
     %   41,
     %   341,
     %   2330,
     %   0.08,
     %   0.48,
     %   5.42,
     %   23,
     %   189,
     %   19,
     %   120,
     %   0.34,
     %   3.68,
     %   51.8,
     %   67.5,
     %   503,
     %   4684,
     %   848,
     %   0.48,
     %   12,
     %   272,
     %   0.34,
     %   3.69,
     %   50.2 /
C
C     TRACE  MeanVelocityError( NbElt )
      CALL TRTABLES( NBK, NBELT, WMVE, 'MeanVelocityError',
     %  NBMETH, METHOD, COULMETH,
     % 'Number of Elements', 'MeanVelocityError %',
     %'Water Cube: MeanVelocityError %= f(Number of Elements & Method)',
     % ' ',
     % ' ' )
C
C     TRACE  MeanPressureError( NbElt )
      CALL TRTABLES( NBK, NBELT, WMPE, 'MeanPressureError',
     % NBMETH, METHOD, COULMETH,
     %'Number of Elements', 'MeanPressureError %',
     %'Water Cube: MeanPressureError %= f(Number of Elements & Method)',
     % ' ',
     % ' ' )
C
C     TRACE CPU( NbElt )
      CALL TRTABLES( NBK, NBELT, WCPU, 'CpuSeconds',
     %  NBMETH, METHOD, COULMETH,
     % 'Number of Elements', 'CPU Seconds',
     % 'Water Cube: CPU Seconds = f( Number of Elements & Method )',
     % ' ',
     % ' ' )
C
      RETURN
      END
