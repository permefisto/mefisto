      SUBROUTINE NS2DTGV
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES RESULTATS DES TESTS SUR LE TOURBILLON DE TAYLOR-GREEN
C ----- OU LA SOLUTION EXACTE EST CONNUE ET POUR LE CAS WATER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  TEXAS A & M UNIVERSITY at QATAR Fevrier 2012
C23456---------------------------------------------------------------012
      PARAMETER (NBMETH=9, NBCAS=28, NBK=19)
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
     %     '# Mefisto  Scheme 1  dt=0.005',
     %     '§ Mefisto  Scheme 2  dt=0.005',
     %     'x Mefisto  Scheme 2  dt=0.01',
     %     '* Mefisto  Scheme 2'' dt=0.01',
     %     '+ Mefisto  Sch3 PISO dt=0.005',
     %     '@ OpenFOAM RAS Model dt=0.001',
     %     '& OpenFOAM Dirichlet dt=0.2 ',
     %     '$ IcoFOAM  LaminFlow dt=0.005',
     %     'O IcoFOAM  Dirichlet dt=0.01'  /
C
C     COULEUR = NUMERO DE LA METHODE
      INTEGER COULMETH(NBCAS)
      DATA    COULMETH /
     %        3, 3, 5, 1, 2, 4, 5,
     %        6, 6, 6,
     %        7, 7, 7,
     %        8, 8, 8,
     %        9, 9, 9,
     %        3, 6, 7,
     %        8, 8, 8,
     %        9, 9, 9  /
C
C     NOMBRE D'ELEMENTS EQUIVALENTS CARRES DU MAILLAGE
      REAL  NBELT(NBCAS)
      DATA  NBELT / 625, 2500,  2500, 10000, 10000, 10000, 10000,
     &              625, 2500, 10000,
     &              625, 2500, 10000,
     &              625, 2500, 10000,
     &              625, 2500, 10000,
     &             2500, 2500,  2500,
     &              625, 2500, 10000,
     &              625, 2500, 10000  /
C
C     MeanVelocityError=  %
      REAL  WMVE(NBCAS)
      DATA  WMVE /
     %            0.692729,
     %            0.324656,
     %            0.288672,
     %            0.143302,
     %            0.143303,
     %            0.143305,
     %            0.203396,
     %            0.221423,
     %            0.0974416,
     %            0.0363806,
     %            0.204664,
     %            0.0438288,
     %            0.0099625,
     %            0.960995,
     %            0.418064,
     %            0.152334,
     %            0.670651,
     %            0.206421,
     %            0.051989,
     %            0.312662,
     %            0.0970982,
     %            0.0438724,
     %            0.960826,
     %            0.417649,
     %            0.151526,
     %            0.669855,
     %            0.205219,
     %            0.0506976 /
C
C     MeanPressureError=  %
      REAL  WMPE(NBCAS)
      DATA  WMPE /
     &            1.33268,
     &            0.234900,
     &            0.170999,
     &            0.0796487,
     &            0.0796488,
     &            0.0796473,
     &            0.276159,
     &            0.322038,
     &            0.12163,
     &            0.0516565,
     &            0.112359,
     &            0.0323194,
     &            0.0106378,
     &            1.38197,
     &            0.530431,
     &            0.220849,
     &            0.331463,
     &            0.113134,
     &            0.0324828,
     &            0.209470,
     &            0.121012,
     &            0.0321864,
     &            1.3818,
     &            0.530004,
     &            0.21913,
     &            0.329874,
     &            0.111332,
     &            0.0307396 /
C
C     CPU on Sony=   seconds  (Uses an equivalent ratio coefficient)
      REAL   WCPU(NBCAS)
      DATA   WCPU /  3.75,  15.83,   20.,   181, 31.67,  77,  84,
     &               6.67,  23.38,   95.,
     &               0.05,   0.12,   0.5,
     &               0.25,   0.67,  4.96,
     &               0.16,   0.36,  1.96,
     &               9.17,  23.38,  0.13,
     &               0.25,   0.67,  4.34,
     &               0.15,   0.36,  2.03 /
C
C     TRACE  MeanVelocityError( NbElt )
      CALL TRTABLES( NBK, NBELT, WMVE, 'MeanVelocityError',
     %  NBMETH, METHOD, COULMETH,
     % 'Number of Elements', 'MeanVelocityError %',
     % 'Water 2d: MeanVelocityError % = f(Number of Elements & Method)',
     % ' ',
     % ' ' )
C
C     TRACE  MeanPressureError( NbElt )
      CALL TRTABLES( NBK, NBELT, WMPE, 'MeanPressureError',
     %  NBMETH, METHOD, COULMETH,
     % 'Number of Elements', 'MeanPressureError %',
     % 'Water 2d: MeanPressureError % = f(Number of Elements & Method)',
     % ' ',
     % ' ' )
C
C     TRACE CPU( NbElt )
      CALL TRTABLES( NBK, NBELT, WCPU, 'CpuSeconds',
     %  NBMETH, METHOD, COULMETH,
     % 'Number of Elements', 'CPU Seconds',
     % 'Water 2d: CPU Seconds = f( Number of Elements & Method )',
     % ' ',
     % ' ' )
C
      RETURN
      END
