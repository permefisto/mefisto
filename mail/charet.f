      SUBROUTINE CHARET( NX, NY, NDIM,   MNXYZS, MNNSEF,
     %                   NUMTQ,  NUTQAD, NUMSO1, NUMSO2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RECHERCHE DE L'ARETE D'UN EF LA PLUS PROCHE DU POINT CLIQUE
C -----     A L'INTERIEUR D'UN EF
C ENTREES:
C --------
C NX     : ABSCISSE ECRAN DU POINT CLIQUE
C NY     : ORDONNEE ECRAN DU POINT CLIQUE
C NDIM   : DIMENSION DE L'ESPACE DE LA TRIANGULATION (2 OU 3)
C MNXYZS : ADRESSE DU TMS XYZSOMMET DE LA SURFACE
C MNNSEF : ADRESSE DU TMS NSEF DE LA SURFACE
C
C SORTIES:
C --------
C NUMTQ  : >0 NUMERO DE L'ELEMENT FINI CONTENANT LE POINT CLIQUE
C          =0 SI POINT CLIQUE A L'EXTERIEUR DU MAILLAGE
C NUTQAD : >0 NUMERO DE L'AUTRE ELEMENT FINI ADJACENT PAR CETTE ARETE
C          =0 SI PAS DE SECOND TRIANGLE CONTENANT L'ARETE (FRONTIERE)
C NUMSO1 : NUMERO DU SOMMET 1 DE L'ARETE COMMUNE
C NUMSO2 : NUMERO DU SOMMET 2 DE L'ARETE COMMUNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY Septembre 2010
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE   (RMCN(1),MCN(1))
      REAL           AXOXYZ(3), AXYZ(3,2)
      REAL           D, DMIN
C
      NUMTQ  = 0
      NUTQAD = 0
      NUMSO1 = 0
      NUMSO2 = 0
C     ADRESSE MCN DU TABLEAU NUSOEF NO DES 4 SOMMETS DES NBTQ EF ACTUELS
      MNSOEL = MNNSEF + WUSOEF
C     NBTQ LE NOMBRE D'EF
      NBTQ   = MCN( MNNSEF + WBEFOB )
C
C     RECHERCHE DE L'EF NUMTQ CONTENANT LE POINT CLIQUE
      MNXYZ = MNXYZS + WYZSOM
      CALL TQPTCLIC( NX, NY, NDIM, MCN(MNXYZ), NBTQ, MCN(MNSOEL), NUMTQ)
      IF( NUMTQ .EQ. 0 ) RETURN
C
C     COORDONNEES REELLES 2D ou AXONOMETRIQUES 3D DU POINT CLIQUE
      AXOXYZ(1) = XOB2PX(NX)
      AXOXYZ(2) = YOB2PX(NY)
C
C     BOUCLE SUR LES ARETES DU TRIANGLE OU QUADRANGLE NUMTQ
      MN = MNSOEL - 5 + 4 * NUMTQ
      IF( MCN(MN+4) .EQ. 0 ) THEN
         NBS = 3
      ELSE
         NBS = 4
      ENDIF
C
      NS2  = MCN(MN+NBS)
      MNS2 = MNXYZ -3 + 3 * NS2
C
C     RECHERCHE DE L'ARETE EN COORDONNEES AXONOMETRIQUES
C     LA PLUS PROCHE DU POINT CLIQUE
      DMIN = 1E28
      DO 10 I=1,NBS
C
C        LES 2 NUMEROS DE SOMMETS DE L'ARETE I-1 DE L'EF NUMTQ
         NS1 = NS2
         NS2 = MCN(MN+I)
C
         MNS1 = MNS2
         MNS2 = MNXYZ -3 + 3 * NS2
C
         IF( NDIM .EQ. 2 ) THEN
            AXYZ(1,1) = RMCN( MNS1   )
            AXYZ(2,1) = RMCN( MNS1+1 )
            AXYZ(3,1) = 0
            AXYZ(1,2) = RMCN( MNS2   )
            AXYZ(2,2) = RMCN( MNS2+1 )
            AXYZ(3,2) = 0
         ELSE
            CALL XYZAXO( RMCN(MNS1), AXYZ(1,1) )
            CALL XYZAXO( RMCN(MNS2), AXYZ(1,2) )
         ENDIF
C
C        DISTANCE DU POINT A L'ARETE I-1
         D = DIPTDRR( AXOXYZ,  AXYZ(1,1),  AXYZ(1,2) )
         IF( D .LT. DMIN ) THEN
            DMIN   = D
            NUMSO1 = NS1
            NUMSO2 = NS2
         ENDIF
C
 10   CONTINUE
C
C     BOUCLE SUR LES EF DU MAILLAGE POUR TROUVER L'EF DE L'AUTRE COTE DE L'ARETE
      DO 30 I=1,NBTQ
C
C          LES NUMEROS DES SOMMETS DE L'EF NUMERO I
           MN = MNSOEL + 4*I - 5
           IF( MCN(MN+4) .EQ. 0 ) THEN
              NBS = 3
           ELSE
              NBS = 4
           ENDIF
C
           NS2 = MCN(MN+NBS)
           DO K=1,NBS
C
C             LE NUMERO DES 2 SOMMETS DE L'ARETE K DE L'EF I
              NS1 = NS2
              NS2 = MCN(MN+K)
C
C             L'ARETE K DE L'EF I EST ELLE L'ARETE NUSOM1-NUSOM2
              IF( NS1 .EQ. NUMSO2 ) THEN
                 IF( NS2 .EQ. NUMSO1 ) THEN
                    IF( I .NE. NUMTQ ) GOTO 8000
                 ENDIF
              ELSE IF( NS1 .EQ. NUMSO1 ) THEN
                 IF( NS2 .EQ. NUMSO2 ) THEN
                    IF( I .NE. NUMTQ ) GOTO 8000
                 ENDIF
              ENDIF
C
           ENDDO
C
 30     CONTINUE
        GOTO 9000
C
C       EF ADJACENT RETROUVE
 8000   NUTQAD = I
C
 9000   RETURN
        END
