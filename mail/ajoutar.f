       SUBROUTINE AJOUTAR( NX, NY, NDIM, MNXYZS, MNNSEF,
     %                     MOXYZS, MONSEF, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DECOUPER UNE ARETE DU MAILLAGE ET LES 2 TRIANGLES ADJACENTS
C -----    en 4 SOUS-TRIANGLES
C ENTREES:
C --------
C NX     : ABSCISSE CLIQUEE POUR L ARETE A PERMUTER
C NY     : ORDONNEE CLIQUEE POUR L ARETE A PERMUTER
C NDIM   : DIMENSION DE L'ESPACE DE LA TRIANGULATION (2 OU 3)
C MNXYZS : ADRESSE DU TMS XYZSOMMET DE LA SURFACE
C MNNSEF : ADRESSE DU TMS NSEF DE LA SURFACE
C
C MODIFIES:
C ---------
C MOXYZS : NOMBRE DE MOTS DU TMS XYZSOMMET DE LA SURFACE FINALE
C MONSEF : NOMBRE DE MOTS DU TMS NSEF      DE LA SURFACE FINALE
C NBSOM  : NOMBRE DE SOMMETS
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: LEVI-OVSIOUK DEA ANALYSE NUMERIQUE UPMC PARIS    JANVIER 2000
C MODIFS : Alain PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 2000
C AJOUTS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      INTEGER           NOSOTQ(1:4)
C
      IERR = 0
C     NBSOM NOMBRE DE SOMMETS
      NBSOM = MCN(MNXYZS+WNBSOM)
C
C     RECHERCHE DE L'ARETE D'UN EF LA PLUS PROCHE DU POINT CLIQUE
C     ===========================================================
      CALL CHARET( NX, NY, NDIM, MNXYZS, MNNSEF,
     %             NUMTQ,  NUTQAD, NUMSO1, NUMSO2 )
      IF( NUMTQ  .LE. 0 ) THEN
         NBLGRC(NRERR)=1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='POINT CLIQUE HORS MAILLAGE'
         ELSE
            KERR(1)='CLICKED POINT OUTSIDE the MESH'
         ENDIF
         CALL LEREUR
         IERR=1
         RETURN
      ENDIF
C
C     MISE A JOUR DU TMS XYZSOMMET
C     ============================
      NBSOM = MCN(MNXYZS+WNBSOM)
      L     = WYZSOM + 3 * NBSOM
      IF ( MOXYZS .LT. L+3 ) THEN
         CALL TNMCAU( 'ENTIER', MOXYZS , L+3 , MOXYZS , MNXYZS)
C        MNXYZS=MN APRES L'AUGMENTATION!
         MOXYZS=L+3
         NBSOM=MCN(MNXYZS+WNBSOM)
      ENDIF

C     LES COORDONNEES DU NOUVEAU POINT MILIEU DE L'ARETE
      MNX1 = MNXYZS + WYZSOM + 3*NUMSO1 -3
      MNX2 = MNXYZS + WYZSOM + 3*NUMSO2 -3
      RMCN(MNXYZS+L  ) = (RMCN( MNX1   )+RMCN( MNX2   ))/2
      RMCN(MNXYZS+L+1) = (RMCN( MNX1+1 )+RMCN( MNX2+1 ))/2
      RMCN(MNXYZS+L+2) = (RMCN( MNX1+2 )+RMCN( MNX2+2 ))/2

      NBSOM = NBSOM + 1
      MCN(MNXYZS+WNBSOM)=NBSOM
      CALL ECDATE(MCN(MNXYZS))
      MCN(MNXYZS+MOTVAR(6))=NONMTD('~>>>XYZSOMMET')

C     MISE A JOUR DU TMS NSEF
C     =======================
      NBTQ = MCN(MNNSEF+WBEFOB)
      L    = WUSOEF + 4 * NBTQ
      IF( MONSEF .LT. L+8 ) THEN
         CALL TNMCAU( 'ENTIER', L , L+8 , L , MNNSEF )
         MONSEF=L+8
      ENDIF

C     ON MODIFIE LES 2 TRIANGLES ou QUADRANGLES AYANT CETTE ARETE COMMUNE
C     POUR EN FORMER 4 AVEC LE MILIEU COMME BARYCENTRE
C
C     PREMIER ELEMENT FINI NUMTQ DIVISE EN 2 SOUS-TRIANGLES ou 1Q+1T
C     --------------------------------------------------------------
      NTEMP = MNNSEF + WUSOEF
      MN    = NTEMP  + 4 * NUMTQ - 4

C     NUMERO DES SOMMETS DU TQ NUMTQ
      NOSOTQ(1) = MCN(MN  )
      NOSOTQ(2) = MCN(MN+1)
      NOSOTQ(3) = MCN(MN+2)
      NOSOTQ(4) = MCN(MN+3)
      IF( NOSOTQ(4) .EQ. 0 ) THEN
         NBS = 3
      ELSE
         NBS = 4
      ENDIF
      NS1 = 0
      DO 10 I=1,NBS
         IF( NOSOTQ(I) .EQ. NUMSO2 ) THEN

C            LE SOMMET NUMSO2 EST ECRASE PAR LE MILIEU DE L'ARETE
             MCN(MN+I-1) = NBSOM

C            RECHERCHE DU SOMMET PRECEDANT OU SUIVANT = NUMSO1
             IF( I .LT. NBS ) THEN
                I1 = I+1
             ELSE
                I1 = 1
             ENDIF

             IF( I .GT. 1 ) THEN
                I2 = I-1
             ELSE
                I2 = NBS
             ENDIF

             IF( NOSOTQ(I1) .EQ. NUMSO1 ) THEN
                NS1 = NOSOTQ( I2 )
             ELSE
                NS1 = NOSOTQ( I1 )
             ENDIF
             GOTO 11

          ENDIF
 10   CONTINUE
C
C     DEUXIEME SOUS TRIANGLE ou 1Q+1T DE NUMTQ DERRIERE LES EF ACTUELS
C     ----------------------------------------------------------------
 11   NTEMP = MNNSEF + WUSOEF + 4*NBTQ
      MCN(NTEMP  ) = NUMSO2
      MCN(NTEMP+1) = NBSOM
      MCN(NTEMP+2) = NS1
      MCN(NTEMP+3) = 0
C
C     TROISIEME ELEMENT FINI DANS LE TRIANGLE ou QUADRANGLE ADJACENT
C     --------------------------------------------------------------
      NS2 = 0
      IF( NUTQAD .NE. 0 ) THEN
         MN = MNNSEF + WUSOEF + 4 * NUTQAD - 4
         NOSOTQ(1) = MCN(MN  )
         NOSOTQ(2) = MCN(MN+1)
         NOSOTQ(3) = MCN(MN+2)
         NOSOTQ(4) = MCN(MN+3)
         IF( NOSOTQ(4) .EQ. 0 ) THEN
            NBS = 3
         ELSE
            NBS = 4
         ENDIF
         DO 20 I=1,NBS
            IF( NOSOTQ(I) .EQ. NUMSO2 ) THEN
C
C              LE SOMMET NUMSO2 EST ECRASE PAR LE MILIEU DE L'ARETE
               MCN(MN+I-1) = NBSOM
C
C              RECHERCHE DU SOMMET PRECEDANT OU SUIVANT = NUMSO1
               IF( I .LT. NBS ) THEN
                  I1 = I+1
               ELSE
                  I1 = 1
               ENDIF
C
               IF( I .GT. 1 ) THEN
                  I2 = I-1
               ELSE
                  I2 = NBS
               ENDIF
C
               IF( NOSOTQ(I1) .EQ. NUMSO1 ) THEN
                  NS2 = NOSOTQ( I2 )
               ELSE
                  NS2 = NOSOTQ( I1 )
               ENDIF
               GOTO 21
C
            ENDIF
 20      CONTINUE
C
C        QUATRIEME ELEMENT FINI DANS L'EF ADJACENT DERRIERE LES EF ACTUELS
 21      NTEMP = MNNSEF + WUSOEF + 4*NBTQ + 4
         MCN(NTEMP  ) = NUMSO2
         MCN(NTEMP+1) = NBSOM
         MCN(NTEMP+2) = NS2
         MCN(NTEMP+3) = 0
      ENDIF
C
C     FIN DE LA MISE A JOUR DU TMS NSEF
      IF ( NUTQAD .NE. 0 ) THEN
        NBTQ = NBTQ + 2
      ELSE
        NBTQ = NBTQ + 1
      ENDIF
      MCN(MNNSEF+WBEFOB)=NBTQ
      MCN(MNNSEF+WUTYOB)=3
      MCN(MNNSEF+WUTYMA)=0
      MCN(MNNSEF+WUTFMA)=0
      MCN(MNNSEF+WBSOEF)=4
      MCN(MNNSEF+WBTGEF)=0
      MCN(MNNSEF+WBEFTG)=0
      MCN(MNNSEF+WBEFAP)=0
C
      CALL ECDATE(MCN(MNNSEF))
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )
C
      RETURN
      END
