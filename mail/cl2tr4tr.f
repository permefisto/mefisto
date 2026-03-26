      SUBROUTINE CL2TR4TR( MNXYZS, MNNSEF, L1ARET, L2ARET, LARETE,
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
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY           Janvier 2021
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
      INTEGER           LARETE(L1ARET,L2ARET)

      IERR = 0

C     CLIQUER SUR UNE ARETE ADJACENTE A 2 TRIANGLES
C     ---------------------------------------------
      CALL SEARCLIC( RMCN(MNXYZS+WYZSOM), L1ARET, L2ARET, LARETE,
     %               NUARPRCL )
      IF( NUARPRCL .LE. 0 ) GOTO 9999

C     LES 2 SOMMETS DE L'ARETE COMMUNE NUARPRCL
      NOST1 = LARETE(1,NUARPRCL)
      NOST2 = LARETE(2,NUARPRCL)

C     LES 2 TRIANGLES ADJACENTS A CETTE ARETE NUARPRCL 
      NOTQ1 = ABS( LARETE(4,NUARPRCL) )
      NOTQ2 = ABS( LARETE(5,NUARPRCL) )

C     AUGMENTATION EVENTUELLE DU TMS XYZSOMMET
C     ----------------------------------------
C     NBSOM NOMBRE DE SOMMETS
      NBSOM = MCN(MNXYZS+WNBSOM)
      L     = WYZSOM + 3 * NBSOM
      IF ( MOXYZS .LT. L+3 ) THEN
         CALL TNMCAU( 'ENTIER', MOXYZS, L+3, MOXYZS, MNXYZS )
C        MNXYZS = MN APRES L'AUGMENTATION!
         MOXYZS = L+3
         NBSOM  = MCN(MNXYZS+WNBSOM)
      ENDIF

C     LES COORDONNEES DU NOUVEAU POINT NBSOM MILIEU DE L'ARETE COMMUNE
C     ----------------------------------------------------------------
      MNX1 = MNXYZS + WYZSOM + 3*NOST1 -3
      MNX2 = MNXYZS + WYZSOM + 3*NOST2 -3
      RMCN(MNXYZS+L  ) = (RMCN( MNX1   )+RMCN( MNX2   ))/2
      RMCN(MNXYZS+L+1) = (RMCN( MNX1+1 )+RMCN( MNX2+1 ))/2
      RMCN(MNXYZS+L+2) = (RMCN( MNX1+2 )+RMCN( MNX2+2 ))/2

      NBSOM = NBSOM + 1
      MCN( MNXYZS+WNBSOM ) = NBSOM
      CALL ECDATE( MCN(MNXYZS) )
      MCN( MNXYZS+MOTVAR(6) ) = NONMTD('~>>>XYZSOMMET')

C      AUGMENTATION EVENTUELLE DU TMS NSEF
C     ------------------------------------
      NBTQ = MCN(MNNSEF+WBEFOB)
      L    = WUSOEF + 4 * NBTQ
      IF( MONSEF .LT. L+8 ) THEN
         CALL TNMCAU( 'ENTIER', L, L+8, L, MNNSEF )
         MONSEF=L+8
      ENDIF

C     ON MODIFIE LES 2 TRIANGLES ou QUADRANGLES AYANT CETTE ARETE COMMUNE
C     POUR EN FORMER 4 AVEC LE MILIEU COMME BARYCENTRE

C     PREMIER ELEMENT FINI NOTQ1 DIVISE EN 2 SOUS-TRIANGLES ou 1Q+1T
C     --------------------------------------------------------------
      NTEMP = MNNSEF + WUSOEF
      MN    = NTEMP  + 4 * NOTQ1 - 4

C     NUMERO DES SOMMETS DU TQ NOTQ1
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
      DO I=1,NBS
         IF( NOSOTQ(I) .EQ. NOST2 ) THEN

C            LE SOMMET NOST2 EST ECRASE PAR LE MILIEU DE L'ARETE
             MCN(MN+I-1) = NBSOM

C            RECHERCHE DU SOMMET PRECEDANT OU SUIVANT = NOST1
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

             IF( NOSOTQ(I1) .EQ. NOST1 ) THEN
                NS1 = NOSOTQ( I2 )
             ELSE
                NS1 = NOSOTQ( I1 )
             ENDIF
             GOTO 10

         ENDIF
      ENDDO

C     DEUXIEME SOUS TRIANGLE ou 1Q+1T DE NOTQ1 DERRIERE LES EF ACTUELS
C     ----------------------------------------------------------------
 10   NTEMP = MNNSEF + WUSOEF + 4*NBTQ
      MCN(NTEMP  ) = NOST2
      MCN(NTEMP+1) = NBSOM
      MCN(NTEMP+2) = NS1
      MCN(NTEMP+3) = 0

C     TROISIEME ELEMENT FINI DANS LE TRIANGLE ou QUADRANGLE ADJACENT
C     --------------------------------------------------------------
      NS2 = 0
      IF( NOTQ2 .GT. 0 ) THEN
         MN = MNNSEF + WUSOEF + 4 * NOTQ2 - 4
         NOSOTQ(1) = MCN(MN  )
         NOSOTQ(2) = MCN(MN+1)
         NOSOTQ(3) = MCN(MN+2)
         NOSOTQ(4) = MCN(MN+3)
         IF( NOSOTQ(4) .EQ. 0 ) THEN
            NBS = 3
         ELSE
            NBS = 4
         ENDIF
         DO I=1,NBS
            IF( NOSOTQ(I) .EQ. NOST2 ) THEN

C              LE SOMMET NOST2 EST ECRASE PAR LE MILIEU DE L'ARETE
               MCN(MN+I-1) = NBSOM

C              RECHERCHE DU SOMMET PRECEDANT OU SUIVANT = NOST1
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

               IF( NOSOTQ(I1) .EQ. NOST1 ) THEN
                  NS2 = NOSOTQ( I2 )
               ELSE
                  NS2 = NOSOTQ( I1 )
               ENDIF
               GOTO 30

            ENDIF
         ENDDO

C        QUATRIEME ELEMENT FINI DANS L'EF ADJACENT DERRIERE LES EF ACTUELS
 30      NTEMP = MNNSEF + WUSOEF + 4*NBTQ + 4
         MCN(NTEMP  ) = NOST2
         MCN(NTEMP+1) = NBSOM
         MCN(NTEMP+2) = NS2
         MCN(NTEMP+3) = 0
      ENDIF

C     FIN DE LA MISE A JOUR DU TMS NSEF
      IF ( NOTQ2 .GT. 0 ) THEN
        NBTQ = NBTQ + 2
      ELSE
        NBTQ = NBTQ + 1
      ENDIF
      MCN(MNNSEF+WBEFOB) = NBTQ
      MCN(MNNSEF+WUTYOB) = 3
      MCN(MNNSEF+WUTYMA) = 0
      MCN(MNNSEF+WUTFMA) = 0
      MCN(MNNSEF+WBSOEF) = 4
      MCN(MNNSEF+WBTGEF) = 0
      MCN(MNNSEF+WBEFTG) = 0
      MCN(MNNSEF+WBEFAP) = 0

      CALL ECDATE(MCN(MNNSEF))
      MCN(MNNSEF+MOTVAR(6)) = NONMTD( '~>>>NSEF' )

 9999 RETURN
      END
