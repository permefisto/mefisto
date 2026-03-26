      SUBROUTINE TTAJPI( NBPTIN, NOPTIN, DISTSO,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                   IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AJOUTER LES POINTS INTERNES UTILISATEUR DANS LE TABLEAU PTXYZD
C -----   METTRE A JOUR NPSOFR NUMERO DU POINT INTERNE

C ENTREES:
C --------
C NBPTIN : NOMBRE DE POINTS INTERNES A AJOUTER
C NOPTIN : NUMEROS DES POINTS INTERNES
C DISTSO : DISTANCE SOUHAITEE AUTOUR DES POINTS INTERNES
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS UTILISABLES DANS PTXYZD
C HEXAPAVE: MIN ET MAX DES COORDONNEES DU PAVAGE
C NBIPAV : NOMBRE D'ARETES DANS LA DIRECTION I DU PAVAGE
C ECHPAV : ECHELLE DANS LA DIRECTION I DU PAVAGE

C ENTREES ET SORTIES:
C -------------------
C NBSOMM : NOMBRE ACTUEL DE POINTS DECLARES DANS PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NPSOFR : 1 000 000  + NUMERO DU  POINT INTERNE DES POINTS AJOUTES
C N1SPAVE: NO DANS PTXYZD DU 1-ER SOMMET DE CHAQUE PAVE
C NOPTSUIV: NO DU POINT SUIVANT DANS LE CHAINAGE DES POINTS DES PAVES

C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 UN POINT INCONNU OU INCORRECT
C          5 SI TABLEAU PTXYZD SATURE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   AVRIL 2008
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))

      INTEGER           NOPTIN(1:NBPTIN),
     &                  NPSOFR(1:MXSOMM)
      REAL              DISTSO(1:NBPTIN)
      DOUBLE PRECISION  PTXYZD(4,MXSOMM), XYZ(3)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3)
      INTEGER           NBIPAV(3), N1SPAVE(*), NOPTSUIV(MXSOMM)

      DO I=1,NBPTIN

C        LE NUMERO DU I-EME POINT IMPOSE PAR L'UTILISATEUR
         NP = NOPTIN(I)

C        LE SOMMET NP EXISTE-T-IL?
         CALL LXNLOU( NTPOIN, NP, NTPO, MNPO )
         IF( NTPO .LE. 0 ) THEN
C           NON
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:10),'(I10)') NP
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'TTAJPI: POINT INCONNU ' // KERR(MXLGER)(1:10)
            ELSE
               KERR(1) = 'TTAJPI: UNKNOWN POINT ' // KERR(MXLGER)(1:10)
            ENDIF
            CALL LEREUR
            IERR = 1
         ELSE
C           OUI
            CALL LXTSOU( NTPO, 'XYZSOMMET', NTSO, MNSO )
            MNSO = MNSO + WYZSOM
            XYZ(1) = RMCN(MNSO)
            XYZ(2) = RMCN(MNSO+1)
            XYZ(3) = RMCN(MNSO+2)

C           AJOUT DU POINT INTERNE IMPOSE I A PTXYZD ET LA GRILLE
            CALL REAJSTPV( HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                     XYZ,DBLE(DISTSO(I)), NBSOMM, MXSOMM, PTXYZD,
     %                     NS )

            IF( NS .EQ. 0 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'ttajpi: TABLEAU PTXYZD SATURE'
               ELSE
                  KERR(1) = 'ttajpi: ARRAY PTXYZD SATURATED'
               ENDIF
               CALL LEREUR
               IERR = 2
               RETURN
            ENDIF

C           LE NUMERO DE POINT DANS LE LEXIQUE DES POINTS
            NPSOFR( NS ) = 1 000 000 + NP

            IF( LANGAG .EQ. 0 ) THEN
              WRITE(IMPRIM,*)'ttajpi: POINT IMPOSE',I,': NO SOMMET=',NS,
     %                       ' NO POINT=',NP
               WRITE(IMPRIM,*)'ttajpi:  XYZ=',XYZ,
     %                     ' LONGUEUR ARETE SOUHAITEE=', DBLE(DISTSO(I))
            ELSE
             WRITE(IMPRIM,*)'ttajpi: IMPOSED POINT',I,': VERTEX NB=',NS,
     %                      ' POINT NB=',NP
               WRITE(IMPRIM,*) 'ttajpi: XYZ=',XYZ,
     %                         ' WISHED EDGE LENGTH=', DBLE(DISTSO(I))
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END
