       SUBROUTINE S3TRIA( NDIM, MNXYZS, MONSEF, MNNSEF, NUTQCR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CLIQUER 3 SOMMETS ET EN CREER UN TRIANGLE
C -----
C ENTREES:
C --------
C NDIM   : DIMENSION (2 ou 3) DE L'ESPACE DE LA SURFACE
C MNXYZS : ADRESSE DU TABLEAU XYZSOMMET   DE LA SURFACE
C MONSEF : NOMBRE DE MOTS DU TMS NSEF     DE LA SURFACE
C MNNSEF : ADRESSE DU TABLEAU NSEF        DE LA SURFACE

C SORTIES:
C --------
C MONSEF : NOMBRE DE MOTS DU TMS NSEF     DE LA SURFACE
C NUTQCR : >0 NUMERO DU TRIANGLE CREE
C          =0 SI PAS CREE
C IERR   : =0 SI PAS D'ERREUR
C          =1 SI POINT CLIQUE EN DEHORS DE LA TRIANGULATION
C          =2 SI 2 SOMMETS SONT IDENTIQUES-> TRIANGLE INCORRECT NON CREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Labo J-L. LIONS  UPMC  PARIS   SEPTEMBRE 2007
C MODIFS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NUMSTO(3)

      IERR = 0

C     NOMBRE DE TRIANGLES ou QUADRANGLES du MAILLAGE
      NBEFOB = MCN( MNNSEF + WBEFOB )

C     NOMBRE DE SOMMETS du MAILLAGE
      NBSOM = MCN( MNXYZS + WNBSOM )

      MNX = MNXYZS + WYZSOM
      DO K=1,3

C        RECHERCHE DU SOMMET LE PLUS PROCHE NUMSTO(K) DU POINT CLIQUE
C        ============================================================
         CALL SESTCLIC( RMCN(MNXYZS+WYZSOM), ITSTCLIC, NSCLIC )
         IF( NSCLIC .LE. 0 ) THEN
C           NSCLIC=0 SI LE POINT CLIQUE EST HORS MAILLAGE ou ABANDON DEMANDE
            IERR   = 1
            NUTQCR = 0
            GOTO 9999
         ENDIF

C        TRACE EN JAUNE DE L'ITEM SOMMET '+NSCLIC'
         NUMSTO(K) = NSCLIC
         MN = MNXYZS + WYZSOM + 3*NSCLIC -3
         PRINT*,'s3tria: SOMMET',K,' CLIQUE: No=',NSCLIC,
     %          ' X=',RMCN(MN),' Y=',RMCN(MN+1),' Z=',RMCN(MN+2)

      ENDDO

C     VERIFICATION SI 2 SOMMETS DES 3 SOMMETS SONT IDENTIQUES
C     =======================================================
      IF( NUMSTO(1) .EQ. NUMSTO(2)  .OR.
     %    NUMSTO(2) .EQ. NUMSTO(3)  .OR.
     %    NUMSTO(3) .EQ. NUMSTO(1)  ) THEN
         PRINT *,'INCORRECT TRIANGLE de SOMMETS=',NUMSTO
         IERR   = 2
         NUTQCR = 0
         GOTO 9999
      ENDIF

C     LE TRIANGLE EST IL DE SURFACE POSITIVE?
C     =======================================
      IF( NDIM .EQ. 2 ) THEN
         MN = MNXYZS + WYZSOM - 3
         SURF = SURTR2( RMCN(MN+3*NUMSTO(1)),
     %                  RMCN(MN+3*NUMSTO(2)),
     %                  RMCN(MN+3*NUMSTO(3)) )
         IF( SURF .LE. 0 ) THEN
C           PAR PERMUTATION ST 2-3 LA SURFACE DEVIENT POSITIVE
            K         = NUMSTO(2)
            NUMSTO(2) = NUMSTO(3)
            NUMSTO(3) = K
         ENDIF
      ENDIF

C     AJOUT DU TRIANGLE NUMSTO(1:3) AU TMS NSEF
C     =========================================
      L = WUSOEF + 4 * NBEFOB
      IF( MONSEF .LT. L+4 ) THEN
         CALL TNMCAU( 'ENTIER', MONSEF, L+64, L, MNNSEF )
         MONSEF=L+64
      ENDIF
      MN = MNNSEF + WUSOEF + 4 * NBEFOB - 1
      DO K=1,3
         MCN(MN+K) = NUMSTO(K)
      ENDDO
      MCN(MN+4) = 0

C     UN TRIANGLE DE PLUS
      NBEFOB = NBEFOB + 1
      MCN(MNNSEF+WBEFOB) = NBEFOB

C     NUMERO DU TRIANGLE CREE
      NUTQCR = NBEFOB

      CALL ECDATE( MCN(MNNSEF) )
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )

 9999 RETURN
      END
