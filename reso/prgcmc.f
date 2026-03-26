      SUBROUTINE PRGCMC( MNTOPO, MNNPEF, MNXYZN, NBDLNO, NCODSA,
     %                   MNLPLI, MNLPCO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE STOCKAGE CONDENSE DE LA MATRICE DU MAILLAGE POUR
C -----    NBDLNO DEGRES DE LIBERTE ( CONSTANT EN TOUT NOEUD ) PAR NOEUD

C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU TABLEAU TOPOLOGIE   DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"XXXX DE L'OBJET
C MNXYZN : ADRESSE  MCN DU TABLEAU XYZNOEUD    DE L'OBJET
C NBDLNO : NOMBRE CONSTANT DE DEGRES DE LIBERTE PAR NOEUD
C NCODSA : CODE DE STOCKAGE DE LA MATRICE PROFIL
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE

C SORTIES:
C --------
C MNLPLI : ADRESSE MCN DU POINTEUR SUR LES COEFFICIENTS DIAGONAUX
C          DE LA MATRICE DU SYSTEME LINEAIRE
C MNLPCO : ADRESSE MCN DE LA LISTE DES NUMEROS DE COLONNES
C          DES COEFFICIENTS NON NULS DE LA MATRICE
C IERR   : 0 SI PAS D'ERREUR
C          1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY       ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1989
C MODIFS : ALAIN  PERRONNET  Saint PIERRE du PERRAY         FEVRIER 2021
C23456---------------------------------------------------------------012
C     RUSE POUR CALCULER 2**31-1 PLUS GRAND ENTIER STOCKABLE SUR 32 BITS
C                                        -2 pour la FIN de BOUCLE
      PARAMETER         ( MAXENTIER=2**30-2+2**30 )
 
      include"./incl/langue.inc"
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/gsmenu.inc"

      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
      INTEGER         MNNPEF(1:*)
      INTEGER, allocatable, dimension(:) :: LPTDVOI, LISTVOI

      IERR = 0
      IF( NCODSA .EQ. 0 ) RETURN

C     TABLEAUX LPTDVOI et LISTVOI NON ALLOUES
      IERALPTDVOI = 1
      IERALISTVOI = 1

C     NOMBRE DE COORDONNEES 2 3 ou 6
      NBCOOR = MCN(MNXYZN+WBCOON)
C     NOMBRE DE NOEUDS DE L'OBJET
      NBNOE = MCN(MNXYZN+WNBNOE)
C     NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET
      NTDL = NBNOE * NBDLNO

C     ALLOCATION DU TABLEAU LPTDVOI DES POINTEURS SUR LE DERNIER VOISIN
      NBNOE1 = 1 + NBNOE
      IF( LANGAG .EQ. 0 ) THEN
        PRINT*,'prgcmc: DEMANDE  ALLOCATION LPTDVOI(',NBNOE1,') ENTIERS'
         ALLOCATE ( LPTDVOI( 0:NBNOE ), STAT=IERALPTDVOI )
         IF( IERALPTDVOI .NE. 0 ) THEN
          PRINT*,'prgcmc: ERREUR ALLOCATION LPTDVOI(',NBNOE1,') ENTIERS'
            IERR = IERALPTDVOI
            GOTO 9990
         ENDIF
        PRINT*,'prgcmc: CORRECTE ALLOCATION LPTDVOI(',NBNOE1,') ENTIERS'
      ELSE
        PRINT*,'prgcmc: ALLOCATION DEMAND  LPTDVOI(',NBNOE1,') INTEGERS'
         ALLOCATE ( LPTDVOI( 0:NBNOE ), STAT=IERALPTDVOI )
         IF( IERALPTDVOI .NE. 0 ) THEN
          PRINT*,'prgcmc: ALLOCATION ERROR LPTDVOI(',NBNOE1,') INTEGERS'
            IERR = IERALPTDVOI
            GOTO 9990
         ENDIF
        PRINT*,'prgcmc: CORRECT ALLOCATION LPTDVOI(',NBNOE1,') INTEGERS'
      ENDIF

C     CREATION DE LA LISTE DU NOMBRE D'ELEMENTS FINIS DE CHAQUE NOEUD
      CALL NBEFNOEUD( MNTOPO, MNNPEF, NBNOE,
     %                LPTDVOI(1), MXNOE1EF, IERR )

C     NOINTE=1 : AXISYMETRIQUE DEGRE 1 <-> NOEUDS=SOMMETS
C            2 : AXISYMETRIQUE DEGRE 2 <-> NOEUDS=SOMMETS+MILIEUX
C            3 : LAGRANGE de DEGRE 1   <-> NOEUDS=SOMMETS
C            4 : LAGRANGE de DEGRE 2   <-> NOEUDS=SOMMETS+MILIEUX

C     ESTIMATION DU NOMBRE MOYEN DE NOEUDS VOISINS DANS CHAQUE EF
      NBNOVO = MXNOE1EF/2

C     ALLOCATION DU TABLEAU LISTVOI DES NOEUDS VOISINS
      NBAGMX = 0

 10   NBAGMX = NBAGMX + 1
      IF( NBAGMX .GT. 6 ) THEN
C        TROP D'AUGMENTATIONS DE LA VALEUR DE NBNOVO
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='prgcmc: TROP d''AUGMENTATIONS de NBNOVO'
            KERR(2)='POUR DECLARER LE TABLEAU DES VOISINS DES NOEUDS'
            KERR(3)='REDUIRE le MAILLAGE'
         ELSE
            KERR(1)='prgcmc: TOO AUGMENTATIONS of NBNOVO'
            KERR(2)='to DECLARE the ARRAY of NODE NEIGHBORS of NODES'
            KERR(3)='REDUCE the MESH'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9990
      ELSE
         NBNOVO = NBNOVO + 1
         PRINT*,'prgcmc: NB VOISINS 1 EF NBNOVO=',NBNOVO,
     %          'MAX NOEUDS 1EF=',MXNOE1EF
      ENDIF

      IF( NBAGMX .GT. 1 ) THEN
C        RECREER LA LISTE DU NOMBRE D'ELEMENTS FINIS DE CHAQUE NOEUD
         CALL NBEFNOEUD( MNTOPO, MNNPEF, NBNOE,
     %                   LPTDVOI(1), MXNOE1EF, IERR )
      ENDIF

C     ESTIMATION DU NOMBRE DE VOISINS DE CHAQUE NOEUD A PARTIR DU
C     NOMBRE D'EF DE CHAQUE NOEUD ET DU NOMBRE MOYEN DE VOISINS
      LPTDVOI(0) = 0
      DO K = 1, NBNOE
C        NOMBRE D'EF DU NOEUD K
         NBEF = LPTDVOI(K)
C        NOMBRE MOYEN DE VOISINS PAR EF
         IF( NBEF .LE. 6 ) THEN
            NBNV = MXNOE1EF - 1
         ELSE
            NBNV = NBNOVO
         ENDIF
C        ESTIMATION DU NOMBRE DE NOEUDS VOISINS DU NOEUD K
         LPTDVOI(K) = LPTDVOI(K-1) + NBNV * NBEF
      ENDDO
      MXLISTV = LPTDVOI(NBNOE)

      IF( MXLISTV .LE. 0 ) THEN
C        MXLISTV > 2**31-1
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='prgcmc: SATURATION DU MAX DES ENTIERS'
            KERR(2)='POUR DECLARER LE TABLEAU DES VOISINS DES NOEUDS'
            KERR(3)='REDUIRE le MAILLAGE'
         ELSE
            KERR(1)='prgcmc: SATURATION of INTEGER MAXIMUM'
            KERR(2)='to DECLARE the ARRAY of NODE NEIGHBORS of NODES'
            KERR(3)='REDUCE the MESH'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9990
      ENDIF

C     DECLARATION DU TABLEAU LISTVOI(MXLISTV) DES NOEUDS VOISINS
C     SI TAILLE<MAXENTIER=2**31-1
      IF( MXLISTV .LE. 0 .OR. MXLISTV .GE. MAXENTIER ) THEN
C         MXLISTV EST TROP GRAND POUR ETRE DECLARE
         PRINT*,'prgcmc: declaration de LISTVOI(',MXLISTV,')'
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: TAILLE DE LA LISTE DES VOISINS avec'
            KERR(2) = 'DEPASSEMENT DU MAX DES ENTIERS 2**31-1'
            KERR(3) = 'REDUIRE LA TAILLE DU MAILLAGE'
         ELSE
            KERR(1) = 'ERROR: NEIGHBORS LIST SIZE'
            KERR(2) = 'OVERFLOW of INTEGER MAX 2**31-1'
            KERR(3) = 'REDUCE THE SIZE OF THE MESH'
         ENDIF
         CALL LERESU
         IERR = 2
         GOTO 9990
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'DEMANDE  ALLOCATION LISTVOI(',MXLISTV,') ENTIERS'
         ALLOCATE ( LISTVOI( 1:MXLISTV ), STAT=IERALISTVOI )
         IF( IERALISTVOI .NE. 0 ) THEN
            PRINT*,'ERREUR ALLOCATION LISTVOI(',MXLISTV,') ENTIERS'
            IERR = IERALISTVOI
            GOTO 9990
         ENDIF
         PRINT*,'CORRECTE ALLOCATION LISTVOI(',MXLISTV,') ENTIERS'
      ELSE
         PRINT*,'ALLOCATION DEMAND  LISTVOI(',MXLISTV,') INTEGERS'
         ALLOCATE ( LISTVOI( 1:MXLISTV ), STAT=IERALISTVOI )
         IF( IERALISTVOI .NE. 0 ) THEN
            PRINT*,'ALLOCATION ERROR LISTVOI(',MXLISTV,') INTEGERS'
            IERR = IERALISTVOI
            GOTO 9990
         ENDIF
         PRINT*,'CORRECT ALLOCATION LISTVOI(',MXLISTV,') INTEGERS'
      ENDIF

      CALL LISTNOVO( MNTOPO, MNNPEF, NBNOE, MXLISTV,
     %               LPTDVOI(0), LISTVOI, MXVOIS, IERR )

      IF( IERR .EQ. 5 ) THEN
C        REALLOCATION DE LISTVOI DEMANDEE
         IF( IERALISTVOI .EQ. 0 ) then
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'prgcmc: LISTVOI est DESALLOUE'
            ELSE
               PRINT*,'prgcmc: LISTVOI is DESALLOCATED'
            ENDIF
            DEALLOCATE( LISTVOI )
            IERALISTVOI = 1
         ENDIF
         GOTO 10
      ENDIF


C     RESERVATION du TABLEAU POINTEUR SUR LES LIGNES LPLIGN
C     =====================================================
      IF( MNLPLI .EQ. 0 ) THEN
C        NTDL LE NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET
         CALL TNMCMX( 'ENTIER', MAXVAR )
         IF( MAXVAR .LT. NTDL+1 ) GOTO 9900
         CALL TNMCDC( 'ENTIER', NTDL+1, MNLPLI )
         IF( MNLPLI .LE. 0 ) GOTO 9900
      ENDIF

C     RESERVATION DU TABLEAU LPCOLO DU NUMERO DES COLONNES DE CHAQUE LIGNE
C     ====================================================================
C     MAJORATION DU NOMBRE DE COLONNES A STOCKER
C     EGAL AU NOMBRE DE COEFFICIENTS NON NUL DE LA MATRICE
      IF( NCODSA .GT. 0 ) THEN
C        MATRICE SYMETRIQUE
         LOLPCI = NBDLNO**2 * ( LPTDVOI(NBNOE) + NBNOE*2 + MXVOIS ) / 2
      ELSE
C        MATRICE NON SYMETRIQUE
         LOLPCI = NBDLNO**2 * ( LPTDVOI(NBNOE) + NBNOE )
      ENDIF
C     LE NO DES COLONNES POUR CHAQUE LIGNE
      PRINT*,'prgcmc: DECLARATION du TABLEAU du NUMERO des COLONNES de',
     %        LOLPCI,' entiers'
      CALL TNMCMX( 'ENTIER', MAXVAR )
      IF( MAXVAR .LT. LOLPCI ) GOTO 9900
      CALL TNMCDC( 'ENTIER', LOLPCI, MNLPCO )
      IF( MNLPCO .LE. 0 ) GOTO 9900


C     FORMATION DES TABLEAUX LPLIGN ET LPCOLO
C     =======================================
      MCN(MNLPLI) = 0
      IA  = 0
      IAN = 0
      IF( NCODSA .GT. 0 ) THEN

C        MATRICE SYMETRIQUE
C        ------------------
         MN1 = 1
         DO 30 I=1,NBNOE
            IANVOI = 0
            MN2 = LPTDVOI(I)

C           CALCUL DU NOMBRE DE VOISINS DU NOEUD I DE NUMERO<I
            NBV = 0
            DO MNV=MN1,MN2
C              NUMERO DU NOEUD VOISIN
               NUMV = LISTVOI(MNV)
               IF( NUMV .LT. I ) THEN
                  NBV = NBV + 1
               ENDIF
            ENDDO

C           PARCOURS DES NBV NOEUDS VOISINS DU NOEUD I
            DO MNV=MN1,MN1+NBV-1
C              NUMERO DU NOEUD VOISIN
ccc               NUMV = MCN(MNV)
               NUMV = LISTVOI(MNV)
C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD NUMV
               NUMIV = ( NUMV - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,NBDLNO
                     IANV = IA + JN + ( ( IN - 1 ) * IN ) / 2
     %                    + NBV * NBDLNO * ( IN - 1 )
     %                    + IANVOI
                     MCN(MNLPCO+IANV-1) = NUMIV + JN
                  ENDDO
               ENDDO
               IANVOI = IANVOI + NBDLNO
            ENDDO

C           LES NBDLNO DEGRES DE LIBERTE DU NOEUD I
            IAN  = IA
            NUMI = ( I - 1 ) * NBDLNO
            DO IN=1,NBDLNO
               DO JN=1,IN
                  IAN = IA + JN + ( ( IN - 1 ) * IN ) / 2
     %                + NBV * NBDLNO *  IN
                  MCN(MNLPCO+IAN-1) = NUMI + JN
               ENDDO
               MCN(MNLPLI+NUMI+IN) = IAN
            ENDDO

            IA  = IAN
            MN1 = MN2 + 1
 30      ENDDO

      ELSE

C        MATRICE NON SYMETRIQUE
C        ----------------------
         MN1 = 1
         DO I=1,NBNOE
            IANVOI = 0
            MN2 = LPTDVOI(I)
            NBV = MN2 - MN1 + 1
            DO MNV=MN1,MN2
               NUMV = LISTVOI(MNV)
C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD NUMV
               NUMIV = ( NUMV - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,NBDLNO
                     IANV = IA + JN + IANVOI
     %               + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     MCN(MNLPCO+IANV-1) = NUMIV + JN
                  ENDDO
               ENDDO
               IANVOI = IANVOI + NBDLNO
            ENDDO
C           LES NBDLNO DEGRES DE LIBERTE DU NOEUD I
            NUMI = ( I - 1 ) * NBDLNO
            DO IN=1,NBDLNO
               KN = 0
               DO JN=1,NBDLNO
                  IF (JN.NE.IN) THEN
                     KN = KN + 1
C                    PEUT ETRE UN PROBLEME SUR MCN(MNLPLI+I) QUI SUIT ???
                     IAN = IA + KN + MCN(MNLPLI+I) * NBDLNO
     %                + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     MCN(MNLPCO+IAN-1) = NUMI + JN
                  ENDIF
               ENDDO
               IAN = IA + ( NBV + 1 ) * NBDLNO *  IN
               MCN(MNLPCO+IAN-1) = NUMI + IN
               MCN(MNLPLI+NUMI+IN) = IAN
            ENDDO
            IA  = IAN
            MN1 = MN2 + 1
         ENDDO
      ENDIF

C     NOMBRE DE COEFFICIENTS DE LA MATRICE CONDENSEE
C     ----------------------------------------------
      LOLPCO = MCN(MNLPLI+NTDL)
      IF( LOLPCI .LT. LOLPCO ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'prgcmc: VALEUR DE LOLPCI A AUGMENTER'
         ELSE
            KERR(1) = 'prgcmc: AUGMENT THE LOLPCI VALUE'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9900
      ENDIF

C     REDUCTION DE LA TAILLE DU TABLEAU DES NUMEROS DE COLONNES
      CALL TNMCRA( 'ENTIER', LOLPCI, LOLPCO, MNLPCO )
      GOTO 9990

C     SATURATION DU SUPER-TABLEAU MCN
C     -------------------------------
 9900 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='prgcmc: PAS ASSEZ DE MEMOIRE MCN MOTMCN=             '
         WRITE(KERR(1)(47:61),'(I15)') MOTMCN
      ELSE
         KERR(1) = 'prgcmc: NOT ENOUGH MEMORY MCN MOTMCN=              '
         WRITE(KERR(1)(44:58),'(I15)') MOTMCN
      ENDIF
      CALL LEREUR
      IERR = 20

C     DESTRUCTION DES TABLEAUX AUXILIAIRES POUR LES VOISINS
C     -----------------------------------------------------
 9990 IF( IERALISTVOI .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'prgcmc: LISTVOI est DESALLOUE'
         ELSE
            PRINT*,'prgcmc: LISTVOI is DESALLOCATED'
         ENDIF
         DEALLOCATE( LISTVOI )
         IERALISTVOI = 1
      ENDIF

      IF( IERALPTDVOI .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'prgcmc: LPTDVOI est DESALLOUE'
         ELSE
            PRINT*,'prgcmc: LPTDVOI is DESALLOCATED'
         ENDIF
         DEALLOCATE( LPTDVOI )
         IERALPTDVOI = 1
      ENDIF

      RETURN
      END
