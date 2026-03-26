      SUBROUTINE CALVVP( NOMPB,  NTLXOB, NTDL,
     %                   METVVP, EIGMIN, EIGMAX, NBROOT,
     %                   NBDLFX, MNDLFX, MNVDLX,
     %                   NCODSK, MNMUKG, KG,
     %                   NCODSM, MNMUMG, MG,
     %                   NTVVPR, MNVVPR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES VALEURS EIGV ET VECTEURS PROPRES VP TELS QUE
C -----    ( K + EIGV M ) VP = 0   ET   tM VP M = Id
C
C ENTREES:
C --------
C NOMPB  : NOM EN 4 CARACTERES DE LA CLASSE DU PROBLEME TRAITE
C          'THER' ou 'ELAS' ou 'FLUI'
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET
C NTDL   : NOMBRE DE LIGNES ET COLONNES DES MATRICES K ET M EN ENTREE
C          NOMBRE TOTAL DE DL BLOQUES+LIBRES
C
C METVVP : NUMERO DE LA METHODE DE CALCUL DES VALEURS ET VECTEURS PROPRES
C          =0 A DEFINIR PAR LECTURE DES DONNEES DE L'UTILISATEUR
C          =1 METHODE DES SOUS-ESPACES
C             EIGMIN et EIGMAX DOIVENT ETRE INITIALISES
C             EIGMIN BORNE INFERIEURE DES VALEURS PROPRES A CALCULER
C             EIGMAX BORNE SUPERIEURE DES VALEURS PROPRES A CALCULER
C          =2 METHODE DE L'ITERATION INVERSE
C             NBROOT et eventuellement EIGMAX DOIVENT ETRE INITIALISES
C             NBROOT>0 NOMBRE DES PLUS PETITES VALEURS PROPRES A CALCULER
C             NBROOT=0 EIGMAX BORNE SUPERIEURE DES VALEURS PROPRES A CALCULER
C
C NBDLFX : NOMBRE DE DEGRES DE LIBERTE BLOQUES
C MNDLFX : ADRESSE MCN DES NUMEROS DE 1 A NTDL DES NODLFX DEGRES FIXES
C          DE LIBERTE BLOQUES
C MNVDLX : ADRESSE MCN DES VALEURS DES NODLFX DEGRES FIXES
C
C NCODSK : CODE DE STOCKAGE DE LA MATRICE PROFIL K
C          =0 SI DIAGONALE
C          >0 SI PROFIL SYMETRIQUE
c          <0 SI PROFIL NON SYMETRIQUE
C MNMUK  : ADRESSE MCN DU POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE K
C KG     : MATRICE PROFIL REELLE DOUBLE PRECISION de RAIDEUR
C
C NCODSM : CODE DE STOCKAGE DE LA MATRICE PROFIL M
C          =0 SI DIAGONALE
C          >0 SI PROFIL SYMETRIQUE
C          <0 SI PROFIL NON SYMETRIQUE
C MNMUMG : ADRESSE MCN DU POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE M
C MG     : MATRICE PROFIL REELLE DOUBLE PRECISION de MASSE
C
C SORTIES:
C --------
C NTVVPR : NUMERO      DU TMS DES VALEURS ET VECTEURS PROPRES
C MNVVPR : ADRESSE MCN DU TMS DES VALEURS ET VECTEURS PROPRES
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET      ANALYSE NUMERIQUE UPMC PARIS    AOUT 1998
C23456---------------------------------------------------------------012
      PARAMETER (MXCN0=10)
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/a___vecteur.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  EIGMIN, EIGMAX
      DOUBLE PRECISION  TETA, PI2, VECT(MXCN0)
      DOUBLE PRECISION  COMXKG, COMXMG, COMXPR
      INTEGER           NUM(MXCN0)
      CHARACTER*(*)     NOMPB
C
      DOUBLE PRECISION  MG(*), KG(*)
      DOUBLE PRECISION, allocatable, dimension(:) :: RG
      INTEGER           IERRGALLOC
      INTRINSIC         ALLOCATED
c
cccc-------------------------------------------------------------------------
cccc     premieres valeurs propres du 3cube unite (k**2+m**2+n**2)
cccc     cf $MEFISTO/prpr/picarre.f
ccc      double precision vp3cub1(35), picarre
ccc      picarre = (atan(1d0) * 4d0) ** 2
ccc      vp3cub1(1 ) =  3 * picarre
ccc      vp3cub1(11) = 12 * picarre
ccc      do 2 i=1,3
ccc         vp3cub1(i+ 1) =  6 * picarre
ccc         vp3cub1(i+ 4) =  9 * picarre
ccc         vp3cub1(i+ 7) = 11 * picarre
ccc         vp3cub1(i+11) = 14 * picarre
ccc         vp3cub1(i+14) = 14 * picarre
ccc         vp3cub1(i+17) = 17 * picarre
ccc         vp3cub1(i+20) = 18 * picarre
ccc         vp3cub1(i+23) = 19 * picarre
ccc         vp3cub1(i+26) = 21 * picarre
ccc         vp3cub1(i+29) = 21 * picarre
ccc         vp3cub1(i+32) = 22 * picarre
ccc 2    continue
cccc-------------------------------------------------------------------------
cccc-------------------------------------------------------------------------
cccc     22 premieres valeurs propres du 6cube unite
ccc      double precision vp6cub1(22), picarre
ccc      picarre = (atan(1d0) * 4d0) ** 2
ccc      vp6cub1(1) = 6 * picarre
ccc      do 2 i=2,7
ccc         vp6cub1(i) = 9 * picarre
ccc 2    continue
ccc      do 4 i=8,22
ccc         vp6cub1(i) = 12 * picarre
ccc 4    continue
cccc------------------------------------------------------------------------
c
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
      MNDLIB = 0
      IERRGALLOC = 1
C
C     PREPARATION DES MATRICES K ET M SELON LES DL FIXES A ZERO
C     ---------------------------------------------------------
C     RENUMEROTATION DES DEGRES DE LIBERTE LIBRES NODLIB
C     NODLIB(I) = NUMERO DU DEGRE DE LIBERTE S'IL EST LIBRE
C                 -INDICE DANS LA LISTE DES DL BLOQUES S'IL EST BLOQUE
C     NBDLIB    = NOMBRE TOTAL DE DEGRES DE LIBERTE LIBRES
      CALL REDLIB( NBDLFX, MCN(MNDLFX), NTDL,  NBDLIB, MNDLIB )
C
C     SUPPRESSION DES LIGNES ET COLONNES DES DL BLOQUES
C     DES MATRICES KG ET MG
      IF( NBDLIB .LT. NTDL ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) NBDLFX, NTDL
         ELSE
            WRITE(IMPRIM,20000) NBDLFX, NTDL
         ENDIF
         NBRDAG1 = MCN( MNMUKG + NTDL )
         CALL PCDLID( NTDL,MCN(MNDLIB),NCODSK,MCN(MNMUKG),KG )
         NBRDAG2 = MCN( MNMUKG + NBDLIB )
         CALL PCDLID( NTDL,MCN(MNDLIB),NCODSM,MCN(MNMUMG),MG )
         NBRDAG3 = MCN( MNMUMG + NBDLIB )
      ENDIF
10000 FORMAT(/'SUPPRESSION DES',I10,' LIGNES et COLONNES des ',
     %' DL FIXES a ZERO sur [K] et [M] parmi les',I10,' DL')
20000 FORMAT(/'SUPPRESSION of',I10,' LINES and COLUMNS of ',
     %' DoF FIXED at NULL on [K] and [M] among the',I10,' DoF')
C
C
C
C ====================================================================
cccC AJOUT POUR Weichung WANG et Tsung-Minh HWANG september 24-th 2006
cccC     MISE SUR FICHIER DES MATRICES PROFIL KG et MG
ccc      CALL FILESKYLINE( 'KG', NBDLIB, MCN(MNMUKG), KG, IERR )
ccc      CALL FILESKYLINE( 'MG', NBDLIB, MCN(MNMUMG), MG, IERR )
C ====================================================================
C
C
C
C     CALCUL DU COEFFICIENT MAXIMAL DES 2 MATRICES PROFILS REDUITES
C     -------------------------------------------------------------
      CALL PRCOMX( NBDLIB, NCODSK, MCN(MNMUKG), KG, COMXKG )
      CALL PRCOMX( NBDLIB, NCODSM, MCN(MNMUMG), MG, COMXMG )
      COMXPR = MAX( COMXKG, COMXMG )
C
cccC     MISE A L'ECHELLE DES 2 MATRICES
cccC     -------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10001) 'K',COMXKG, 'M',COMXMG
ccc         WRITE(IMPRIM,10002) COMXPR
ccc         WRITE(IMPRIM,10001) 'K',COMXKG / COMXPR, 'M',COMXMG / COMXPR
      ELSE
         WRITE(IMPRIM,20001) 'K',COMXKG, 'M',COMXMG
ccc         WRITE(IMPRIM,20002) COMXPR
ccc         WRITE(IMPRIM,20001) 'K',COMXKG / COMXPR, 'M',COMXMG / COMXPR
      ENDIF
10001 FORMAT('COEFFICIENT MAXIMAL de la MATRICE [',A,'] =',G15.7)
20001 FORMAT('MAXIMUM COEFFICIENT of the MATRIX [',A,'] =',G15.7)
ccc10002 FORMAT('DIVISION des 2 MATRICES [K]&[M] par ',G15.7, '=>' )
ccc20002 FORMAT('DIVISION of 2 MATRICES  [K]&[M] by ',G15.7, '=>' )
C
ccc      D = 1D0 / COMXPR     ccc le 20 mars 2009
ccc      CALL CTEPRO( D, NBDLIB, NCODSK, MCN(MNMUKG), KG )
ccc      CALL CTEPRO( D, NBDLIB, NCODSM, MCN(MNMUMG), MG )
ccc      COMXMG = COMXMG / COMXPR
ccc      COMXKG = COMXKG / COMXPR
CCCC
CCCC     COEFFICIENT DE TRANSLATION TETA TEL QUE K <-- K + TETA M
CCCC     --------------------------------------------------------
CCC      CALL INVITE( 96 )
CCC      NCVALS = 0
CCC      CALL LIRRDP( NCVALS, TETA )
CCC      IF( NCVALS .LT. 0 ) GOTO 9900
CCCC
CCCC     EVENTUELLE TRANSLATION DE K EN K+TETA*M
CCCC     ---------------------------------------
CCC      IF( TETA .NE. 0D0 ) THEN
CCC         CALL MUA2PD( NBDLIB, 1D0, NCODSK, MCN(MNMUKG), KG,
CCC     %                       TETA, NCODSM, MCN(MNMUMG), MG,
CCC     %                                     MCN(MNMUKG), KG )
CCC      ENDIF
      TETA = 0D0
C
C     SAUVEGARDE DANS LE TABLEAU RG DE LA MATRICE KG
C     APRES PRISE EN COMPTE DES DL FIXES A ZERO
C     ----------------------------------------------
      NBRDRG = MCN( MNMUKG + NBDLIB )
C
C     ALLOCATION DYNAMIQUE EN FORTRAN 90 DE LA MATRICE PROFIL RG
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'ALLOCATION DEMAND of',NBRDRG,
     %                ' DOUBLE PRECISION of the [RG] MATRIX'
      ALLOCATE ( RG(1:NBRDRG), STAT=IERRGALLOC )
      IF( IERRGALLOC .NE. 0 ) THEN
       WRITE(IMPRIM,*)'ALLOCATION ERROR  of',NBRDRG,
     %                ' DOUBLE PRECISION of the [MG] MATRIX'
         IERR = IERRGALLOC
         GOTO 9900
      ENDIF
      WRITE(IMPRIM,*) 'ALLOCATION DONE   of',NBRDRG,
     %                ' DOUBLE PRECISION of the [RG] MATRIX'
C
C     COPIE DE SAUVEGARDE DE KG DANS RG
      MOTSRG = NBRDRG * MOREE2
      CALL TRTATA( KG, RG, MOTSRG )
C
C     CHOIX DE LA METHODE DE CALCUL DES MODES PROPRES
C     -----------------------------------------------
      IF( METVVP .LE. 0 ) THEN
         INIDONN = 0
         CALL LIMTCL( 'methvvpr' , METVVP )
      ELSE
         INIDONN = 1
      ENDIF
      IF( METVVP .LE. 0 ) GOTO 9900
      IF( METVVP .EQ. 1 ) THEN
C
C        METHODE DES ITERATIONS PAR SOUS-ESPACES
C        ***************************************
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10032)
         ELSE
            WRITE(IMPRIM,20032)
         ENDIF
10032    FORMAT(/'METHODE DES SOUS-ESPACES')
20032    FORMAT(/'SUB-SPACE METHOD')
         CALL SSPACE( NOMPB,  NTLXOB, NTDL, INIDONN, EIGMIN, EIGMAX,
     %                NBDLFX, MNDLFX, MNVDLX, NBDLIB, MNDLIB,
     %                NCODSK, MNMUKG, KG,
     %                NCODSM, MNMUMG, MG, COMXMG,
     %                MOTSRG, RG,
     %                NTVVPR, MNVVPR, IERR )
C
         IF( IERR .EQ. 3 ) THEN
         NBLGRC(NRERR) = 5
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MAXIMUM ATTEINT des ITERATIONS'
            KERR(2) = 'SORTIE avec des VALEURS PROPRES NON CONVERGEES'
            KERR(3) = 'VOIR LES VALEURS PROPRES CALCULEES'
            KERR(4) = 'CHOISIR UNE BORNE PLUS PROCHE (mais PLUS PETITE)
     %DU MINIMUM DES VALEURS PROPRES'
            KERR(5) = 'ET      UNE BORNE PLUS GRANDE DU MAXIMUM DES VALE
     %URS PROPRES'
         ELSE
            KERR(1) = 'MAXIMUM ITERATIONS REACHED'
            KERR(2) = 'EXIT with NOT CONVERGED EIGENVALUES'
            KERR(3) = 'SEE the COMPUTED EIGENVALUES'
            KERR(4) = 'CHOOSE A NEARER (but SMALLER) BOUND of EIGENVALUE
     % MINIMUM'
            KERR(5) = 'CHOOSE A GREATER BOUND of EIGENVALUE MAXIMUM'
         ENDIF
         CALL LEREUR
         ENDIF
C
      ELSE
C
C        METHODE DES ITERATIONS INVERSES AVEC TRANSLATION
C        ************************************************
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10033)
         ELSE
            WRITE(IMPRIM,20033)
         ENDIF
10033    FORMAT(/'METHODE de l''ITERATION INVERSE avec TRANSLATIONS')
20033    FORMAT(/'INVERSE ITERATION METHOD with TRANSLATIONS')
         CALL ITEINV( NOMPB,  NTLXOB, NTDL, INIDONN, NBROOT, EIGMAX,
     %                NBDLFX, MNVDLX, NBDLIB, MNDLIB,
     %                NCODSK, MNMUKG, KG,
     %                NCODSM, MNMUMG, MG,
     %                TETA,   MOTSRG, RG,
     %                NTVVPR, MNVVPR, IERR )
C
      ENDIF
C
C     GESTION DE L'ERREUR
      IF( IERR .NE. 0 .AND. IERR .NE. 3 ) GOTO 9900
C
C     AFFICHAGE DES VALEURS PROPRES ET DES COMPOSANTES DES VECTEURS PROPRES
C     =====================================================================
C     FREQUENCE = SQRT( VALEUR PROPRE ) / (2 PI) EN ELASTICITE
C     LE NOMBRE DE VALEURS PROPRES ET DE DL
      NBROOT = MCN(MNVVPR+WBVECT)
      NTDL   = MCN(MNVVPR+WBCOVE)
      IF( IERR .NE. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,13000) NBROOT
         ELSE
            WRITE(IMPRIM,23000) NBROOT
         ENDIF
13000 FORMAT(/' Les ',I6,' VALEURS PROPRES ou FREQUENCES PROPRES(=SQRT(V
     %P)/2Pi):',/1X,80(1H=))
23000 FORMAT(/' The ',I6,' EIGENVALUES or FREQUENCIES(=SQRT(E)/2Pi :',
     %/1X,80(1H=))
C
         MNR = MNVVPR + WECTEU + MOREE2*NTDL*NBROOT - 1
C        2 Pi
         PI2 = ATAN(1D0) * 8D0
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,13001) (I,RMCN(MNR+I),
     %           SQRT(ABS(RMCN(MNR+I)))/PI2,I=1,NBROOT)
         ELSE
            WRITE(IMPRIM,23001) (I,RMCN(MNR+I),
     %           SQRT(ABS(RMCN(MNR+I)))/PI2,I=1,NBROOT)
         ENDIF
13001 FORMAT(' VALEUR PROPRE',I5,' =',G15.7,'  FREQUENCE =',G15.7,' Hz')
23001 FORMAT(' EIGENVALUE',I5,' =',G15.7,'  FREQUENCY =',G15.7,' Hz')
C
C
cccC        AFFICHAGE AVEC DES VALEURS PROPRES CONNUES CF Goong CHEN
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            WRITE(IMPRIM,13002) (I, vp3cub1(i), RMCN(MNR+I),
ccc     %           100*abs(vp3cub1(i)-RMCN(MNR+I))/vp3cub1(i),I=1,NBROOT )
ccc         ELSE
ccc            WRITE(IMPRIM,23002) (I, vp3cub1(i), RMCN(MNR+I),
ccc     %           100*abs(vp3cub1(i)-RMCN(MNR+I))/vp3cub1(i),I=1,NBROOT )
ccc         ENDIF
ccc13002    FORMAT(' VALEUR PROPRE EXACTE',I5,'=',G15.6,' VP CALCULEE=',
ccc     %G15.6,  ' RELATIVE ERROR=',G12.3,' %')
ccc
ccc23002 FORMAT(' EXACT EIGENVALUE',I5,' =',G15.6,' COMPUTED EIGENVALUE=',
ccc     %G15.6, ' RELATIVE ERROR=',G12.3,' %')
      ENDIF
C
cccccC     K/c - VAL PROPRE M/c = 0  AVEC tV M/c V = 1
cccccC     VV = V / racine carree(c) POUR AVOIR
cccccC     tVV M VV = tV/sqrt(c) M V/sqrt(c) = 1
ccccc      D = SQRT( 1D0 / COMXPR )
cccccC
cccccC     MISE DANS VVPR DES VECTEURS PROPRES VV = V / racine carree(c)
ccccc
C
C     MISE DANS VVPR DES VECTEURS PROPRES
      MND = ( MNVVPR + WECTEU - 1 ) / MOREE2
      MXCVEC = MIN(MXCN0,NTDL)
      DO 3010 I=1,NBROOT
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,13009) MXCVEC,NTDL,I
         ELSE
            WRITE(IMPRIM,23009) MXCVEC,NTDL,I
         ENDIF
         M = 0
         DO 3005 K=1,NTDL
cccC           REMISE A L'ECHELLE DES COMPOSANTES DES VECTEURS PROPRES
ccc            DMCN(MND+K) = DMCN(MND+K) * D
C           REMPLISSAGE DES MXCVEC PREMIERES COMPOSANTES NON NULLES
            IF( M .LT. MXCVEC .AND. DMCN(MND+K) .NE. 0D0 ) THEN
               M       = M + 1
               NUM(M)  = K
               VECT(M) = DMCN(MND+K)
            ENDIF
 3005    CONTINUE
C        AFFICHAGE EFFECTIF
         WRITE(IMPRIM,13010) (NUM(K),VECT(K),K=1,MXCVEC)
         MND = MND + NTDL
 3010 CONTINUE
13009 FORMAT(/I6,' PREMIERES COMPOSANTES NON NULLES DES',I8,
     %' COMPOSANTES DU VECTEUR PROPRE',I4)
23009 FORMAT(/I6,' FIRST NO NULL COMPONENTS of',I8,
     %' COMPONENTS of the EIGENVECTOR',I4)
13010 FORMAT(5(I6,':',G14.6))
C
C     SUPPRESSION DES TABLEAUX MC DEVENUS INUTILES
 9900 IF(MNDLIB .GT. 0) CALL TNMCDS( 'ENTIER', NTDL, MNDLIB )
      IF(IERRGALLOC .EQ. 0) DEALLOCATE( RG )
C
      RETURN
      END
