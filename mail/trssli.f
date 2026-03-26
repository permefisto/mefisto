      SUBROUTINE TRSSLI( KNOM, NULIGN )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES SOMMETS SIMPLES D'UNE LIGNE NON FERMEE
C -----    ET NON STRUCTUREE
C
C ENTREES:
C --------
C KNOM   : LE NOM DE LA LIGNE A TRAITER
C
C SORTIE :
C --------
C NULIGN : NUMERO DE LA LIGNE DANS LE LEXIQUE DES LIGNES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        MARS 1997
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOM
      CHARACTER*10      NMSOMM
C
C     LE LEXIQUE DE LA LIGNE
      CALL NUOBNM( 'LIGNE', KNOM, NULIGN )
      CALL LXNLOU( NTLIGN, NULIGN, NTLXOB , MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LIGNE INCONNUE:' // KNOM
         ELSE
            KERR(1) = 'UNKNOWN LINE:' // KNOM
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU 'SOAR' DE LA LIGNE
      CALL LXTSOU( NTLXOB , 'NSEF' , NTTSMA , MNTSMA )
      IF( NTTSMA .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LIGNE NON MAILLEE ' // KNOM
         ELSE
            KERR(1) = 'NOT MESHED LINE' // KNOM
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU 'XYZSOMMET' DE LA LIGNE
      CALL LXTSOU( NTLXOB , 'XYZSOMMET' , NTSOMM , MNSOMM )
      IF( NTSOMM .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LIGNE SANS XYZSOMMET '// KNOM
         CALL LEREUR
         RETURN
      ENDIF
      NBSOM = MCN( MNSOMM + WNBSOM )
C
C     LE TYPE DE FERMETURE DU MAILLAGE DE LA LIGNE
      NUTFMA = MCN( MNTSMA + WUTFMA )
      IF( NUTFMA .EQ. 1 ) THEN
C        LA LIGNE A DEJA ETE TESTEE COMME ETANT FERMEE
C        => PAS D'ARETES SIMPLES
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PAS DE SOMMET SIMPLE'
         ELSE
            KERR(1) = 'NO SIMPLE VERTEX'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TYPE DU MAILLAGE
      NUTYMA = MCN( MNTSMA + WUTYMA )
      IF( NUTYMA .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SI LE MAILLAGE EST STRUCTURE =>'
            KERR(2) = 'PAS DE TRACE DES SOMMETS SIMPLES'
         ELSE
            KERR(1) = 'STRUCTURED MESHED =>'
            KERR(2) = 'NO DRAWING of SIMPLE VERTICES'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     NOMBRE DE SOAR
      NBEFOB = MCN( MNTSMA + WBEFOB )
C     NOMBRE DE SOMMETS ET TANGENTES PAR EF
      NBSOEF = MCN( MNTSMA + WBSOEF )
C     ADRESSE DU NUMERO DU 1-ER SOMMET DE LA 1-ERE FACE
      MNSS   = MNTSMA + WUSOEF
C
C     CONSTRUCTION DU TABLEAU DES SOMMETS DES ARETES DE LA LIGNE
      MNSOAR = 0
      CALL GESOAR( RMCN(MNSOMM+WYZSOM),
     %             NBSOEF , NBEFOB , MCN(MNSS) , 2 ,
     %             L1SOAR , L2SOAR , MNSOAR , IERR )
      IF( IERR .LT. 0 ) RETURN
C
C     LA PILE DES SOMMETS DES SOMMETS SIMPLES
      MNNUST = 0
      CALL TNMCDC( 'ENTIER', NBSOM, MNNUST )
      CALL AZEROI( NBSOM, MCN(MNNUST) )
      MNNUST = MNNUST - 1
C
C     SI LE TABLEAU TRACE DE CET OBJET EXISTE ALORS LES OPTIONS DE TRACE
C     SONT INITIALISEES AVEC SES VALEURS
      CALL COUOBJ( 2, KNOM, MNTRAC )
C
C     LA VISEE DE L'OBJET SELON LES PARAMETRES DU COMMON / TRVARI /
      CALL VISEE1
C
C     RECHERCHE DES SOMMETS APPARTENANT A UNE SEULE ARETE
      MN = MNSOAR
      DO 90 JS = 1, L2SOAR
         IF( MCN(MN) .NE. 0 ) THEN
C           LE SOMMET EST INITIALISE
            IF( MCN(MN+3) .EQ. 0 ) THEN
C
C              LE SOMMET APPARTIENT A UNE SEULE ARETE <=> IL EST SIMPLE
               IERR = IERR + 1
C
C              LE NUMERO DU SOMMET
               NOS  = MCN(MN)
C
C              AFFICHAGE DU SOMMET SIMPLE
               MNSS = MNSOMM + WYZSOM + 3 * NOS - 3
               WRITE(IMPRIM,10090) NOS,RMCN(MNSS),RMCN(MNSS+1),
     %                             RMCN(MNSS+2)
C              LE SOMMET EST EMPILE
               MCN(MNNUST+IERR) = NOS
C
C              TRACE DU SOMMET SIMPLE ET DE SON NUMERO
               IF( NDIMLI .EQ. 2 ) THEN
                  WRITE( NMSOMM, '(I9)' ) NOS
                  L = NUDCNB( NMSOMM )
                  NMSOMM = '*' // NMSOMM(1:L)
                  CALL SANSBL( NMSOMM, L )
                  CALL SYMBOLE2D( NCBLAN, RMCN(MNSS), RMCN(MNSS+1),
     %                            NMSOMM(1:L) )
               ELSE
                  WRITE( NMSOMM, '(I9)' ) NOS
                  CALL SANSBL( NMSOMM, L )
                  CALL SYMBOLE3D( NCBLAN, RMCN(MNSS), NMSOMM(1:L) )
               ENDIF
C
            ENDIF
         ENDIF
         MN = MN + L1SOAR
 90   CONTINUE
C
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'LIGNE AVEC ',IERR,' SOMMETS SIMPLES'
         WRITE(IMPRIM,*) 'CES SOMMETS SONT PEUT ETRE IDENTIFIABLES'
         WRITE(IMPRIM,*) 'CF OPTION 20 DU MENU DEBUT'
      ELSE
         WRITE(IMPRIM,*) 'LINE WITH ',IERR,' SIMPLE VERTICES'
         WRITE(IMPRIM,*) 'THESE VERTICES MAY BE IDENTIFIED'
         WRITE(IMPRIM,*) 'CF OPTION 20 of DEBUT MENU'
      ENDIF
C
C     RECHERCHE POUR CHAQUE SOMMET DU SOMMET LE PLUS PROCHE
      DO 140 J=1,IERR
         IMIN = 0
         DMIN = 1E28
         NS1  = MCN(MNNUST+J)
         MNS1 = MNSOMM + WYZSOM + 3 * NS1 - 3
         DO 130 I=1,IERR
            IF( I .NE. J ) THEN
C              CALCUL DE LA DISTANCE
               MNSS = MNSOMM + WYZSOM + 3 * MCN(MNNUST+I) - 3
               D    = ( RMCN(MNS1  ) - RMCN(MNSS  ) ) ** 2
     %              + ( RMCN(MNS1+1) - RMCN(MNSS+1) ) ** 2
     %              + ( RMCN(MNS1+2) - RMCN(MNSS+2) ) ** 2
               IF( D .LT. DMIN ) THEN
                  IMIN = I
                  DMIN = D
               ENDIF
            ENDIF
 130     CONTINUE
C        IMIN EST LE SOMMET LE PLUS PROCHE DU SOMMET J
         NOS  = MCN(MNNUST+IMIN)
         MNSS = MNSOMM + WYZSOM + 3 * NOS - 3
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10140) NS1, NOS, SQRT(DMIN)
         ELSE
            WRITE(IMPRIM,20140) NS1, NOS, SQRT(DMIN)
         ENDIF
10140 FORMAT(/' DISTANCE entre le SOMMET ',I6,
     %        ' et son PLUS PROCHE SOMMET ',I6,' = ',G14.6)
20140 FORMAT(/' DISTANCE between VERTEX ',I6,
     %        ' and its NEAREST VERTEX ',I6,' = ',G14.6)
         WRITE(IMPRIM,10090) NS1,RMCN(MNS1),RMCN(MNS1+1),RMCN(MNS1+2)
         WRITE(IMPRIM,10090) NOS,RMCN(MNSS),RMCN(MNSS+1),RMCN(MNSS+2)
 140  CONTINUE
10090 FORMAT(' SIMPLE SOMMET',I6,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
C
C     DESTRUCTION DE LA PILE
      MNNUST = MNNUST + 1
      CALL TNMCDS( 'ENTIER', NBSOM, MNNUST )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'INFORMATION: La LIGNE: ', KNOM
      ELSE
         WRITE(IMPRIM,*) 'INFORMATION: The LINE: ', KNOM
      ENDIF
      IF( IERR .EQ. 0 ) THEN
C        LIGNE FERMEE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)
     %' est une LIGNE FERMEE de R**3  .................................'
         ELSE
            WRITE(IMPRIM,*)
     %' is a CLOSED LINE of R**3  .....................................'
         ENDIF
         MCN( MNTSMA + WUTFMA ) = 1
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)
     %' N''est PAS une ligne FERMEE de R**3  ..........................'
         ELSE
            WRITE(IMPRIM,*)
     %' is NOT a CLOSED line of R**3  .................................'
         ENDIF
C        LIGNE NON-FERMEE
         MCN( MNTSMA + WUTFMA ) = 0
      ENDIF
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
      CALL TNMCDS( 'ENTIER' , L1SOAR * L2SOAR , MNSOAR )
      RETURN
      END
