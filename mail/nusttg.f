      SUBROUTINE NUSTTG( MNNSEF, MNXYZS,  MONUST, MNNUST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE TABLEAU DU NUMERO DE SOMMET DE CHAQUE TANGENTE
C -----    D'UN MAILLAGE DEFINI PAR LES 2 TMS XYZSOMMET ET NSEF
C
C ENTREES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES EF
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES EF
C          CF ~/td/d/a___nsef'
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DES XYZ DES SOMMETS
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DES XYZ DES SOMMETS
C          CF '~/td/d/a___xyzsommet'
C
C SORTIES:
C --------
C MONUST : NOMBRE D'ENTIERS DU TABLEAU NUST
C MNNUST : ADRESSE MCN DU TABLEAU NUST DES NUMEROS DE SOMMET DE
C          CHAQUE TANGENTE ( AVEC MEME LA TANGENTE 0 )
C          ATTENTION : LE TABLEAU EST A DETRUIRE ENSUITE DANS MCN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS         MAI 1996
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
      INTEGER          NOSOEL(64)
C
C     LE NOMBRE DE SOMMETS ET TANGENTES DU MAILLAGE
      NBSOM = MCN( MNXYZS + WNBSOM )
      NBTGS = MCN( MNXYZS + WNBTGS )
C
C     RESERVATION DU TABLEAU NUST DANS MCN
      MNNUST = 0
      MONUST = NBTGS + 1
      CALL TNMCDC( 'ENTIER', MONUST, MNNUST )
      CALL AZEROI( MONUST, MCN(MNNUST) )
C
C     LES CARACTERISTIQUES DU MAILLAGE NSEF
      IF( MNNSEF .LE. 0 ) RETURN
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     PARCOURS DES EF DU MAILLAGE
C     NOMBRE DE TG PAR SOMMET
      NBTGST = NBTGEF / NBSOEF
      DO 30 N=1,NBEFOB
C
C        LE NUMERO DES NBSOEF SOMMETS ET DES NBTGEF TGS DE L'EF N
         CALL NSEFNS( N, NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C
C        EN SORTIE NBSOEF NUMEROS DE SOMMETS SONT INITIALISES
C                  LES DERNIERS PEUVENT ETRE NULS
C                  PAR EXEMPLE LE QUATRIEME D'UN TRIANGLE
C                              LE CINQUIEME D'UN TETRAEDRE ...
C                  DE MEME POUR LES NUMEROS DES TANGENTES
C
C        BOUCLE SUR LES TANGENTES EN CHAQUE SOMMET DE L'EF
         L = 0
         DO 20 I=1,NBSOEF
C           LE NUMERO GLOBAL DU SOMMET I DE L'EF N
            NST = NOSOEL( I )
            IF( NST .GT. 0 ) THEN
               DO 10 K=1,NBTGST
C                 LE NUMERO DE LA TANGENTE K DU SOMMET NST
                  L   = L + 1
                  NTG = ABS( NOSOEL( NBSOEF + L ) )
C                 STOCKAGE DU NO DE SOMMET DE LA TANGENTE
                  MCN( MNNUST + NTG ) = NST
 10            CONTINUE
            ENDIF
 20      CONTINUE
 30   CONTINUE
C
C     LA TANGENTE NULLE N'A PAS DE NUMERO DE SOMMET
      MCN( MNNUST ) = 0
      END
