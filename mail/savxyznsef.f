      SUBROUTINE SAVXYZNSEF( NUTYOB, KNMPLS, MNXYZS, MNNSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXPORTER SUR UN FICHIER de CARACTERES ASCII
C -----    les COORDONNEES DES SOMMETS et TANGENTES
C          les NUMEROS DES SOMMETS DU MAILLAGE D'UN PLSV
C          C-A-D les TMS XYZSOMMET et NSEF Mefisto

C ENTREES:
C --------
C NUTYOB : NUMERO du TYPE de PLSV=1:Point 2:Ligne 3:Surface 4:Volume
C KNMPLS : NOM DU PLSV DANS LE LEXIQUE DES PLSV
C MNXYZS : >0 ADRESSE MCN du TMS 'XYZSOMMET'
C MNNSEF : =0 ADRESSE MCN du TMS 'NSEF' pour un POINT
C          >0 ADRESSE MCN du TMS 'NSEF' pour une LIGNE ou SURFACE ou VOLUME

C SORTIE :
C --------
C IERR   : =0  le FICHIER xyznsef.PLSV.KNMPLS est ECRIT dans
C              le REPERTOIRE du PROJET
C          >0  ERREUR RENCONTREE dans le TMS 'NSEF'. PAS de FICHIER CREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1997
C MODIFS : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*10       NMTYOB
      CHARACTER*24       KNMPLS
      CHARACTER*34       KNMFIC
      CHARACTER*1        KSUFIX(4)
      INTEGER            NOSOEF(64)
      LOGICAL            LEXIST,LOPEN
      DATA  KSUFIX / 'p', 'l', 's', 'v' /


      IF( MNXYZS .LE. 0 .OR. ( NUTYOB.GE.2 .AND. MNNSEF.LE.0 )  ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *, 'Le fichier xyznsef N''EST PAS CREE'
         ELSE
            PRINT *, 'The file xyznsef IS NOT CREATED'
         ENDIF
         IERR = 1
         RETURN
      ENDIF

C     EXISTENCE OU NON DE TANGENTES
      NBTGS = MCN( MNXYZS + WNBTGS )
      IF( NBTGS .GT. 0 ) THEN
         NOTG = 1
      ELSE
         NOTG = 0
      ENDIF

C     DECLARATION OUVERTURE DU FICHIER xyznsef.PLSV.NOM_du_PLSV
C     =========================================================
C     LE NOM DU FICHIER a ECRIRE dans le REPERTOIRE du PROJET
      CALL MINUSC( KNMPLS )
      I = NUDCNB( KNMPLS )
      KNMFIC = 'xyznsef.' // KSUFIX(NUTYOB) // '.' // KNMPLS(1:I)

C     SI LE FICHIER KNMFIC EXISTE ALORS IL EST DETRUIT PUIS RECONSTRUIT
      INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( LEXIST ) THEN
C        LE FICHIER KNMFIC EXISTE
         IF( .NOT. LOPEN ) THEN
C           OUVERTURE DU FICHIER KNMFIC
            CALL TRUNIT( NF )
            OPEN( FILE=KNMFIC, UNIT=NF, STATUS='OLD' )
         ENDIF
C        DESTRUCTION DU FICHIER
         CLOSE( NF, STATUS='DELETE' )
      ENDIF

C     CREATION DU FICHIER KNMFIC
C     --------------------------
      CALL TRUNIT( NF )
      OPEN( UNIT=NF, ERR=9999, STATUS='NEW',
     %      FILE=KNMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )

C     ECRITURE DU NO et NOM du TYPE D'OBJET en PREMIERE LIGNE du FICHIER
C     ------------------------------------------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10050) NUTYOB, NMTYOB(NUTYOB)
10050    FORMAT('{ ',I2,'; ',A10,'; { Numero et Type du PLSV } }')
      ELSE
         WRITE(NF,10051) NUTYOB, NMTYOB(NUTYOB)
10051    FORMAT('{ ',I2,'; ',A10,'; { Number and PLSV Type } }')
      ENDIF

C     LE NOMBRE de COORDONNEES DES SOMMETS DU MAILLAGE
      NBCOOR = MCN( MNXYZS + WBCOOR )
C     LE NOMBRE DE SOMMETS DU MAILLAGE
      NBSOM = MCN( MNXYZS + WNBSOM )
C     LE NOMBRE DE TANGENTES DU MAILLAGE
      IF( NOTG .GT. 0 ) THEN
         NBTGS = MCN( MNXYZS + WNBTGS )
      ELSE
         NBTGS = 0
      ENDIF

C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS = MNXYZS + WYZSOM - 3

C     ECRITURE DU NOM DU PLSV en SECONDE LIGNE du FICHIER
C     ---------------------------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10059) KNMPLS
10059    FORMAT('{ ',A24,'; { NOM du PLSV } }')
      ELSE
         WRITE(NF,10060) KNMPLS
10060    FORMAT('{ ',A24,'; { PLSV''s NAME } }')
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN

C        ECRITURE DU NOMBRE DE SOMMETS DU PLSV en TROISIEME LIGNE du FICHIER
C        -------------------------------------------------------------------
         WRITE(NF,10061) NBSOM
10061    FORMAT( I8, ';   {NBSOM  NOMBRE de SOMMETS DU PLSV}' )
C        ECRITURE DE NBTGS
         WRITE(NF,10062) NBTGS
10062    FORMAT( I8, ';   {NBTGS  NOMBRE de TANGENTES DU PLSV}' )
C        ECRITURE DE NBCOOR
         WRITE(NF,10063) NBCOOR
10063    FORMAT( I8, ';   {NBCOOR NOMBRE de COORDONNEES D UN SOMMET}')
C
C        ECRITURE DE XYZ DES NBSOM SOMMETS
10065    FORMAT( 6(E15.7,'; ') )
C        ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
         MNS = MNXYZS + WYZSOM
         DO I=1,NBSOM
            DO K=0,NBCOOR-1
C              MISE A ZERO DE LA VALEUR ABSOLUE de XYZ < 1E-20
               IF( ABS( RMCN(MNS+K) ) .LE. 1E-20 ) RMCN(MNS+K)=0.
            ENDDO
            WRITE(NF,10065) (RMCN(MNS+K),K=0,NBCOOR-1)
            MNS = MNS + NBCOOR
         ENDDO

C        ECRITURE DE XYZ DES NBTGS TANGENTES
C        MNS ADRESSE DE LA 1-ERE COORDONNEE DE LA 1-ERE TG DU TMS 'XYZSOMMET'
         IF( NBTGS .GT. 0 ) THEN
            WRITE(NF,10066) (RMCN(MNS+K),K=0,2)
10066       FORMAT( 3(E15.7,'; '), '; { XYZ de la PREMIERE TANGENTE }' )
            MNS = MNS + 3
            DO I=2,NBTGS
               WRITE(NF,10065) (RMCN(MNS+K),K=0,2)
               MNS = MNS + 3
            ENDDO
         ENDIF

      ELSE

C        ECRITURE DU NOMBRE DE SOMMETS DU PLSV en TROISIEME LIGNE du FICHIER
C        -------------------------------------------------------------------
         WRITE(NF,10081) NBSOM
10081    FORMAT( I8, ';   {NBSOM PLSV''s VERTICES NUMBER}' )
C        ECRITURE DE NBTGS
         WRITE(NF,10082) NBTGS
10082    FORMAT( I8, ';   {NBTGS PLSV''s TANGENT VECTORS NUMBER}' )
         WRITE(NF,10083) NBCOOR
10083    FORMAT( I8, ';   {NBCOOR PLSV''s VERTICE COORDINATES NUMBER}' )

C        ECRITURE DE XYZ DES NBSOM SOMMETS
C        ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
         MNS = MNXYZS + WYZSOM
         DO I=1,NBSOM
            DO K=0,NBCOOR-1
C              MISE A ZERO DE LA VALEUR ABSOLUE de XYZ < 1E-20
               IF( ABS( RMCN(MNS+K) ) .LE. 1E-20 ) RMCN(MNS+K)=0.
            ENDDO
            WRITE(NF,10065) (RMCN(MNS+K),K=0,NBCOOR-1)
            MNS = MNS + NBCOOR
         ENDDO

C        ECRITURE DE XYZ DES NBTGS TANGENTES
C        MNS ADRESSE DE LA 1-ERE COORDONNEE DE LA 1-ERE TG DU TMS 'XYZSOMMET'
         IF( NBTGS .GT. 0 ) THEN
            WRITE(NF,10086) (RMCN(MNS+K),K=0,2)
10086       FORMAT( 3(E15.7,';'),' { XYZ of FIRST TANGENT VECTOR }')
            MNS = MNS + 3
            DO I=2,NBTGS
              WRITE(NF,10065) (RMCN(MNS+K),K=0,2)
              MNS = MNS + 3
           ENDDO
         ENDIF
      ENDIF
      IF( NUTYOB .EQ. 1 ) GOTO 9900


C     LA TOPOLOGIE : TMS NSEF NON STRUCTURE
C     ============== ----------------------
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN

C     ECRITURE DES GENERALITES DU TABLEAU NSEF
      IF( NOTG .GT. 0 ) THEN
C        LES TANGENTES SONT A TRITER
         NBTGEF = MCN(MNNSEF+WBTGEF)
         NBEFTG = MCN(MNNSEF+WBEFTG)
         NBEFAP = MCN(MNNSEF+WBEFAP)
      ELSE
         NBTGEF = 0
         NBEFTG = 0
         NBEFAP = 0
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
      WRITE(NF,10091) MCN(MNNSEF+WUTYOB)
10091 FORMAT( I8, '; {NUTYOB  type du PLSV(1 2 3 4)}')
      WRITE(NF,10092) MCN(MNNSEF+WUTFMA)
10092 FORMAT( I8, '; {NUTFMA  type FERME(1) ou NON(0) ou INCONNU(-1)}')
      WRITE(NF,10093) NBSOEF
10093 FORMAT( I8, '; {NBSOEF  nombre de sommets par EF}')
      WRITE(NF,10094) NBTGEF
10094 FORMAT( I8, '; {NBTGEF  nombre de tangentes par EF}')
      WRITE(NF,10095) NBEFOB
10095 FORMAT( I8, '; {NBEFOB  nombre des EF du PLSV}')
      WRITE(NF,10096) NBEFTG
10096 FORMAT( I8, '; {NBEFTG  nombre des EF avec TG}')
      WRITE(NF,10097) NBEFAP
10097 FORMAT( I8, '; {NBEFAP  nombre des EF avec POINTEUR sur EF a TG}')
      WRITE(NF,10098) MCN(MNNSEF+WUTYMA)
10098 FORMAT( I8, '; {NUTYMA  type structure(1,...,7) ou non(0)}')
      ELSE
      WRITE(NF,11091) MCN(MNNSEF+WUTYOB)
11091 FORMAT( I8, '; {NUTYOB  type of the PLSV(1 2 3 4)}')
      WRITE(NF,11092) MCN(MNNSEF+WUTFMA)
11092 FORMAT(I8, '; {NUTFMA  type: CLOSED(1) or NOT(0) or UNKNOWN(-1)}')
      WRITE(NF,11093) NBSOEF
11093 FORMAT( I8, '; {NBSOEF  number of vertices for 1 FE}')
      WRITE(NF,11094) NBTGEF
11094 FORMAT( I8, '; {NBTGEF  number of tangent vectors for 1 FE}')
      WRITE(NF,11095) NBEFOB
11095 FORMAT( I8, '; {NBEFOB  number of FE of the PLSV}')
      WRITE(NF,11096) NBEFTG
11096 FORMAT( I8, '; {NBEFTG  number of FE with tangent vectors}')
      WRITE(NF,11097) NBEFAP
11097 FORMAT( I8, '; {NBEFAP  number of FE with POINTER on FE with TG}')
      WRITE(NF,11098) MCN(MNNSEF+WUTYMA)
11098 FORMAT( I8, '; {NUTYMA  structured type(1,...,7) or NOT(0)}')
      ENDIF

C     LA BOUCLE SUR LES EF DU MAILLAGE
C     ================================
      DO 100 NEF = 1, NBEFOB

C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( NEF   , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEF, IERR )

C        ECRITURE DES NUMEROS DES SOMMETS DE L'EF
         WRITE(NF,11108) (NOSOEF(K),K=1,NBSOEF)

 100  ENDDO
11108 FORMAT(  8(I8,'; ') )
11110 FORMAT( 10(I8,'; ') )
C
C     ADRESSE DERRIERE LES NUMEROS DES SOMMETS DANS NSEF
      IF( NBTGS .GT. 0 ) THEN
         MNS = MNNSEF + LDAPEF - 1
         IF( NBEFAP .GT. 0 ) THEN
C           ECRITURE DES numero>0 de l'EF a TG sinon 0
            WRITE(NF,11110) (MCN(MNS+K),K=1,NBEFAP)
         ENDIF
         MNS = MNNSEF + LDNGEF - 1
         IF( NBEFTG .GT. 0 ) THEN
C           ECRITURE du numero geometrique des EF a TG
C          ( 0 : 'C1 degre 3' ,
C            1 : 'CERCLE',   2 : 'ELLIPSE' ,     3 : 'COURBE B-SPLINE' ,
C           11 : 'SPHERE',  12 : 'ELLIPSOIDE',  13 : 'CYLINDRE', 14 : 'CONE',
C           15 : 'TORE'  ,  16 : 'TRANSFINI' ,  17 : 'SURFACE B-SPLINE') ;
            WRITE(NF,11110) (MCN(MNS+K),K=1,NBEFTG)
            IF( NBTGEF .GT. 0 ) THEN
               MNS = MNNSEF + LDTGEF - 1
               DO I=1,NBEFTG
C                 ECRITURE du +-no des TANGENTES de l'EF a TG
                  WRITE(NF,11110) (MCN(MNS+K),K=1,NBTGEF)
                  MNS = MNS + NBTGEF
               ENDDO
            ENDIF
         ENDIF
      ENDIF

 9900 CLOSE( NF )

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'Le fichier ',KNMFIC,' EST ECRIT'
      ELSE
         PRINT *, 'The file ',KNMFIC,' is WRITTEN'
      ENDIF
      PRINT*

      RETURN

C     PB A L'OUVERTURE DU FICHIER
 9999 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'PB a l''OUVERTURE DU FICHIER xyznsef'
      ELSE
         KERR(1) = 'FILE xyznsef CAN NOT BE OPENED'
      ENDIF
      CALL LEREUR

      RETURN
      END
