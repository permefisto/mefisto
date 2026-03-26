      SUBROUTINE NSEFPA( LTNSEF,
     %                   NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   NBEFOB, NX, NY, NZ,
     %                   IERR   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LE NOMBRE TOTAL DE SOMMETS DU MAILLAGE,
C -----     DE SOMMETS DE CHAQUE EF, DU NOMBRE D'EF,
C           LE NOMBRE D'ARETES PAR DIRECTION SI MAILLAGE STRUCTURE,
C           LES DECALAGES POUR ARRIVER SUR LES TABLEAUX LPEFTG, NGEFTG, NUTGEF
C
C ENTREES:
C --------
C LTNSEF : LE TABLEAU NSEF
C
C SORTIES:
C --------
C NUTYMA :  NUMERO DE TYPE DU MAILLAGE
C           -1 : 'NON STRUCTURE AVEC TG', 0 : 'NON STRUCTURE SANS TG '
C            2 : 'SEGMENT   STRUCTURE'  ,
C            3 : 'TRIANGLE  STRUCTURE' , 4 : 'QUADRANGLE STRUCTURE' ,
C            5 : 'TETRAEDRE STRUCTURE' , 6 : 'PENTAEDRE  STRUCTURE' ,
C            7 : 'HEXAEDRE  STRUCTURE'
C NBSOEL :  NOMBRE REEL DE SOMMETS DE CHAQUE EF SI MAILLAGE STRUCTURE
C           0 SI MAILLAGE NON STRUCTURE
C NBSOEF :  NOMBRE DE SOMMETS   STOCKES PAR EF (0 POUR COMPLETER)
C NBTGEF :  NOMBRE DE TANGENTES STOCKES PAR EF (0 POUR COMPLETER)
C           ( TRIANGLE STRUCTURE SANS TG => NBSOEL=3 NBSOEF=4 NBTGEF=0
C             TRIANGLE STRUCTURE AVEC TG => NBSOEL=3 NBSOEF=4 NBTGEF=8 ... )
C
C LDAPEF :  NOMBRE DE MOTS DE DECALAGE POUR POINTER SUR LDEFAP(1)
C           ADRESSE MCN = MNNSEF + LDAPEF
C LDNGEF :  NOMBRE DE MOTS DE DECALAGE POUR POINTER SUR NGEFTG(1)
C           ADRESSE MCN = MNNSEF + LDNGEF
C LDTGEF :  NOMBRE DE MOTS DE DECALAGE POUR POINTER SUR NUTGEF(1,1)
C           ADRESSE MCN = MNNSEF + LDTGEF
C
C NBEFOB :  NOMBRE D'ELEMENTS FINIS DU MAILLAGE
C NX,NY,NZ: LE NOMBRE D'ARETES DANS LES DIRECTION X Y Z
C           CF LE TMS ~/td/d/a___nsef
C IERR   :  =0 SI PAS DE PROBLEME
C           >0 SI UNE ERREUR EST RENCONTREE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS       MARS 1989
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS       MAI  1996
C2345X7..............................................................012
      IMPLICIT INTEGER  (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/a___nsef.inc"
      INTEGER           LTNSEF(0:*)
C
C     INITIALISATIONS
      IERR   = 0
      NBEFOB = 0
      NX     = 0
      NY     = 0
      NZ     = 0
C
C     LE TYPE DU MAILLAGE
      NUTYMA = LTNSEF( WUTYMA )
C
C     LE NOMBRE D'EF DE CE MAILLAGE
      NBEFOB = LTNSEF( WBEFOB )
C
C     LE NOMBRE DE NUMEROS DE SOMMETS STOCKES POUR LES EF DE CE MAILLAGE
C     DIFFERENT DE NBSOEL = NOMBRE DE SOMMETS REELS DE CET EF
C     (TRIANGLE NBSOEF=4 ET NBSOEL=3 ...)
      NBSOEF = LTNSEF( WBSOEF )
C
C     LE NOMBRE DE TANGENTES PAR EF
      NBTGEF = LTNSEF( WBTGEF )
C
C     LE NOMBRE D'EF AVEC TANGENTES
      NBEFTG = LTNSEF( WBEFTG )
C
C     LE NOMBRE DE MOTS DU TABLEAU POINTEUR SUR LES EF A TG
      NBEFAP = LTNSEF( WBEFAP )
C
      IF( NUTYMA .EQ. 0 ) THEN
C
C        MAILLAGE NON STRUCTURE AVEC OU SANS TANGENTES
C        =============================================
         NBSOEL = 0
         LDAPEF = WUSOEF + NBSOEF * NBEFOB
C
      ELSE IF ( NUTYMA .GT. 0 ) THEN
C
C        MAILLAGE STRUCTURE
         GOTO ( 10, 20, 30, 40, 50, 60, 70, 80 ) , NUTYMA
C
C        NOEUDSOMMET STRUCTURE
C        =====================
  10     NX     = LTNSEF( WBARNS )
         NBSOEL = 1
         LDAPEF = WBARNS + 1
         GOTO 100
C
C        SEGMENT STRUCTURE
C        =================
  20     NX     = LTNSEF( WBARSE )
         NBSOEL = 2
         LDAPEF = WBARSE + 1
         GOTO 100
C
C        TRIANGLE STRUCTURE
C        ==================
  30     NX     = LTNSEF( WBARTR )
         NBSOEL = 3
         LDAPEF = WBARTR + 1
         GOTO 100
C
C        QUADRANGLE STRUCTURE
C        ====================
  40     NX     = LTNSEF( WBARXQ )
         NY     = LTNSEF( WBARYQ )
         NBSOEL = 4
         LDAPEF = WBARYQ + 1
         GOTO 100
C
C        TETRAEDRE STRUCTURE
C        ===================
  50     NX     = LTNSEF( WBARTE )
         NBSOEL = 4
         LDAPEF = WBARTE + 1
         GOTO 100
C
C        PENTAEDRE STRUCTURE
C        ===================
  60     NX     = LTNSEF( WBARTP )
         NY     = LTNSEF( WBARZP )
         NBSOEL = 6
         LDAPEF = WBARZP + 1
         GOTO 100
C
C        HEXAEDRE STRUCTURE
C        ==================
  70     NX     = LTNSEF( WBARXH )
         NY     = LTNSEF( WBARYH )
         NZ     = LTNSEF( WBARZH )
         NBSOEL = 8
         LDAPEF = WBARZH + 1
         GOTO 100
C
C        6-CUBE STRUCTURE
C        ================
C        NOMBRE D'ARETES DANS UNE DIRECTION
  80     NX     = LTNSEF( WANCUB )
C        NY     = DIMENSION n du 6-CUBE
         NY     = LTNSEF( WDNCUB )
         NBSOEL = 2 ** NY
         LDAPEF = WDNCUB + 1
         GOTO 100
C
      ELSE
C
C        ERREUR DE TYPE
         NBLGRC(NRERR) = 1
         KERR(1) = 'NSEFPA:NUTYMA TYPE INCORRECT'
     %          // KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
      ENDIF
C
C     LES DECALAGES POUR ARRIVER SUR LPEFAP ET NGEFTG
C     -----------------------------------------------
 100  LDNGEF = LDAPEF + NBEFAP
      LDTGEF = LDNGEF + NBEFTG
C
      IF( NBEFOB .LE. 0 ) THEN
         IERR = 1
      ENDIF
      END
