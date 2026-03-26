      SUBROUTINE STNOST( NTLXOB,
     %                   NTNSE1, MNNSE1, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DESTRUCTURATION D'UN PLSV STRUCTURE EN NON STRUCTURE
C -----
C
C ENTREES:
C --------
C NTLXOB : NUMERO DU TABLEAU TMS DU LEXIQUE DU PLSV
C
C SORTIES:
C --------
C NTNSE1 : NUMERO      DU TMS 'NSEF' DU PLSV NON STRUCTURE
C MNNSE1 : ADRESSE MCN DU TMS 'NSEF' DU PLSV NON STRUCTURE
C          CF ~td/d/a___nsef
C IERR   : 0 SI PAS D'ERREUR RENCONTREE
C          1 SI       ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/a___nsef.inc"
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
C
      IERR   = 0
      NTNSE1 = 0
      MNNSE1 = 0
C
C     LE PLSV A DESTRUCTURER
C     ======================
C     LE TABLEAU 'NSEF' DE CE PLSV
      CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )
      IF( NTNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PLSV SANS LE MAILLAGE NSEF'
         ELSE
            KERR(1) = 'PLSV WITHOUT the MESH NSEF'
         ENDIF
         CALL LEREUR
         IERR   = 1
         RETURN
      ENDIF
C
C     LE PLSV EST IL STRUCTURE ?
      NUTYMA = MCN( MNNSEF + WUTYMA )
      IF( NUTYMA .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'PLSV DEJA NON STRUCTURE'
         ELSE
            KERR(1) = 'PLSV ALREADY NOT STRUCTURED'
         ENDIF
         CALL LERESU
         NTNSE1 = NTNSEF
         MNNSE1 = MNNSEF
         RETURN
      ENDIF
C
C     LES CARACTERISTIQUES DU MAILLAGE STRUCTURE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF,
     %             NBEFOB, NX, NY, NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE DE MOTS DU TMS NSEF NON STRUCTURE
      NBTGEF = MCN( MNNSEF + WBTGEF )
      NBEFTG = MCN( MNNSEF + WBEFTG )
      NBEFAP = MCN( MNNSEF + WBEFAP )
      IF( NBEFTG .GT. 0 ) THEN
         MOTS2 = NBEFAP + NBEFTG * (1+NBTGEF)
      ELSE
C        POUR CONTENIR LE TABLEAU DES NO DES SOMMETS ET TGS DU DERNIER EF
         MOTS2 = 0
      ENDIF
C
C     64 MOTS DE PLUS POUR CONTENIR LE TABLEAU DES NO DES SOMMETS ET TGS DU DERN
      MOTS3 = 64
C
C     NOMBRE DE MOTS DE LA PREMIERE PARTIE DU TMS NSEF
      MOTS1  = WUSOEF + NBEFOB * NBSOEF
      MNNONS = 0
      CALL TNMCDC( 'MOTS', MOTS1+MOTS2+MOTS3, MNNONS )
C
C     COPIE DE LA PREMIERE PARTIE DU TMS NSEF
      CALL TRTATA( MCN(MNNSEF), MCN(MNNONS), WUSOEF )
C
C     LE TYPE DU MAILLAGE EST MAINTENANT NON STRUCTURE
      MCN( MNNONS + WUTYMA ) = 0
C
C     LE PARCOURS DES EF DU MAILLAGE STRUCTURE
      MN = MNNONS + WUSOEF
C
C     REMPLISSAGE DES NUMEROS DES NBSOEF SOMMETS DES EF
      DO 10 NUELEM = 1, NBEFOB
C        LE NUMERO DES SOMMETS DE L'ELEMENT NUELEM
         CALL NSEFNS( NUELEM, NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG,
     %                MCN(MN), IERR )
         IF( IERR .NE. 0 ) GOTO 9900
         MN = MN + NBSOEF
 10   CONTINUE
C
C     COPIE DE LA DERNIERE PARTIE DU TMS NSEF
      IF( MOTS2 .GT. 0 ) THEN
         IF( NUTYMA .EQ. 1 ) THEN
            MN = MNNSEF + WBARNS
         ELSE IF( NUTYMA .EQ. 2 ) THEN
            MN = MNNSEF + WBARSE
         ELSE IF( NUTYMA .EQ. 3 ) THEN
            MN = MNNSEF + WBARTR
         ELSE IF( NUTYMA .EQ. 4 ) THEN
            MN = MNNSEF + WBARYQ
         ELSE IF( NUTYMA .EQ. 5 ) THEN
            MN = MNNSEF + WBARTE
         ELSE IF( NUTYMA .EQ. 6 ) THEN
            MN = MNNSEF + WBARZP
         ELSE IF( NUTYMA .EQ. 7 ) THEN
            MN = MNNSEF + WBARZH
         ELSE IF( NUTYMA .EQ. 8 ) THEN
            MN = MNNSEF + WDNCUB
         ENDIF
         CALL TRTATA( MCN(MN+1),
     %                MCN(MNNONS+WUSOEF+NBEFOB*NBSOEF), MOTS2 )
      ENDIF
C
C     DESTRUCTION PUIS RECREATION DU PLSV NSEF NON STRUCTURE
C     ======================================================
      CALL LXTSDS( NTLXOB, 'NSEF' )
      CALL LXTNDC( NTLXOB, 'NSEF', 'ENTIER', MOTS1+MOTS2 )
      CALL LXTSOU( NTLXOB, 'NSEF',  NTNSE1 , MNNSE1 )
C
C     COPIE DU TMC DANS LE TMS
      CALL TRTATA( MCN(MNNONS), MCN(MNNSE1), MOTS1+MOTS2 )
C
C     DESTRUCTION DU TMC NSEF NON STRUCTURE
 9900 CALL TNMCDS( 'MOTS', MOTS1+MOTS2+MOTS3, MNNONS )
      RETURN
      END
