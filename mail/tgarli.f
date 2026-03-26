      SUBROUTINE TGARLI( MNARLI, MNSOLI,
     %                   NBARLI, MNXYTG, MNNTGL, MNCGEF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : STOCKAGE DES 3 COMPOSANTES DES 2 TANGENTES DES ARETES D'UNE LIGNE
C -----
C
C ENTREES:
C --------
C MNARLI : ADRESSE MCN DU TABLEAU 'nsef'      DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TABLEAU 'xyzsommet' DE LA LIGNE
C
C SORTIES:
C --------
C NBARLI : NOMBRE D'ARETES DE LA LIGNE  (0 SI PROBLEME RENCONTRE)
C MNXYTG : ADRESSE MCN DU TABLEAU DES 3 COMPOSANTES DES TANGENTES
C MNNTGL : ADRESSE MCN DU TABLEAU DES 2 NUMEROS DES TANGENTES DES
C          NBARLI ARETES DE LA LIGNE    (0 SI PROBLEME RENCONTRE)
C          TABLEAU A DETRUIRE EN FIN D'UTILISATION
C MNCGEF : ADRESSE MCN DU TABLEAU DU CODE GEOMETRIQUE DES NBARLI ARETES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           NOSOEL(4)
C
      MNXYTG = 0
      MNNTGL = 0
      MNCGEF = 0
C
C     LE NOMBRE DE TANGENTES STOCKEES POUR CETTE LIGNE
      IF( MCN( MNSOLI + WNBTGS ) .LE. 0 ) GOTO 9999
C
C     LE NOMBRE D'ARETES DE LA LIGNE
      NBARLI = MCN( MNARLI + WBEFOB )
      IF( NBARLI .LE. 0 ) GOTO 9999
C
C     RESERVATION DU TABLEAU NTGL
      CALL TNMCDC( 'ENTIER', 3*NBARLI, MNNTGL )
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNARLI) ,
     %             NUTYMA , NBSOEL , NBSOEF , NBTGEF,
     %             LDAPEF , LDNGEF , LDTGEF , NBARLI ,
     %             NX     , NY     , NZ     ,
     %             IERR   )
      IF( NBTGEF .LE. 0 .OR. IERR .NE. 0 ) GOTO 9990
C
C     PARCOURS DES ARETES DE LA LIGNE
      MNXYTG = MNSOLI + WYZSOM + 3 * MCN(MNSOLI+WNBSOM)
      MNNTG  = MNNTGL - 1
      MNCGEF = MNNTGL + 2 * NBARLI
      MNNCG  = MNCGEF - 1
C
      DO 50 NA=1,NBARLI
C        LE NUMERO DES NBSOEF SOMMETS DE L'EF NA
         CALL NSEFNS( NA     , NUTYMA , NBSOEF , NBTGEF,
     %                LDAPEF , LDNGEF , LDTGEF,
     %                MNARLI , NX , NY , NZ ,
     %                NCOGEL , NUGEEF , NUEFTG, NOSOEL , IERR )
C
C        LES NUMEROS DES 2 TANGENTES OU 0
         IF( NUEFTG .EQ. 0 ) THEN
C           LES 2 TANGENTES SONT NULLES
            DO 10 I=1,2
               MCN( MNNTG + I ) = 0
 10         CONTINUE
C           CODE GEOMETRIQUE
            MCN( MNNCG + 1 ) = 0
         ELSE
            DO 20 I=1,2
               MCN( MNNTG + I ) = NOSOEL( 2 + I )
 20         CONTINUE
C           CODE GEOMETRIQUE
            MCN( MNNCG + 1 ) = NUGEEF
         ENDIF
         MNNTG = MNNTG + 2
         MNNCG = MNNCG + 1
 50   CONTINUE
      RETURN
C
C     DESTRUCTION DU TABLEAU NTGL EN CAS D'ERREUR
 9990 IF( MNNTGL .GT. 0 ) CALL TNMCDS( 'ENTIER', 3*NBARLI, MNNTGL )
 9999 NBARLI = 0
      RETURN
      END
