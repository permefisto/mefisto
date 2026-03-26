      SUBROUTINE FAEFTG( NTNSEF, MNNSEF, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER DU TMS NSEF TOUS LES FAUX EF A TANGENTES
C -----
C
C ENTREES:
C --------
C NTNSEF : NO TMS 'NSEF' DU MAILLAGE DU PLSV
C MNNSEF : ADRESSE MCN DU TMS 'NSEF'
C
C SORTIES:
C --------
C IERR   : =0 SI PAS D'ERREUR RENCONTREE
C          >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
C
C     LE TABLEAU 'NSEF' DU PLSV
      IF( NTNSEF .LE. 0 .OR. MNNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FAEFTG: PLSV SANS NSEF'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DU PLSV
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX,     NY,     NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE D'EF A TG DU PLSV
      NBEFTG = MCN( MNNSEF + WBEFTG )
      IF( NBEFTG .LE. 0 ) RETURN
C
C     MISE A JOUR DU NUMERO DES TANGENTES DES EF A TG
      NBFAEF = 0
      MNTG0  = MNNSEF + LDTGEF - 1
      DO 20 N = 1, NBEFTG
C        +- LE NUMERO DES NBTGEF TG
         K = 0
         DO 10 NT=1,NBTGEF
C           +- NUMERO DE LA TG
C           DETECTION DES FAUX EF A TG
            K = K + ABS( MCN(MNTG0+NT) )
 10      CONTINUE
         IF( K .EQ. 0 ) NBFAEF = NBFAEF + 1
         MNTG0 = MNTG0 + NBTGEF
 20   CONTINUE
C
      IF( NBFAEF .EQ. 0 ) RETURN
C
C     SUPPRESSION DES FAUX EF A TG (TOUS LES NO DE TANGENTES SONT NULS)
      MNCG1   = 0
      NBEFTG1 = NBEFTG - NBFAEF
      IF( NBEFTG1 .EQ. 0 ) THEN
C
C        IL NE RESTE PLUS D'EF A TG
         NBEFAP = 0
         GOTO 9000
      ENDIF
C
C     IL RESTE AU MOINS UN EF A TG
      NBEFAP = NBEFOB
      CALL TNMCDC( 'ENTIER', NBEFTG1*(1+NBTGEF), MNCG1 )
      MNCG1 = MNCG1 - 1
      MNTG1 = MNCG1  + NBEFTG1
C
      MNAP0 = MNNSEF + LDAPEF - 1
      MNCG0 = MNAP0  + NBEFOB
      MNTG0 = MNCG0  + NBEFTG
C
      NBEFTG = 0
      DO 50 N = 1, NBEFOB
C
C        NUMERO D'EF A TG
         NUEFTG = MCN( MNAP0 + N )
         IF( NUEFTG .EQ. 0 ) GOTO 50
C
C        DETECTION DES FAUX EF A TG
         MNTG = MNTG0 + NBTGEF * NUEFTG - NBTGEF
         K = 0
         DO 30 NT=1,NBTGEF
C           +- LE NUMERO DES NBTGEF TG
            K = K + ABS( MCN(MNTG+NT) )
 30      CONTINUE
C
         IF( K .EQ. 0 ) THEN
C
C           FAUX EF A TG
            MCN( MNAP0 + N ) = 0
C
         ELSE
C
C           VRAI EF A TG
            NBEFTG = NBEFTG + 1
C           LE POINTEUR SUR L'EF A TG
            MCN( MNAP0 + N ) = NBEFTG
C           LE CODE GEOMETRIQUE
            MCN( MNCG1 + NBEFTG ) = MCN( MNCG0 + N )
C           LES NBTGEF NUMEROS DE TG
            DO 40 NT=1,NBTGEF
               MCN( MNTG1 + NT ) = MCN( MNTG + NT )
 40         CONTINUE
            MNTG1 = MNTG1 + NBTGEF
C
         ENDIF
 50   CONTINUE
C
C     RECOPIE DANS LE TMS 'NSEF'
C     LES CODES GEOMETRIQUES
      MNCG0 = MNCG0 + 1
      MNCG1 = MNCG1 + 1
      CALL TRTATA( MCN(MNCG1), MCN(MNCG0), NBEFTG1 )
C     LES NUMEROS DES TANGENTES
      MNTG0 = MNCG0 + NBEFTG1
      MNTG1 = MNCG1 + NBEFTG1
      CALL TRTATA( MCN(MNTG1), MCN(MNTG0), NBEFTG1*NBTGEF )
C
C     MISE A JOUR
 9000 MCN( MNNSEF + WBEFTG ) = NBEFTG1
      MCN( MNNSEF + WBEFAP ) = NBEFAP
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNNSEF) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE
      IF( MNCG1 .GT. 0 )CALL TNMCDS('ENTIER', NBEFTG1*(1+NBTGEF), MNCG1)
C
C     LE TMS EST RACCOURCI
      CALL TAMSRA( NTNSEF, LDAPEF + NBEFAP + NBEFTG1 * (1+NBTGEF) )
      END
