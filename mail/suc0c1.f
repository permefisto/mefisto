      SUBROUTINE SUC0C1( NMSURF, NTLXSF,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    POUR LE MAILLAGE D'UNE SURFACE, TRANSFORMER LES ARETES C0-C1
C -----    EN ARETES C1 ET MODIFIER EN CONSEQUENCE LES EF SURFACIQUES
C          DU MAILLAGE
C
C ENTREES:
C --------
C NMSURF : NOM DE LA SURFACE
C NTLXSF : NO TMS DU LEXIQUE DE LA SURFACE DE MAILLAGE A TRAITER
C
C SORTIES:
C --------
C NTNSEF : NO TMS 'NSEF'      DU MAILLAGE DE LA SURFACE
C MNNSEF : ADRESSE MCN DU TMS 'NSEF'
C NTXYZS : NO TMS 'XYZSOMMET' DU MAILLAGE DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET'
C IERR   : =0 SI PAS D'ERREUR RENCONTREE
C          >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMSURF
C
      INTEGER           NOSOEL(12)
C
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'NSEF', NTNSEF, MNNSEF )
      IF( NTNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUC0C1: SURFACE SANS NSEF'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DE LA SURFACE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE D'EF A TG DE LA SURFACE
      NBEFTG = MCN( MNNSEF + WBEFTG )
      IF( NBEFTG .LE. 0 ) RETURN
C
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUC0C1: SURFACE SANS XYZSOMMET'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LE NOMBRE DE TANGENTES DU MAILLAGE DE LA SURFACE
      NBTGSF = MCN( MNXYZS + WNBTGS )
      IF( NBTGSF .LE. 0 ) RETURN
C     LE NOMBRE DE SOMMETS DU MAILLAGE DE LA SURFACE
      NBSOSF = MCN( MNXYZS + WNBSOM )
C
C     CONSTRUCTION DU TABLEAU DES ARETES DE LA SURFACE
      CALL HAARSU( MCN(MNXYZS+WYZSOM), NMSURF, 4, NTARSU, MNARSU,
     %             IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE D'ENTIERS PAR ARETE
      MOARET = MCN( MNARSU + WOARET )
C     LA MAJORATION DU NOMBRE D'ARETES
      MXARET = MCN( MNARSU + WXARET )
C     LE NOMBRE D'ARETES A TG
      NBARTG = MCN( MNARSU + WBARTG )
      IF( NBARTG .LE. 0 ) RETURN
C
C     ADRESSE MCN DU TABLEAU LARETE
      MNARET = MNARSU + W1LGFR
      IF( MCN(MNARSU+WUMXLF) .GT. 0 ) THEN
         MNARET = MNARET + MCN(MNARSU+WUMXLF) - MCN(MNARSU+WUMILF) + 1
      ENDIF
C
C     ADRESSE MCN DU TABLEAU NUTGAR DU TMS A___ARETE
      MNTGAR = MNARET + MOARET*MXARET - 2
C
C     RESERVATION D'UN TABLEAU DES NOUVEAUX EF A TG
      MNNETG = 0
      CALL TNMCDC( 'ENTIER', NBARTG*NBTGEF, MNNETG )
      NBNETG = 0
C
C     PARCOURS DU TABLEAU ARETES
C     ==========================
      MNA = MNARET - MOARET - 1
      DO 100 NUAR = 1,MXARET
C
         MNA = MNA + MOARET
C        L'EVENTUEL NUMERO D'ARETE A TG
         NUARTG = MCN(MNA+7)
         IF( NUARTG .GT. 0 ) THEN
C
C           ARETE A TG : NUMERO DES 2 TANGENTES DE L'ARETE
            MN   = MNTGAR + 2 * NUARTG
            NTA1 = MCN( MN   )
            NTA2 = MCN( MN+1 )
            IF( NTA1 .EQ. 0 .AND. NTA2 .EQ. 0 ) GOTO 100
C
C           IL EXISTE AU MOINS UNE TANGENTE A CETTE ARETE
            IF( MCN(MNA+5) .EQ. 0 ) GOTO 100
C
C           L'ARETE EST ADJACENTE A 2 EF
C           LE NUMERO DES 2 SOMMETS DE L'ARETE NUAR
            NS1 = MCN(MNA+1)
            NS2 = MCN(MNA+2)
C
            DO 80 K=1,2
C
C              NUMERO DE L'EF CONTENANT CETTE ARETE
               NEF = MCN(MNA+3+K)
C
C              LE NUMERO DES SOMMETS DE L'EF NEF DE LA SURFACE
               CALL NSEFNS( ABS(NEF), NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNNSEF, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C              PAS DE TEST SUR IERR CAR DEJA FAIT DANS CALL HAARSU
C
C              RECHERCHE DU NO DE L'ARETE NS1 NS2 DANS NEF
               IF( NEF .GT. 0 ) THEN
C                 ARETE NS1 NS2 VUE DANS LE SENS DIRECT
                  DO 10 I1=1,NCOGEL
                     IF( NOSOEL(I1) .EQ. NS1 ) THEN
                        IF( I1 .NE. NCOGEL ) THEN
                           I2 = I1 + 1
                        ELSE
                           I2 = 1
                        ENDIF
                        IF( NOSOEL(I2) .EQ. NS2 ) GOTO 30
                     ENDIF
 10               CONTINUE
               ELSE
C                 ARETE NS1 NS2 VUE DANS LE SENS INDIRECT
                  DO 20 I1=1,NCOGEL
                     IF( NOSOEL(I1) .EQ. NS2 ) THEN
                        IF( I1 .NE. NCOGEL ) THEN
                           I2 = I1 + 1
                        ELSE
                           I2 = 1
                        ENDIF
                        IF( NOSOEL(I2) .EQ. NS1 ) GOTO 30
                     ENDIF
 20               CONTINUE
               ENDIF
C
C              I1 EST LE NUMERO DE L'ARETE NS1 NS2 DANS NEF
C              POSITION DES 2 TANGENTES DE CETTE ARETE DANS NUEFTG
 30            NTG1 = 2 * I1 - 1
               IF( I1 .NE. NCOGEL ) THEN
                  NTG2 = NTG1 + 3
               ELSE
                  NTG2 = 2
               ENDIF
C              ADRESSE DANS NUTGEF
               IF( NUEFTG .GT. 0 ) THEN
C                 EF A TG
                  MN = MNNSEF + LDTGEF + NBTGEF * (NUEFTG-1)
C                 LES 2 TANGENTES DE L'ARETE SONT IMPOSEES A L'ARETE DE L'EF
C                 PUISQUE LES 2 TG DE L'ARETE SONT LA SYNTHESE DE CELLES
C                 DES 2 EF ADJACENTS PAR CETTE ARETE
               ELSE
C                 EF SANS TG AU DEPART
                  MNNS = MNNSEF + LDAPEF + ABS(NEF)-1
                  IF( MCN(MNNS) .EQ. 0 ) THEN
C                    EF ENCORE SANS TG => IL DEVIENT UN EF A TG
                     NBNETG = NBNETG + 1
C                    LE SIGNE - POUR INDIQUER UN NOUVEL EF A TG
                     MCN( MNNS ) = -(NBEFTG+NBNETG)
C                    TOUS LES NUMEROS DE TANGENTES SONT MIS A ZERO
                     MN = MNNETG + NBTGEF * (NBNETG-1)
                     CALL AZEROI( NBTGEF, MCN(MN) )
                  ELSE
C                    EF CREE AVEC DEJA AU MOINS UNE ARETE A TG
                     MN = MNNETG + NBTGEF * (-MCN(MNNS)-NBEFTG-1)
                  ENDIF
               ENDIF
C              LES 2 TANGENTES DE L'ARETE SONT AJOUTEES A PARTIR DE MN
               IF( NEF .GT. 0 ) THEN
C                 PROTECTION DE NTA1 NTA2 POUR L'EF K SUIVANT
                  NT1 = NTA1
                  NT2 = NTA2
                ELSE
                  NT1 = NTA2
                  NT2 = NTA1
               ENDIF
               MN = MN - 1
               MCN( MN+NTG1 ) = NT1
               MCN( MN+NTG2 ) = NT2
 80         CONTINUE
         ENDIF
 100  CONTINUE
C
C     ALLONGEMENT DU TMS 'NSEF' ET MISE A JOUR
C     ========================================
      IF( NBNETG .GT. 0 ) THEN
         NBEFAP = MCN(MNNSEF+WBEFAP)
         L      = LDAPEF + NBEFAP + (NBEFTG+NBNETG) * (1+NBTGEF)
         CALL TAMSAU( NTNSEF, L )
         CALL TAMSOU( NTNSEF, MNNSEF )
C
C        LE POINTEUR SUR LES EF A TG REDEVIENT POSITIF
         MN = MNNSEF + LDAPEF
         DO 190 I = MN, MN+NBEFAP-1
            MCN(I) = ABS( MCN(I) )
 190     CONTINUE
C
C        DECALAGE DE NBNETG DU TABLEAU DES NUMEROS DE TG DES EF A TG
         MN = MN + NBEFAP + NBEFTG
         DO 200 I = MN+NBTGEF*NBEFTG-1, MN, -1
            MCN(I+NBNETG) = MCN(I)
 200     CONTINUE
C
C        MISE A ZERO DU CODE GEOMETRIQUE (ANCIEN EF C0)
         DO 210 I = MN, MN+NBNETG-1
            MCN(I) = 0
 210     CONTINUE
C
C        COPIE DES TG DES NOUVEAUX EF A TG A LA SUITE DES ANCIENS
         MN = MN + NBNETG + NBEFTG * NBTGEF
         CALL TRTATA( MCN(MNNETG), MCN(MN), NBNETG*NBTGEF )
C
C        LE NOMBRE FINAL D'EF A TG
         MCN( MNNSEF + WBEFTG ) = NBEFTG + NBNETG
      ENDIF
C
C     DESTRUCTION DU TABLEAU NETG
      CALL TNMCDS( 'ENTIER', NBARTG*NBTGEF, MNNETG )
      END
