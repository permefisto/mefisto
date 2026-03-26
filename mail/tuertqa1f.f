      SUBROUTINE TUERTQA1F( MINAR1F, L1ARFA, L2ARFA, MNARFA,
     %                       NBAR1F, MNAR1F, MNNSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE LES FACES (TRIANGLES ou QUADRANGLES) AYANT
C ------   AU MOINS 2 ou 1 UNE ARETE SIMPLE (APPARENANT A UN SEUL TQ)

C ENTREES:
C --------
C MINAR1F: MINIMUM D'ARETES SIMPLES DANS 1 SEUL TQ POUR DETRUIRE SES TQ
C          =1 AU MOINS 1 ARETE  SIMPLE  POUR DETRUIRE SES TQ
C          =2 AU MOINS 2 ARETES SIMPLES POUR DETRUIRE SES TQ
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE MAXIMAL D'ARETES DECLARABLES DANS LE TABLEAU MCN(MNARFA)
C MNARFA : ADRESSE DANS MCN DU TABLEAU NARFA DES QTANGLES DU MAILLAGE
C          NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR LE QTANGLE SUIVANT
C          NARFA(4:L1ARFA,I)= NO NOSTQT DU QTANGLE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE L1ARFA-3 QTANGLES, 
C          LE DERNIER NUMERO DE QTANGLE EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES QTANGLES EST INCOMPLETE
C NBAR1F : NOMBRE D'ARETES APPARTENANT A UNE SEULE FACE
C MNAR1F : >0 ADRESSE MCN DU TMC NAR1F NO NARFA DES ARETES SIMPLES
C          =0 SI NBAR1F=0
C MNNSEF : ADRESSE DU TABLEAU NSEF DE LA SURFACE

C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          2 UNE ARETE NON RETROUVEE DANS LE TABLEAU NARFA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C2345X7..............................................................012
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NOSOAR(2)

      IERR = 0

      IF( NBAR1F .LE. 0 ) GOTO 9999

C     SUPPRESSION EVENTUELLE DES FACES DES NBAR1F ARETES SIMPLES
C     ----------------------------------------------------------
      NBFSUP   = 0
      MNNOSOEF = MNNSEF + WUSOEF
      MNAR1F1  = MNAR1F - 1
      DO 10 N = 1, NBAR1F

C        LA N-EME ARETE SIMPLE APPARTENANT A UNE SEULE FACE
         MNAR1 = MNAR1F1 + N
         NAR1F = MCN( MNAR1 )
         IF( NAR1F .LE. 0 ) GOTO 10

C        ADRESSE MCN DE L'ARETE DANS LE TABLEAU NARFA
         MNAR = MNARFA - L1ARFA + L1ARFA * NAR1F

C        NUMERO DE LA FACE TQ APPARTENANT A UN SEUL TQ
         N1TQ = ABS( MCN( MNAR+3 ) )

C        SI MINAR1F=1 ALORS
C           TOUTE FACE TQ AYANT AU MOINS UNE ARETE SIMPLE EST DETRUITE
C        SINON SI MINAR1F=2 ALORS
C           SEULEMENT SI N1TQ A AU MOINS 2 ARETES SIMPLES EST DETRUITE

         IF( MINAR1F .EQ. 2 ) THEN

C           SI N1TQ A AU MOINS 2 ARETES SIMPLES IL SERA DETRUIT
C           ---------------------------------------------------
            NB = 0
            DO 5 NA = 1, 3

C              LE NO DES 2 SOMMETS DU N1TQ
               NS1 = MCN( MNNOSOEF + 4 * N1TQ -5 + NA )
               IF( NA .EQ. 3 ) THEN
                  NA1 = 1
               ELSE
                  NA1 = NA+1
               ENDIF
               NS2 =  MCN( MNNOSOEF + 4 * N1TQ -5 + NA1 )

C              QUEL EST SON NO D'ARETE?
C              RECHERCHE DU QTANGLE AJACENT PAR L'ARETE NS1 NS2
               IF( NS1 .LT. NS2 ) THEN
                  NOSOAR(1) = NS1
                  NOSOAR(2) = NS2
               ELSE
                  NOSOAR(1) = NS2
                  NOSOAR(2) = NS1
               ENDIF

C              NOARET LE NUMERO DE L'ARETE NS1 NS2 DANS LARETE
               CALL HACHAR(2, NOSOAR,L1ARFA,L2ARFA,MCN(MNARFA),3,NOARET)

               IF( NOARET .LE. 0 ) THEN
C                  NOARET =0 SI LE TABLEAU NOSOAR N'A PAS ETE RETROUVE
C                         >0 SI LE TABLEAU NOSOAR   A     ETE RETROUVE
                  IERR = 2
                  GOTO 9999
               ENDIF

C              L'ARETE NS1 NS2 EST ELLE UNE ARETE SIMPLE
C              C-A-D APPARTIENT ELLE AU TABLEAU AR1F?
               DO K = 1, NBAR1F
                  IF( NOARET .EQ. ABS(MCN(MNAR1F1+K)) ) THEN
C                    NOARET EST UNE ARETE SIMPLE
                     NB = NB + 1
                     GOTO 5
                  ENDIF
               ENDDO
C              NON: NOARET N'APPARTIENT PAS AU TABLEAU AR1F

 5          ENDDO

            IF( NB .LE. 1 ) GOTO 10
C           N1TQ A AU PLUS 1 ARETE SIMPLE -> PAS DE DESTRUCTION

         ENDIF

C        DESTRUCTION DE CETTE FACE N1TQ
         NBFSUP = NBFSUP + 1
         MN1  = MNNOSOEF + 4 * N1TQ - 5
         DO L=1,4
            MCN( MN1+L ) = 0
         ENDDO

C        NUMERO DE L'ARETE SIMPLE N1 EST RENDU NEGATIF
         MCN( MNAR1 ) = -NAR1F

 10   ENDDO


C     COMPRESSION PAR SUPPRESSION DES TQ VIDES
C     ----------------------------------------
      MNTQ0 = MNNOSOEF -1
      MNTQ1 = MNNOSOEF -1
      NBTQ1 = 0
      DO N = 1, MCN(MNNSEF+WBEFOB)

         IF( MCN( MNTQ0+1 ) .GT. 0 ) THEN
            NBTQ1 = NBTQ1 + 1
            DO K=1,4
               MCN( MNTQ1 + K ) = MCN( MNTQ0 + K )
            ENDDO
            MNTQ1 = MNTQ1 + 4
         ENDIF

         MNTQ0 = MNTQ0 + 4

      ENDDO

C     MISE A JOUR DU NOMBRE DE TQ
      MCN(MNNSEF+WBEFOB) = NBTQ1

C     MISE A JOUR DU TMS NSEF
      CALL ECDATE( MCN(MNNSEF) )
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )

 9999 PRINT*,'tuertqa1f:',NBFSUP,
     %       ' FACES TQ AYANT UNE ARETE SIMPLE ont ete SUPPRIMEES'
      RETURN
      END
