      SUBROUTINE CREETQ1FC( L1ARFA, MNARFA, NBAR1F, MNAR1F,
     %                      MNNSEF, MONSEF, NBTRAJ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREER TOUT Q-TRIANGLE avec 2 ARETES d'1 QT et un SOMMET COMMUN
C ------   C-A-D 2 ARETES CHACUNES APPARTENANT a 1 QT et CES 2 ARETES
C          ONT UN SOMMET COMMUN

C ENTREES:
C --------
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C MNARFA : ADRESSE DANS MCN DU TABLEAU NARFA DES QTANGLES DU MAILLAGE
C          EN SORTIE NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C                    NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C                    NARFA(3,I)= CHAINAGE HACHAGE SUR LE QTANGLE SUIVANT
C                    NARFA(4:L1ARFA,I)= NO NOSTQT DU QTANGLE
C                                       CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE L1ARFA-3 QTANGLES, 
C          LE DERNIER NUMERO DE QTANGLE EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES QTANGLES EST INCOMPLETE
C NBAR1F : NOMBRE D'ARETES APPARTENANT A UNE SEULE FACE
C MNAR1F : >0 ADRESSE MCN DU TMC NAR1F NO NARFA DES ARETES TRIPLES
C          =0 SI NBAR1F=0

C MODIFIES:
C ---------
C MNNSEF : ADRESSE DU TABLEAU NSEF DES FACES QT DE LA SURFACE
C          + le TMS NSEF de la SURFACE a ete AUGMENTE de NBTRAJ TRIANGLES
C MONSEF : NOMBRE DE MOTS DU TMS NSEF DE LA SURFACE

C SORTIE :
C --------
C NBTRAJ : NOMBRE DE TRIANGLES AJOUTES EN POSITION NBEFOB-NBTRAJ+1:NBEFOB
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C2345X7..............................................................012
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      INTEGER        NUMST(3)

      IF( NBAR1F .LE. 0 ) RETURN

C     RECHERCHE DES ARETES SIMPLES AR1F AYANT UN SOMMET COMMUN
C     --------------------------------------------------------
      MNAR1F1 = MNAR1F - 1
      NBTRAJ  = 0
      DO 20 N1 = 1, NBAR1F

C        LA N1-EME ARETE SIMPLE DU TABLEAU AR1F(APPARTENANT A UNE SEULE FACE)
C        NUMERO NARFA DE L'ARETE N1
         NAR1  = MCN( MNAR1F1 + N1 )
         IF( NAR1 .LE. 0 ) GOTO 20

C        ADRESSE MCN DE L'ARETE N1 DANS LE TABLEAU NARFA
         MNA1  = MNARFA - L1ARFA + L1ARFA * NAR1

C        LE NUMERO DES 2 SOMMETS DE L'ARETE SIMPLE N1
         NA1S1 = MCN( MNA1     )
         NA1S2 = MCN( MNA1 + 1 )

C        LE NUMERO DE FACE TQ DE L'ARETE SIMPLE N1
         NF1 = MCN( MNA1 + 4 )

C        RECHERCHE D'UNE AUTRE ARETE SIMPLE DE SOMMET NS1 ou NS2
         MNAR2 = MNAR1F - 1
         DO 10 N2 = 1, NBAR1F

C           LA N2-EME ARETE SIMPLE APPARTENANT A UNE SEULE FACE
            IF( N2 .EQ. N1 ) GOTO 10

C           NUMERO NARFA DE L'ARETE N2
            NAR2 = MCN( MNAR1F1 + N2 )
            IF( NAR2 .LE. 0 ) GOTO 10

C           ADRESSE MCN DE L'ARETE N2 DANS LE TABLEAU NARFA
            MNA2 = MNARFA - L1ARFA + L1ARFA * NAR2

C           LE NUMERO DES 2 SOMMETS DE L'ARETE SIMPLE N2
            NA2S1 = MCN( MNA2     )
            NA2S2 = MCN( MNA2 + 1 )

C           LE NUMERO DE FACE TQ DE L'ARETE SIMPLE N2
            NF2 = MCN( MNA2 + 4 )

            IF( NF1 .EQ. NF2 ) GOTO 10

            NBCPLE = 0
            IF( NA1S1 .EQ. NA2S1 .OR. NA1S2 .EQ. NA2S1 ) THEN
               NBCPLE = 1
               NUMST(3) = NA2S2
            ENDIF

            IF( NA1S1 .EQ. NA2S2 .OR. NA1S2 .EQ. NA2S2) THEN
               NBCPLE = 1
               NUMST(3) = NA2S1
            ENDIF

            IF( NBCPLE .EQ. 1 ) THEN

C              2 ARETES SIMPLES AVEC UN SOMMET COMMUN -> 1 TRIANGLE
C              AJOUT DU TRIANGLE NUMST(1:3) FORME DES 2 ARETES AU TMS NSEF
C              -----------------------------------------------------------
               NBTRAJ = NBTRAJ + 1
               NUMST(1) = NA1S1
               NUMST(2) = NA1S2

               NBEFOB = MCN(MNNSEF+WBEFOB)
               MONSEF = WUSOEF + 4 * NBEFOB
               IF( MONSEF .LT. L+4 ) THEN
                  CALL TNMCAU( 'ENTIER', L, L+64, L, MNNSEF )
                  MONSEF=L+64
               ENDIF
               MN = MNNSEF + WUSOEF + 4 * NBEFOB - 1
               DO K=1,3
                  MCN(MN+K) = NUMST(K)
               ENDDO
               MCN(MN+4) = 0

C              UN TRIANGLE DE PLUS DANS LE MAILLAGE NSEF
               NBEFOB = NBEFOB + 1
               MCN(MNNSEF+WBEFOB) = NBEFOB

               CALL ECDATE( MCN(MNNSEF) )
               MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )

C              LES 2 ARETES NAR1 et NAR2 SONT TRAITEES
               MCN( MNAR1F1 + N1 ) = -NAR1
               MCN( MNAR1F1 + N2 ) = -NAR2

               PRINT*,'creetq1fc: CREATION du TRIANGLE',NBEFOB,' de St:'
     %               ,NUMST

               GOTO 20

            ENDIF

 10      ENDDO

 20   ENDDO

      PRINT*,'creetq1fc: CREATION de',NBTRAJ,' TRIANGLES'

C     REMISE POSITIVE DES NUMEROS D'ARETES SIMPLES
      DO N1 = 1, NBAR1F
C        LA N1-EME ARETE SIMPLE APPARTENANT A UNE SEULE FACE
C        NUMERO NARFA DE L'ARETE N1
         MNA1 = MNAR1F1 + N1
         NAR1 = MCN( MNA1 )
         IF( NAR1 .LT. 0 ) MCN( MNA1 ) = -NAR1
      ENDDO

      RETURN
      END
