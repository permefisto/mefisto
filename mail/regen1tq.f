       SUBROUTINE REGEN1TQ( MNNSEF, NOSOTQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    REGENERER le DERNIER QUADRANGLE ou TRIANGLE DETRUIT
C -----    DE NUMERO DE SOMMETS NOSOTQ
C -----
C ENTREES:
C --------
C MNNSEF : ADRESSE DU TABLEAU NSEF DE LA SURFACE
C NOSOTQ : NUMERO DES 4 SOMMETS DU TQ DETRUIT A REGENERER

C SORTIES:
C --------
C NOSOTQ : NUMERO DES 4 SOMMETS DU TQ DETRUIT A REGENERER
C          LA REGENERATION EST EFFECTUEE SEULEMENT SI NOSOTQ(1)>0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NOSOTQ(4)

      IERR = 0

      IF( NOSOTQ(1) .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
           PRINT*,'regen1tq: le TQ NOSOTQ',NOSOTQ,' N''EST PAS REGENERE'
         ELSE
           PRINT*,'regen1tq: the TQ NOSOTQ',NOSOTQ,' is NOT REGENERATED'
         ENDIF
         GOTO 9999
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
           PRINT*,'regen1tq: le TQ NOSOTQ St:',NOSOTQ,' est REGENERE'
         ELSE
           PRINT*,'regen1tq: the TQ NOSOTQ St',NOSOTQ,' is REGENERATED'
         ENDIF

      ENDIF

C     LE TRIANGLE ou QUADRANGLE EST AJOUTE EN DERNIERE POSITION
C     ---------------------------------------------------------
      NBTQ = MCN( MNNSEF+WBEFOB )
      NBTQ = NBTQ + 1
      MN   = MNNSEF + WUSOEF + 4 * NBTQ  - 5

      DO K = 1, 4
         MCN( MN + K ) = NOSOTQ( K )
      ENDDO

C     UN EF TQ DE PLUS
      MCN( MNNSEF+WBEFOB ) = NBTQ

      CALL ECDATE( MCN(MNNSEF) )
      MCN( MNNSEF+MOTVAR(6) )=NONMTD( '~>>>NSEF' )

 9999 RETURN
      END
