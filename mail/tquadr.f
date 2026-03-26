      SUBROUTINE TQUADR( NTLXSU, NBS1, NBS2, NBS4, NBTGS,
     %                   NTFASU, MNFASU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE TABLEAU NSEF D'UNE QUADRANGULATION-TRIANGULATION
C -----    D'UN QUADRANGLE COURBE EN PRESENCE OU NON DE TANGENTES
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TAMS LEXIQUE DU QUADRANGLE COURBE
C NBS1   : NOMBRE DE SOMMETS DU COTE 1 ET 3 DU QUADRANGLE
C NBS2   : NOMBRE DE SOMMETS DU COTE 2 DU QUADRANGLE   NBS2 > NBS4
C NBS4   : NOMBRE DE SOMMETS DU COTE 4 DU QUADRANGLE
C NBTGS  : >0 LES EF ONT DES TANGENTES, 0 SINON
C
C SORTIES:
C --------
C NBTGS  : NOMBRE DE TANGENTES STOCKEES POUR LE MAILLAGE DU QUADRANGLE COURBE
C NTFASU : NUMERO DU TAMS DU TABLEAU NSEF
C MNFASU : ADRESSE MCN DU TABLEAU NSEF ( cf ./td/d/a___nsef )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a___nsef.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     QUELQUES INITIALISATIONS
      NBA1 = NBS1 - 1
      NBA4 = NBS4 - 1
      N24  = NBS2 - NBS4
      NDSQ = NBS1 - N24
C
C     LE NOMBRE DE QUADRANGLES
      NBQ  = (NDSQ-1) * NBA4
C
C     LE NOMBRE DE TRIANGLES
      I = 1
      DO 10 J=2,N24
        I = I + 2 * J - 1
 10   CONTINUE
      NBT = (NBS1-NDSQ)*2*NBA4 + I
C
C     LE NOMBRE D'EF DU MAILLAGE
      NBTQ = NBT + NBQ
C
C     LE NOMBRE DE TG PAR EF
C     LE NOMBRE DE TANGENTES STOCKEES POUR LE MAILLAGE DU QUADRANGLE COURBE
C     LE NOMBRE D'EF SUPPORTS DE TANGENTES
C     LE NOMBRE D'ENTIERS J POUR LES NUMEROS DES TANGENTES DES EF
      IF( NBTGS .GT. 0 ) THEN
         NBTGEF = 8
         NBTGS  = 6 * NBT + 8 * NBQ
         NBEFTG = NBTQ
         J      = NBTQ + NBTQ + 8 * NBTQ
      ELSE
         NBTGEF = 0
         NBTGS  = 0
         NBEFTG = 0
         J      = 0
      ENDIF
C
C     LE NOMBRE DE SOMMETS DE LA TRIANGULATION-QUADRANGULATION
      CALL LXTNDC( NTLXSU , 'NSEF' , 'ENTIER', WUSOEF+4*NBTQ+J )
      CALL LXTSOU( NTLXSU , 'NSEF' ,  NTFASU , MNFASU )
C
C     TYPE DE L'OBJET : SURFACE
      MCN ( MNFASU + WUTYOB ) = 3
C
C     SURFACE NON FERMEE
      MCN ( MNFASU + WUTFMA ) = 0
C
C     NUMERO DU TYPE DE MAILLAGE : NON STRUCTURE
      MCN ( MNFASU + WUTYMA ) = 0
C
C     NBSOEF  NOMBRE DE SOMMETS PAR NSEF
      MCN ( MNFASU + WBSOEF ) = 4
C
C     NBEFOB  NOMBRE D EF DE L'OBJET
      MCN ( MNFASU + WBEFOB ) = NBTQ
C
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C1
      MCN ( MNFASU + WBTGEF ) = NBTGEF
C
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNFASU + WBEFOB ) = NBTQ
C
C     LE NOMBRE D'EF AVEC TANGENTES DE LA SURFACE
      MCN ( MNFASU + WBEFTG ) = NBEFTG
C
C     LE NOMBRE D'EF AVEC POINTEUR SUR LES EF  A TG DE LA SURFACE
      MCN ( MNFASU + WBEFAP ) = NBEFTG
C
C     GENERATION DES QUADRANGLES ET TRIANGLES DES COUCHES REGULIERES
C     LE DEBUT DU TABLEAU NUSOEF ET DES POINTEURS SUR LES TANGENTES
      MN   = MNFASU + WUSOEF
C
      DO 40 J=1,NBA4
         DO 20 I=1,NDSQ-1
C           LES QUADRANGLES
            MCN( MN   ) = NUSOTQ( NBS1, NBS2, NBS4, I  , J   )
            MCN( MN+1 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J   )
            MCN( MN+2 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J+1 )
            MCN( MN+3 ) = NUSOTQ( NBS1, NBS2, NBS4, I  , J+1 )
            MN = MN + 4
 20      CONTINUE
         DO 30 I=NDSQ,NBA1
C           LES TRIANGLES
            MCN( MN   ) = NUSOTQ( NBS1, NBS2, NBS4, I  , J   )
            MCN( MN+1 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J+1 )
            MCN( MN+2 ) = NUSOTQ( NBS1, NBS2, NBS4, I  , J+1 )
            MCN( MN+3 ) = 0
            MN = MN + 4
            MCN( MN   ) = NUSOTQ( NBS1, NBS2, NBS4, I  , J   )
            MCN( MN+1 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J   )
            MCN( MN+2 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J+1 )
            MCN( MN+3 ) = 0
            MN = MN + 4
 30      CONTINUE
 40   CONTINUE
C
C     LES TRIANGLES DU COIN SUPERIEUR DROIT
      J = NBS4
      DO 90 I=NDSQ,NBA1
C        LE PREMIER TRIANGLE DE LA LIGNE
         MCN( MN   ) = NUSOTQ( NBS1, NBS2, NBS4, I  , J   )
         MCN( MN+1 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J   )
         MCN( MN+2 ) = NUSOTQ( NBS1, NBS2, NBS4, I+1, J+1 )
         MCN( MN+3 ) = 0
         MN = MN + 4
         DO 80 II=I+1,NBA1
C           LES TRIANGLES RESTANTS PAR COUPLES
            MCN( MN   ) = NUSOTQ( NBS1, NBS2, NBS4, II  , J   )
            MCN( MN+1 ) = NUSOTQ( NBS1, NBS2, NBS4, II+1, J+1 )
            MCN( MN+2 ) = NUSOTQ( NBS1, NBS2, NBS4, II  , J+1 )
            MCN( MN+3 ) = 0
            MN = MN + 4
            MCN( MN   ) = NUSOTQ( NBS1, NBS2, NBS4, II  , J   )
            MCN( MN+1 ) = NUSOTQ( NBS1, NBS2, NBS4, II+1, J   )
            MCN( MN+2 ) = NUSOTQ( NBS1, NBS2, NBS4, II+1, J+1 )
            MCN( MN+3 ) = 0
            MN = MN + 4
 80      CONTINUE
         J = J + 1
 90   CONTINUE
C
C     LE TABLEAU DES POINTEURS ET DES CODES GEOMETRIQUES
      IF( NBTGS .GT. 0 ) THEN
         MN = MNFASU + WUSOEF + 4 * NBTQ - 1
         DO 100 N=1,NBEFTG
            MN = MN + 1
            MCN(MN) = N
C           CODE GEOMETRIQUE 'INTERPOLATION TRANSFINIE'
            MCN(MN+NBEFTG) = 16
 100     CONTINUE
         MN = MN + NBEFTG
C
C        LES NUMEROS DES TANGENTES DES EF SUIVENT LA NUMEROTATION DES EF
         N = 0
         DO 140 J=1,NBA4
            DO 120 I=1,NDSQ-1
C              LE QUADRANGLE
               DO 115 L=1,8
                  N = N + 1
                  MCN(MN+L) = N
 115           CONTINUE
               MN = MN + 8
 120        CONTINUE
            DO 130 I=NDSQ,NBA1
C              LES 2 TRIANGLES
               DO 128 K=1,2
                  DO 125 L=1,6
                     N = N + 1
                     MCN(MN+L) = N
 125              CONTINUE
                  MCN(MN+7) = 0
                  MCN(MN+8) = 0
                  MN = MN + 8
 128           CONTINUE
 130        CONTINUE
 140     CONTINUE
C
         DO 190 I=NDSQ,NBA1
C           LE PREMIER TRIANGLE DE LA LIGNE
            DO 150 L=1,6
               N = N + 1
               MCN(MN+L) = N
 150        CONTINUE
            MCN(MN+7) = 0
            MCN(MN+8) = 0
            MN = MN + 8
            DO 180 II=I+1,NBA1
C              LES 2 TRIANGLES RESTANTS PAR COUPLES
               DO 170 K=1,2
                  DO 160 L=1,6
                     N = N + 1
                     MCN(MN+L) = N
 160              CONTINUE
                  MCN(MN+7) = 0
                  MCN(MN+8) = 0
                  MN = MN + 8
 170           CONTINUE
 180        CONTINUE
 190     CONTINUE
      ENDIF
C
C     LA DATE DE CREATION
      CALL ECDATE ( MCN( MNFASU ) )
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN ( MNFASU + MOTVAR(6) ) = NONMTD ( '~>>>NSEF' )
      END
