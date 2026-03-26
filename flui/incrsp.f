      SUBROUTINE INCRSP( NTDL,   NTDLFX, NDIM, NBNOVI, NDDLNO,
     %                   LPLIGN, LPCOLO, AG,
     %                   LPLIGC, LPCOLC, AGC,
     %                   IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DE LA MATRICE DE PRECONDITIONNEMENT SYMETRIQUE
C -----    AGC A PARTIR DE LA MATRICE MORSE AG
C          TOUT COEFFICIENT DE AG SUPPORT D'UN DL FIXE EST MIS A ZERO
C          DANS AGGC ET LE COEFFICIENT DIAGONAL DE AGC EST MIS A 1
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE AG (ET AGC)
C NTDLFX : NTDLFX(N) =0 SI LE DL N EST LIBRE, 1 SI LE DL N EST FIXE
C          NO TEMOIN D'UN DEGRE DE LIBERTE FIXE'
C LPLIGN : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AG
C LPCOLO : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AG
C AG     : MATRICE INITIALE NON FACTORISEE
C LPLIGC : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AGC
C LPCOLC : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AGC
C
C SORTIES:
C --------
C AGC    : MATRICE FACTORISEE INCOMPLETE <=> MATRICE DE PRECONDITIONNEMENT
C IERR   : =0 SI PAS D'ERREUR DETECTEE
C          =1 SI UNE COLONNE DE AG N'EST PAS RETROUVEE DANS AG
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  MARS    2009
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NTDLFX(NTDL), NDDLNO(0:NBNOVI)
      INTEGER           LPLIGN(NTDL+1), LPCOLO(*)
      INTEGER           LPLIGC(NTDL+1), LPCOLC(*)
      DOUBLE PRECISION  AG(*), AGC(*)
      DOUBLE PRECISION  S, PIVMAX, PIVMIN
C
C     INITIALISATION
C     --------------
      IERR = 0
C     RECHERCHE DU COEFF MAXIMAL EN VALEUR ABSOLUE SAUF SUR DL FIXES
C     A supprimer apres mise au point ...
      PIVMIN = 1D100
      PIVMAX =-PIVMIN
      DO 30 K=1,NTDL
         IF( NTDLFX(K) .EQ. 0 ) THEN
C           DL LIBRE, NON FIXE
            S = AG( LPLIGN(1+K) )
            IF( S .GT. PIVMAX ) PIVMAX = S
            IF( S .LT. PIVMIN ) PIVMIN = S
            IF( S .LT. 0D0 ) print *,'AG(',K,',',K,')=',S,' <0 !...'
         ENDIF
 30   CONTINUE
      print *,'PIVOT MINIMUM=',PIVMIN,'  PIVOT MAXIMUM=',PIVMAX
      PIVMAX = MAX( ABS(PIVMIN), ABS(PIVMAX) )
C
C     CONSTRUCTION DE AGC LA MATRICE DE PRECONDITIONNEMENT
C     ----------------------------------------------------
C     MISE A ZERO DE LA MATRICE DE PRECONDITIONNEMENT
      DO 40 K=1,LPLIGC(NTDL+1)
         AGC(K) = 0.D0
 40   CONTINUE
C
C     INITIALISATION DE LA MATRICE DE PRECONDITIONNEMENT
      ID0 = 0
      DO 100 N=1,NBNOVI
         ID1 = NDDLNO(N)
C
         IF( (ID1-ID0) .EQ. (NDIM+1) ) THEN
            IF( NTDLFX( ID1 ) .EQ. 0 ) THEN
C              SI LE DL EST LIBRE ET DE PRESSION IL EST MARQUE EN NEGATIF
               NTDLFX(ID1) = -ID1
               print *,'DL PRESSION NODL=',ID1
            ENDIF
         ENDIF
C
         DO 90 I=ID0+1, ID1
C
C        INDICE DU COEFFICIENT DIAGONAL I DE AG
         NDIAG = LPLIGN(I+1)
C
         IF( NTDLFX(I) .EQ. 0 ) THEN
C
C           LIGNE DU DEGRE DE LIBERTE I LIBRE
            IC1 = LPLIGC(I)
            IC2 = LPLIGC(I+1)
            IC0 = IC1
            DO 101 KT = LPLIGN(I)+1, NDIAG-1
               IC1  = IC0
               NCOL = LPCOLO(KT)
               DO 102 KC = IC1+1, IC2-1
                  IF( NCOL .EQ. LPCOLC(KC) ) THEN
C                    NUMERO DE COLONNE IDENTIQUE
                     IF( NTDLFX(NCOL) .EQ. 0 ) THEN
C                       DL LIBRE
                        AGC(KC) = AG(KT)
                     ELSE
C                       DL FIXE: COEF DE LA COLONNE ANNULE
                        AGC(KC) = 0D0
                     ENDIF
                     IC0 = KC
                     GOTO 101
                  ENDIF
 102           CONTINUE
               WRITE (IMPRIM,10102) I, NCOL
               IERR = 1
10102          FORMAT(' ERREUR INCRCO : LIGNE ',I12,
     %         ' INDICE DE COLONNE',I12,' DE AG NON RETROUVE DANS AGC')
               RETURN
 101        CONTINUE
C
C           LE COEFFICIENT DIAGONAL
            AGC( IC2 ) = AG( NDIAG )
C
         ELSE
C
C           LIGNE D'UN DEGRE DE LIBERTE FIXE =>
C           LIGNE RAMENEE A 1 COEFFICIENT DIAGONAL DE VALEUR 1
C           LES AUTRES COEFFICIENTS SONT NULS
            NDIAGC = LPLIGC(I+1)
            DO 104 KC = LPLIGC(I)+1, NDIAGC-1
               AGC( KC ) = 0D0
 104        CONTINUE
            AGC( NDIAGC ) = 1D0
C
         ENDIF
 90   CONTINUE
      ID0 = ID1
C
 100  CONTINUE
C
C     REMISE A ZERO DES DL PRESSION LIBRES
      DO 200 I=1,NTDL
         IF( NTDLFX(I) .LT. 0 ) NTDLFX(I) = 0
 200  CONTINUE
C
c
      print *
      print *,'MATRICE AGC en sortie de INCRSP'
      DO i=1,ntdl
         print 12345,(i,LPCOLC(m),agc(m),m=LPLIGC(i)+1, LPLIGC(i+1))
      enddo
12345 format(5('  agc(',i3,',',i3,')=',D15.6))
C
      RETURN
      END
