      SUBROUTINE BLDLGC( NTDL,   NDSM,   NBDLFX, NODLFX, VADLFX,
     %                   NCODSA, LPLIGN, LPCOLO, AG,     BG  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT : PRISE EN COMPTE DES C.L. DIRICHLET SUR AG MATRICE MORSE
C  ----- ET BG LES NDSM SECONDS MEMBRES
C        LES COEFFICIENTS DES LIGNES DES DL FIXES SONT ANNULES SAUF
C        LE COEFFICIENT DIAGONAL QUI DEVIENT 1
C
C       | AG 11  AG 12 | | XG 1 |   | BG 1 |
C                        |      | =            DEVIENT
C                        | XD 2 |
C
C       | AG 11    0   | | XG 1 |   | BG 1 - AG 12 x XD 2 |
C       |              | |      | = |                     |
C       |  0  Identite | | XG 2 |   | XD 2                |
C
C ENTREES:
C  -------
C  NTDL   : NOMBRE D'INCONNUES DU SYSTEME
C  NDSM   : NOMBRE DE SECONDS MEMBRES
C  NBDLFX : NOMBRE DE  DEGRES DE LIBERTE FIXES
C  NODLFX : NUMERO DES DEGRES DE LIBERTE FIXES
C  VADLFX : VALEURS IMPOSEES
C  NCODSA : CODE DE STOCKAGE DE LA MATRICE A
C           1     SYMETRIQUE MORSE
C          -1 NON SYMETRIQUE MORSE
C  LPLIGN : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE
C  LPCOLO : NUMERO DE COLONNE DE CHACUN DES COEFFICIENTS NON NULS
C  AG     : MATRICE MORSE   AG
C  BG     : SECONDS MEMBRES BG(NTDL,NDSM)
C
C SORTIES :
C  --------
C  AG     : EST MODIFIE EN SORTIE
C  BG     : EST MODIFIE EN SORTIE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS         OCTOBRE 1989
C MODIFS : ALAIN PERRONNET LABORATOIRE J-L.LIONS UPMC PARIS OCTOBRE 2007
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION  VADLFX(NDSM,NBDLFX), AG(*), BG(NTDL,NDSM)
      INTEGER           LPLIGN(NTDL+1), LPCOLO(*), NODLFX(NBDLFX)
C
      IF ( NBDLFX .LE. 0 ) RETURN
C
      IF ( NCODSA .EQ. 1 ) THEN
C
C        -----------------------------------
C        ------   MATRICE SYMETRIQUE  ------
C        -----------------------------------
C        PLUS GRANDE DIFFERENCE ENTRE NO DE COLONNES SUR UNE LIGNE
         MXDIFV = LARVOI( NTDL, LPLIGN, LPCOLO )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'MAX DES DIFFERENCES ENTRE COLONNES d''une L
     %IGNE de la MATRICE=',MXDIFV
         ELSE
            WRITE(IMPRIM,*) 'MAX DIFFERENCES between COLUMNS of a LINE o
     %f the MATRIX=',MXDIFV
         ENDIF
C
         DO 1 N=1,NBDLFX
            NODL = NODLFX(N)
            K1   = LPLIGN(NODL) + 1
            KK2  = LPLIGN(NODL+1)
            K2   = KK2 - 1
            DO K=K1,K2
C     ---   CONTRIBUTION DE V AUX SECONDS MEMBRES   ---
               DO L=1,NDSM
                  BG(LPCOLO(K),L) = BG(LPCOLO(K),L)
     %                            - AG(K) * VADLFX(L,N)
               ENDDO
C     ---   MISE A ZERO DE A(NODL,J) J=1,...,NODL   ---
               AG(K) = 0.D0
            ENDDO
C     ---   MISE A ZERO DE A(J,NODL) J=NODL+1,...   ---
            DO 4 I=NODL+1, MIN(NTDL,NODL+MXDIFV)

               K1 = LPLIGN(I) + 1
               K2 = LPLIGN(I+1)
               DO 5 K=K1,K2

ccc               IF(LPCOLO(K)-NODL) 5,6,4

                  IF( LPCOLO(K) .GT. NODL ) GOTO 4
                  IF( LPCOLO(K) .EQ. NODL ) THEN
                     DO L=1,NDSM
                        BG(I,L) = BG(I,L) - AG(K) * VADLFX(L,N)
                     ENDDO
                     AG(K) = 0.D0
                     GOTO 4
                  ENDIF

 5             ENDDO

 4          ENDDO
C     ---   LE COEFFICIENT DIAGONAL   ---
            DO L=1,NDSM
               BG(NODL,L) = VADLFX(L,N)
            ENDDO
            AG(KK2) = 1.D0
 1       ENDDO
C
      ELSE
C
C        ----------------------------------------
C        ------   MATRICE NON SYMETRIQUE   ------
C        ----------------------------------------
         DO N=1,NBDLFX
            NODL = NODLFX(N)
            K1   = LPLIGN(NODL) + 1
            KK2  = LPLIGN(NODL+1)
            K2   = KK2 - 1
C     ---   MISE A ZERO DE A(NODL,J) J=1,...,NTDL   ---
            DO K=K1,K2
               AG(K) = 0.D0
            ENDDO
C     ---   LE SECOND MEMBRE   ---
            DO L=1,NDSM
               BG(NODL,L) = VADLFX(L,N)
            ENDDO
C     ---   LE COEFFICIENT DIAGONAL   ---
            AG(KK2) = 1.D0
         ENDDO
      ENDIF
C
      RETURN
      END
