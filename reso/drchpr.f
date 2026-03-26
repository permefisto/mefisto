      SUBROUTINE DRCHPR( NTDL, NCODSA, LPDIAG, A, B0, NIVEAU,  B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: DESCENTE ET/OU REMONTEE D UN SYSTEME LINEAIRE FACTORISE SOUS
C ---- LA FORME  A * B = L * TL * B = B0  DITE DE CHOLESKY
C      L EST UNE MATRICE TRIANGULAIRE INFERIEURE
C      L EST RANGEE SOUS FORME PROFIL DANS LE TABLEAU A
C      ATTENTION: cholpr.f DOIT ETRE EXECUTE AVANT
C
C ENTREES:
C --------
C NTDL   : ORDRE DE LA MATRICE A
C NCODSA : =0 SI MATRICE DIAGONALE, >0 SI MATRICE SYMETRIQUE NON DIAGONALE
C LPDIAG : SI NCODSA>0 ALORS POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL DE A
C          LPDIAG(1) = 0, LPDIAG(I+1) = ADRESSE DANS A DU I-EME COEFFICIENT
C                                       DIAGONAL SI A NON DIAGONALE
C A      : MATRICE L RESULTAT DE LA FACTORISATION DE CHOLESKY ISSUE DE CHOLPR
C B0     : TABLEAU B0(NTDL) DU SECOND MEMBRE
C NIVEAU : =0  TL*B=B0   REMONTEES SEULEMENT VALABLE SI B=B0 SINON ERREUR
C          =1  L *B=B0   DESCENTES SEULEMENT
C          >1  L*TL*B=B0 DESCENTES ET REMONTEES
C
C SORTIES:
C --------
C B      : TABLEAU B(NTDL) DES SOLUTIONS
C          SI NIVEAU > 1 B0 ET B PEUVENT ETRE IDENTIQUES A L APPEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    AVRIL 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY    AVRIL 2013
C23456---------------------------------------------------------------012
C$    use OMP_LIB
      include"./incl/threads.inc"
      INTEGER           LPDIAG(0:NTDL)
      DOUBLE PRECISION  A(1:*), B(NTDL), B0(NTDL), S, PROSCD

      IF( NTDL .LE. MININDSTHR .OR. NTDL .LE. NBTHREADS*MININD1THR )THEN

C        -----------------------------------
C        EXECUTION SEQUENTIELLE SUR 1 THREAD
C        -----------------------------------
         IF( NCODSA .GT. 0 ) THEN

C           ================================
C           MATRICE NON DIAGONALE SYMETRIQUE
C           ================================
C
C           DESCENTE ET/OU REMONTEE?
            IF( NIVEAU .EQ. 0 ) GOTO 10
C
C           LA DESCENTE
C           -----------
            LPDM1 = 0
            DO I = 1, NTDL
C              NO DU COEFFICIENT DIAGONAL
               LPD = LPDIAG(I)
C              NOMBRE DE COEFFICIENTS DE LA LIGNE I
               IH  = LPD - LPDM1
               IB  = I - IH
               S   = 0.D0
               DO K = 1, IH-1
                  S = S + A(LPDM1+K) * B(IB+K)
               ENDDO
               B(I)  = ( B0(I) - S ) / A(LPD)
C              PASSAGE A LA LIGNE SUIVANTE DE LA MATRICE
               LPDM1 = LPD
            ENDDO
C
            IF( NIVEAU .EQ. 1 ) RETURN
C
C           LA REMONTEE
C           -----------
 10         LPD = LPDIAG(NTDL)
            DO I = NTDL, 1, -1
               LPDM1 = LPDIAG(I-1)
               IH    = LPD - LPDM1 - 1
C              LA SOLUTION I
               B(I)  = B(I) / A(LPD)
               DO K = 1, IH
                  IK = I - K
                  B(IK) = B(IK) - A(LPD-K) * B(I)
               ENDDO
               LPD = LPDM1
            ENDDO
C
            RETURN

         ELSE
C
C           =================
C           MATRICE DIAGONALE
C           =================
C
C           LA DESCENTE
C           -----------
            IF( NIVEAU .EQ. 0 ) GOTO 20
C
            DO I = 1, NTDL
               B(I) = B0(I) / A(I)
            ENDDO
C
C           LA REMONTEE
C           -----------
 20         IF( NIVEAU .EQ. 1 ) RETURN
C
            DO I = 1, NTDL
               B(I) = B(I) / A(I)
            ENDDO
C
            RETURN
         ENDIF

      ELSE

C        ---------------------------------
C        EXECUTION PARALLELE SUR NBTHREADS
C        ---------------------------------
         IF( NCODSA .GT. 0 ) THEN

C           ================================
C           MATRICE NON DIAGONALE SYMETRIQUE
C           ================================
C
C           DESCENTE ET/OU REMONTEE?
            IF( NIVEAU .EQ. 0 ) GOTO 100
C
C           LA DESCENTE
C           -----------
            LPDM1 = 0
            DO I = 1, NTDL
C              NO DU COEFFICIENT DIAGONAL
               LPD = LPDIAG(I)
C              NOMBRE DE COEFFICIENTS DE LA LIGNE I
               IH = LPD - LPDM1
               IB = I - IH
               S  = PROSCD( A(LPDM1+1), B(IB+1), IH-1 )
               B(I) = ( B0(I) - S ) / A(LPD)
C              PASSAGE A LA LIGNE SUIVANTE DE LA MATRICE
               LPDM1 = LPD
            ENDDO
C
            IF( NIVEAU .EQ. 1 ) RETURN
C
C           LA REMONTEE
C           -----------
 100        LPD = LPDIAG(NTDL)
            DO I = NTDL, 1, -1
               LPDM1 = LPDIAG(I-1)
               IH    = LPD - LPDM1 - 1
C              LA SOLUTION I
               B(I)  = B(I) / A(LPD)
               DO K = 1, IH
                  IK = I - K
                  B(IK) = B(IK) - A(LPD-K) * B(I)
               ENDDO
               LPD = LPDM1
            ENDDO

         ELSE

C           =================
C           MATRICE DIAGONALE
C           =================

C           LA DESCENTE
C           -----------
            IF( NIVEAU .EQ. 0 ) GOTO 200
C
C///////////////////////////////////////////////////////////////////////
C$OMP PARALLEL PRIVATE(I)
C$OMP DO SCHEDULE( STATIC, NTDL/NBTHREADS )
            DO I = 1, NTDL
               B(I) = B0(I) / A(I)
            ENDDO
C$OMP END DO
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////

C           LA REMONTEE
C           -----------
 200        IF( NIVEAU .EQ. 1 ) RETURN

C///////////////////////////////////////////////////////////////////////
C$OMP PARALLEL PRIVATE(I)
C$OMP DO SCHEDULE( STATIC, NTDL/NBTHREADS )
            DO I = 1, NTDL
               B(I) = B(I) / A(I)
            ENDDO
C$OMP END DO
C$OMP END PARALLEL
C///////////////////////////////////////////////////////////////////////

         ENDIF

      ENDIF

      RETURN
      END
