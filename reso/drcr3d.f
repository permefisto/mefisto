      SUBROUTINE DRCR3D( MUDL, A, B0, NTDL, NDSM, NIVEAU,   B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DESCENTES ET/OU REMONTEES D'UN SYSTEME DONT LA MATRICE EST
C ----- FACTORISEE SOUS LA FORME A = L * U DITE DE GAUSS
C       L EST UNE MATRICE TRIANGULAIRE INFERIEURE A DIAGONALE UNITE
C       U EST UNE MATRICE TRIANGULAIRE SUPERIEURE
C       L ET U SONT RANGEES DANS A SOUS FORME PROFIL
C       A(I,J) NON NUL => A(J,I) NON NUL  (PROFIL SYMETRIQUE)
C       A NON SYMETRIQUE MUDL POINTE SUR LES COEFFICIENTS DIAGONAUX
C       LE SP CRMC3D A DU ETRE EXECUTE AUPARAVANT
C       VERSION REELLE DOUBLE PRECISION

C PARAMETRES D'ENTREE:
C --------------------
C MUDL  : MUDL(1)=0,MUDL(I+1)= ADRESSE DU I-EME COEF DIAGONAL
C A     : MATRICE FACTORISEE A=L*U
C B0    : LES NDSM SECONDS MEMBRES  B0(NDSM,NTDC)
C NTDL  : ORDRE DE LA MATRICE A
C NDSM  : NBRE DE SECONDS MEMBRES
C NIVEAU: 0  U*B=B0    REMONTEES SEULES EXACTES SEULEMENT SI B = B0
C         1  L*B=B0    DESCENTES SEULES
C         2  L*U*B=B0  DESCENTES,REMONTEES

C PARAMETRES  DE  SORTIE:
C -----------------------
C B     :  TABLEAU DES RESULTATS  B (NDSM,NTDL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET  LAN 189 ET IRIA PARIS               OCTOBRE 1977
C MODIFS : A.PERRONNET  Saint Pierre du Perray             NOVEMBRE 2020
C.......................................................................
      DOUBLE PRECISION A(*), B(*), B0(*), S
      INTEGER          MUDL(*)

C     COHERENCE DU PARAMETRE NIVEAU
C     -----------------------------
      IF( NIVEAU.GE.0  .AND.  NIVEAU.LE.2 ) GOTO 1

      PRINT 2,NIVEAU
    2 FORMAT(1X,130(1H%)/' ERREUR:DANS DRCR3D (A NON SYMETRIQUE) NIVEAU=
     %',I12,' VALEUR INTERDITE'/1X,130(1H%))

      STOP

C     LES DESCENTES
C     -------------
    1 IF( NIVEAU .EQ. 0 ) GOTO 10

      IA=0
      MUDLI=0
      DO I=1,NTDL

         I1=I-1
         MUDLI1=MUDL(I+1)
         IH=(MUDLI1-MUDLI)/2
         IMI= I-IH
         IA0= IA
         IB0=(IMI-1) *NDSM

         DO NC=1,NDSM
            S=0.D0
            IA=IA0
            IB=IB0+NC

            DO K=1,IH
               S=S+A(IA+K)*B(IB)
               IB=IB+NDSM
            ENDDO

            B(IB)=B0(IB)-S
         ENDDO

         IA=IA+MUDLI1-MUDLI
         MUDLI=MUDLI1

      ENDDO

      IF( NIVEAU .EQ. 1 ) GOTO 9999

C     LES REMONTEES
C     -------------
 10   NTDL1=NTDL+1
      IB=NTDL1
      MUDLI1=MUDL(NTDL1)
      IA=MUDLI1+1
      DO II=1,NTDL

         I=NTDL1-II
         MUDLI=MUDL(I)
         IH=(MUDLI1-MUDLI)/2
         IA=IA-1
         IB=IB-1
         IBN=IB*NDSM
         S=1.D0/A(IA)

         DO NC=1,NDSM
            B(IBN) = B(IBN)*S
            IBN = IBN-1
         ENDDO

         IF( IH .LE. 0 ) GOTO 20
         IB0= IBN+NDSM+1

         DO K=1,IH
            IA=IA-1
            IBB=IB0
            S=A(IA)

            DO NC=1,NDSM
               B(IBN)=B(IBN)-B(IBB-NC)*S
               IBN=IBN-1
            ENDDO

         ENDDO

         IA=IA-IH
 20      MUDLI1=MUDLI

      ENDDO

 9999 RETURN
      END
