      SUBROUTINE CRMC3D( MUDL, A0, NTDL, EPS, NENTRE, IMPRE, A, NRETOU )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:   FACTORISATION  A=L.U POUR UNE MATRICE INVERSIBLE NON SYMETRIQUE
C ----   MAIS  A PROFIL SYMETRIQUE

C PARAMETRES D'ENTREE:
C --------------------
C MUDL  : MUDL(1)=0
C         MUDL(I+1)=ADRESSE DU I-EME COEFFICIENT DIAGONAL
C A0    : MATRICE AVANT FACTORISATION
C NTDL  : ORDRE DE LA MATRICE A0 ET A
C EPS   : PRECISION DE LA FACTORISATION
C NENTRE:=0 RETOUR AU PROGRAMME DES DETECTION D'UN PIVOT < EPS
C        =1 LES CALCULS SE POURSUIVENT  DANS TOUS LES CAS

C PARAMETRES DE SORTIE:
C ---------------------
C A     : MATRICE L ET  U  APRES FACTORISATION
C NRETOU: =0 SI AUCUN PIVOT  N'EST < EPS EN RELATIF
C         =1 SINON

C REMARQUES:
C ----------
C PIVOT =(A0(IDIAGONAL)-SA)/A0(IDIAGONAL)
C A0 ET A PEUVENT  ETRE CONFONDUES A L APPEL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR: A.PERRONNET  LAN189 ET IRIA  PARIS           OCTOBRE 1977
C MODIFS     : A.PERRONNET  Saint Pierre du Perray         NOVEMBRE 2020
C ......................................................................
      DOUBLE PRECISION  A0(*), A(*)
      INTEGER           MUDL(*)
      DOUBLE PRECISION  AA, S, S1, SAU, DABS, SA
      COMMON /UNITES/ LECTEU,IMPRIM,FILUNI(30)

   15    FORMAT(131('%')/' RESULTATS A VERIFIER: LE',I7,' COEFFICIENT
     +   DIAGONAL EST EN RELATIF AU DESSOUS DE LA PRECISION '/
     +   'EPS=',G15.7,' SAU=',G15.7,' A0= ',G15.7/131('%'))

      NRETOU=0
      MUDLI=0
      DO 1 I=1,NTDL
C        I1=I-1    OPTIM

C        TRAITEMENT DE LA I-EME LIGNE DE L ET COLONNE DE U     (A=L*U)

         MUDLI1=MUDL(I+1)
         IH=(MUDLI1-MUDLI)/2
C        IH=NOMBRE DE COEFFICIENTS DE LA I-EME LIGNE (DIAGONALE EXCLUE)
         SA=0.D0
         IF( IH .LE. 0 ) GOTO 2
         IMI=I-IH
C        IMI=N0 DE LA 1-ERE COLONNE  NON NULLE DE LA I-EME LIGNE
         IA0=MUDLI
         IA1=MUDLI+IH
         MUDLJ=MUDL(IMI)
         IMI1=IMI-1

C        BOUCLE SUR LES COLONNES  DE LA I-EME  LIGNE DE L
C        BOUCLE SUR LES LIGNES    DE LA I-EME  COLONNE DE U
         DO 4 JJ=1,IH
            J =IMI1+JJ
C           J1=J-1   OPTIM
            MUDLJ1=MUDL(J+1)
            JH=(MUDLJ1-MUDLJ)/2
            JMI=J-JH
            JA0=MUDLJ
            JA1=MUDLJ+JH
C           JA0 POINTE EN TETE DE LA J-EME LIGNE
C           JA1 POINTE EN TETE DE LA J-EME COLONNE
            IS=IMI-JMI
            IF( IS .GT. 0 ) GOTO 6

            MA=JMI
            IA=IA0-IS
            JA=JA0
            IAA=IA1-IS
            JAA=JA1
            GOTO 7

 6          MA=IMI
            IA=IA0
            JA=JA0+IS
            IAA=IA1
            JAA=JA1+IS

 7          S=0.D0
            S1=0.D0

            K1=J-MA
            IF( K1 .LE. 0 ) GOTO 9
C
            DO 10 K=1,K1
C              L(I,J)
               S=S+A(IA+K)*A(JAA+K)
C              U(J,I)
               S1=S1+A(JA+K)*A(IAA+K)
 10         ENDDO
C
 9          IA=IA+K1+1
            IAA=IAA+K1+1
            A(IA)=(A0(IA)-S)/A(MUDLJ1)
            A(IAA)=A0(IAA)-S1
            MUDLJ=MUDLJ1
 4       ENDDO
C
C        TRAITEMENT DU I-EME COEFFICIENT DIAGONAL
         DO K=1,IH
            SA=SA+A(IA0+K)* A(IA1+K)
         ENDDO

    2    AA=A0(MUDLI1)
         SAU=AA-SA

C        TEST SUR LA PRECISION DU PIVOT
         IF( DABS(SAU) .GT. EPS*DABS(AA) ) GOTO 14

C        AU DESSOUS DE LA PRECISION
         WRITE (IMPRIM,15) I,EPS,SAU,AA
         NRETOU=1
         IF( NENTRE .EQ. 0 ) GOTO 16

 14      A(MUDLI1)=SAU
         MUDLI=MUDLI1
 1    ENDDO

 16   NL=NTDL
ccc   IF(IABS(IMPRE)-6) 9999,18,19
      IF( IABS(IMPRE) .GT. 6 ) GOTO 19
      IF( IABS(IMPRE) .LT. 6 ) GOTO 9999
      NL=10
 19   CALL AFPROF( NL, MUDL, A )
C
 9999 RETURN
      END
