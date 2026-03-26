C     FICHIER ~/incl/pp.inc
C**************************************************************************************
C     FICHIER DEPENDANT MACHINE: VERSION STANDARD PC 32 ou 64 Bits
C
C     DEFINITION DU NOMBRE DE MOTS DU COMMON BLANC
C     MOTMCN = NOMBRE TOTAL DE MOTS DISPONIBLES DANS LE
C                SUPER  TABLEAU NUMERIQUE MCN DU COMMON BLANC
C     DANS LEQUEL TOUS LES TABLEAUX de MEFISTO SONT STOCKES
C     COMMON MCN( MOTMCN ) est utilise ensuite dans Mefisto
C
C     Pour - MC le nombre d'octets de la memoire centrale du PC
C            - SWAP le nombre d'octets de la zone de SWAP sur disque
C     il faut calculer MOTMCN = ((MC+SWAP)-16 000 000) / 4
C     ou ce qui revient au meme
C              MC + SWAP = STACK = 4 * MOTMCN + 16 000 000
C     avec 2 seuils
C     MOTMCN <=   532 000 000 SUR PC 32bits pour compilateur gfortran
C     MOTMCN <= 2 130 000 000 SUR PC 64bits pour compilateur gfortran m64
C     IMPOSSIBLE d'aller au dela sans modifier tout le code pour que les 
C     adresses memoire de Mefisto stockees dans des entiers puissent depasser
C     2**31=2 147 483 256  ce qui N'EST PAS LE CAS ACTUELLEMENT
C
C     QUELQUES VALEURS ANCIENNES:
C     PARAMETER (MOTMCN =   28 000 000) VERSION MC+SWAP= 128 Mega
C     PARAMETER (MOTMCN =  128 000 000) VERSION MC+SWAP= 528 Mega
C     PARAMETER (MOTMCN =  232 000 000) VERSION MC+SWAP=1024 Mega
C     PARAMETER (MOTMCN =  532 000 000) VERSION MC+SWAP>=2144 Mega
C     PARAMETER (MOTMCN =   38 000 000) CYGWIN sous WINDOWS 512Mega
C     PARAMETER (MOTMCN =  246 000 000) CYGWIN sous WINDOWS 1  Giga
C     PARAMETER (MOTMCN =  232 000 000) PC=AP        1GO MC + 0.98G SWAP
C     PARAMETER (MOTMCN =  532 000 000) PC=ALPE    0.9GO MC + 1.7G  SWAP
C     PARAMETER (MOTMCN =  532 000 000) PC=ALAIN   3.1GO MC + 1.9G  SWAP
C     PARAMETER (MOTMCN = 1532 000 000) PC=AP         8GO MC + 8Giga SWAP
C     PARAMETER (MOTMCN = 1532 000 000) PC=AP        32GO MC + 1Giga SWAP

      INTEGER MOTMCN
      PARAMETER (MOTMCN = 1 532 000 000)
