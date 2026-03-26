C     FICHIER ~/incl/ppmck.inc
C**************************************************************************************
C     FICHIER DEPENDANT MACHINE: VERSION STANDARD PC 4 CARACTERES PAR MOT
C
C     DEFINITION DU COMMON / MSMCK  /  ZONE DES CARACTERES DE MEFISTO
C
      PARAMETER  (NC1MCK=48 ,  MOTMCK=1024)
      CHARACTER*(NC1MCK)  MCK(MOTMCK)
      COMMON / MSMCK  /   MCK
