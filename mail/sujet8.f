      SUBROUTINE SUJET8(PARM,NBEFOB,NUSSOS,NVSOB,NVNUM)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :
C -----       1)  TETRAEDRISATION D'UN MAILLAGE
C                          PENTAEDRIQUE-TETRAEDRIQUE
C             2)  PENTAEDRISATION D'UN MAILLAGE
C                              HEXAEDRIQUE-PENTAEDRIQUE
C             3)  TETRAEDRISATION  "      "          "
C ENTREES:
C --------
C PARM   : = 1 --> 1)
C          = 2 --> 2)
C          = 3 --> 3)
C NBEFOB : NOMBRE DE S-OBJS DU MAILLAGE
C NUSSOS : NUSSOS(I,J) : NUMERO DU I EME SOMMET DU J EME S-OBJ
C
C SORTIES:
C --------
C NVSOB : NOMBRE DE NOUVEAUX S-OBJS
C NVNUM : NOUVEAU TABLEAU DE NUMEROS
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : DACHRAOUI TAREK
C ===================================================================
      INTEGER PARM,NBEFOB
     %      ,NUSSOS(8,NBEFOB)
     %      ,NVSOB,NVNUM(8,NVSOB)
     +      , NVSOB1,NVNUM1(8,600)
      IF (PARM.EQ.1) THEN
        CALL TETPEN(NBEFOB,NUSSOS,NVSOB,NVNUM)
      ELSEIF (PARM.EQ.2) THEN
        CALL PENHEX(NBEFOB,NUSSOS,NVSOB,NVNUM)
      ELSEIF (PARM.EQ.3) THEN
        CALL PENHEX(NBEFOB,NUSSOS,NVSOB1,NVNUM1)
        CALL TETPEN(NVSOB1,NVNUM1,NVSOB,NVNUM)
      ELSE
        WRITE(*,*) ' PARAMETRE ILLEGAL :'
      ENDIF
      END
