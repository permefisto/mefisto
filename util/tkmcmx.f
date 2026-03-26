      SUBROUTINE TKMCMX( MAXVAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOMBRE MAXIMAL DE CARACTERES DE LA PLUS GRANDE
C ----- ZONE LIBRE DU SUPER-TABLEAU MCK
C
C SORTIE :
C --------
C MAXVAR : NOMBRE DE CARACTERES CONTIGUS DANS LA PLUS
C          GRANDE ZONE LIBRE DU SUPER-TABLEAU MCK
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  NOVEMBRE 1983
C.......................................................................
      include"./incl/motmcg.inc"
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
C
C     RECHERCHE DU MAXIMUM
C     VARIABLES CARACTERES
      CALL TAMCMX( 'CARACTERE' , MCG(MGZLMK) , MAXVAR )
      RETURN
      END
