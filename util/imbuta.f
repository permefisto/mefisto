      SUBROUTINE IMBUTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPRIMER LE NOMBRE ET LES DESCRITEURS DE BUFFERS DU
C ----- REPERTOIRE RETAMS DES TABLEAUX MS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  DECEMBRE 1983
C.......................................................................
      include"./incl/motmcg.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
C
      IF( MGBUTA .GT. 0 ) THEN
         WRITE(IMPRIM,10000) NBBUTA,(MCG(MGBUTA+I),I=0,3*NBBUTA-1)
      ENDIF
10000 FORMAT(' NOMBRE DE BUFFERS DES TABLEAUX MS=',I7/
     %(' ADRESSE MCG DU BUFFER=',I10,' NO PAGE DANS LE BUFFER=',I7,
     %' NO DERNIER APPEL A CE BUFFER=',I12))
      RETURN
      END
