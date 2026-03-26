      SUBROUTINE MSDS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETRUIRE LA MEMOIRE SECONDAIRE
C ----- DETRUIRE TOUS LES FICHIERS NUMERIQUES ET LE SUPER-FICHIER
C
C ATTENTION : LA MEMOIRE SECONDAIRE DOIT AUPARAVANT ETRE OUVERTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS    DECEMBRE 1983
C.......................................................................
      include"./incl/ppmck.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/motmcg.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
C
C     LA MEMOIRE SECONDAIRE EST ELLE OUVERTE ?
      IF( MGFIMS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MSDS: LA MS N EST PAS OUVERTE'
         ELSE
            KERR(1) = 'MSDS: The MS IS NOT OPEN'
         ENDIF
         CALL LEREUR
         CALL ARRET(-1)
      ENDIF
C
C        LES FICHIERS DECLARES SONT DETRUITS
C        ===================================
         IG = MGFIMS
         DO 40 I=1,M2FIMS
            IF( MCG(IG) .GT. 0 ) THEN
C              LE FICHIER I EST DETRUIT
               CLOSE( UNIT=MCG(IG) , ERR=20 , STATUS='DELETE' ,
     %                IOSTAT=IOERR )
               GOTO 30
C
C              TRAITEMENT DE L'ERREUR DETECTEE PAR CLOSE
 20            NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') MCG(IG)
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'LE FICHIER '//KERR(MXLGER)(1:4)//
     %                   ' NE PEUT ETRE DETRUIT'
               ELSE
                  KERR(1) = 'The FILE '//KERR(MXLGER)(1:4)//
     %                   ' CAN NOT BE DELETED'
               ENDIF
               CALL LEREUR
               CALL IMFIMS
            ENDIF
C
 30         IG = IG + M1FIMS
 40      CONTINUE
C
C     LE SUPER-FICHIER EST DETRUIT
      CLOSE( UNIT=NOFISF , ERR=9900 , STATUS='DELETE' , IOSTAT=IOERR )
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'DESTRUCTION TOTALE DE LA MEMOIRE SECONDAIRE'
      ELSE
         KERR(1) = 'TOTAL DELETION of SECONDARY MEMORY'
      ENDIF
      CALL LEREUR
      GOTO 9999
C
C     NON FERMETURE DU SUPER-FICHIER
 9900 NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NOFISF
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'SUPER FICHIER NON FERMABLE '//KERR(MXLGER)(1:4)
      ELSE
         KERR(1) = 'SUPER FILE CAN NOT BE CLOSED '//KERR(MXLGER)(1:4)
      ENDIF
      CALL LEREUR
C
C     LA MS EST DETRUITE => ARRET NORMAL DU TRAVAIL
9999  STOP
      END
