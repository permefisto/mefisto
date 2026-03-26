      SUBROUTINE FINDEL( NOMTAB, IA, NFPOBA, MOPAGE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: RECHERCHE DANS LE COMMON ELFINI QUI CONTIENT LE PREMIER
C ---  ENREGISTREMENT DU FICHIER ACCES DIRECT NFPOBA DU TABLEAU
C      DE NOM NOMTAB ET RESTAURATION DANS M A L ADRESSE IA
C
C PARAMETRES D ENTREE :
C ----------------------
C NOMTAB : NOM DU TABLEAU
C NFPOBA : NO DU FICHIER D ACCES DIRECT SUPPORT DES
C          VALEURS DES POLYNOMES DE BASE AUX POINTS D INTEGRATION
C MOPAGE : NOMBRE DE MOTS D'UNE PAGE DU FICHIER NFPOBA
C
C PARAMETRE  RESULTAT :
C ---------------------
C IA     : ADRESSE DANS M DU TABLEAU NOMTAB APRES RESTAURATION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN  PERRONNET ANALYSE NUMERIQUE UPMC PARIS  MAI 1989
C ......................................................................
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON            M(1)
      COMMON / ELFINI / NELEXI,NOPAGE,NELFI(510)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*4       CHARX,CAR4
C
10001 FORMAT(1X,130(1H%)/' ERREUR FINDEL: Le TABLEAU DE NOM:',A,
     &' N EST PAS DECRIT dans le COMMON ELFINI'/1X,130(1H%)/
     &' NO FICHIER POBA =',I5/
     &' COMMON ELFINI: NB TABLEAUX=',I7,' NO 1-ERE PAGE LIBRE=',I7/
     &5(1X,A4,I7,2I5,4X))
10003 FORMAT(' Les',I6,' MOTS DU TABLEAU ',A,' SONT RESTAURES DANS MCN')
20003 FORMAT(' The',I6,' WORDS of the ARRAY ',A,' are RESTORED in MCN')
C
C     RECHERCHE DU NOM DU TABLEAU
C     ---------------------------
      NB=4*(NELEXI-1)+1
      DO 1 I=1,NB,4
         IF( NELFI(I) .EQ. NOMTAB ) THEN
             IA1=I
             GO TO 2
         ENDIF
    1 CONTINUE
C
C     LE TABLEAU N A PAS ETE RETROUVE
C     -------------------------------
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TABLEAU '//CHARX(NOMTAB)//' INCONNU'
      ELSE
         KERR(1) = 'ARRAY '//CHARX(NOMTAB)//' UNKNOWN'
      ENDIF
      CALL LEREUR
      WRITE(IMPRIM,10001) CHARX(NOMTAB),NFPOBA,NELEXI,NOPAGE,NELFI
      RETURN
C
C     LE TABLEAU EST RETROUVE. ADRESSAGE
C     ----------------------------------
    2 IA    = NELFI(IA1+1)
      L     = NELFI(IA1+2)
      IF( IA .LE. 0 ) THEN
C        LE TABLEAU NON PRESENT DANS MCN EST RESTAURE
C        A PARTIR DE NFPOBA
         CALL TNMCDC( 'REEL2' , L/2 + 1 , IA )
         NPAG         = NELFI(IA1+3)
         NOPAG        = NPAG-1
         NELFI(IA1+1) = IA
C        LECTURE DU TABLEAU NOMTAB
         CALL ESTASF(IA,L,1,NFPOBA,MOPAGE,NOPAG)
C
C        IMPRESSION DE LA RESTAURATION
C        -----------------------------
         CAR4 = CHARX( NOMTAB )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10003) L,CAR4
         ELSE
            WRITE(IMPRIM,20003) L,CAR4
         ENDIF
      ENDIF
      END
