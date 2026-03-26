      SUBROUTINE XYZDSPT( RXYZ, NUMPT, XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER LE POINT DE COORDONNEES DXYZ PARMI LES POINTS ACTUELS
C -----
C ENTREE :
C --------
C RXYZ   : 3 COORDONNEES DU POINT A IDENTIFIER
C
C SORTIES:
C --------
C NUMPT  : 0 SI DXYZ N'EST PAS IDENTIFIE A L'UN DES POINTS
C         >0 NUMERO DU POINT IDENTIFIE
C XYZ    : 3 COORDONNES DU POINT NUMPT (NON MODIFIE SI NUMPT=0)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET St PIERRE DU PERRAY & LJLL UPMC  FEVRIER 2010
C23456---------------------------------------------------------------012
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1))
      REAL             RXYZ(3), XYZ(3)
C
C     RECHERCHE DE L'ADRESSE MCN DU LEXIQUE DES POINTS
      CALL TAMSOU( NTPOIN , MNPOIN )
C
C     LE DEBUT DU CHAINAGE DES POINTS OCCUPES DANS LE LEXIQUE
      NUMPT = MCN( MNPOIN + 5 )
C
C     NOMBRE D'ENTIERS POUR STOCKER UN NOM DE POINT
      NBENNM = MCN( MNPOIN + 2 )
C
C     LA BOUCLE SUR LES POINTS OCCUPES
C     ================================
 10   IF( NUMPT .GT. 0 ) THEN
C
C        ADRESSE MCN DU DEBUT DU POINT DANS LE LEXIQUE
         MNPT = MNPOIN + MCN(MNPOIN) * NUMPT
C
C        LE LEXIQUE DE CE POINT EXISTE-T-IL ?
         NTPT = MCN( MNPT + NBENNM + 2 )
C
         IF( NTPT .GT. 0 ) THEN
C
C           CE POINT EXISTE : RECHERCHE DU TMS 'XYZSOMMET'
            CALL LXTSOU( NTPT, 'XYZSOMMET', NTXYZ, MNXYZ )
            IF( NTXYZ .GT. 0 ) THEN
C
C              ADRESSE DE X DU POINT NUMPT
               MNP = MNXYZ + WYZSOM
C
C              IDENTIFICATION?
               CALL XYZIDE( RMCN(MNP), RXYZ, IDENTQ )
               IF( IDENTQ .EQ. 1 ) THEN
C
C                 POINT IDENTIFIE
                  XYZ(1) = RMCN(MNP  )
                  XYZ(2) = RMCN(MNP+1)
                  XYZ(3) = RMCN(MNP+2)
                  RETURN
C
               ENDIF
            ENDIF
         ENDIF
C
C        PASSAGE AU POINT SUIVANT
         NUMPT = MCN( MNPT + NBENNM )
         GOTO 10
C
      ENDIF
C
      RETURN
      END
