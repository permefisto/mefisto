      SUBROUTINE SDDEF3( KNOMOB , KNOMTY , LISTE , KPBMCL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LECTURE AU NIVEAU DE L'OBJET DES DONNEES AUX LIMITES
C -----    DANS LE CAS D'UNE RESOLUTION PAR SOUS-DOMAINES
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET ASSOCIE AUX DONNEES
C KNOMTY : TYPE DE L'OBJET ASSOCIE AUX DONNEES
C LISTE : LISTE DES MOTS CLES
C KPBMCL : NOM DU TABLEAU DES C.L.
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C2345X7..............................................................012
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/gsmenu.inc"
      CHARACTER*(*)     KNOMOB,KNOMTY
      CHARACTER*(*)     LISTE(1:*),KPBMCL
      CHARACTER*24      KNOMTD
      CHARACTER*80      NOMTS
C
      IERR = 0
C
C     RECHERCHE DU 1-ER BLANC DANS LES NOMS
      L1 = INDEX( KNOMOB , ' ' )
      IF( L1 .GT. 0 ) THEN
         L1 = L1 - 1
      ELSE
         L1 = LEN( KNOMOB )
      ENDIF
      L2 = INDEX( KNOMTY , ' ' )
      IF( L2 .GT. 0 ) THEN
         L2 = L2 - 1
      ELSE
         L2 = LEN( KNOMTY )
      ENDIF
C
C     IMPRESSION D'UN MESSAGE
110   NBLGRC(NRHIST) = 1
      KHIST(1) = KNOMTY // KNOMOB
      CALL LHISTO
C
C     LE NOM DU TABLEAU A ENTRER
      CALL LIMTCL( KPBMCL , NMTCL )
      IF( NMTCL .LE. 0 ) RETURN
C
C     LE NOM DU TABLEAU TMS A REMPLIR
      NOMTS = '~>' // KNOMTY(1:L2) // '>' // KNOMOB(1:L1) // '>'
     %            // LISTE(NMTCL)
      L3 = INDEX( NOMTS , ' ' )
      IF( L3 .GT. 0 ) THEN
         L3 = L3 - 1
      ELSE
         L3 = LEN( NOMTS )
      ENDIF
      KNOMTD = '~>>>' // LISTE(NMTCL)
      CALL MOTSTD( KNOMTD , NOMTS(1:L3) , IERR )
      IF( IERR .EQ. 0 ) CALL AFTSTD( NOMTS(1:L3) )
      GOTO 110
C
      END
