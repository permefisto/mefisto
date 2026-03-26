      SUBROUTINE COMMOR( NTLXOB , KNOMS1 , KNOMS2 , MNINSD )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : COMPRIMER UNE STRUCTURE MORSE EN NE GARDANT COMME COLONNES
C ----- QUE LES D.L. LIES AUX NOEUDS DE L'INTERFACE
C ENTREES :
C ---------
C NTLXOB : NUMERO DU LEXIQUE DE L'OBJET
C KNOMS1 : NOM DE LA STRUCTURE A COMPRIMER
C MNINSD : ADRESSE MCN DU TABLEAU DES D.L. DE L'INTERFACE
C
C SORTIES :
C ---------
C KNOMS2  : NOM DE LA STRUCTURE COMPRIMEE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___morse.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*12      KNOMS1,KNOMS2
      CHARACTER*24      KNOMS
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     OUVERTURE DE LA MATRICE MORSE
C     -----------------------------
C
      KNOMS = 'MORSE"'//KNOMS1
      L = INDEX( KNOMS , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMS )
      ENDIF
      CALL LXTSOU( NTLXOB , KNOMS(1:L) , NTMORS , MNMORS )
C
C     LES TMS 'MORSE"KNOMST'
C     ----------------------
C
      NTDL = MCN(MNMORS+WBLIMO)
      LPMORS = MCN( MNMORS + WPMORS )
      MNLPLI = MNMORS + WPLIGN
      MNLPCO = MNLPLI + NTDL + 1
      MNAG   = MNMORS + LPMORS
C     INITIALISATION DES TABLEAUX
C     COMPRESSION DE LA STRUCTURE MORSE
      NCODSA = 1
      MNINS1 = MNINSD+1
      CALL COMMO1( NCODSA, MCN(MNLPLI), MCN(MNLPCO), MCN(MNAG) )
C
C     DECLARATION DE LA NOUVELLE STRUCTURE MORSE
C     ------------------------------------------
C
      KNOMS = 'MORSE"'//KNOMS2
      L = INDEX( KNOMS , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMS )
      ENDIF
      CALL LXTSOU( NTLXOB , KNOMS(1:L) , NTMORS , MNMORS )
      IF( NTMORS .GT. 0 ) THEN
C        LA MATRICE MORSE EST DETRUITE POUR ETRE REDECLAREE
         CALL LXTSDS( NTLXOB , KNOMS(1:L) )
      ENDIF
      LOLPCO = MCN(MNLPLI+NTDL)
      LPMORS = WPLIGN + 1 + NTDL + LOLPCO
      IF( MOD(LPMORS,MOREE2) .EQ. 1 ) LPMORS = LPMORS + 1
      LO = LPMORS / MOREE2 + LOLPCO
      CALL LXTNDC( NTLXOB , KNOMS(1:L) , 'REEL2' , LO )
      CALL LXTSOU( NTLXOB , KNOMS(1:L) , NTMORS , MNMORS )
C     COPIE DU TABLEAU LPLIGN
      MNLPLS = MNMORS + WPLIGN
      CALL TRTATA( MCN(MNLPLI) , MCN(MNLPLS) , NTDL+1 )
C     COPIE DU TABLEAU LPCOLO
      MNLPCS = MNLPLS + NTDL + 1
      CALL TRTATA( MCN(MNLPCO) , MCN(MNLPCS) , LOLPCO )
C     COPIE DE LA MATRICE
      MNAGS  = MNMORS + LPMORS
      CALL TRTATA( MCN(MNAG) , MCN(MNAGS) , LOLPCO*MOREE2 )
C
C     MISE A JOUR DES TMS 'MORSE"KNOMST'
C     ----------------------------------
C
      MCN( MNMORS + WBLIMO ) = NTDL
      MCN( MNMORS + WPMORS ) = LPMORS
C     LA DATE
      CALL ECDATE( MCN(MNMORS) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNMORS + MOREE2 ) = NONMTD( '~>>>MORSE' )
C
C     DESTRUCTION DE LA STRUCTURE INITIALE
C     ------------------------------------
C
      KNOMS = 'MORSE"'//KNOMS1
      L = INDEX( KNOMS , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMS )
      ENDIF
      CALL LXTSOU( NTLXOB , KNOMS(1:L) , NTMORS , MNMORS )
      CALL LXTSDS( NTLXOB , KNOMS(1:L) )
C
      RETURN
      END
