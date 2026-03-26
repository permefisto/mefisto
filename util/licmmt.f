      SUBROUTINE LICMMT( NL0, NC0 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     EFFACER DE KLG LES COMMENTAIRES ENCADRES DE { { }} {}
C -----
C
C ENTREES ET SORTIES :
C --------------------
C NL0 NC0 : NUMERO DE LIGNE ET COLONNE   DU 1-ER CARACTERE A EXAMINER
C           NUMERO DE LIGNE ET COLONNE+1 DU } FINAL   EN SORTIE
C NL0     :  0 EN SORTIE SI UN PROBLEME A ETE RENCONTRE
C           >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     JUILLET 2002
C MODIFS : PERRONNET ALAIN LJLL UPMC & St PIERRE DU PERRAY  AVRIL   2013
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      NL1    = NL0
      NC1    = NC0

 1    NBPOUV = 0
      NBPFER = 0
      NLPO = 0
      NCPO = 0
      NLPF = 0
      NCPF = 0

C     NOMBRE DE } et {  DE LA LIGNE NL1 A PARTIR DE NC1
 5    DO NC = NC1, NBCALI

         IF( KLG(NL1)(NC:NC) .EQ. '{' ) THEN
C           RECHERCHE POSITION PREMIERE {
            IF( NBPOUV .EQ. 0 ) THEN
               NLPO = NL1
               NCPO = NC
            ENDIF
            NBPOUV = NBPOUV+1
         ENDIF

         IF( KLG(NL1)(NC:NC) .EQ. '}' ) THEN
            NBPFER = NBPFER+1
            NLPF   = NL1
            NCPF   = NC
         ENDIF

         IF( NBPOUV .GT. 0 . AND. NBPOUV .EQ. NBPFER ) GOTO 50

         IF( NBPOUV .LT. NBPFER ) GOTO 9000

      ENDDO
C
      IF( NBPOUV .EQ. 0 ) RETURN
C
C     IL EXISTE { { ... }  et } MANQUANTE
      IF( INTERA .GE. 3 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER),'(I3)') NBPOUV-NBPFER+1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='LU: il RESTE '//KERR(MXLGER)(1:3)
     %          //'{ COMMENTAIRE(S) a TERMINER par }'
         ELSE
            KERR(1) ='LU: there are'//KERR(MXLGER)(1:3)
     %          //'{ COMMENTS to finish by }'
         ENDIF
         CALL LERESU
      ENDIF
C
C     LECTURE D'UNE LIGNE SUPPLEMENTAIRE DE DONNEES DANS KLG(LHKLG)
      CALL LIRLIG( I )
      IF( I .GT.  0 ) GOTO 9010
C     SI LE CARACTERE D'ABANDON @ EXISTE DANS CETTE LIGNE
C     IL EST IGNORE CAR CONSIDERE COMME UN CARACTERE QUELCONQUE
      NL1 = LHKLG
      NC1 = 1
      GOTO 5

C     IL EXISTE UNE PAIRE { } FINISSANTE QUI EST EFFACEE
 50   NCF = 0
      DO NL = NLPO, NLPF
         IF( NL .EQ. NLPO ) THEN
            NCD = NCPO
         ELSE
            NCD = 1
         ENDIF
         IF( NL .EQ. NLPF ) THEN
            NCF = NCPF
         ELSE
            NCF = NBCALI
         ENDIF
         KLG(NL)(NCD:NCF) = ' '
      ENDDO

C     SUPPRESSION DES LIGNES BLANCHES DE KLG
      IF( NL0 .GT. 0 .AND. NBPOUV .GT. 0 ) THEN
         NBLD = 0
         DO NL = NLPO, NLPF
            IF( NUDCNB( KLG(NL) ) .EQ. 0 ) THEN
C              NOMBRE DE LIGNES DE KLG A DECALER
               NBLD = NBLD + 1
            ELSE
               IF( NBLD .GT. 0 ) THEN
C                 ECRASEMENT DE LA LIGNE
                  KLG(NL-NBLD) = KLG(NL)
               ENDIF
            ENDIF
         ENDDO
         LHKLG = MAX( 1, LHKLG - NBLD )
      ENDIF

C     PRESENCE IMMEDIATE de {?
C     RECHERCHE DU PROCHAIN CARACTERE NON BLANC
      NL1 = LHKLG
      NC1 = NCF
      CALL CARPNB( NL1, NC1 )
      IF( NL1 .GT. 0 ) THEN
         IF( KLG(NL1)(NC1:NC1) .EQ. '{' ) GOTO 1
      ENDIF
      RETURN
C
C     ERREUR
 9000 NBLGRC(NRERR) = 2 + LHKLG - NL0
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='UNE } SANS UNE { QUI LA PRECEDE dans LES LIGNES'
         DO NL=NL0,LHKLG
            KERR(2+NL-NL0) = KLG(NL)
         ENDDO
      ELSE
         KERR(1) ='ONE } WITHOUT a PREVIOUS { in LINES'
         DO NL=NL0,LHKLG
            KERR(2+NL-NL0) = KLG(NL)
         ENDDO
      ENDIF
      GOTO 9090
C
C     DEBORDEMENT DU TABLEAU KLG EN ATTENTE DE }
 9010 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='PARENTHESE OUVRANTE { NON TERMINEE PAR  } dans '
      ELSE
         KERR(1)='OPENED { NOT FINISHED BY  } in'
      ENDIF
      KERR(2) = KLG(LHKLG)
C
C     AFFICHAGE DE L'ERREUR
 9090 CALL LEREUR
      NL0 = 0

      RETURN
      END
