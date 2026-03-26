      SUBROUTINE MAJXYZNSE( NBSOM,  XYZSOM, NEWXYZ,
     %                      L1SOEF, NBEFOB, NOSOEF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MISE A JOUR DU TABLEAU XYZSOM ET NOSOEF EN
C -----    IDENTIFIANT LES SOMMETS PROCHES
C          RENUMEROTANT LES SOMMETS
C          ELIMINANT LES EF DESACTIVES ou AYANT 2 SOMMETS DE MEME NUMERO
C                    ou TRIANGLES REDUITS a une ARETE

C ENTREE :
C --------
C L1SOEF : NOMBRE DE MOTS STOCKES POUR CHAQUE EF DANS LE TABLEAU NOSOEF

C MODIFIES:
C ---------
C NBSOM  : NOMBRE DE SOMMETS AVANT ET APRES
C XYZSOM : 3 COORDONNEES DES NBSOM  SOMMETS
C NEWXYZ : NOUVEAU NUMERO DES SOMMETS APRES RECENSEMENT
C          TABLEAU AUXILIAIRE
C NBEFOB : NOMBRE D'EF AVANT ET APRES
C NOSOEF : NUMERO DES NBSOEF SOMMETS DES NBEFOB EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LJLL UPMC et St PIERRE DU PERRAY     Octobre 2015
C2345X7..............................................................012
      PARAMETER (MXCOSTID=1024)
      INTEGER    NOSOEF(L1SOEF,NBEFOB), NEWXYZ(0:NBSOM),
     %           NOCOSTID(2,MXCOSTID)
      REAL       XYZSOM(3,NBSOM), VECNOR(3)

      NBSOM0  = NBSOM
      NBEFOB0 = NBEFOB

C     NOUVEAU NUMERO DES SOMMETS
      DO NS = 0, NBSOM 
         NEWXYZ( NS ) = 0
      ENDDO

C     RECENSEMENT DES SOMMETS ACTIFS et IDENTIFICATION selon EPZERO EPSXYZ
C     --------------------------------------------------------------------
      NBCOSTID = 0
      DO 30 NEF = 1, NBEFOB

         IF( NOSOEF(1,NEF) .GT. 0 ) THEN

C           NBSTEF: NOMBRE DE SOMMETS DE L'EF NEF
            DO NBSTEF = L1SOEF, 1, -1
               IF( NOSOEF( NBSTEF, NEF ) .NE. 0 ) GOTO 5
            ENDDO
            GOTO 30

C           ESSAI D'IDENTIFICATION DES SOMMETS de l'EF NEF
 5          DO K = 1, NBSTEF-1

               NS = NOSOEF( K, NEF )
               IF( NS .GT. 0 ) THEN

                  DO 20 I = K+1, NBSTEF

                     NS2 = NOSOEF( I, NEF )
                     IF( NS2 .LE. 0 .OR. NS2 .EQ. NS ) GOTO 20

                     CALL XYZIDE( XYZSOM(1,NS), XYZSOM(1,NS2), IDENTQ )
C                    IDENTQ : 1 SI LES 2 POINTS SONT JUGES IDENTIQUES
C                             0 SINON
                     IF( IDENTQ .EQ. 1 ) THEN

C                       LE COUPLE DE SOMMETS NS-NS2 EST IDENTIFIE
C                       LE PLUS GRAND NO EST REMPLACE PAR LE PLUS PETIT
                        IF( NS .LE. NS2 ) THEN
                           NSP = NS
                           NSG = NS2
                        ELSE
                           NSP = NS2
                           NSG = NS
                        ENDIF
C                       CE COUPLE EXISTE T IL DEJA?
                        DO N = 1, NBCOSTID
                           IF( NSP .EQ. NOCOSTID(1,N) .AND.
     %                         NSG .EQ. NOCOSTID(2,N) ) GOTO 20
                        ENDDO

C                       NON: IL EST AJOUTE
                        IF( NBCOSTID .GE. MXCOSTID ) THEN
                         PRINT*,'majxyzse: AUGMENTER MXCOSTID=',MXCOSTID
                           GOTO 20
                        ENDIF

C                       UN COUPLE DE SOMMETS IDENTIFIE DE PLUS
                        NBCOSTID = NBCOSTID + 1
C                       LE PLUS GRAND NO SERA REMPLACE PAR LE PLUS PETIT
                        NOCOSTID(1,NBCOSTID) = NSP
                        NOCOSTID(2,NBCOSTID) = NSG

                     ENDIF
 20               ENDDO
               ENDIF

            ENDDO

C           AJOUT DES NO DES SOMMETS DE L'EF DANS NEWXYZ
            DO K= 1, NBSTEF
               NS = NOSOEF( K, NEF )
               NEWXYZ( NS ) = NS
            ENDDO

         ENDIF
 30   ENDDO

C     SUPPRESSION DES SOMMETS IDENTIFIES   (LE PLUS GRAND NO DES 2)
C     ----------------------------------
      IF( NBCOSTID .GT. 0 ) THEN
         PRINT*
         PRINT*,'majxyznse:',NBCOSTID,' COUPLES DE SOMMETS IDENTIFIES'
         DO N = 1,NBCOSTID
            PRINT*,'majxyznse: (', NOCOSTID(1,N),',',NOCOSTID(2,N),')'
         ENDDO
         DO N = 1,NBCOSTID
            NEWXYZ( NOCOSTID(2,N) ) = -NOCOSTID(1,N)
         ENDDO
      ENDIF

C     RENUMEROTATION DES SOMMETS ACTIFS
C     ---------------------------------
      NBS = 0
      DO NS=1,NBSOM
         IF( NEWXYZ( NS ) .GT. 0 ) THEN
C           UN SOMMET ACTIF DE PLUS AVEC L'IDENTIFICATION
            NBS = NBS + 1
            NEWXYZ( NS ) = NBS
         ENDIF
      ENDDO

C     MISE A JOUR DU NO DES SOMMETS DU TABLEAU NOSOEF
C     -----------------------------------------------
      NBEF = 0
      DO 50 NEF = 1, NBEFOB

         IF( NOSOEF(1,NEF) .GT. 0 ) THEN

C           NBSTEF: NOMBRE DE SOMMETS DE L'EF NEF
            DO NBSTEF = L1SOEF, 1, -1
               IF( NOSOEF( NBSTEF, NEF ) .NE. 0 ) GOTO 35
            ENDDO
            GOTO 50

C           LE NOUVEAU NO DES SOMMETS DE L'EF NBEF
 35         NBEF = NBEF + 1
            DO K = 1, NBSTEF
C              NUMERO DU SOMMET AVANT RENUMEROTATION
               NS = ABS( NOSOEF( K, NEF ) )
C              NUMERO DU SOMMET APRES RENUMEROTATION
               NSNEW = NEWXYZ( NS )
 40            IF( NSNEW .LT. 0 ) THEN
C                 SOMMET IDENTIFIE
                  NSNEW = NEWXYZ( -NSNEW )
C                 SOMMET PEUT ETRE LUI MEME IDENTIFIE
                  GOTO 40
               ENDIF
               NOSOEF( K, NBEF ) = NSNEW
            ENDDO

            IF( NBSTEF .EQ. 3 ) THEN
C              NBEF EST UN TRIANGLE. EST IL DEGENERE?
               CALL NORFA3( XYZSOM( 1, NOSOEF(1,NBEF) ),
     %                      XYZSOM( 1, NOSOEF(2,NBEF) ),
     %                      XYZSOM( 1, NOSOEF(3,NBEF) ),
     %                      VECNOR, IER )
               IF( IER .NE. 0 ) THEN
                  PRINT*
                  PRINT*,'majxyznse: le TRIANGLE',NEF,' de SOMMETS',
     %                   (NOSOEF(L,NEF),L=1,NBSTEF),
     %' REDUIT a UNE ARETE est SUPPRIME ce QUI PEUT OUVRIR le MAILLAGE!'
                  PRINT*
                  NBEF = NBEF - 1
                  GOTO 50
               ENDIF
            ENDIF

C           MISE A ZERO AU DELA DE NBSTEF
            DO K = NBSTEF+1, L1SOEF
               NOSOEF( K, NBEF ) = 0
            ENDDO

C           RECHERCHE DE SOMMET NUL ou 2 SOMMETS DE MEME NUMERO
            DO M = 1, NBSTEF

C              LE SOMMET M DE L'EF NBEF
               NS = NOSOEF( M, NBEF )
               IF( NS .LE. 0 ) THEN
C                 UN SOMMET DE NUMERO NUL -> NEF SUPPRIME
                  GOTO 45
               ENDIF

               DO L = M+1, NBSTEF
                  NS2 = NOSOEF( L, NBEF )
                  IF( NS .EQ. NS2 ) THEN
C                    NEF AVEC 2 SOMMETS IDENTIQUES -> NEF SUPPRIME
                     GOTO 45
                  ENDIF
               ENDDO

            ENDDO

            GOTO 50

C           SUPPRESSION DE L'EF NEF DEVENU NBEF MOMENTANEMENT
 45         PRINT*
            PRINT*,'majxyznse: l''EF',NEF,' de SOMMETS',
     %             (NOSOEF(L,NEF),L=1,L1SOEF),' est DEVENU'
            PRINT*,'majxyznse: l''EF',NBEF,' de SOMMETS',
     %             (NOSOEF(L,NBEF),L=1,L1SOEF),' est SUPPRIME'
            NBEF = NBEF - 1
            GOTO 50

         ENDIF

 50   ENDDO

C     COMPRESSION DES XYZ DANS XYZSOM
C     -------------------------------
C     SUPPRESSION DES SOMMETS IDENTIFIES
      DO N = 1, NBCOSTID
         NEWXYZ( ABS( NOCOSTID(2,N) ) ) = 0
      ENDDO

      NBS = 0
      DO N = 1, NBSOM
         NS = NEWXYZ( N )    
         IF( NS .GT. 0 ) THEN
            NBS = NBS + 1
            DO K= 1, 3
               XYZSOM( K, NS ) = XYZSOM( K, N )
            ENDDO
         ENDIF
      ENDDO

C     NOMBRE FINAL DE SOMMETS ACTIFS
      NBSOM = NBS
      PRINT *,'majxyznse:',NBSOM0,
     %        ' SOMMETS INITIAUX avant IDENTIFICATION'
      PRINT *,'majxyznse:',NBSOM,
     %        ' SOMMETS FINAUX   apres IDENTIFICATION'


C     NOMBRE FINAL DES EF ACTIFS VERIFIES NON DEGENERES
      NBEFOB = NBEF

C     QUALITE MINIMALE DES NBEFOB EF FINAUX
C     -------------------------------------
      NBPO = 0
      NBAR = 0
      NBTR = 0
      NBQU = 0
      NBTE = 0
      NBPE = 0
      NBHE = 0

      NBEFQM = 0
      QUAMIN = 2.
      DO 100 NEF = 1, NBEFOB

C        NBSTEF: NOMBRE DE SOMMETS DE L'EF NEF
         DO NBSTEF = L1SOEF, 1, -1
            IF( NOSOEF( NBSTEF, NEF ) .NE. 0 ) GOTO 60
         ENDDO

C        LA QUALITE DE L'EF
 60      GOTO( 61, 62, 64, 64, 68, 68, 68, 68 ), NBSTEF

 61      NCOGEL = 1
         NBPO   = NBPO + 1
         GOTO 70

 62      NCOGEL = 2
         NBAR   = NBAR + 1
         GOTO 70

 64      IF( NOSOEF(4,NEF) .EQ. 0 ) THEN
            NCOGEL = 3
            NBTR   = NBTR + 1
         ELSE
            NCOGEL = 4
            NBQU   = NBQU + 1
         ENDIF
         GOTO 70

 68      IF( NOSOEF(5,NEF) .EQ. 0 ) THEN
            NCOGEL = 5
            NBTE   = NBTE + 1
         ELSE IF( NOSOEF(6,NEF) .EQ. 0 ) THEN
            NCOGEL = 6
            NBPE   = NBPE + 1
         ELSE
            NCOGEL = 8
            NBHE   = NBHE + 1
         ENDIF

 70      CALL QUALEF( NCOGEL,   NOSOEF(1,NEF), NBSOM, XYZSOM,
     %                SURFVOLU, QUALIT, IERR )

         IF( IERR .NE. 0 ) THEN
            PRINT*,'majxyznse: PROBLEME: EF ',NEF,' de SOMMETS:',
     %             (NOSOEF(I,NEF),I=1,NCOGEL),' de QUALITE INCALCULABLE'
C           L'EF NEF EST MIS A ZERO
            DO I=1,L1SOEF
               NOSOEF(I,NEF) = 0
            ENDDO
            GOTO 100
         ENDIF

         IF( QUALIT .LT. QUAMIN ) THEN
            QUAMIN=QUALIT
         ENDIF

         IF( QUALIT .LT. 0.03 ) THEN
            NBEFQM = NBEFQM + 1
            PRINT*,'majxyznse: Attention: EF',NEF,' de SOMMETS:',
     %             (NOSOEF(I,NEF),I=1,NCOGEL),
     %             ' de QUALITE=',QUALIT,'<0.03'
         ENDIF

 100  ENDDO

      PRINT *,'majxyznse:',NBEFOB0,' EF INITIAUX'
      PRINT *,'majxyznse:',NBEFOB ,' EF FINAUX de QUALITE MINIMALE',
     %         QUAMIN,' et',NBEFQM,' EF de QUALITE <0.03'

      IF( NBPO .GT. 0 ) PRINT*,'majxyznse:',NBPO,' POINTS'
      IF( NBAR .GT. 0 ) PRINT*,'majxyznse:',NBAR,' ARETES'
      IF( NBTR .GT. 0 ) PRINT*,'majxyznse:',NBTR,' TRIANGLES'
      IF( NBQU .GT. 0 ) PRINT*,'majxyznse:',NBQU,' QUADRANGLES'
      IF( NBTE .GT. 0 ) PRINT*,'majxyznse:',NBTE,' TETRAEDRES'
      IF( NBPE .GT. 0 ) PRINT*,'majxyznse:',NBPE,' PENTAEDRES'
      IF( NBHE .GT. 0 ) PRINT*,'majxyznse:',NBHE,' HEXAEDRES'

      RETURN
      END
