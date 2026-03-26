      SUBROUTINE SU3ANARF( L1ARFA, L2ARFA, NARFA, NOQT, NOSTQT, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER DU HACHAGE LES ARETES DU QTANGLE NOQT AU TABLEAU
C -----    DES ARETES DES QTANGLES
C          UNE ARETE APPARTIENT A AU PLUS MXQTAR=L1ARFA-3 QTANGLES
C ENTREES:
C --------
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE DE QTANGLES DU TABLEAU NARFA
C NOQT   : NUMERO NOSTQT DU QTANGLE A SUPPRIMER DANS NARFA
C NOSTQT : LES 4 NUMEROS DES SOMMETS DES QTANGLES DE LA SURFACE

C MODIFIES:
C ---------
C NARFA  : TABLEAU NARFA DES ARETES DES QTANGLES DU MAILLAGE
C          NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR LE QTANGLE SUIVANT
C          NARFA(4:3+MXQTAR,I)= +-NO NOSTQT DU QTANGLE CONTENANT L'ARETE
C                               SI ARETE DANS LE SENS DIRECT ou INDIRECT
C IERR   : =2 SI UNE ARETE DE NOQT N'EST PAS STOCKEE DANS LE TABLEAU NARFA
C          =0 SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Veulettes sur mer               Octobre 2019
C....................................................................012
      include"./incl/langue.inc"
      INTEGER  NOSTQT(4,*), NARFA(L1ARFA,L2ARFA), NS(2)

C     BOUCLE SUR LES ARETES DU QTANGLE NOQT
C     -------------------------------------
      IERR = 0

C     CALCUL DU NOMBRE D'ARETES DU QTANGLE NOQT
      DO NBAR = 4, 1, -1
         IF( NOSTQT( NBAR, NOQT ) .GT. 0 ) GOTO 10
      ENDDO
      GOTO 9999

 10   J1  = NBAR
      NS1 = NOSTQT( J1, NOQT )
      DO 50 J2 = 1, NBAR

C        L'ARETE J2 DU QTANGLE NOQT
         NS2 = NOSTQT( J2, NOQT )
         IF( NS1 .LE. NS2 ) THEN
            NS(1) = NS1
            NS(2) = NS2
         ELSE
            NS(1) = NS2
            NS(2) = NS1
         ENDIF
         J1  = J2
         NS1 = NS2

C        RECHERCHE DE L'ARETE DE SOMMET NS DANS NARFA
         CALL HACHAR( 2, NS, L1ARFA, L2ARFA, NARFA, 3, NAR )

         IF( NAR .LE. 0 ) THEN
C           ARETE NS NON RETROUVEE DANS NARFA
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'su3anarf: ARETE NON RETROUVEE de SOMMETS',NS,
     %                ' du QTANGLE',NOQT,' St:',(NOSTQT(kk,NOQT),kk=1,4)
     %               ,' NAR=',NAR
            ELSE
               PRINT*,'su3anarf: NOT RECOVERED EDGE of VERTICES',NS,
     %          ' of QTANGLE',NOQT,' Vertices:',(NOSTQT(kk,NOQT),kk=1,4)
     %         ,' NAR=',NAR
            ENDIF
            IERR = 2
            GOTO 9999
         ENDIF

C        RECHERCHE DU QTANGLE NOQT SUR L'ARETE NAR DE NARFA
C        --------------------------------------------------
         DO K = 4, L1ARFA

            IF( ABS( NARFA( K, NAR ) ) .EQ. NOQT ) THEN

C              MISE A ZERO DU QTANGLE NOQT DE L'ARETE NAR
C              ------------------------------------------
               NARFA( K, NAR ) = 0

C              DECALAGE POUR SUPPRIMER LE ZERO EN POSITION K DE NAR
C              COMPRESSION DES NO DE QTANGLES DE L'ARETE NAR
               NBZERO = 0
               DO KK = 4, L1ARFA
                  NQT = NARFA( KK, NAR )
                  IF( NQT .EQ. 0 ) THEN
                     NBZERO = NBZERO + 1
                  ELSE
                     NARFA( KK-NBZERO, NAR ) = NQT
                  ENDIF
               ENDDO
               DO KK = L1ARFA+1-NBZERO, L1ARFA
                  NARFA( KK, NAR ) = 0
               ENDDO

               IF( NBZERO .LT. L1ARFA-3 ) THEN
C                 L'ARETE NAR EST CONSERVEE CAR ELLE APPARTIENT
C                 A AU MOINS UN QTANGLE
                  GOTO 50
               ENDIF

C              NAR APPARTIENT A AUCUN QTANGLE =>
C              ELLE EST SUPPRIMEE DE NARFA
C              ---------------------------------
C              RECHERCHE DU CHAINAGE PRECEDANT DE L'ARETE NAR
               NAR0 = 0
               NAR1 = NS(1) + NS(2)
               NAR1 = MOD( ABS(NAR1), L2ARFA )
               IF( NAR1 .EQ. 0 ) NAR1 = L2ARFA

 15            IF( NAR1 .NE. NAR ) THEN
C                 PASSAGE A L'ARETE SUIVANTE
                  NAR0 = NAR1
                  NAR1 = NARFA( 3, NAR1 )
                  GOTO 15
               ENDIF

C              NAR0 PRECEDE NAR DANS LE CHAINAGE DE NARFA ET SUIVIE PAR NAR1
               NAR1 = NARFA( 3, NAR )
               IF( NAR0 .GT. 0 ) THEN

C                 NAR EST RETIRE DU CHAINAGE
                  NARFA( 3, NAR0 ) = NAR1

C                 L'ARETE NAR EST REINITIALISEE A ZERO
                  DO KK = 1, L1ARFA
                     NARFA( KK, NAR ) = 0
                  ENDDO

               ELSE

                  IF( NAR1 .GT. 0 ) THEN

C                    L'ARETE NAR1 EST TRANSFEREE EN L'ARETE NAR
                     DO KK = 1, L1ARFA
                        NARFA( KK, NAR  ) = NARFA( KK, NAR1 )
                        NARFA( KK, NAR1 ) = 0
                     ENDDO

                  ELSE

C                    PAS D'ARETE SUIVANTE DE NAR
C                    L'ARETE NAR EST REINITIALISEE A ZERO
                     DO KK = 1, L1ARFA
                        NARFA( KK,NAR ) = 0
                     ENDDO

                  ENDIF

               ENDIF

               GOTO 50

            ENDIF
 
         ENDDO

 50   ENDDO

 9999 RETURN
      END
