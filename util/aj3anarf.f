      SUBROUTINE AJ3ANARF( LIBREA, L1ARFA, L2ARFA, NARFA,
     %                     NOQT,   NOSTQT, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER PAR HACHAGE LES ARETES DU QTANGLE NOQT AU TABLEAU
C -----    DES ARETES DES QTANGLES
C          UNE ARETE APPARTIENT A AU PLUS MXQTAR=L1ARFA-3 QTANGLES
C ENTREES:
C --------
C LIBREA : 1-ERE ARETE LIBRE A PARTIR DES DERNIERES PAR VALEUR DECROISSANTE
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE DE QTANGLES DU TABLEAU NARFA
C NOQT   : NUMERO NOSTQT DU QTANGLE A AJOUTER DANS NARFA
C NOSTQT : LES 4 NUMEROS DES SOMMETS DES QTANGLES DE LA SURFACE

C MODIFIES:
C ---------
C NARFA  : TABLEAU NARFA DES ARETES DES QTANGLES DU MAILLAGE
C          NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR LE QTANGLE SUIVANT
C          NARFA(4:3+MXQTAR,I)= +-NO NOSTQT DU QTANGLE CONTENANT L'ARETE
C                               SI ARETE DANS LE SENS DIRECT ou INDIRECT
C IERR   : -1 000 000 SI L2ARFA INSUFFISANT DU TABLEAU ARFA DES ARETES
C           1 SI MXQTAR TROP PETIT
C           0 SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Veulettes sur mer                Octobre 2019
C....................................................................012
      include"./incl/langue.inc"
      INTEGER  NOSTQT(4,*), NARFA(L1ARFA,L2ARFA), NS(2)

C     BOUCLE SUR LES ARETES DU QTANGLE NOQT
C     -------------------------------------
      IERR   = 0

C     MXQTAR=NOMBRE MAXIMUM DE QTANGLES POUR UNE ARETE (STOCKES DANS NARFA)
      MXQTAR = L1ARFA - 3

C     CALCUL DU NOMBRE D'ARETES DU QTANGLE NOQT
      DO NBAR = 4, 1, -1
         IF( NOSTQT(NBAR,NOQT) .GT. 0 ) GOTO 10
      ENDDO
      GOTO 9999

 10   J1  = NBAR
      NS1 = NOSTQT( J1, NOQT )
      DO 80 J2 = 1, NBAR

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

C        ADJONCTION DE L'ARETE SI ELLE N'EXISTE PAS DEJA
         CALL HACHAG( 2, NS, L1ARFA, L2ARFA, NARFA, 3,  LIBREA, NAR )

C        AJOUT DU NUMERO DU QTANGLE DE CETTE ARETE
C        NAR > 0 SI ARETE RETROUVEE DANS NARFA
C        NAR < 0 SI ARETE AJOUTEE   DANS NARFA
C        NAR = 0 SI LE TABLEAU NARFA EST SATURE
         IF( NAR .EQ. 0 ) THEN

            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'aj3anarf: SATURATION TABLEAU ARETES des QTANGLES'
            ELSE
               PRINT*,'aj3anarf: EDGES ARRAY of SURFACE IS SATURATED'
            ENDIF
            IERR = -1 000 000
            RETURN

         ENDIF

C        LE K-EME QTANGLE DE L'ARETE NAR EST STOCKE
         IF( NAR .LT. 0 ) NAR = -NAR
         DO K = 4, L1ARFA

            IF( ABS( NARFA( K, NAR ) ) .EQ. ABS( NOQT ) ) THEN
C              NOQT EST DEJA STOCKE
               GOTO 80
            ENDIF

            IF( NARFA( K, NAR ) .EQ. 0 ) THEN
C              STOCKAGE DU QT NOQT EN POSITION K DE L'ARETE NAR
               NARFA( K, NAR ) = NOQT
               GOTO 80
            ENDIF

         ENDDO

C        ICI  PROBLEME: >MXQTAR QTANGLES POUR CETTE ARETE NAR
C        ----------------------------------------------------
C        VALEUR IMPOSEE NEGATIVE DU NO DE LA FACE EN DERNIERE POSITION
         NARFA( L1ARFA, NAR ) = -NOQT

C        AFFICHAGE DE L'ARETE APPARTENANT A >MXQTAR QTANGLES
         PRINT*
         IF( LANGAG .EQ. 0 ) THEN

           PRINT*,'aj3anarf: ARETE COMMUNE A PLUS DE',MXQTAR,' QTANGLES'
            PRINT*,'aj3anarf: face',NOQT,' arete',J2,' Sommets:',NS,
     %           ' des QTANGLES',(NARFA(K,NAR),K=4,L1ARFA)
            DO K=1,MXQTAR
               NOQT2 = ABS( NARFA(3+K,NAR) )
               PRINT*,'QTANGLE',NOQT2,
     %              ' Sommets:',(NOSTQT(L,NOQT2),L=1,4)
            ENDDO

         ELSE

            PRINT*,'aj3anarf: COMMON EDGE IN MORE of',MXQTAR,' QTANGLES'
            PRINT*,'aj3anarf: face',NOQT,' edge',J2,' Vertices:',NS,
     %           ' of QTANGLES',(NARFA(K,NAR),K=4,L1ARFA)
            DO K=4,L1ARFA
               NOQT2 = ABS( NARFA(K,NAR) )
               PRINT*,'QTANGLE',NOQT2,
     %              ' Vertices:',(NOSTQT(L,NOQT2),L=1,4)
            ENDDO

         ENDIF

         IERR = 1

 80   ENDDO

 9999 RETURN
      END
