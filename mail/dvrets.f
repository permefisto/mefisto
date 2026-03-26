      SUBROUTINE DVRETS( NS,  NOTRIA, NOTRSO,
     %                   NBS, MXTRNS, NOTRNS, NOAROP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EMPILER AU PLUS MXTRNS DES TRIANGLES DE SOMMET NS SUPPOSE
C ------   INTERNE
C
C          ATTENTION : SI NS SUR LE BORD OU ANOMALIE => NBS=0
C
C ENTREES:
C --------
C NS     : NUMERO DU SOMMET A TRAITER
C NOTRIA : LISTE CHAINEE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                           ADJACENT PAR L'ARETE I
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C MXTRNS : NOMBRE MAXIMAL DE TRIANGLES EMPILABLES DANS NOTRNS (ET NOAROP)
C
C SORTIES:
C --------
C NBS    : NOMBRE DE TRIANGLES DE SOMMET NS
C          0 SI PROBLEME ( SOMMET FRONTALIER, OU PILE NOTRNS SATUREE... )
C NOTRNS : NUMERO DES NBS TRIANGLES DE SOMMET NS
C NOAROP : NUMERO 1 A 3 DE L'ARETE DANS LE TRIANGLE OPPOSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1994
C....................................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NOTRIA(1:6,*),
     %                  NOTRSO(*)
      INTEGER           NOTRNS(MXTRNS),NOAROP(MXTRNS)
C
C     NT UN TRIANGLE CONTENANT LE SOMMET NS
      NT  = NOTRSO(NS)
      NBS = 0
      IF( NT .GT. 0 )  THEN
         IF( NOTRIA(1,NT) .LE. 0 ) THEN
C           ANOMALIE A CORRIGER . OUBLI DE METTRE NOTRSO(NS) A ZERO
            WRITE(IMPRIM,*) 'DVRETS: OUBLI DE METTRE NOTRSO(NS)=0 AVANT'
            NBS = 0
            RETURN
         ENDIF
C
C        LE POINT EST INTERNE ET PEUT ETRE UN POINT IMPOSE PAR L'UTILISATEUR
C        PARCOURS DES TRIANGLES DE SOMMET NS DANS UN SEUL SENS
C        -------------------------------------------------------------------
C        LE NUMERO DES SOMMETS DE NT
         DO 10 I=1,3
            IF( NOTRIA(I,NT) .EQ. NS ) GOTO 20
   10    CONTINUE
         WRITE(IMPRIM,*)'DVRETS: SOMMET ',NS,' NON DANS TRIANGLE0 ',NT
         RETURN
C
   20    J1 = I + 1
         IF( J1 .EQ. 4 ) J1 = 1
         NS1 = NOTRIA(J1 ,NT)
C
C        LE TRIANGLE DE L'AUTRE COTE DE L'ARETE NS-NS1
         NT1 = NOTRIA(I+3,NT)
C
C        LA CONTRIBUTION DE NS1 A LA PILE DES TRIANGLES DE SOMMET NS
         NBS  = 1
         NOTRNS( 1 ) = NT
         NOAROP( 1 ) = J1
C
C        PARCOURS DES TRIANGLES DE SOMMET NS PAR L'ARETE NS NS1
C        ------------------------------------------------------
   25    IF( NT1 .NE. NT ) THEN
C
C           RECHERCHE DE L'ARETE NS-NS1
            DO 30 J=1,3
               IF( NOTRIA(J,NT1) .EQ. NS1 ) GOTO 40
   30       CONTINUE
C
C           L'ARETE SUIVANTE DE SOMMET NS
   40       IF( J .LT. 3 ) THEN
               J1 = J + 1
            ELSE
               J1 = 1
            ENDIF
            IF( J1 .LT. 3 ) THEN
               J2 = J1 + 1
            ELSE
               J2 = 1
            ENDIF
            NS2 = NOTRIA(J2,NT1)
C
            IF( NBS .GE. MXTRNS ) THEN
C              DEBORDEMENT DE LA PILE => PROBLEME DANS LE CHAINAGE
               WRITE(IMPRIM,*) 'DVRETS: SOMMET ',NS,' TROP DE TRIANGLES'
               NBS = 0
               RETURN
            ENDIF
C
C           LE TRIANGLE NT1 PEUT ETRE EMPILE
            NBS  = NBS + 1
            NOTRNS( NBS ) = NT1
            NOAROP( NBS ) = J2
C
C           PASSAGE AU TRIANGLE SUIVANT
            NT1 = NOTRIA(J1+3,NT1)
            IF( NT1 .LE. 0 ) THEN
               WRITE(IMPRIM,*) 'DVRETS:',NS,' SOMMET NON INTERNE'
               NBS = 0
               RETURN
            ENDIF
C
C           PASSAGE AU TRIANGLE DE SOMMET NS SUIVANT
            NS1 = NS2
            GOTO 25
         ENDIF
      ENDIF
      END
