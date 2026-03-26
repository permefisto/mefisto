      SUBROUTINE NBARXFA( L1ARFA, L2ARFA, NARFA,
ccc                          NBSOEF, NOSOEF,
     %                    NBARXF, NOARPB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    COMPTER LE NOMBRE D'ARETES DE NARFA APPARTENANT
C -----    A 1, 2, >=L1ARFA-3 FACES (QUADRANGLE ou TRIANGLE)

C ENTREES:
C --------
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE DE FACES DU TABLEAU NARFA
C NARFA  : NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR L'ARETE SUIVANTE
C          NARFA(4:L1ARFA,I)= NO DANS NOSOFA DE LA FACE CONTENANT L'ARETE I
C          SI UNE ARETE APPARTIENT A PLUS DE L1ARFA-3 FACES, 
C          LE NUMERO DE FACE L1ARFA EST RENDU NEGATIF POUR INDIQUER
C          QUE LA LISTE DES FACES EST INCOMPLETE

cccC NBSOEF : NOMBRE DE SOMMETS ET TANGENTES PAR FACE
cccC          (4 SANS TG ET 12 AVEC TG )
cccC NOSOEF : LES 4 NUMEROS DES SOMMETS DES FACES DE LA SURFACE

C SORTIES:
C --------
C NBARXF : NBARXF(n) NOMBRE D'ARETES APPARTENANT A n  FACES(ou QT) n=1,2
C          NBARXF(6) NOMBRE D'ARETES APPARTENANT A >5 FACES
C NOARPB : >0 NO NARFA DE LA DERNIERE ARETE APPARTENANT A AU MOINS 6 QT
C IERR   : =1 SI L1ARFA-3 < 2
C          =0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2015
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray             Avril 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
ccc      INTEGER        NOSOEF(NBSOEF,*)
      INTEGER        NARFA(L1ARFA,L2ARFA),
     %               NBARXF( 6 )

      IERR = 0
      DO NA = 1, 6
         NBARXF( NA ) = 0
      ENDDO

      NOARPB = 0
      MXFAAR = L1ARFA - 3
      IF( MXFAAR .LT. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'nbarxfa: MXFAAR=',MXFAAR,' TROP PETIT'
         ELSE
            PRINT*,'nbarxfa: MXFAAR=',MXFAAR,' TOO SHORT'
         ENDIF
         IERR = 1
         RETURN
      ENDIF

C     RECHERCHE DES ARETES APPARTENANT A UNE SEULE FACE=QT
      NBARET = 0
      DO NA = 1, L2ARFA

         IF( NARFA(1,NA) .NE. 0 ) THEN

C           L'ARETE EST INITIALISEE
            NBARET = NBARET + 1

C           QUEL EST LE NOMBRE DE QT CONTENANT CETTE ARETE?
            NBF = 0
            DO K = 4, L1ARFA
               IF( NARFA(K,NA) .NE. 0 ) NBF = NBF + 1
            ENDDO

C           L'ARETE APPARTIENT A NBF QT
            IF( NBF .GE. 3 ) THEN
ccc               PRINT*,'nbarxfa: ATTENTION l''ARETE',NA,
ccc     %                ' St1',NARFA(1,NA),' St2',NARFA(2,NA),
ccc     %                ' APPARTIENT aux Faces',(NARFA(3+K,NA),K=1,NBF)
               NOARPB = NA
            ENDIF

ccc            IF( NBF .NE. 2 ) THEN
cccC              PLUS DE 2 QT CONTENANT L'ARETE NA
ccc               DO K = 4, L1ARFA
ccc                  NOQT = ABS( NARFA(K,NA) )
ccc                  IF( NOQT .GT. 0 ) THEN
cccC                    LE NOMBRE DE SOMMETS du QT
ccc                     DO NBS=NBSOEF,1,-1
ccc                        IF( NOSOEF(NBS,NOQT) .GT. 0 ) GOTO 10
ccc                     ENDDO
ccc 10                 PRINT*,'nbarxfa: ARETE(',NA,')=',(NARFA(L,NA),L=1,2)
ccc     %                 ,' de l''EF',NOQT,' St:',(NOSOEF(L,NOQT),L=1,NBS)
ccc                  ENDIF
ccc               ENDDO
ccc            ENDIF

C           NOMBRE de QT CONTENANT L'ARETE NA de NARFA
            IF( NBF .GT. 6 ) NBF = 6
            NBARXF( NBF ) = NBARXF( NBF ) + 1

         ENDIF

      ENDDO

C     BILAN du NOMBRE de Faces PAR ARETE
C     ----------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10000, NBARET
      ELSE
         PRINT 20000, NBARET
      ENDIF
10000 FORMAT(' nbarxfa:',I9,' ARETES du MAILLAGE' )
20000 FORMAT(' nbarxfa:',I9,' MESH EDGES' )

      DO K=1,5
         IF( NBARXF(K) .GT. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT 10001, NBARXF(K), '  ', K
            ELSE
               PRINT 20001, NBARXF(K), '  ', K
            ENDIF
         ENDIF
      ENDDO
      IF( NBARXF(6) .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10001, NBARXF(6), ' >', 5
         ELSE
            PRINT 20001, NBARXF(6), ' >', 5
         ENDIF
      ENDIF

10001 FORMAT(' nbarxfa:',I9,' ARETES APPARTENANT a',A2,I1,' Faces')
20001 FORMAT(' nbarxfa:',I9,' EDGES BELONGING to',A2,I1,' Faces')

      RETURN
      END
