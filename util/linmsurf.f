      SUBROUTINE LINMSURF( MNDFOB,  NBNOSU, MNNOSU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LIRE LES NOMS D'UNE LISTE DE SURFACES D'UN OBJET
C -----
C
C ENTREES:
C --------
C MNDFOB : ADRESSE MCN DU TABLEAU DE DEFINITION DE L'OBJET
C NBNOSU : <0 DEMANDE TOUTES LES SURFACES DE L'OBJET SANS DEMANDER
C             LE NOMBRE ET LE NOM DES SURFACES
C
C SORTIES:
C --------
C NBNOSU : >0 NOMBRE DE SURFACES DEMANDEES, LUES ET RETROUVEES
C          =0 PAS DE LECTURE DE SURFACES OU ABANDON
C MNNOSU : NUMERO DES SURFACES DANS LE LEXIQUE DES SURFACES
C          DE NOM LU ET RETROUVE PARMI LES SURFACES DE L'OBJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY SEPTEMBRE 2010
C23456---------------------------------------------------------------012
      include"./incl/a_objet__definition.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
C
C     CALCUL DU NOMBRE DE SURFACES DE L'OBJET
C     ---------------------------------------
C     NBDOBJ = NOMBRE DE PLSV de L'OBJET
      NBDOBJ = MCN( MNDFOB + WBDOBJ )
      MNTYOB = MNDFOB + WTYOBJ - 2
      NBSUR = 0
      DO J = 1, NBDOBJ
         MN = MNTYOB + 2 * J
         IF( MCN(MN) .EQ. 3 ) THEN
C           LE TYPE EST CELUI D'UNE SURFACE
            NBSUR = NBSUR + 1
         ENDIF
      ENDDO
C
C     LECTURE DU NOMBRE DE SURFACES DE TRACE DE LA SOLUTION
C     -----------------------------------------------------
      MNNOSU = 0
      IF( NBNOSU .LT. 0 ) THEN
C        TOUTES LES SURFACES DE L'OBJET SANS DEMANDER LEURS NOMS
         NBNOSU = NBSUR
      ELSE
C        LES SURFACES DEMANDEES
         NBNOSU = 0
         CALL INVITE( 37 )
         CALL LIRENT( NCVALS, NBNOSU )
         IF( NCVALS .LE. 0 ) RETURN
         IF( NBNOSU .LE. 0 ) RETURN
         IF( NBNOSU .GT. NBSUR ) THEN
C           => TOUTES LES SURFACES DE L'OBJET SANS DEMANDER LEURS NOMS
            NBNOSU = NBSUR
         ENDIF
      ENDIF
      CALL TNMCDC( 'ENTIER', NBNOSU, MNNOSU )
C
      IF( NBNOSU .LT. NBSUR ) THEN
C
C        LECTURE DU NOM DES NBNOSU SURFACES
C        ----------------------------------
         DO K = 1, NBNOSU
C
C           INVITE POUR ENTRER LE NOM DE LA SURFACE
 10         CALL INVITE( 43 )
            CALL OBJENU( 'SURFACE' , NUSUTR )
            IF( NUSUTR .LE. 0 ) THEN
             IF( MNNOSU .GT. 0 ) CALL TNMCDS( 'ENTIER', NBNOSU, MNNOSU )
             NBNOSU = 0
             RETURN
            ENDIF
C
C           VERIFICATION: CETTE SURFACE EST ELLE UNE SURFACE DE L'OBJET?
            N = 0
            DO J = 1, NBDOBJ
               MN = MNTYOB + 2 * J
               IF( MCN(MN) .EQ. 3 ) THEN
C                 LE TYPE EST CELUI D'UNE SURFACE
                  N = N + 1
                  IF( MCN(MN+1) .EQ. NUSUTR ) GOTO 20
               ENDIF
            ENDDO
C
C           SURFACE NON RETROUVEE DANS L'OBJET
            NBLGRC(NRERR) = 2 + N
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOM DE SURFACE NON RECONNUE PARMI'
               KERR(2) = 'LES NOMS DE SURFACE DE L''OBJET'
            ELSE
               KERR(1) = 'SURFACE NAME NOT KNOWN AMONG'
               KERR(2) = 'THE SURFACE NAMES of THE OBJECT'
            ENDIF
            N = 2
            DO J = 1, NBDOBJ
               MN = MNTYOB + 2 * J
               IF( MCN(MN) .EQ. 3 ) THEN
C                 LE TYPE EST CELUI D'UNE SURFACE
                  N = N + 1
C                 LE NOM DE LA SURFACE
                  CALL NMOBNU( 'SURFACE', MCN(MN+1), KERR(N) )
               ENDIF
            ENDDO
            CALL LEREUR
            GOTO 10
C
C           SURFACE DEJA DONNEE?
 20         DO J = 1, K-1
               IF( MCN( MNNOSU-1+J ) .EQ.  NUSUTR ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'NOM DE SURFACE DEJA DONNE'
                  ELSE
                     KERR(1) = 'SURFACE NAME ALREADY GIVEN'
                  ENDIF
                  CALL LEREUR
                  GOTO 10
               ENDIF
            ENDDO
C
C           SURFACE DE L'OBJET RETROUVEE
            MCN( MNNOSU - 1 + K ) = NUSUTR
C
         ENDDO
C
      ELSE
C
C        TOUTES LES SURFACES DE L'OBJET SONT LISTEES
C        -------------------------------------------
         NBSUR = 0
         DO J = 1, NBDOBJ
            MN = MNTYOB + 2 * J
            IF( MCN(MN) .EQ. 3 ) THEN
C              LE TYPE EST CELUI D'UNE SURFACE
               NBSUR = NBSUR + 1
C              LA SURFACE DE L'OBJET EST AJOUTEE
               MCN( MNNOSU - 1 + NBSUR ) = MCN(MN+1)
            ENDIF
         ENDDO
C
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TOUTES LES SURFACES DE L''OBJET SONT LISTEES'
         ELSE
            KERR(1) = 'ALL OBJECT SURFACES ARE USED'
         ENDIF
         CALL LERESU
C
      ENDIF
C
      RETURN
      END
