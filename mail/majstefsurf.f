      SUBROUTINE MAJSTEFSURF( NBSOM, XYZSOM, NEWNST, NBEFQT, NOSOEF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TRANSFORMER TOUT QUADRANGLE AYANT 2 SOMMETS IDENTIQUES EN 1 TRIANGLE
C -----   SUPPRIMER LES TRIANGLES AYANT 2 SOMMETS IDENTIQUES
C         RENUMEROTER LES SOMMETS DE 1 A NBSOM NOUVEAU

C MODIFIES:
C ---------
C NBSOM  : NOMBRE DE SOMMETS DE XYZSOM AVANT ET APRES IDENTIFICATIONS
C NEWNST : EN ENTREE NOUVEAU NUMERO DE CHACUN DES NBSOM SOMMETS INITIAUX
C          EN SORTIE NOUVEAU NUMERO DE CHACUN DES NBSOM SOMMETS FINAUX
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE DE LA SURFACE
C NBEFQT : NOMBRE  D'EF TRIANGLES ou QUADRANGLES AVANT et APRES
C NOSOEF : NUMERO DES 4 SOMMETS DES NBEFQT QUAD-TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY NOVEMBRE 2015
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NOSOEF(4,NBEFQT), NEWNST(0:NBSOM)
    
      PRINT *,'majstefsurf: AVANT', NBSOM,' SOMMETS',
     %         NBEFQT,' TRIANGLES ou QUADRANGLES'

C     IDENTIFICATION DU NUMERO DES SOMMETS DANS LES QT
      NBNEWEF = 0
      DO 50 K = 1, NBEFQT

         DO L=1,4
C           LE NOUVEAU NUMERO DU SOMMET
            NOSOEF(L,K) = NEWNST( NOSOEF(L,K) )
         ENDDO

C        NOMBRE DE SOMMETS DE NUMERO NUL CREES
         IF( NOSOEF(4,K) .EQ. 0 ) THEN
C           TRIANGLE
            NB0 = 1
            NBS = 3
         ELSE
C           QUADRANGLE
            NB0 = 0
            NBS = 4
         ENDIF

C        RECHERCHE DE LA MULTIPLICITE DES SOMMETS DU QT
         DO L=1,NBS-1

C           SOMMET L DU QT K
            NS0 = NOSOEF(L,K)

C           MULTIPLICITE DU SOMMET L
 10         NB = 1
            DO M = L+1, NBS
               IF( NOSOEF(M,K) .EQ. NS0 ) THEN
                  NB = NB + 1
                  M2 = M
               ENDIF
            ENDDO

            IF( NB .GE. 2 ) THEN
C              NOSOEF(L,K) = NOSOEF(M2,K)
               NOSOEF(M2,K) = 0
               NB0 = NB0 + 1
               GOTO 10
            ENDIF

         ENDDO

         IF( NB0 .EQ. 1 ) THEN

C           LE ZERO EST IL EN POSITION 4?
            IF( NOSOEF(4,K) .NE. 0 ) THEN
C              NON: MISE EN 4-EME POSITION
               DO M=1,3
                  IF( NOSOEF(M,K) .EQ. 0 ) THEN
                     DO N=M+1,4
                        NOSOEF(N-1,K)=NOSOEF(N,K)
                     ENDDO
                     NOSOEF(4,K) = 0
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF( NB0 .GE. 2 ) THEN

C           EF A SUPPRIMER
            GOTO 50

         ENDIF

C        MISE EN POSITION DU QT
         NBNEWEF = NBNEWEF + 1
         DO L=1,4
            NOSOEF(L,NBNEWEF) = NOSOEF(L,K)
         ENDDO

 50   ENDDO

C     NOMBRE FINAL DE QT
      NBEFQT = NBNEWEF

C     COMPRESSION XYZSOM DES SOMMETS ELIMINES
C     INITIALISATION DU NOUVEAU NUMERO DES SOMMETS
      DO N=0,NBSOM
         NEWNST( N ) = 0
      ENDDO

      DO K = 1, NBEFQT
         DO L = 1, 4
            N = NOSOEF(L,K)
            IF( N .GT. 0 ) THEN
               NEWNST( N ) = 1
            ENDIF
         ENDDO
      ENDDO

C     RENUMEROTATION DES SOMMETS RESTANTS
      NBS = 0
      DO N=1,NBSOM
         IF( NEWNST( N ) .GT. 0 ) THEN
            NBS = NBS + 1
            NEWNST( N ) = NBS
         ENDIF
      ENDDO

      DO K = 1, NBEFQT
         DO L=1,4
C           LE NOUVEAU NUMERO DU SOMMET
            NOSOEF(L,K) = NEWNST( NOSOEF(L,K) )
         ENDDO
      ENDDO

C     COMPRESSION DE XYZSOM
      DO N = 1, NBSOM
         K = NEWNST( N )
         IF( K .GT. 0 ) THEN
            DO L=1,3
               XYZSOM(L,K) = XYZSOM(L,N)
            ENDDO
         ENDIF
      ENDDO

C     NOMBRE FINAL DE SOMMETS
      NBSOM = NBS

      PRINT *,'majstefsurf: APRES', NBSOM,' SOMMETS',
     %         NBEFQT,' TRIANGLES ou QUADRANGLES'
      RETURN
      END
