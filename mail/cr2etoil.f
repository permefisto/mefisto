      SUBROUTINE CR2ETOIL( KTITRE, PTXYZD, NPSOFR, NFLPER,
     %                     MXFACO, LEFACO, NO0FAR, NBTRCF, NOTRCF,
     %                     NOTETR, N1TETS,
     %                     MXTECF, NBTECF, NOTECF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DE L'ETOILE DES NBTRCF FACES PERDUES
C -----    DE LA FRONTIERE PAR AJOUT DES TETRAEDRES DES SOMMETS
C          SUPPRIMABLES PROCHES DES FACES PERDUES

C ENTREES:
C --------
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0  SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE

C NFLPER : NUMERO LEFACO DE LA FACE PERDUE A TRAITER
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          45:  NO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          678: NO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          9: ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C             => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C             LEFACO(9,*) -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          10: HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C              LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C              NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C              SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C              NF = LEFACO( 9, NF )  ...
C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON
cccC          12: = NO FACEOC DE 1 A NBFACES D'OC
C NO0FAR : NUMERO DES 3 SOMMETS DES FACES AJOUTEES AU CF

C NBTRCF : NOMBRE DE FACES DE NOTRCF
C NOTRCF : >0 NUMERO DANS LEFACO DES TRIANGLES PERDUS  DU CF
C          <0 NUMERO DANS NO0FAR DES TRIANGLES AJOUTES AU CF

C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS

C SORTIES:
C --------
C NBTECF : NOMBRE DE TETRAEDRES DE L'ETOILE
C NOTECF : NUMERO NOTETR DES NBTECF TETRAEDRES DE L'ETOILE
C IERR   : =0  SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY             Mars 2018
C2345X7..............................................................012
      PARAMETER        (MXSTSU=256)
      CHARACTER*(*)     KTITRE

      DOUBLE PRECISION  PTXYZD(4,*), D, ARLG,
     %                  CBTR(3), PTPROJ(3)

      INTEGER           NOTETR(8,*), N1TETS(*), NFETOI(3),
     %                  LEFACO(11,0:MXFACO), NO0FAR(3,*), NPSOFR(*),
     %                  NOTRCF(NBTRCF), NOTECF(MXTECF), NOSTSU(MXSTSU)

      INTEGER           NOSOTR(3), NOSOARTE(2,6), NOSOFATE(3,4)
      DATA              NOSOARTE / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      IERR    = 0
      NBTECF0 = NBTECF

C     -----------------------------------------------------------------
C     EXISTE T IL UN SOMMET DE L'ETOILE, SUPPRIMABLE, TRES PROCHE DES
C     FACES PERDUES DE L'ETOILE?
C     (RISQUE D'UNE ETOILE APLATIE EN CE SOMMET => ABANDON FACE PERDUE)
C     AJOUT DES TETRAEDRES DE CES SOMMETS SUPPRIMABLES
C     -----------------------------------------------------------------
C     CONSTRUCTION DE LA LISTE DES SOMMETS SUPPRIMABLES DES TETRAEDRES
C     DE L'ETOILE POUR NE PAS FAIRE PLUSIEURS FOIS LE CALCUL DE PROXIMITE
      NBSTSU = 0
      NBTECF = 0
      DO N = 1, NBTECF0
         NTE = NOTECF( N )
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN
            NBTECF = NBTECF + 1
            NOTECF( NBTECF ) = NTE
            DO 10 K=1,4
               NS = NOTETR( K, NTE )
               IF( NPSOFR( NS ) .GT. 0 ) GOTO 10
C              NS EST SUPPRIMABLE. DEJA STOCKE?
               DO M=1,NBSTSU
                  IF( NS .EQ. NOSTSU(M) ) GOTO 10
               ENDDO
C              NON STOCKE: NS EST UN SOMMET SUPPRIMABLE DE PLUS
               IF( NBSTSU .GE. MXSTSU ) THEN
                  PRINT*,'cr2etoil: AUGMENTER MXSTSU=',MXSTSU
                  GOTO 20
               ENDIF
               NBSTSU = NBSTSU + 1
               NOSTSU( NBSTSU ) = NS
 10         ENDDO
         ENDIF
      ENDDO
      NBTECF0 = NBTECF

      print*,'cr2etoil: 5) Tableau nostsu des',NBSTSU,
     %       ' sommets supprimables :',(nostsu(n),n=1,nbstsu)

 20   DO 40 N = 1, NBSTSU

C        LE SOMMET SUPPRIMABLE
         NSU = NOSTSU( N )

C        CALCUL DE LA DISTANCE DE NSU AUX FACES PERDUES DU CF
         NBTEC0 = NBTECF
         DO 35 L=1,NBTRCF

            NTR = NOTRCF( L )
            IF( NTR .GT. 0 ) THEN
               NOSOTR(1) = LEFACO(1,NTR)
               NOSOTR(2) = LEFACO(2,NTR)
               NOSOTR(3) = LEFACO(3,NTR)
            ELSE
               NOSOTR(1) = NO0FAR(1,-NTR)
               NOSOTR(2) = NO0FAR(2,-NTR)
               NOSOTR(3) = NO0FAR(3,-NTR)
            ENDIF

C           POINT DE PROJECTION DE NSU SUR LE PLAN DE LA FACE PERDUE NTR
            CALL PRPTPLD( PTXYZD(1,NSU),
     %                    PTXYZD(1,NOSOTR(1)),
     %                    PTXYZD(1,NOSOTR(2)),
     %                    PTXYZD(1,NOSOTR(3)),
     %                    PTPROJ, IER )
C           PTPROJ(3) : XYZ DU POINT PROJETE SUR LE PLAN NOSOTR(1:3)
            IF( IER .NE. 0 ) GOTO 35

C           3 COORDONNEES BARYCENTRIQUES DU POINT PT DANS LE TRIANGLE NTR
            CALL CBPTTR( PTXYZD(1,NOSOTR(1)),
     %                   PTXYZD(1,NOSOTR(2)),
     %                   PTXYZD(1,NOSOTR(3)),
     %                   PTPROJ,  CBTR )

C           PTPROJ EST IL INTERNE OU PROCHE DU TRIANGLE?
            IF( ABS(CBTR(1))+ABS(CBTR(2))+ABS(CBTR(3)) .LE. 1.1D0 ) THEN

C              OUI: CALCUL DE SA DISTANCE NSU-PTPROJ
               D = SQRT( ( PTXYZD(1,NSU) - PTPROJ(1) ) **2
     %                 + ( PTXYZD(2,NSU) - PTPROJ(2) ) **2
     %                 + ( PTXYZD(3,NSU) - PTPROJ(3) ) **2 )

C              PERIMETRE DES 3 ARETES DU TRIANGLE NTR
               ARLG = 0D0
               NS1  = NOSOTR(3)
               DO M=1,3
                  NS2 = NOSOTR(M)
                  ARLG = ARLG + SQRT(
     %                 ( PTXYZD(1,NS2) - PTXYZD(1,NS1) ) ** 2
     %               + ( PTXYZD(2,NS2) - PTXYZD(2,NS1) ) ** 2
     %               + ( PTXYZD(3,NS2) - PTXYZD(3,NS1) ) ** 2 )
                  NS1 = NS2
               ENDDO

               IF( D .LT. ARLG * 0.1D0 ) THEN

C                 LE SOMMET NSU EST PROCHE DE LA FACE NTR
C                 LES TETRAEDRES DE CE SOMMET SONT AJOUTES A L'ETOILE
                  CALL TETR1S( NSU,   N1TETS, NOTETR,
     %                         NBTNS, MXTECF-NBTECF, NOTECF(NBTECF+1),
     %                         IERR )
                  NBTECF = NBTECF + NBTNS

C                 UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
                  CALL UNITABL( NOTECF, NBTECF )

                  IF( IERR .NE. 0 ) THEN
                 PRINT*,'cr2etoil: PAS de PROBLEME dans TETR1S NSU=',NSU
                     GOTO 38
                  ENDIF

                  IF( NBTECF .GT. NBTEC0 ) THEN
                     PRINT*,'cr2etoil: FACE PERDUE',NFLPER,
     %                      ' SOMMET PROCHE',NSU,' SUPPRIME et',
     %                       NBTECF-NBTEC0,' TETRAEDRES AJOUTES'

                     DO M=NBTEC0+1,NBTECF
                        NT = NOTECF( M )
                        PRINT*,'cr2etoil: Ajout du NOTETR(',NT,')=',
     %                         (NOTETR(J,NT),J=1,8)
                     ENDDO

                     NBTEC0 = NBTECF
                  ENDIF

C                 PASSAGE AU SOMMET SUIVANT
                  GOTO 40

               ENDIF

            ENDIF

 35      ENDDO

C        LE SOMMET SUPPRIMABLE NSU N'A PAS ETE SUPPRIME
 38      NOSTSU( N ) = 0

 40   ENDDO

C     UNICITE DES NO DES SOMMETS SUPPRIMES
      CALL UNITABL( NOSTSU, NBSTSU )
      print*,'cr2etoil: 5) Tableau nostsu des',NBSTSU,
     %       ' sommets SUPPRIMES :',(nostsu(n),n=1,nbstsu)

C     UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
      CALL UNITABL( NOTECF, NBTECF )

      IF( NBTECF0 .NE. NBTECF ) THEN
      KTITRE='cr2etoil: 5)       TETRAEDRES AJOUTES POUR SUPPRIMER
     %     SOMMETS'
         WRITE(KTITRE(14:18),'(I5)') NBTECF
         WRITE(KTITRE(54:58),'(I5)') NBSTSU
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRFETO4( KTITRE(1:L), PTXYZD, 0, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
      ENDIF

      RETURN
      END
