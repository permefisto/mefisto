      SUBROUTINE AJTRAUCF( KTITRE, NFP,    PTXYZD, NOTETR,
     %                     NBTECF, NOTECF,
     %                     MXFACO, LEFACO, NO0FAR, NB0FAR,
     %                     NBSTCF, NOSTCF, MXTRCF, NBTRCF, NOTRCF,
     %                     NBCF,   N1ARCF, NOARCF, N1FEOC, NFETOI,
     %                     MXFAN2, NBFAN2, NOFAN2, NDFAN2,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     AJOUTER AUX TRIANGLES NOTRCF DU CF DE LA FACE D'ANGLE
C -----     MAXIMAL / TRIANGLE DU CF DU TETRAEDRE S'ENROULANT
C           AUTOUR DE L'ARETE N'APPARTENANT PAS A 2 ET SEULEMENT 2 FACES
C           DE L'ETOILE

C ENTREES:
C --------
C KTITRE : TITRE DU TRACE COMPLETE PAR LE NOMBRE DE TETRAEDRES
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 1
C          1: NUMERO DU TETRAEDRE DANS NOTETR
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3  =0
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C          EN VERSION 2 (NFETOI(3,.)>0)
C          1: NUMERO NOTETR DU TETRAEDRE AU DELA DE LA FACE
C             =0 SI INCONNU
C          2: NUMERO PTXYZD DU 1-ER SOMMET DE LA FACE
C          3: NUMERO PTXYZD DU 2-ME SOMMET DE LA FACE
C          4: NUMERO PTXYZD DU 3-ME SOMMET DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE
C             L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C NBSTCF : NOMBRE DE SOMMETS DU TABLEAU NOSTCF
C NOSTCF : NUMERO PTXYZD DES SOMMETS DES TRIANGLES DU CF

C MXTRCF : NOMBRE MAXIMAL DE TRIANGLES DU TABLEAU NOTRCF
C NBTRCF : NOMBRE         DE TRIANGLES DU TABLEAU NOTRCF
C          C-A-D DU POLYGONE ENCORE DIT ENSUITE ETOILE
C NOTRCF : SI NOTRCF(*)>0 NUMERO LEFACO DU TRIANGLE
C                      <0 NUMERO NO0FAR DU TRIANGLE AJOUTE AU CF

C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C          NORMALE VERS L'INTERIEUR DU TETRAEDRE LA CONTENANT
C NB0FAR : NOMBRE DE FACES AJOUTEES AU CF

C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C N1ARCF : NUMERO DU 1-ER SOMMET OU 1-ERE ARETE DES LIGNES DU CONTOUR FERME
C          N1ARCF(0) POINTEUR SUR LA PREMIERE ARETE VIDE
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE

C NBTECF : NOMBRE DE TETRAEDRES DE LA LISTE DES TETRAEDRES DE L'ETOILE
C NOTECF : NUMERO DANS NOTETR DES TETRAEDRES DE L'ETOILE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C MXFAN2 : NOMBRE MAXIMAL D'ENTIERS DU TABLEAU NOFAN2
C NBFAN2 : NOMBRE D'ARETES DU CF N'APPARTENANT PAS A 2 ET SEULEMENT 2 FACES
C NOFAN2 : NOFAN2(NDFAN2+1)=NOMBRE NBF DE FACES NFETOI DE L'ARETE
C          NOFAN2(NDFAN2+2)=NUMERO PTXYZD DU SOMMET 1  DE L'ARETE
C          NOFAN2(NDFAN2+3)=NUMERO PTXYZD DU SOMMET 2  DE L'ARETE
C          NOFAN2(NDFAN2+3+1...NBF)=NUMERO NF DANS NFETOI DE LA FACE
C NDFAN2 : ADRESSE-1 DANS NOFAN2 DE LA PREMIERE DONNEE DE L'ARETE A TRAITER

C SORTIES:
C --------
C NBTRCF : NOMBRE DE TRIANGLES DU TABLEAU NOTRCF Y COMPRIS CEUX AJOUTES
C NOTRCF : SI NOTRCF(*)>0 NUMERO LEFACO DU TRIANGLE
C                      <0 NUMERO NO0FAR DU TRIANGLE AJOUTE AU CF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2017
C2345X7..............................................................012
      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*), COS2FA, COSMIN
      INTEGER           NOTETR(8,*), NOTECF(NBTECF),
     %                  LEFACO(11,0:MXFACO), NO0FAR(3,*),NOSTCF(NBSTCF),
     %                  NOTRCF(MXTRCF), N1ARCF(0:*), NOARCF(1:3,*),
     %                  NFETOI(5,*),    NOFAN2(MXFAN2)
      INTEGER           NOSOTR(3), NOSOTR2(3), NFA(2)
      INTEGER           NOSOAR(2,6),
     %                  NOSOFATE(3,4),
     %                  NOFAARTE(2,6)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
      DATA              NOSOAR   / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE
C     NUMERO DES 2 FACES DU TETRAEDRE AYANT UNE ARETE COMMUNE K
      DATA              NOFAARTE / 1,4,  1,2,  1,3,  3,4,  2,4,  2,3 /

      KTITRE='ajtraucf: FACE PERDUE         AJOUT DE TRIANGLES AU CF'
      WRITE( KTITRE(23:28), '(I6)' ) NFP
      CALL SANSDBL( KTITRE, L )
      PRINT *, KTITRE(1:L)
      CALL TRFETO11( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %               NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %               NBTECF, NOTECF, NOTETR,
     %               NBFAN2, NOFAN2 )

C     PARCOURS DES ARETES N'APPARTENANT PAS A 2 ET SEULEMENT 2 FACES
C     --------------------------------------------------------------
      NBTRCF0 = NBTRCF

C     NOMBRE DE FACES NFETOI DE L'ARETE NAS12 NS1-NS2
      NBF = NOFAN2( NDFAN2 + 1 )
      IF( NBF .NE. 0 ) THEN
         IERR=9
         RETURN
      ENDIF

C     NUMERO DES 2 SOMMETS DE L'ARETE DANS AUCUNE FACE SIMPLE DE L'ETOILE
      NS1 = NOFAN2( NDFAN2 + 2 )
      NS2 = NOFAN2( NDFAN2 + 3 )
      NCF = NOFAN2( NDFAN2 + 4 )
      PRINT*,'ajtraucf: TRAITEMENT ARETE du CF de St',NS1,NS2,
     %       ' dans',NBF,' FACE SIMPLE de l''ETOILE'

C     -------------------------------------------------------------
C     0 FACE DE L'ETOILE POUR CETTE ARETE NS1-NS2 DU CF
C     AJOUT AUX TRIANGLES DU CF DE LA FACE DU TETRAEDRE S'ENROULANT
C     AUTOUR DE L'ARETE NS1-NS2 D'ANGLE MAXIMAL / TRIANGLE DU CF
C     -------------------------------------------------------------
C     RETROUVER LE TRIANGLE NTR DU CF D'ARETE NS1-NS2
      NS3 = 0
      DO K=NBTRCF,1,-1

C        NUMERO DU TRIANGLE K DU CF
         NTR = NOTRCF(K)
         IF( NTR .GT. 0 ) THEN

C           NTR EST UNE FACE DE LEFACO
            DO L=1,3
               NAS1 = LEFACO(L,NTR)
               IF( L .EQ. 3 ) THEN
                  LL = 1
               ELSE
                  LL = L+1
               ENDIF
               NAS2 = LEFACO(LL,NTR)
               IF( NAS1 .EQ. NS1 .AND. NAS2 .EQ. NS2 .OR.
     %             NAS2 .EQ. NS1 .AND. NAS1 .EQ. NS2 ) THEN
C                  L'ARETE L DU TRIANGLE NTR EST L'ARETE NA1 DU CF
C                  NS3 NUMERO DU SOMMET DU TRIANGLE NTR NON NS1 ou NS2
                  IF( LL .EQ. 3 ) THEN
                     LLL = 1
                  ELSE
                     LLL = LL+1
                  ENDIF
                  NS3 = LEFACO(LLL,NTR)
                  GOTO 10
               ENDIF
            ENDDO

         ELSE

C           NTR EST UNE FACE NO0FAR AJOUTEE AU CF
            NTRA = -NTR
            DO L=1,3
               NAS1 = NO0FAR(L,NTRA)
               IF( L .EQ. 3 ) THEN
                  LL = 1
               ELSE
                  LL = L+1
               ENDIF
               NAS2 = NO0FAR(LL,NTRA)
               IF( NAS1 .EQ. NS1 .AND. NAS2 .EQ. NS2 .OR.
     %             NAS2 .EQ. NS1 .AND. NAS1 .EQ. NS2 ) THEN
C                  L'ARETE L DU TRIANGLE NTRA EST L'ARETE NA1 DU CF
C                  NS3 NUMERO DU SOMMET DU TRIANGLE NTRA NON NS1 ou NS2
                  IF( LL .EQ. 3 ) THEN
                     LLL = 1
                  ELSE
                     LLL = LL+1
                  ENDIF
                  NS3 = NO0FAR(LLL,NTRA)
                  GOTO 10
               ENDIF
            ENDDO

         ENDIF

      ENDDO


C     LE TRIANGLE NTR DE LEFACO OU -NTR DE NO0FAR A
C     L'ARETE NS1-NS2 DU CF
C     RECHERCHE DES TETRAEDRES DE L'ETOILE D'ARETE NS1-NS2
C     ET DE CELUI D'ANGLE MAXIMAL AVEC LE TRIANGLE NTR
C     ----------------------------------------------------
 10   COSMIN = 2D0
      NFAMIN = 0
      NTEMIN = 0
      NMIN   = 0
      DO N = 1, NBTECF

C        NUMERO NOTETR DU TETRAEDRE N DE L'ETOILE
         NTE = NOTECF( N )

         DO K=1,6

C           ARETE K DU TETRAEDRE NTE
            NAS1 = NOTETR( NOSOAR(1,K), NTE )
            NAS2 = NOTETR( NOSOAR(2,K), NTE )

            IF( NAS1 .EQ. NS1 .AND. NAS2 .EQ. NS2 .OR.
     %          NAS2 .EQ. NS1 .AND. NAS1 .EQ. NS2 ) THEN

C              L'ARETE K DU TETRAEDRE NTE EST L'ARETE NS1-NS2 DU CF
C              NOFAARTE(1:2,K) NUMEROS DES 2 FACES AYANT CETTE ARETE K
               NFA(1) = NOFAARTE(1,K)
               NFA(2) = NOFAARTE(2,K)

               DO MM=1,2

C                 NS4 NUMERO DU SOMMET DE LA FACE NFA(MM) DE NTE
C                 DIFFERENT DE NAS1 ET NAS2
                  DO KK=1,3
                     NS4 = NOTETR( NOSOFATE(KK,NFA(MM)), NTE )
                     IF( NS4 .NE. NAS1 .AND. NS4 .NE. NAS2)THEN
                        GOTO 15
                     ENDIF
                  ENDDO

C                 COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 FACES
C                 CELLES DU CF NAS1-NAS2-NS3 et du TETRAEDRE
C                 NAS1-NAS2-NS4 ORIENTEES SELON LE SENS NAS1-NAS2
 15               CALL COS2TD( PTXYZD(1,NAS1),
     %                         PTXYZD(1,NAS2),
     %                         PTXYZD(1,NS3),
     %                         PTXYZD(1,NAS1),
     %                         PTXYZD(1,NAS2),
     %                         PTXYZD(1,NS4),
     %                         COS2FA, IERR1, IERR2 )
                  IF( IERR1.EQ.0 .AND. IERR2.EQ.0 .AND.
     %                COS2FA .LT. COSMIN ) THEN
                     COSMIN = COS2FA
                     NMIN   = N
C                    NUMERO NOTETR DU TETRAEDRE
                     NTEMIN = NTE
C                    NUMERO LOCAL DE SA FACE
                     NFAMIN = NFA(MM)
C                    AUTRE SOMMET QUE NS1-NS2
                     NS4MIN = NS4
                  ENDIF

               ENDDO

            ENDIF
         ENDDO

      ENDDO

C     LA FACE NFAMIN DU TETRAEDRE NTEMIN DE COSINUS MINIMAL EST NMIN
      IF( NMIN .EQ. 0 ) THEN
C        PAS DE FACE DE TETRAEDRE D'ARETE NS1-NS2
C        => ABANDON DE L'AJOUT DE FACE AU CF
         PRINT*,'ajtraucf: arete CF',ns1,ns2,
     %          'dans AUCUNE FACE des TETRAEDRES'
C        ESSAI DE REPRISE=5
         IERR = 5
         GOTO 8000
      ENDIF

cccC              LE TETRAEDRE NTEMIN EST SUPPRIME DES TETRAEDRES DE L'ETOILE
ccc               PRINT *,'ajtraucf: 10) arete',ns1,ns2,'dans',nbfs12,' faces',
ccc     %         ' => retrait tetraedre  NOTECF(NMIN)=',NOTECF(NMIN),
ccc     %         ' st:',(notetr(kl,NOTECF(NMIN)),kl=1,8)
ccc               NOTECF( NMIN ) = 0

ccc               PRINT*,'ajtraucf: 10) arete CF',ns1,ns2,'dans',nbfs12,' faces'
ccc     %         ,' => ajout dans le CF de la face',NFAMIN,' du tetra',
ccc     %         NTEMIN,' St:',(NOTETR(KL,NTEMIN),KL=1,8)


C     LA FACE NFAMIN DU TETRAEDRE NTEMIN EST AJOUTEE
C     DANS NOTRCF AUX TRIANGLES DU CF ET NOARCF ARETES DU CF
C     ------------------------------------------------------
      IF( NBTRCF .GE. MXTRCF ) THEN
         PRINT*,'ajtraucf: TROP DE FACES DANS LE CF. AUGMENTER MXTRCF'
         IERR = 3
         GOTO 8000
      ENDIF

C     LA FACE NFAMIN DE NTEMIN EST ELLE UNE FACE LEFACO?
      CALL NULEFT( NFAMIN, NTEMIN, NOTETR, MXFACO, LEFACO, NF1 )
      IF( NF1 .GT. 0 ) THEN

C        LA FACE NFAMIN EST LA FACE NF1 DANS LEFACO
         NBTRCF = NBTRCF + 1
         NOTRCF( NBTRCF ) = NF1
         NOLFAR = NF1

      ELSE

C        LA FACE N'EST PAS DANS LEFACO.  EST ELLE UNE FACE NFETOI?
         DO L=1,3
C           NUMERO DES SOMMETS DE LA FACE
            NOSOTR(L) = NOTETR( NOSOFATE(L,NFAMIN), NTEMIN )
         ENDDO
         CALL TRI3NO( NOSOTR, NOSOTR )

         NF1 = N1FEOC
 18      IF( NF1 .GT. 0 ) THEN

            IF( NFETOI(3,NF1) .EQ. 0 ) THEN
C              SI NFETOI EN VERSION 1
C              NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
               NTE1 = ABS( NFETOI(1,NF1) )
C              NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE1
               I1 = ABS( NFETOI(2,NF1) )
C              LA FACE EST REINITIALISEE NON TRAITEE
               NFETOI(2,NF1) = I1
C              LE TETRAEDRE OPPOSE PAR LA FACE I1 DE NTE1
               NTEOP = NOTETR( 4+I1, NTE1 )
C              LES 3 SOMMETS DE LA FACE SIMPLE DE L'ETOILE
C              NORMALE VERS L'INTERIEUR DU TETRAEDRE NTE1 (DONC DE L'ETOILE)
C              CAR PERMUTATION DES SOMMETS 2 ET 3 DE NOSOFATE
               NOSOTR2(1) = NOTETR( NOSOFATE(1,I1), NTE1 )
               NOSOTR2(2) = NOTETR( NOSOFATE(3,I1), NTE1 )
               NOSOTR2(3) = NOTETR( NOSOFATE(2,I1), NTE1 )
            ELSE
C              NFETOI EN VERSION 2
C              LE TETRAEDRE OPPOSE PAR LA FACE I1 DE NTE1
               NTEOP = NFETOI(1,NF1)
               NOSOTR2(1) = NFETOI(2,NF1)
               NOSOTR2(2) = NFETOI(3,NF1)
               NOSOTR2(3) = NFETOI(4,NF1)
            ENDIF
            CALL TRI3NO( NOSOTR2, NOSOTR2 )
            IF( NOSOTR(1) .EQ. NOSOTR2(1) .AND.
     %          NOSOTR(2) .EQ. NOSOTR2(2) .AND.
     %          NOSOTR(3) .EQ. NOSOTR2(3) ) THEN
C              OUI: ELLE EST DEJA DANS NFETOI
C              AJOUT EVENTUEL DU TETRAEDRE OPPOSE DANS L'ETOILE
               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP).GT. 0 ) THEN
                  DO K=1,NBTECF
                     IF( NOTECF(K) .EQ. NTEOP ) GOTO 9000
                  ENDDO
C                 LE TETRAEDRE OPPOSE EST AJOUTE A L'ETOILE
                  NBTECF = NBTECF + 1
                  NOTECF( NBTECF ) = NTEOP
               ENDIF
               GOTO 9000
            ENDIF

C           PASSAGE A LA FACE SUIVANTE DE L'ETOILE
            NF1 = NFETOI(5,NF1)
            GOTO 18

         ENDIF

C        LA FACE N'EST PAS DANS LEFACO et NFETOI
C        EST ELLE DEJA DANS NO0FAR?
         DO K=1,NB0FAR
            DO L=1,3
               NOSOTR2(L) = NO0FAR(L,K)
            ENDDO
            CALL TRI3NO( NOSOTR2, NOSOTR2 )
            IF( NOSOTR(1) .EQ. NOSOTR2(1) .AND.
     %          NOSOTR(2) .EQ. NOSOTR2(2) .AND.
     %          NOSOTR(3) .EQ. NOSOTR2(3) ) THEN
C              OUI: ELLE EST DEJA DANS NO0FAR
               GOTO 9000
            ENDIF
         ENDDO

C        NON. LA FACE EST AJOUTEE A NO0FAR
C        .................................
         NB0FAR = NB0FAR + 1
C        LE SIGNE - DANS NOTRCF POUR INDIQUER 3 SOMMETS RANGES
C        DANS NO0FAR

         NOLFAR = -NB0FAR
         IF( NBTRCF .GE. MXTRCF ) THEN
            PRINT*,'ajtraucf: TROP DE FACES CF. AUGMENTER MXTRCF'
            IERR = 3
            GOTO 8000
         ENDIF

         NBTRCF = NBTRCF + 1
         NOTRCF( NBTRCF ) = NOLFAR

C        NUMERO DES SOMMETS DE LA FACE, NORMALE VERS
C        L'INTERIEUR DU TETRAEDRE
         NO0FAR(1,NB0FAR)=NOTETR(NOSOFATE(1,NFAMIN), NTEMIN)
         NO0FAR(2,NB0FAR)=NOTETR(NOSOFATE(3,NFAMIN), NTEMIN)
         NO0FAR(3,NB0FAR)=NOTETR(NOSOFATE(2,NFAMIN), NTEMIN)

         PRINT*,'ajtraucf: AJOUT NO0FAR(',NB0FAR,'):',
     %          (NO0FAR(L,NB0FAR),L=1,3)

      ENDIF

C     APRES L'AJOUT DE LA FACE NBTRCF DU TABLEAU NOTRCF C-A-D
C     DE LEFACO(NOTRCF(NBTRCF)) OU NO0FAR(-NOTRCF(NBTRCF)
C     MISE A JOUR DES ARETES SIMPLES, DES SOMMETS,
C     DES SOMMETS ISOLES DU CONTOUR FERME CF
C     --------------------------------------------------------
      CALL MJDUCF( PTXYZD, MXFACO, LEFACO, NO0FAR,
     %             MXTRCF, NBTRCF, NOTRCF,
     %             MXARCF, NBCF,   N1ARCF, NOARCF,
     %             MXSTCF, NBSTCF, NOSTCF,
     %             MXSTIS, NBSTIS, NOSTIS,  IERR )


 8000 KTITRE='ajtraucf: FACE PERDUE         FIN AJOUT.        TRIANGLES 
     %FORMENT le CF'
      WRITE( KTITRE(23:28), '(I6)' ) NFP
      WRITE( KTITRE(42:47), '(I6)' ) NBTRCF
      CALL SANSDBL( KTITRE, L )
      PRINT *, KTITRE(1:L)
      CALL TRFETO11( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %               NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %               NBTECF, NOTECF, NOTETR,
     %               NBFAN2, NOFAN2 )

C     LISTE DES TRIANGLES DU CF ET AJOUTES
 9000 DO K=1,NBTRCF
         NF = NOTRCF( K )
         IF( NF .GT. 0 ) THEN
            PRINT*,'ajtraucf: TRIANGLE du CF:', NF,' dans LEFACO St:',
     %             (LEFACO(L,NF),L=1,3)
         ELSE
            PRINT*,'ajtraucf: TRIANGLE du CF:', NF,' dans NO0FAR St:',
     %             (NO0FAR(L,-NF),L=1,3)
         ENDIF
      ENDDO

      RETURN
      END
