      SUBROUTINE MJDUCF( PTXYZD, MXFACO, LEFACO, NO0FAR,
     %                   MXTRCF, NBTRCF, NOTRCF,
     %                   MXARCF, NBCF,   N1ARCF, NOARCF,
     %                   MXSTCF, NBSTCF, NOSTCF,
     %                   MXSTIS, NBSTIS, NOSTIS,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    APRES L'AJOUT DE LA FACE NBTRCF DU TABLEAU NOTRCF C-A-D
C -----    DE LEFACO(NOTRCF(NBTRCF)) OU NO0FAR(-NOTRCF(NBTRCF)
C          MISE A JOUR DES ARETES SIMPLES, DES SOMMETS,
C          DES SOMMETS ISOLES DU CONTOUR FERME CF
C
C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C MXFACO : NOMBRE MAXIMAL DE TRIANGLES PERMIS POUR LA TRIANGULATION
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C MXTRCF : NOMBRE MAXIMAL DE TRIANGLES DANS NOTRCF
C MXSTCF : NOMBRE MAXIMAL DE SOMMETS DU CF DECLARABLES DANS NOSTCF
C MXSTIS : NOMBRE MAXIMAL DE SOMMETS ISOLES DU CF DECLARABLES DANS NOSTIS
C
C ENTREES ET SORTIES:
C -------------------
C LEFACO : ENSEMBLE DES FACES=TRIANGLES DE LA PEAU OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          123: NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          45 : NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          678: NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C               ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C                          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          9  : LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C               HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C               LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C               NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C               SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C               NF  = LEFACO( 9, NF )  ...
C          10 : LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C          11 : LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
cccC          12 : LEFACO(12,.) = NO FACEOC DE 1 A NBFACES D'OC

C SORTIES :
C ---------
C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C NBTRCF : NOMBRE DE TRIANGLES DU TABLEAU NOTRCF SUR LEFACO
C          C-A-D DU POLYGONE ENCORE DIT ENSUITE ETOILE
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DE L'ETOILE
C N1ARCF : NUMERO DU 1-ER SOMMET OU 1-ERE ARETE DES LIGNES DU CONTOUR FERME
C          N1ARCF(0) POINTEUR SUR LA PREMIERE ARETE VIDE
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE
C NBSTCF : NOMBRE DE SOMMETS DU CF
C NOSTCF : NUMERO PTXYZD DES NBSTCF SOMMETS PERIPHERIQUES DU CF
C NBSTIS : NOMBRE DE SOMMETS ISOLES DANS L'ETOILE
C NOSTIS : NUMERO DES SOMMETS ISOLES N'APPARTENANT PAS AU CONTOUR
C IERR   : 0 SI PAS D'ERREUR
C          1 ANOMALIE 2 TETRAEDRES NON REELLEMENT OPPOSES PAR UNE FACE
C          3 NBTRCF>MXTRCF => TROP DE TRIANGLES PERDUS COAGULES
C          4 NBSTIS>MXSTIS => TROP DE SOMMETS ISOLES
C          5 NBSTCF>MXSTCF => TROP DE SOMMETS DU CF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  St PIERRE DU PERRAY             Fevrier 2018
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE
      LOGICAL                  TRACTE0
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(MXTRCF), NOSTCF(MXSTCF), NOSTIS(MXSTIS),
     %                  N1ARCF(0:*), NOARCF(3,MXARCF),
     %                  NOSOTR(3), NAR(3), NAR0(3)
      CHARACTER*128     KTITRE

      TRACTE0 = TRACTE

C     LES 3 SOMMETS DE LA NOUVELLE FACE NBTRCF AJOUTEE AU CF
C     ------------------------------------------------------
      NFLFAR = 0
      NFL = NOTRCF( NBTRCF )

      IF( NFL .GT. 0 ) THEN
         NOSOTR(1) = LEFACO( 1, NFL )
         NOSOTR(2) = LEFACO( 2, NFL )
         NOSOTR(3) = LEFACO( 3, NFL )
      ELSE
         NFLFAR = -NFL
         NOSOTR(1) = NO0FAR( 1, NFLFAR )
         NOSOTR(2) = NO0FAR( 2, NFLFAR )
         NOSOTR(3) = NO0FAR( 3, NFLFAR )
      ENDIF

C     AJOUT OU RETRAIT DES 3 ARETES DE LA NOUVELLE FACE NOSOTR DU CF
C     DANS LE TABLEAU NOARCF. RECHERCHE DES ARETES SIMPLES DU CF
C     DU TRIANGLE NOSOTR
C     --------------------------------------------------------------
      NBAREP = 0
      NBAREN = 0
      NBARE  = 0
      DO K=1,3

         IF( K .NE. 3 ) THEN
            KK = K+1
         ELSE
            KK = 1
         ENDIF

C        NOSOTR(K)-NOSOTR(K+1) EST ELLE UNE ARETE DU CF ACTUEL?
         NSK1 = NOSOTR(K)
         NSK2 = NOSOTR(KK)
         CALL ARARCF( NSK1, NSK2, NBCF, N1ARCF, NOARCF,
     %                NLFCF, NAR0(K), NAR(K) )
C        NAR(K)>0   NUMERO NOARCF DE L'ARETE NSK1->NSK2 DANS NOARCF
C              <0  -NUMERO NOARCF DE L'ARETE NSK2->NSK1 DANS NOARCF
C              =0   ARETE NSK1-NSK2 NON DANS NOARCF C-A-D LE CF
C        NAR0(K)>0  NUMERO NOARCF DE L'ARETE QUI PRECEDE NAR(K) DANS NOARCF
C               =0  ARETE NSK1-NSK2 NON DANS NOARCF C-A-D LES ARETES SIMPLES CF

         IF( NAR(K) .NE. 0 ) THEN
C           UNE ARETE SIMPLE DE PLUS DU CF POUR CETTE FACE NBTRCF
            NBARE = NBARE + 1
            IF( NAR(K) .GT. 0 ) THEN
               NBAREP = NBAREP + 1
            ELSE
               NBAREN = NBAREN + 1
            ENDIF
         ENDIF

      ENDDO

ccc         IF( NBAREP .GT. 0 .AND. NBAREN .GT. 0 ) THEN
cccC           PROBLEME UNE ARETE DE NOSOTR DANS LE SENS INVERSE
cccC           ET UNE DANS LE MEME SENS
ccc            PRINT *,'PB mjducf: 10) NS1=',NS1,' NS2=',NS2,
ccc     %              ' NAR0=',NAR0,' NAR=',NAR,
ccc     %              ' AVEC DES SENS INCORRECTS'
ccc         ENDIF

      IF( NBAREP .GT. 0 ) THEN

C        NOSOTR N'EST PAS DANS LE BON SENS PAR RAPPORT
C        AUX ARETES SIMPLES PERIPHERIQUES DU CF
C        ---------------------------------------------
C        PERMUTATION DES SOMMETS 2<->3
         KKK         = NOSOTR( 2 )
         NOSOTR( 2 ) = NOSOTR( 3 )
         NOSOTR( 3 ) = KKK

C        PERMUTATION DES ARETES 1<->3
         KKK     = NAR0(1)
         NAR0(1) = NAR0(3)
         NAR0(3) = KKK

         KKK    = NAR(1)
         NAR(1) = NAR(3)
         NAR(3) = KKK

         IF( NFL .GT. 0 ) THEN
C           LEFACO(NFL) N'EST PAS DANS LE BON SENS ET EST RECTIFIEE
            LEFACO(1,NFL) = NOSOTR( 1 )
            LEFACO(2,NFL) = NOSOTR( 2 )
            LEFACO(3,NFL) = NOSOTR( 3 )
C           PERMUTATION DES ARETES 1<->3
            KKK           = LEFACO(6,NFL)
            LEFACO(6,NFL) = LEFACO(8,NFL) 
            LEFACO(8,NFL) = KKK
         ELSE
C           NO0FAR(NFLFAR) N'EST PAS DANS LE BON SENS ET EST RECTIFIEE
            NO0FAR(1,NFLFAR) = NOSOTR( 1 )
            NO0FAR(2,NFLFAR) = NOSOTR( 2 )
            NO0FAR(3,NFLFAR) = NOSOTR( 3 )
         ENDIF

      ELSE

C        LE NUMERO D'ARETE REDEVIENT POSITIF
         DO K=1,3
            NAR0(K) = ABS( NAR0(K) )
            NAR (K) = ABS( NAR (K) )
         ENDDO

      ENDIF

C     MISE A JOUR DES ARETES ET SOMMETS DU CF
C     AJOUT ET SUPPRESSION D'ARETES DU CF DANS NOARCF
C     -----------------------------------------------
      IF( NBARE .EQ. 1 ) THEN

C        TRIANGLE NOSOTR AVEC UNE ARETE DU CF A SUPPRIMER
C        ET LES 2 AUTRES ARETES A AJOUTER
C        ------------------------------------------------
         DO K=1,3

            IF( NAR(K) .NE. 0 ) THEN

C              L'ARETE K DE NOSOTR = NOSOTR(K)-NOSOTR(KK)
C              EST A SUPPRIMER DU CF
C              LES 2 AUTRES ARETES DU CF SONT A AJOUTER AU CF
C              L'ARETE DU CF EST NAR(K) AVEC LE CHAINAGE
C              NAR0(K)->NAR(K)->NOARCF(2,NAR(K))

               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF

               IF( KK .NE. 3 ) THEN
                  KKK = KK+1
               ELSE
                  KKK = 1
               ENDIF

C              ACTUELLEMENT CES 3 ARETES SE SUIVENT
               NA0 = NAR0(K)
               NA1 = NAR(K)
               NA2 = NOARCF( 2, NA1 )
               PRINT*
               PRINT*,'mjducf: AVANT'
               print*,'NOARCF(',NA0,')=',(NOARCF(J,NA0),J=1,3)
               print*,'NOARCF(',NA1,')=',(NOARCF(J,NA1),J=1,3)
               print*,'NOARCF(',NA2,')=',(NOARCF(J,NA2),J=1,3)

C              RECHERCHE D'UNE ARETE VIDE NA DANS NOARCF
               NA = N1ARCF( 0 )
               IF( NA .LE. 0 ) THEN
                  PRINT*,'mjducf: AUGMENTER MXARCF TROP FAIBLE'
                  IERR = 3
                  GOTO 8000
               ENDIF
C              MISE A JOUR DE LA 1-ERE ARETE VIDE DANS NOARCF
               N1ARCF( 0 ) = NOARCF( 2, NA )
               print*,'Nouvelle arete NA=',NA,' N1ARCF(0)=',N1ARCF(0),
     %                ' NOSOTR=',NOSOTR

C              LA NOUVELLE ARETE EST INTERCALEE ENTRE NA1 et NA2
ccc            NOARCF( 2, NA0 ) = NA1   DEJA VRAI
               NOARCF( 2, NA1 ) = NA
               NOARCF( 2, NA  ) = NA2

C              L'ARETE NA1 EST NOSOTR(K)->NOSOTR(KK) et
C              EST REMPLACEE PAR L'ARETE NOSOTR(K)->NOSOTR(KKK)
C              SUIVIE DE L'ARETE NOSOTR(KKK)->NOSOTR(KK)
ccc            NOARCF( 1, NA1 ) = NOSOTR( KK  )   DEJA VRAI
               NOARCF( 1, NA  ) = NOSOTR( KKK )
ccc            NOARCF( 1, NA2 ) = NOSOTR( K   )   DEJA VRAI

C              LE NUMERO DU TRIANGLE DE L'AUTRE COTE DE CETTE ARETE
               NOARCF( 3, NA1 ) = NFL
               NOARCF( 3, NA  ) = NFL
               PRINT*,'mjducf: APRES'
               print*,'NOARCF(',NA0,')=',(NOARCF(J,NA0),J=1,3)
               print*,'NOARCF(',NA1,')=',(NOARCF(J,NA1),J=1,3)
               print*,'NOARCF(',NA, ')=',(NOARCF(J,NA ),J=1,3)
               print*,'NOARCF(',NA2,')=',(NOARCF(J,NA2),J=1,3)

C              NOSOTR(KKK) EST AJOUTE AUX SOMMETS DU CF
               NBSTCF = NBSTCF + 1
               NOSTCF(NBSTCF) = NOSOTR(KKK)

C              REDEPART A PARTIR DE L'ARETE NA1
               N1ARCF( NLFCF ) = NA1

               IERR = 0
               GOTO 8000

            ENDIF

         ENDDO


      ELSE IF( NBARE .EQ. 2 ) THEN

C        TRIANGLE NOSOTR AVEC 2 ARETES DU CF A SUPPRIMER
C        ET L'AUTRE ARETE NON CF A AJOUTER
         DO K=1,3

            IF( NAR(K) .EQ. 0 ) THEN

C              L'ARETE K DE NOSOTR N'APPARTIENT PAS AU CF EST A AJOUTER
C              LES 2 AUTRES ARETES KK et KKK APPARTIENNENT AU CF
C              ET ELLES SONT A SUPPRIMER

C              LA PREMIERE ARETE DU CF DE NOSOTR EST NAR(KKK)
C              CAR LE CF EST EN SENS INVERSE DE NOSOTR PAR SA
C              CONSTRUCTION
C              AVEC LE CHAINAGE NAR0(KKK)->NAR(KKK)->NAR(KK)

               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF

               IF( KK .NE. 3 ) THEN
                  KKK = KK+1
               ELSE
                  KKK = 1
               ENDIF

C              L'ARETE NOSOTR(K)-NOSOTR(KK) REMPLACE L'ARETE NAR(KK) DU CF
               NOARCF( 2, NAR0(KKK) ) = NAR(KK)

C              LE 1-ER SOMMET DE L'ARETE
               NOARCF( 1, NAR(KK) ) = NOSOTR(K)

C              L'ARETE NOARCF SUIVANTE
ccc               NOARCF( 2, NAR(KK) ) = NOARCF( 2, NAR(KK) ) deja fait

C              LE NUMERO DU TRIANGLE DE L'AUTRE COTE DE CETTE ARETE
               NOARCF( 3, NAR(KK) ) = NFL

C              L'ARETE NAR(KKK) DU CF EST SUPPRIMEE DU CF
               IF( NAR(KKK) .EQ. N1ARCF(NLFCF) ) THEN
                  N1ARCF( NLFCF ) = NAR0(KKK)
               ENDIF

C              NAR(KKK) DEVIENT UNE ARETE VIDE DU CF
               NOARCF( 2, NAR(KKK) ) = N1ARCF(0)
               N1ARCF(0) = NAR(KKK)

C              SUPPRESSION DU SOMMET NOSOTR(KKK) DES SOMMETS DU CF
               DO L=1,NBSTCF
                  IF( NOSTCF(L) .EQ. NOSOTR(KKK) ) GOTO 24
               ENDDO
 24            DO LL=L+1,NBSTCF
                  NOSTCF(LL-1) = NOSTCF(LL)
               ENDDO
               NBSTCF = NBSTCF - 1

C              REDEPART A PARTIR DE LA DERNIERE ARETE AJOUTEE AU CF
               N1ARCF( NLFCF ) = NAR(KK)

               IERR = 0
               GOTO 8000

            ENDIF

         ENDDO


      ELSE IF( NBARE .EQ. 3 ) THEN

C        TRIANGLE NOSOTR AVEC 3 ARETES DU CF
C        LES 3 ARETES DU CF SONT A SUPPRIMER DU CF
C        REPERAGE DE L'ARETE AVANT LA PREMIERE ARETE DU CF
         DO K=1,3
            NAR(K) = ABS( NAR(K) )
         ENDDO

         DO 26 K=1,3
            NAE = NAR0(K)
            DO L=1,3
               IF( NAE .EQ. NAR(L) ) GOTO 26
            ENDDO
            GOTO 28
 26      ENDDO

C        REPERAGE DE L'ARETE SUIVANTE DES 3 ARETES DU TRIANGLE
 28      DO 21 K=1,3
            NAS = NOARCF( 2, NAR(K) )
            DO L=1,3
               IF( NAS .EQ. NAR(L) ) GOTO 21
            ENDDO
            GOTO 25
 21      ENDDO

C         SUPPRESSION DES 3 ARETES DU TRIANGLE DANS NOARCF          
 25       DO K=1,3

C            L'ARETE NAR(K) DU CF EST SUPPRIMEE DU CF
             NA  = NAR(K)
             NA3 = NOARCF( 2, NA )

C            SUPPRESSION DE NA DANS LE CHAINAGE DU CF
             NOARCF( 2, NAR0(K) ) = NA3
             IF( NA .EQ. N1ARCF(NLFCF) ) THEN
                N1ARCF( NLFCF ) = NA3
             ENDIF

C            NA DEVIENT LA PREMIERE ARETE VIDE
             NOARCF( 2, NA ) = N1ARCF(0)
             N1ARCF(0) = NA

          ENDDO

C         CHAINAGE DE L'ARETE AVANT ENTREE DU TRIANGLE SUR
C         L'ARETE DE SORTIE DU TRIANGLE
          NOARCF( 2, NAE ) = NAS
          N1ARCF( NLFCF ) = NAE

C         CONSTRUCTION DU TABLEAU DES NBSTCF NO DES SOMMETS DU CF
          CALL CRSTCF( NBCF,   N1ARCF, NOARCF, NBARCF,
     %                 MXSTCF, NBSTCF, NOSTCF, IERR )

C         REDEPART A L'ARETE INITIALE DU CF
          IERR = 0
          GOTO 8000


       ELSE IF( NBARE .EQ. 0 ) THEN

C         SITUATION NORMALEMENT IMPOSSIBLE
          PRINT*,'PB mjducf: TRIANGLE',NOSOTR,
     %           ' SANS ARETE DU CF',
     %            NS1,NS2,' CAS IMPOSSIBLE!...'
          IERR = 1
          GOTO 9000

       ENDIF

C     AFFICHAGE ET TRACE FINAL
 8000 KTITRE='mjducf:        FACES TRCF        SOMMETS CF         SOMMET
     %S ISOLES'
      WRITE( KTITRE( 9:14),'(I6)' ) NBTRCF
      WRITE( KTITRE(27:32),'(I6)' ) NBSTCF
      WRITE( KTITRE(44:49),'(I6)' ) NBSTIS
      CALL SANSDBL( KTITRE, L )
      PRINT *, KTITRE(1:L)
      CALL TRDUCF( KTITRE, PTXYZD, NBTRCF, NOTRCF, NBSTIS, NOSTIS,
     %             LEFACO, NO0FAR )

C     LISTE DES TRIANGLES DU CF ET AJOUTES
 9000 DO K=1,NBTRCF
         NFL = NOTRCF( K )
         IF( NFL .GT. 0 ) THEN
            PRINT*,'mjducf: TRIANGLE du CF:', NFL,' dans LEFACO St:',
     %             (LEFACO(L,NFL),L=1,3)
         ELSE
            PRINT*,'mjducf: TRIANGLE du CF:', NFL,' dans NO0FAR St:',
     %             (NO0FAR(L,-NFL),L=1,3)
         ENDIF
      ENDDO

C     AFFICHAGE DES ARETES DES CONTOURS FERMES DU TABLEAU NOARCF
      CALL AFNOARCF( NBCF, N1ARCF, NOARCF )

      RETURN
      END
