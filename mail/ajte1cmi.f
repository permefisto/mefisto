      SUBROUTINE AJTE1CMI( KTITRE, PTXYZD, NOTETR, N1TETS,
     %                     MXFACO, LEFACO, NO0FAR, NBTRCF, NOTRCF, 
     %                     NBSTIS, NOSTIS, NBSTCF, NOSTCF, 
     %                     NBCF,   MXARCF, N1ARCF, NOARCF,
     %                     MXPTIN, NBPTIN, PTINTERS,
     %                     MXTECF, NBTECF, NOTECF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SI UNE ARETE SIMPLE DU CF N'A PAS 2 TETRAEDRES ADJACENTS DANS
C -----    L'ETOILE NOTECF ALORS
C          AJOUT A L'ETOILE NOTECF DES TETRAEDRES DE CETTE ARETE
C          RETRAIT DU TETRAEDRE DE FACE DE COSINUS MINIMAL et DE SES
C          TETRAEDRES OPPOSES AUX 2 FACES ADJACENTES A CETTE ARETE

C ENTREES:
C --------
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS

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
C NBSTIS : NOMBRE DE SOMMETS ISOLES DU CF
C NOSTIS : NUMERO PTXYZD DES NBSTIS SOMMETS ISOLES
C NBSTCF : NOMBRE DE SOMMETS DES ARETES PERIPHERIQUES DU CF
C NOSTCF : NUMERO PTXYZD DES NBSTCF SOMMETS DES ARETES PERIPHERIQUES

C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C MXARCF : MAXIMUM D'ARETES DECLARABLES DANS N1ARCF et NOARCF
C N1ARCF : POINTE SUR LE DEBUT DES ARETES DE CHAQUE LIGNE FERMEE DU CF
C          0 POUR LES PLACES VIDES
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE

C MXPTIN : NOMBRE MAXIMAL DE POINTS D'INTERSECTION ARETE TRIANGLE
C NBPTIN : NOMBRE DE POINTS D'INTERSECTION
C PTINTERS : XYZ DES NBPTIN POINTS D'INTERSECTION

C SORTIES:
C --------
C NBTECF : NOMBRE DE TETRAEDRES DE L'ETOILE
C NOTECF : NUMERO NOTETR DES NBTECF TETRAEDRES DE L'ETOILE
C IERR   : =0 SI PAS D'ERREUR DETECTEE
C          =1 UNE ANORMALE ARETE du CF APPARTIENT A AUCUN TETRAEDRE
C             du MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY            Avril 2018
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0
      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*), PTINTERS(3,MXPTIN)
      INTEGER           NOTETR(8,*), N1TETS(*),
     %                  LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(NBTRCF), NOSTCF(NBSTCF), NOSTIS(NBSTIS),
     %                  N1ARCF(0:*), NOARCF(1:3,1:*), NOTECF(MXTECF)
ccc      INTEGER           NOSOTR(3)

      TRACTE0 = TRACTE
ccc      tracte  = .true.
      IERR    = 0

C     ---------------------------------------------------------------
C     SI UNE ARETE DU CF N'A PAS 2 TETRAEDRES ADJACENTS DANS L'ETOILE
C     ALORS RETRAIT DU TETRAEDRE DE FACE DE COSINUS MINIMAL
C     ---------------------------------------------------------------
      NBARPB  = 0
      NBTECF0 = NBTECF
      DO NCF=1, NBCF

C        LA PREMIERE ARETE DU CF NCF
         NA0 = N1ARCF(NCF)
         IF( NA0 .LE. 0 ) THEN
            print*,'ajte1cmi: ERREUR N1ARCF=',NA0,' VALEUR INCORRECTE'
            GOTO 9999
         ENDIF
         NA1 = NA0

C        NUMERO NOARCF DE L'ARETE SUIVANTE DE L'ARETE NA1 DU CF
 40      NA2 = NOARCF(2,NA1)

C        LES 2 SOMMETS DE L'ARETE NA1->NA2 DU CF
         NSA1 = NOARCF(1,NA1)
         NSA2 = NOARCF(1,NA2)
         NBTEAR = 0

         DO 70 N=1,NBTECF

            NTE = NOTECF( N )

C           NTE TETRAEDRE DE SOMMETS NS1 NS2?
            DO N1=1,4
               NS1 = NOTETR(N1,NTE)
               IF( NS1 .EQ. NSA1 ) GOTO 50
            ENDDO
            GOTO 70

 50         DO N2=1,4
               NS2 = NOTETR(N2,NTE)
               IF( NS2 .EQ. NSA2 ) GOTO 60
            ENDDO
            GOTO 70

C           L'ARETE NS1-NS2 DE NTE EST L'ARETE NA1 DU CF
 60         NBTEAR = NBTEAR + 1
            IF( NBTEAR .GE. 2 ) GOTO 80

 70      ENDDO

C        L'ARETE NSA1-NSA2 DU CF N'EST PAS UNE ARETE DOUBLE
C        DES TETRAEDRES NOTECF.
C        ELLE PEUT APPARTENIR A 1 ou AUCUN TETRAEDRE DE NOTECF
C        -----------------------------------------------------
         PRINT*
         PRINT*,'ajte1cmi: ARETE',NA1,' du CF',NSA1,NSA2,
     %          ' NON ARETE DOUBLE des TETRAEDRES NOTECF MULTPLICITE=',
     %           NBTEAR

cccC        LE NUMERO DU TRIANGLE OPPOSE AU TRIANGLE DU CF AU DELA DE CETTE ARETE
ccc         NTROPCF = NOARCF( 3, NA1 )
ccc         IF( NTROPCF .GT. 0 ) THEN
cccC           NTROPCF EST UNE FACE LEFACO QUI APPARTIENT A UN TETRAEDRE
ccc            NTEOP = LEFACO( 11, NTROPCF )
ccc            PRINT*,'ajte1cmi: la FACE LEFACO(',NTROPCF,')=',
ccc     %              (LEFACO(kk,NTROPCF),kk=1,11)
cccC           EST ELLE UNE FACE DU TETRAEDRE NTEOP?
ccc            DO K=1,3
ccc               NOSOTR( K ) = LEFACO( K, NTROPCF )
ccc            ENDDO
ccc            CALL TRI3NO( NOSOTR, NOSOTR )
ccc            CALL NO1F1T( NOSOTR, NOTETR(1,NTEOP), NF )
ccc            IF( NF .EQ. 0 ) THEN
ccc               PRINT*,'ajte1cmi: N''EST PAS UNE FACE du TETRAEDRE',
ccc     %                 NTEOP,':',(NOTETR(kk,NTEOP),kk=1,8)
ccc            ELSE
ccc               PRINT*,'ajte1cmi: EST UNE FACE du TETRAEDRE',
ccc     %                 NTEOP,':',(NOTETR(kk,NTEOP),kk=1,8)
ccc            ENDIF
ccc         ENDIF

C        AJOUT DE TOUS LES TETRAEDRES D'ARETE NSA1 NSA2 A NOTECF
         CALL TETR1A( NSA1,  NSA2,  N1TETS, NOTETR,
     %                NBTEA, MXTECF-NBTECF, NOTECF(NBTECF+1), IERR )

         IF( NBTEA .LE. 0 ) THEN

            PRINT*,'ajte1cmi: PB l''ARETE NOARCF',NA1,' du CF:',
     %              NSA1,NSA2,' EST dans AUCUN TETRAEDRE'
            IERR = 1
C           TRACE DES FACES TRIANGULAIRES DU CF
            tracte = .true.
            KTITRE='ajte1cmi: ANORMALE ARETE                  DANS AUCUN
     % TETRAEDRE du MAILLAGE'
            WRITE(KTITRE(26:32),'(I7)') NSA1
            WRITE(KTITRE(34:40),'(I7)') NSA2
            CALL TRCFFAPE( KTITRE, PTXYZD, MXFACO, LEFACO, NO0FAR,
     %                     NOTETR, NBTECF, NOTECF, NBTRCF, NOTRCF,
     %                     MXARCF, NBCF,   N1ARCF, NOARCF,
     %                     NBSTIS, NOSTIS, NBPTIN, PTINTERS )
            GOTO 9999

         ENDIF
         NBTECF = NBTECF + NBTEA

C        UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
         CALL UNITABL( NOTECF, NBTECF )

C        RECHERCHE DU TRIANGLE DU CF D'ARETE NSA1-NSA2
         DO N=1,NBTRCF
            NTR = NOTRCF( N )
            DO NA = 1, 3

               IF( NA .EQ. 3 ) THEN
                  N2 = 1
               ELSE
                  N2 = NA+1
               ENDIF

               IF( NTR .GT. 0 ) THEN
                  NS1 = LEFACO(NA,NTR)
                  NS2 = LEFACO(N2,NTR)
               ELSE
                  NS1 = NO0FAR(NA,-NTR)
                  NS2 = NO0FAR(N2,-NTR)
               ENDIF

C              NS1-NS2 EST ELLE L'ARETE NSA1-NSA2 du CF?
               IF( NS1 .EQ. NSA1 .AND. NS2 .EQ. NSA2 .OR.
     %             NS1 .EQ. NSA2 .AND. NS2 .EQ. NSA1 ) THEN
C                 OUI: L'ARETE NS1-NS2 EST L'ARETE NA1 DU CF
C                 RECHERCHE DU 3-EME SOMMET DE NTR
                  IF( N2 .EQ. 3 ) THEN
                     N3 = 1
                  ELSE
                     N3 = N2+1
                  ENDIF
                  IF( NTR .GT. 0 ) THEN
                     NS3 = LEFACO(N3,NTR)
                  ELSE
                     NS3 = NO0FAR(N3,-NTR)
                  ENDIF

C                 RETIRER DE NOTECF LE TETRAEDRE NTCMIN DE SOMMETS
C                 NS1 et NS2 ET DONT UNE DES 2 FACES FORMENT UN ANGLE
C                 DE COSINUS MINIMAL AVEC LE TRIANGLE NS1 NS2 NS3
C                 ---------------------------------------------------------
                  CALL TETR1RCMI( NSA1, NSA2, NS3, NOTETR, PTXYZD,
     %                            NBSTCF, NOSTCF, NBTECF, NOTECF,
     %                            NTCMIN )

                  IF( NTCMIN .EQ. 0 ) THEN
C                    PAS DE RETRAIT DE TETRAEDRE DANS NOTECF
                     PRINT*,'ajte1cmi: PROBLEME PAS de RETRAIT DE TETRAE
     %DRE SUR L''ARETE',NSA1,NSA2,' du CF'
                     NBARPB = NBARPB + 1
                  ENDIF

C                 UNICITE DES NO>0 DES TETRAEDRES DU TABLEAU NOTECF
                  CALL UNITABL( NOTECF, NBTECF )

                  GOTO 80
               ENDIF
            ENDDO
         ENDDO

 80      IF( NA2 .NE. NA0 .AND. NA2 .NE. 0 ) THEN
C           PASSAGE A L'ARETE SUIVANTE DU CF
            NA1 = NA2
            GOTO 40
         ENDIF

      ENDDO

C     TRACE DES TETRAEDRES, TRIANGLES, ARETES DU CF et POINTS D'INTERSECTION
      IF( NBARPB .GT. 0 ) THEN
         tracte = .true.
         KTITRE='ajte1cmi:        TETRAEDRES avec NBARPB=      '
         WRITE(KTITRE(11:15),'(I5)') NBTECF
         WRITE(KTITRE(41:45),'(I5)') NBARPB
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRCFFAPE( KTITRE, PTXYZD, MXFACO, LEFACO, NO0FAR, NOTETR,
     %                  NBTECF, NOTECF, NBTRCF, NOTRCF,
     %                  MXARCF, NBCF,   N1ARCF, NOARCF,
     %                  NBSTIS, NOSTIS, NBPTIN, PTINTERS )
      ENDIF

 9999 TRACTE = TRACTE0
      RETURN
      END
