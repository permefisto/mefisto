      SUBROUTINE LISTNOVO( MNTOPO,  MNNPEF,   NBNOE,  MXLISTV,
     %                     LPVOIS,  LISTVOI,  MXVOIS, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   GENERER LA LISTE DES VOISINS DE CHACUN DES NOEUDS DU MAILLAGE
C ----- ( LE NOEUD N'EST PAS PRIS EN COMPTE DANS SA LISTE )
C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU  TABLEAU  TOPOLOGIE DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"     DE L'OBJET
C NBNOE  : NOMBRE DE NOEUDS DE L'OBJET
C MAXVOI : NOMBRE MAXIMAL DE VOISINS D'UN NOEUD
C MXLISTV: NOMBRE D'ENTIERS DU TABLEAU LISTVOI MXLISTV = NBNOE * MAXVOI
C LPVOIS : POINTEUR SUR LE DERNIER NO DE VOISIN DE CHAQUE NOEUD
C          LPVOIS(0)=0
C SORTIES:
C --------
C LPVOIS : POINTEUR SUR LE DERNIER NUMERO DES VOISINS
C          DE CHACUN DES NOEUDS
C          LPVOIS(I)=INDICE DANS LISTVOI DU DERNIER NO
C                    DE VOISIN DU NOEUD I
C          LPVOIS(0)=0
C LISTVOI: TABLEAU LISTE DES VOISINS DES NOEUDS
C          ATTENTION: LE NOEUD NE FAIT PAS PARTIE DES VOISINS
C MXVOIS : NOMBRE MAXIMAL DE VOISINS D'UN NOEUD DU MAILLAGE
C IERR   : =0 SI PAS D'ERREUR
C          =1 TYPE D'EF NON RETROUVE
C          =5 REALLOCATION DE LISTVOI DEMANDEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE UPMC PARIS    Octobre 1989
C MODIFS : ALAIN  PERRONNET Saint PIERRE du PERRAY             Mars 2021
C MODIFS : ALAIN  PERRONNET Saint PIERRE du PERRAY            Avril 2023
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)

      INTEGER  MNNPEF(1:*), LISTVOI(MXLISTV), LPVOIS(0:NBNOE),
     %         NONOEF(MXNOEL)

      PRINT*
      PRINT*,'listnovo: CONSTRUCTION de la LISTE des NOEUDS VOISINS de C
     %HAQUE NOEUD'

      IERR   = 0
      MXVOIS = 0

C     MISE A ZERO DU TABLEAU LISTVOI DES NOEUDS VOISINS
      CALL AZEROI( MXLISTV, LISTVOI )

C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS DU MAILLAGE
C     ====================================================
      NBTYEL = MCN( MNTOPO + WBTYEL )
      DO NOTYEL=1,NBTYEL

C        MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MNNPEF( NOTYEL )
         IF( MNELE .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: TYPE INCONNU D''ELEMENT FINI'
            ELSE
               KERR(1) = 'ERROR: UNKNOW TYPE OF FINITE ELEMENT'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF

C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NOTYEL
C        ==================================================
         NBELEM = MCN( MNELE + WBELEM )
         NBNDEL = MCN( MNELE + WBNDEL )
         MNNDEL = MNELE + WUNDEL - 1

         DO NUELEM=1,NBELEM

C           LES NUMEROS DES NOEUDS DE L'ELEMENT FINI
C           ----------------------------------------
            MN = MNNDEL + NUELEM
            DO I=1,NBNDEL
               NONOEF(I) = MCN(MN)
               MN = MN + NBELEM
            ENDDO

C           RECHERCHE DES NOEUDS DANS LA LISTE DES VOISINS
C           ----------------------------------------------
            DO I=1,NBNDEL

C              LE NO DU NOEUD I DE L'EF
               NOEI = NONOEF(I)

               DO 10 J=1,NBNDEL

                  IF ( I .EQ. J ) GOTO 10
                  NOEJ = NONOEF(J)

                  DO K = LPVOIS( NOEI-1 )+1, LPVOIS( NOEI )
                     NOEK = LISTVOI(K)
                     IF( NOEK .EQ. 0 ) THEN
C                       ON AJOUTE CE NOEUD NOEJ A LA LISTE
                        LISTVOI(K) = NOEJ
                        GOTO 10
                     ENDIF
                     IF( NOEK .EQ. NOEJ ) THEN
C                       LE NOEUD NOEJ EST DEJA DANS LA LISTE
                        GOTO 10
                     ENDIF
                  ENDDO

C                 IL N'Y A PLUS DE PLACE DANS LISTVOI et
C                 NOEJ N'EST PAS DANS LA LISTE des VOISINS de NOEI!
C                 .................................................
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='listnovo: LISTVOI DECLARE TROP PETIT'
                     KERR(2)='AUGMENTATION DEMANDEE et REALLOCATION'
                  ELSE
                    KERR(1)='listnovo: LISTVOI DECLARATION TOO SMALL'
                    KERR(2)='DEMANDED AUGMENTATION and REALLOCATION'
                  ENDIF
                  CALL LERESU
                  PRINT*,'listnovo: NO VOISINS de ',NOEI,' de',
     %                    LPVOIS( NOEI-1 )+1, ' a',LPVOIS( NOEI ),
     %                    ' INSUFFISANT'
C                 REALLOCATION DE LISTVOI DEMANDEE
                  IERR = 5
                  GOTO 9999

 10            ENDDO

            ENDDO

         ENDDO

      ENDDO

ccc      DO NOE=1,NBNOE
ccc         PRINT*,'listnovo: LISTVOI du NOEUD',NOE,':',
ccc     %           (LISTVOI(K),K=LPVOIS(NOE-1)+1,LPVOIS(NOE))
ccc      ENDDO

C     CLASSEMENT DES NOEUDS VOISINS PAR ORDRE CROISSANT
C     COMPRESSION DU TABLEAU LISTVOI, CREATION DU POINTEUR LPVOIS
C     ===========================================================
      LPVOIS(0) = 0
      LPK0      = 0
      LV1       = 0

      DO NOE=1,NBNOE

C        CALCUL DU NOMBRE DE VOISINS DU NOEUD NOE
         NBV = 0
         DO K1 = 1, LPVOIS(NOE) - LPK0
C           NO DU NOEUD VOISIN K1 DU NOEUD NOE
            NOVOK1 = LISTVOI( LPK0 + K1 )
            IF( NOVOK1 .EQ. 0 ) GOTO 20
C           NBV NOMBRE DE VOISINS DU NOEUD NOE EST AUGMENTE
            NBV = NBV + 1
         ENDDO

C       TRI CROISSANT DES NUMEROS DES NOEUDS VOISINS DE NOE
 20     CALL TRIENT( NBV, LISTVOI(LPK0+1) )

C       COMPRESSION DE LA LIGNE NOE DU TABLEAU LISTVOI SUR LUI-MEME
        DO K1=1,NBV
            LV1 = LV1 + 1
            LISTVOI( LV1 ) = LISTVOI( LPK0 + K1 )
         ENDDO

C        PROTECTION DE LPVOI(NOE-1) PROCHAIN
         LPK0 = LPVOIS( NOE )

C        POINTEUR SUR LE DERNIER NO DE VOISIN DU NOEUD NOE
         LPVOIS( NOE ) = LV1

C        NBV NOMBRE DE VOISINS DU NOEUD NOE
         IF( NBV .GT. MXVOIS ) MXVOIS=NBV

C        VERIFICATION A SUPPRIMER APRES MISE AU POINT
         DO K = LPVOIS(NOE-1)+1, LPVOIS(NOE)-1
            IF( LISTVOI( K ) .GE. LISTVOI( K+1 ) ) THEN
C             LES NO DES NOEUDS VOISINS NE SONT PAS CROISSANTS
              PRINT*,'listnovo: ANOMALIE ',LISTVOI(K),' >=',LISTVOI(K+1)
              PRINT*,'LISTVOI(',NOE,'):',
     %               (LISTVOI(K1),K1=LPVOIS(NOE-1)+1, LPVOIS(NOE))
              PRINT*,'listnovo.f A CORRIGER...'
            ENDIF
         ENDDO

      ENDDO

ccc      PRINT*
ccc      DO NOE=1,NBNOE
ccc         PRINT*,'listnovo: LISTVOI du NOEUD',NOE,':',
ccc     %           (LISTVOI(K),K=LPVOIS(NOE-1)+1,LPVOIS(NOE))
ccc      ENDDO

C     GESTION DES TABLEAUX
C     --------------------
      NBV = LV1 / NBNOE
      IF( LANGAG .EQ. 0 ) THEN
        PRINT*, 'FIN listnovo AVEC'
        PRINT*, 'NOMBRE MAXIMAL DE VOISINS D''UN NOEUD=',MXVOIS
        PRINT*, 'NOMBRE MOYEN   DE VOISINS D''UN NOEUD=',NBV
        PRINT*, 'NOMBRE TOTAL   DE VOISINS DES NOEUDS=',LV1
      ELSE
        PRINT*, 'END listnovo WITH'
        PRINT*, 'NEIGHBORS NUMBER MAXIMUM of a NODE=',MXVOIS
        PRINT*, 'NEIGHBORS NUMBER AVERAGE of a NODE=',NBV
        PRINT*, 'NEIGHBORS TOTAL  NUMBER  of  NODES=',LV1
      ENDIF

C     FIN NORMALE ou PAS ASSEZ DE PLACE EN MEMOIRE
C     --------------------------------------------
 9999 RETURN
      END
