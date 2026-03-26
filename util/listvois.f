      SUBROUTINE LISTVOIS( MNTOPO,  MNNPEF,   NBNOE,  MAXVOI, MXLISTV,
     %                     LPVOIS,  LISTVOI,  MXVOIS, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   GENERER LA LISTE DES VOISINS DE CHACUN DES NOEUDS DU MAILLAGE
C ----- ( LE NOEUD N'EST PAS PRIS EN COMPTE DANS SA LISTE )
C
C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU  TABLEAU  TOPOLOGIE DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"     DE L'OBJET
C NBNOE  : NOMBRE DE NOEUDS DE L'OBJET
C MAXVOI : NOMBRE MAXIMAL DE VOISINS D'UN NOEUD
C MXLISTV: NOMBRE D'ENTIERS DU TABLEAU LISTVOI MXLISTV = NBNOE * MAXVOI

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
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C MODIFS : ALAIN  PERRONNET Saint PIERRE du PERRAY          JANVIER 2021
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           MNNPEF(1:*), LISTVOI(MXLISTV), LPVOIS(0:NBNOE),
     %                  NONOEF(MXNOEL)

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
            GOTO 9000
         ENDIF

C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NOTYEL
C        ==================================================
         NBELEM = MCN( MNELE + WBELEM )
         NBNDEL = MCN( MNELE + WBNDEL )
         MNNDEL = MNELE + WUNDEL - 1

         DO NUELEM=1,NBELEM
            MNNDE = MNNDEL + NUELEM

C           LES NUMEROS DES NOEUDS DE L'ELEMENT FINI
C           ----------------------------------------
            MN = MNNDE
            DO I=0,NBNDEL-1
               NONOEF(I+1) = MCN(MN)
               MN = MN + NBELEM
            ENDDO

C           RECHERCHE DES NOEUDS DANS LA LISTE DES VOISINS
C           ----------------------------------------------
            DO I=0,NBNDEL-1

C              LE NO DU NOEUD I+1 DE L'EF
               NUMI = NONOEF(I+1)
               IAVI = MAXVOI * ( NUMI - 1 ) - 1

               IF( IAVI .LE. 0 .OR. IAVI.GT.MXLISTV ) THEN
C                 DEBORDEMENT DU TABLEAU LISTVOI
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='LISTVOI(            ) REQUIS'
                     KERR(2)='ESSAI de REALLOCATION'
                  ELSE
                     KERR(1)='LISTVOI(            ) REQUESTED'
                     KERR(2)='TRY to REALLOCATE'
                  ENDIF
                  WRITE(KERR(1)(9:20),'(I12)') IAVI
                  CALL LERESU
C                 REALLOCATION DE LISTVOI DEMANDEE
                  IERR = 5
                  GOTO 9000
               ENDIF

               DO 15 J=0,NBNDEL-1
                  IF ( I .EQ. J ) GOTO 15
                  NUMJ = NONOEF(J+1)
                  DO K=1,MAXVOI
                     IAK = 1 + IAVI + K
                     KV  = LISTVOI(IAK)
                     IF (KV.EQ.0) THEN
C                        ON AJOUTE CE NOEUD A LA LISTE
                         LISTVOI(IAK) = NUMJ
                         GOTO 15
                     ELSE  IF (KV.EQ.NUMJ) THEN
C                        LE NOEUD EST DEJA DANS LA LISTE
                         GOTO 15
                     ENDIF
                  ENDDO

C                 IL N'Y A PLUS DE PLACE DANS LISTVOI et
C                 NUMJ N'EST PAS DANS LA LISTE !
C                 ......................................
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:15),'(I15)') MAXVOI
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='MAXVOI TROP PETIT '//KERR(MXLGER)(1:15)
                     KERR(2)='AUGMENTATION DEMANDEE et REALLOCATION'
                  ELSE
                    KERR(1)='MAXVOI TOO SMALL '//KERR(MXLGER)(1:15)
                    KERR(2)='DEMANDED AUGMENTATION and REALLOCATION'
                  ENDIF
                  CALL LERESU
C                 REALLOCATION DE LISTVOI DEMANDEE
                  IERR = 5
                  GOTO 9000
 15            ENDDO

            ENDDO

         ENDDO

      ENDDO

C     CLASSEMENT DES NOEUDS VOISINS PAR ORDRE CROISSANT
C     COMPRESSION DU TABLEAU LISTVOI, CREATION DU POINTEUR LPVOIS
C     ===========================================================
      LPVOIS(0)=0
      IC1 = 0
      DO I=1,NBNOE
         IA = ( I - 1 ) * MAXVOI

C        CLASSEMENT CROISSANT DES VOISINS
         NBV = 0
         DO 21 K1=1,MAXVOI
            IA1 = IA + K1
            IF( LISTVOI(IA1) .EQ. 0 ) GOTO 23
C           NBV NOMBRE DE VOISINS DU NOEUD I
            NBV = NBV + 1
            DO K2=K1+1,MAXVOI
               IA2 = IA + K2
               IF( LISTVOI(IA2) .EQ. 0 ) GOTO 21
               IF( LISTVOI(IA2) .LT. LISTVOI(IA1) ) THEN
                  M           =LISTVOI(IA1)
                  LISTVOI(IA1)=LISTVOI(IA2)
                  LISTVOI(IA2)=M
               ENDIF
            ENDDO
 21      ENDDO

C        COMPRESSION DE LA LIGNE I DU TABLEAU LISTVOI SUR LUI-MEME
 23      DO K1=1,NBV
            IC1 = IC1 + 1
            LISTVOI(IC1) = LISTVOI(IA+K1)
         ENDDO
         LPVOIS(I) = IC1

C        NBV NOMBRE DE VOISINS DU NOEUD I
         IF( NBV .GT. MXVOIS ) MXVOIS=NBV
      ENDDO

C     GESTION DES TABLEAUX
C     --------------------
      NBV = IC1 / NBNOE
      IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,*) 'FIN listvois AVEC'
        WRITE(IMPRIM,*) 'NOMBRE MAXIMAL DE VOISINS D''UN NOEUD=',MXVOIS
        WRITE(IMPRIM,*) 'NOMBRE MOYEN   DE VOISINS D''UN NOEUD=',NBV
        WRITE(IMPRIM,*) 'NOMBRE TOTAL   DE VOISINS DES NOEUDS=',IC1
      ELSE
        WRITE(IMPRIM,*) 'END listvois WITH'
        WRITE(IMPRIM,*) 'NEIGHBORS NUMBER MAXIMUM of a NODE=',MXVOIS
        WRITE(IMPRIM,*) 'NEIGHBORS NUMBER AVERAGE of a NODE=',NBV
        WRITE(IMPRIM,*) 'NEIGHBORS TOTAL  NUMBER  of  NODES=',IC1
      ENDIF

C     FIN NORMALE ou PAS ASSEZ DE PLACE EN MEMOIRE
C     --------------------------------------------
 9000 RETURN
      END
