      SUBROUTINE NBEFNOEUD( MNTOPO, MNNPEF,  NBNOE,
     %                      NBEFNO, MXNOE1EF, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LA LISTE DU NOMBRE D'ELEMENTS FINIS DE CHAQUE NOEUD
C -----
C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU  TABLEAU  TOPOLOGIE DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"     DE L'OBJET
C NBNOE  : NOMBRE DE NOEUDS DE L'OBJET

C SORTIES:
C --------
C NBEFNO  : NBEFNO(N)=NOMBRE d'EF DU NOEUD N
C MXNOE1EF: NOMBRE MAXIMUM DE NOEUD D'UN EF DU MAILLAGE
C IERR    : =0 SI PAS D'ERREUR
C           =1 TYPE D'EF NON RETROUVE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN  PERRONNET Saint PIERRE du PERRAY             Mars 2021
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           MNNPEF(1:*), NBEFNO(NBNOE)

      IERR = 0
      MXNOE1EF = 0

C     MISE A ZERO DU TABLEAU NBEFNO NOMBRE DES EF CONTENANT CHAQUE NOEUD
      CALL AZEROI( NBNOE, NBEFNO )

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
C        NOMBRE DE NOEUD DE L'EF DE TYPE NOTYEL
         NBNDEL   = MCN( MNELE + WBNDEL )
         MXNOE1EF = MAX( MXNOE1EF, NBNDEL )

C        NOMBRE D'EF DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
         MNNDEL = MNELE + WUNDEL - 1

         DO NUELEM=1,NBELEM

            MN = MNNDEL + NUELEM

            DO I = 1, NBNDEL
C              NUMERO DE NOEUD DE L'EF NUELEM
               NOEUD = MCN( MN )
C              LE NOMBRE D'EF CONTENANT LE NOEUD
               NBEFNO( NOEUD ) = NBEFNO( NOEUD ) + 1
               MN = MN + NBELEM
            ENDDO

         ENDDO

      ENDDO

 9999 RETURN
      END
