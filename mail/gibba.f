      SUBROUTINE GIBBA( MNTOPO, MNNPEF, MNRENU,  IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RENUMEROTER LES NOEUDS AVEC LA NOUVELLE NUMEROTATION
C -----
C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU  TABLEAU TOPOLOGIE  DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"XXXX DE L'OBJET
C MNRENU : ADRESSE  MCN DU  TABLEAU DE LA NOUVELLE NUMEROTATION

C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : DEFAIX THIERRY  ANALYSE NUMERIQUE  UPMC  PARIS  JANVIER  1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      INTEGER         MNNPEF(1:*)

C     INITIALISATIONS
C     ===============
      IERR = 0
      MNRENU1 = MNRENU - 1

C     La BOUCLE sur les TYPES D'ELEMENTS FINIS du MAILLAGE
C     ====================================================
C     NBTYEL nb de types d'npef finis
      NBTYEL = MCN( MNTOPO + WBTYEL )
      DO NUTYEL=1,NBTYEL

C        MNELE adresse du tableau NPEF des elements finis de ce TYPE
         MNELE = MNNPEF( NUTYEL )
         IF( MNELE .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'gibba: les TMS NPEF sont INCONNUS'
            CALL LEREUR
            IERR = 1
            RETURN
         ENDIF

C        La BOUCLE sur les ELEMENTS FINIS de ce TYPE
C        -------------------------------------------
C        NBELEM nb d'ef de ce type
         NBELEM = MCN( MNELE + WBELEM )
C        NBNDEL nb de noeuds d'un tel element fini
         NBNDEL = MCN( MNELE + WBNDEL )
C        MNNDEL adresse du tableau des noeuds des EF de ce TYPE
         MNNDEL = MNELE + WUNDEL - 1

         DO NUELEM = 1, NBELEM
C           MN adresse du premier noeud de l'element fini NUELEM
            MN = MNNDEL + NUELEM
            DO I = 0, NBNDEL-1
C              NEW NS = RENU( OLD NS )
               MCN(MN) = MCN( MNRENU1+MCN(MN) )
               MN      = MN + NBELEM
            ENDDO
         ENDDO

      ENDDO

      RETURN
      END
