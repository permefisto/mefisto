      SUBROUTINE COLORE( NBELEM , NBCOLO , MNTOPO , MNELEM , MNNOEU ,
     %                   MNCOEL , MNLPCO , MNELCO , IERR)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : COLORIER LES ELEMENTS D'UN MAILLAGE EN REGROUPANT LES
C ----- ELEMENTS NON-VOISINS DANS UN MEME SOUS-ENSEMBLE
C
C ENTREES :
C ---------
C NBELEM : NOMBRE D'ELEMENTS
C NBCOLO : NOMBRE MAXIMUM DE COULEURS ACCEPTE
C MNTOPO : ADRESSE  MCN DU TABLEAU TOPOLOGIE  DE L'OBJET
C MNELEM : ADRESSES MCN DES TABLEAUX ELEMENTS DE L'OBJET
C MNNOEU : ADRESSE  MCN DU TABLEAU NOEUDS     DE L'OBJET
C
C SORTIES :
C ---------
C MNCOEL : ADRESSE MCN DU TABLEAU COULEUR D'UN ELEMENT
C MNLPCO : ADRESSE MCN DU POINTEUR SUR LA LISTE DES ELEMENTS
C          DE MEME COULEUR
C MNELCO : ADRESSE MCN DE LA LISTE DES ELEMENTS RANGES PAR COULEUR
C IERR   : CODE D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS  JUIN 1990
C23456---------------------------------------------------------------012
      IMPLICIT        INTEGER (W)
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      INTEGER     MNELEM(1:*)
C
      IERR = 0
C
C     ADRESSAGE DES TABLEAUX LPELCO ET ELMCOL
C     ---------------------------------------
      MNLPCO = 0
      CALL TNMCDC( 'ENTIER'  , NBCOLO+1 , MNLPCO )
      MNELCO = 0
      CALL TNMCDC( 'ENTIER'  , NBELEM   , MNELCO )
      MNCOEL = 0
      CALL TNMCDC( 'ENTIER'  , NBELEM   , MNCOEL )
      MNELLO = 0
      CALL TNMCDC( 'LOGIQUE' , NBELEM   , MNELLO )
      MNCLLO = 0
      CALL TNMCDC( 'LOGIQUE' , NBELEM   , MNCLLO )
      MNRNEL = 0
      CALL TNMCDC( 'ENTIER'  , NBELEM   , MNRNEL )
C     MISE A ZERO
      CALL AZEROI( NBELEM    , MCN(MNRNEL) )
C
C     FORMATION DE LA LISTE DES VOISINS
C     =================================
      NBEL = NBELEM
      CALL ELEVOI( MNTOPO , MNELEM , MNNOEU ,
     &             MNLPVO , MNLIVO , NBEL   , IERR )
C
C     FORMATION DES TABLEAUX LPELCO ET ELMCOL
C     =======================================
      CALL COLOR1( NBEL        , NBCOLO      ,
     &             MCN(MNLPVO) , MCN(MNLIVO) ,
     &             MCN(MNLPCO) , MCN(MNELCO) ,
     &             MCN(MNCOEL) , MCN(MNELLO) ,
     &             MCN(MNCLLO) , MCN(MNRNEL) , IERR )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
C     ------------------------------------
      LOLIVO = MCN(MNLPVO+NBEL)
      CALL TNMCDS( 'ENTIER'  , LOLIVO , MNLIVO )
      CALL TNMCDS( 'ENTIER'  , NBEL+1 , MNLPVO )
      CALL TNMCDS( 'LOGIQUE' , NBELEM , MNELLO )
      CALL TNMCDS( 'LOGIQUE' , NBELEM , MNCLLO )
      CALL TNMCDS( 'ENTIER'  , NBELEM , MNRNEL )
C
      RETURN
      END
