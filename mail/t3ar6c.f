      SUBROUTINE T3AR6C( NM6CUB )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES DES 6-CUBES DU MAILLAGE DE L'OBJET NM6CUB
C -----    EN PROJECTION SELON X4 X5 X6 I.E.
C          1 SOMMET = TOUS LES SOMMETS AYANT MEME X1 X2 X3 QUAND
C                     X4 X5 X6 VARIENT
C ENTREES:
C --------
C NM6CUB : NOM DE L'OBJET MAILLE DE 6-CUBES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A M UNIVERSITY & LJLL UPMC Octobre 2005
C MODIFS : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      REAL              XYZ(3)
      CHARACTER*(*)     NM6CUB
      CHARACTER*24      NMSOMM
      include"./incl/nusc1c6.inc"
C
      MNS1 = 0
      MNS2 = 0
C
C     L'OBJET MAILLE EN 6-CUBES
C     =========================
C     LE TABLEAU LEXIQUE DE CE 6-CUBES
      CALL LXLXOU( NTOBJE, NM6CUB, NTLXOB, MN )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = '6-CUBES INCONNU :' // NM6CUB
         CALL LEREUR
         RETURN
      ENDIF
C
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ADRESSAGE DES ADRESSES DES TABLEAUX ELEMENTS DE CET OBJET
      MNELEM = 0
      MXTYEL = 8
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
      MNTELE = MNELEM + MXTYEL
      CALL NDPGEL( NTLXOB,
     %             NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             NBTYEL, MCN(MNTELE), MCN(MNELEM), IERR )
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTXYZP : NUMERO      DU TMS 'XYZPOINT'  DE L'OBJET
C     MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'  DE L'OBJET
C     NTXYZN : NUMERO      DU TMS 'XYZNOEUD'  DE L'OBJET
C     MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'  DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS FINIS DU MAILLAGE
C     MNTELE : NUMERO      DU TMS DES NBTYEL NPEF"...
C     MNELEM : ADRESSE MCN DU TMS DES NBTYEL NPEF"...
C     NDPGST=0   CODE TRAITEMENT DES XYZ DES SOMMETS POINTS NOEUDS DU MAILLAGE
C            0 : NOEUDS=POINTS=SOMMETS
C            1 : NOEUDS=POINTS#SOMMETS
C            2 : NOEUDS#POINTS=SOMMETS
C            3 : NOEUDS#POINTS#SOMMETS
      IF( IERR .NE. 0 ) RETURN
C
C     DIMENSION DE L'ESPACE ou NOMBRE DE COORDONNEES ICI=6
      NBCOOR = MCN( MNXYZP + WBCOOP )
      IF( NBCOOR .NE. 6 ) RETURN
C
C     PARCOURS DU TMS NPEF"6Q1C (NUTYEL=2 ACTUELLEMENT)
      MNST = MNXYZP + WYZPOI - NBCOOR
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF"6Q1C
      MNELE = MCN( MNELEM )
C
C     LE NOMBRE DE TELS ELEMENTS FINIS 6Q1C
      NBELEM = MCN( MNELE + WBELEM )
C
C     LA BOUCLE SUR LES 192 ARETES DES EF 6-CUBES DU MAILLAGE
C     -------------------------------------------------------
      MNEF = MNELE + WUNDEL - 1
      DO 100 N=1,NBELEM
C
C        LE NO DES 64 SOMMETS DE L'EF 6-CUBE N
         MN = MNEF + N - NBELEM
C
         DO 40 K = 1, NB1C6C
C
C           LE NO DU SOMMET INITIAL DE L'ARETE K
            NS1  = MCN( MN + NUS1C6C(1,K) * NBELEM )
            MNS1 = MNST + NBCOOR * NS1
C           LE NO DU SOMMET FINAL DE L'ARETE K
            NS2  = MCN( MN + NUS1C6C(2,K) * NBELEM )
            MNS2 = MNST + NBCOOR * NS2
C           EF SANS TG => TRACE DES ARETES P1
            CALL TRAIT6D( NCOUAF, RMCN(MNS1), RMCN(MNS2) )
C
C           TRACE EVENTUEL DU NO DES SOMMETS=NOEUDS=POINTS
            IF( IAVNSO .NE. 0 ) THEN
               WRITE( NMSOMM , '(I8)' ) NS1
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, RMCN(MNS1), NMSOMM(1:L) )
               WRITE( NMSOMM , '(I8)' ) NS2
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, RMCN(MNS2), NMSOMM(1:L) )
            ENDIF
C
 40      CONTINUE
C
 100  CONTINUE
C
C     LE TRACE DE LA POIGNEE DU 6-CUBES AU MILIEU DE LA DERNIERE ARETE TRACEE
      XYZ(1) = ( RMCN(MNS1  ) + RMCN(MNS2  ) ) * 0.5
      XYZ(2) = ( RMCN(MNS1+1) + RMCN(MNS2+1) ) * 0.5
      XYZ(3) = ( RMCN(MNS1+2) + RMCN(MNS2+2) ) * 0.5
      CALL NUOBNM( 'OBJET', NM6CUB, K )
      CALL ITEMV3(  XYZ,    NM6CUB, K )
      CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNELEM )
      END
