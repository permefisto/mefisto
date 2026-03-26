      SUBROUTINE NDPGEL( NTLXOB,
     %                   NTTOPO, MNTOPO,
     %                   NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %                   NBTYEL, NTNPEF, MNNPEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE NUMERO ET ADRESSE DES TMS
C -----    XYZPOINT XYZNOEUD NPEF D'UNE TOPOLOGIE D'UN OBJET
C
C ENTREES :
C ---------
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET
C
C SORTIES :
C ---------
C NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C NTXYZP : NUMERO      DU TMS 'XYZPOINT'  DE L'OBJET
C MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'  DE L'OBJET
C NTXYZN : NUMERO      DU TMS 'XYZNOEUD'  DE L'OBJET
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'  DE L'OBJET
C NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE
C NTNPEF : TABLEAU DU NUMERO        DU TMS 'NPEF' DES NBTYEL TYPES D'EF
C MNNPEF : TABLEAU DE L ADRESSE MCN DU TMS 'NPEF' DES NBTYEL TYPES D'EF
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*4       NOM4,CHARX
      INTEGER           NTNPEF(1:*),MNNPEF(1:*)
C
C     LE TABLEAU TOPOLOGIE
      IERR = 0
      CALL LXTSOU( NTLXOB, 'TOPOLOGIE', NTTOPO, MNTOPO )
      IF( NTTOPO .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET SANS INTERPOLATION'
            KERR(2) = 'RE-EXECUTER MAILLER POUR LE FAIRE'
         ELSE
            KERR(1) = 'ERROR: OBJECT WITHOUT INTERPOLATION'
            KERR(2) = 'EXECUTE AGAIN MAILLER to DO THAT'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE CODE DE TRAITEMENT DES NOEUDS POINTS SOMMETS
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     RECUPERATION DU TABLEAU XYZSOMMET
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZS, MNXYZS )
C
C     RECUPERATION DU TABLEAU XYZNOEUD
      CALL LXTSOU( NTLXOB, 'XYZNOEUD' , NTXYZN, MNXYZN )
C
C     RECUPERATION DU TABLEAU XYZPOINT
      CALL LXTSOU( NTLXOB, 'XYZPOINT' , NTXYZP, MNXYZP )
C
C     AFFECTATION SELON LE TYPE
C     NDPGST : CODE GENERAL DE TRAITEMENT DES
C              NOEUDS POINTS GEOMETRIQUES ET SOMMETS
C               0 : NOEUDS=POINTS=SOMMETS
C                   LE  TABLEAU  SOMMETS EXISTE
C               1 : NOEUDS=POINTS#SOMMETS
C                   LES TABLEAUX NOEUDS  ET SOMMETS EXISTENT
C               2 : NOEUDS#POINTS=SOMMETS
C                   LES TABLEAUX NOEUDS  ET SOMMETS EXISTENT
C               3 : NOEUDS#POINTS#SOMMETS
C                   LES TABLEAUX NOEUDS POINTS SOMMETS EXISTENT
C
      IF( NTXYZP .GT. 0 .AND. NTXYZN .LE. 0 ) THEN
         NTXYZN = NTXYZP
         MNXYZN = MNXYZP
      ENDIF
      IF( NTXYZN .LE. 0 ) THEN
         NTXYZN = NTXYZS
         MNXYZN = MNXYZS
      ENDIF
      IF( NTXYZP .LE. 0 ) THEN
         IF( NDPGST .EQ. 0 .OR. NDPGST .EQ. 2 ) THEN
            NTXYZP = NTXYZS
            MNXYZP = MNXYZS
         ELSE
            NTXYZP = NTXYZN
            MNXYZP = MNXYZN
         ENDIF
      ENDIF
C
C     VERIFICATION
      IF( NTXYZN .LE. 0 .OR. NTXYZP .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR NDPGEL: TABLEAU POINTS OU NOEUDS ABSENT'
         ELSE
            KERR(1) = 'ERROR NDPGEL: TMS POINTS or NODES UNKNOWN'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LES TYPES D'ELEMENTS FINIS SONT RETROUVES
      NBTYEL = MCN( MNTOPO + WBTYEL )
      IF( NBTYEL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBTYEL
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR NDPGEL:' // KERR(MXLGER)(1:4)
     %             //' TYPES d''ELEMENTS FINIS'
         ELSE
            KERR(1) = 'ERROR NDPGEL:' // KERR(MXLGER)(1:4)
     %             //' TYPES of FINITE ELEMENTS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     BOUCLE SUR LES TYPES D'EF
      DO  10 I=1,NBTYEL
C
C        LE NOM DU TABLEAU NPEF"
         NOM4 = CHARX( MCN(MNTOPO+WMTYEL-1+I) )
         CALL LXTSOU( NTLXOB, 'NPEF"'//NOM4,  NTNPEF(I), MNNPEF(I) )
         IF( NTNPEF(I) .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR NDPGEL: NPEF"'//NOM4//' ABSENT'
            ELSE
               KERR(1) = 'ERROR NDPGEL: NPEF"'//NOM4//' UNKNOWN'
            ENDIF
            CALL LEREUR
            IERR = 4
            RETURN
         ENDIF
C
cccC        NUMERO DU TYPE DE L'EF I
ccc         NUTYEL = MCN( MNNPEF(I) + WUTYEL )
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            WRITE(IMPRIM,10010) I,NOM4,NUTYEL
ccc         ELSE
ccc            WRITE(IMPRIM,20010) I,NOM4,NUTYEL
ccc         ENDIF
C
 10   CONTINUE
C
ccc10010 FORMAT(' Type EF',I1,': Nom ', A4,' No=',I2)
ccc20010 FORMAT(' FE Type',I1,': Name ',A4,' No=',I2)
      RETURN
      END
