       SUBROUTINE AJOUTPT( NX, NY, NDIM, MNXYZS, MNNSEF,
     %                     MOXYZS, MONSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER UN POINT DANS LE MAILLAGE A L'INTERIEUR
C -----    D'UN TRIANGLE 2D OU 3D
C ENTREES:
C --------
C NX     : ABSCISSE DU POINT A RAJOUTER
C NY     : ORDONNEE DU POINT A RAJOUTER
C NDIM   : DIMENSION (2 ou 3) DE L'ESPACE DE LA SURFACE
C MNXYZS : ADRESSE DU TABLEAU XYZSOMMET
C MNNSEF : ADRESSE DU TABLEAU NSEF
C MOXYZS : NOMBRE DE MOTS DECLARES DU TMS XYZSOMMET
C MONSEF : NOMBRE DE MOTS DECLARES DU TMS NSEF
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR, =1 SI POINT CLIQUE EXTERIEUR AU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      REAL              XYZPT(3)
      INTEGER           MOXYZS, MONSEF
      INTEGER           NOSOEF(1:4)
      INTEGER           NOUSOM
C
      IERR = 0
C     ADRESSE MCN DU TABLEAU NUSOEF NO DES 4 SOMMETS DES NBTQ EF ACTUELS
      MNSOEL = MNNSEF + WUSOEF

C     NBTQ LE NOMBRE D'EF
      NBTQ   = MCN( MNNSEF + WBEFOB )

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZS+WNBSOM)
C
      IF( NDIM .LE. 2 ) THEN
C
C        COORDONNEES OBJET DU POINT CLIQUE
         XYZPT(1)=XOB2PX(NX)
         XYZPT(2)=YOB2PX(NY)
         XYZPT(3)=0
C
C        POINT CLIQUE DANS UN TRIANGLE/QUADRANGLE 2D ?
         CALL TQPTCLIC( NX,NY,NDIM,RMCN(MNXYZS+WYZSOM),NBTQ,MCN(MNSOEL),
     %                  NUMTQ )
C
      ELSE
C
C        NUMERO DU SOMMET NUMST VISIBLE LE PLUS PROCHE DU POINT CLIQUE
         CALL NOSTCLIC( NX, NY, NDIM, NBSOM, RMCN(MNXYZS+WYZSOM),
     %                  NBTQ,  MCN(MNSOEL),
     %                  NUMTQ, XYZPT, NUMST )
C
      ENDIF
C
      IF( NUMTQ .EQ. 0 ) THEN
         NBLGRC(NRERR)=1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='POINT CLIQUE HORS MAILLAGE'
         ELSE
            KERR(1)='CLICKED POINT OUTSIDE THE MESH'
         ENDIF
         CALL LEREUR
         IERR=1
         RETURN
      ENDIF
C
C     MISE A JOUR DU TABLEAU XYZSOMMET
C     ================================
      NBSOM = MCN(MNXYZS+WNBSOM)
      L     = WYZSOM + 3 * NBSOM
C     ON AJUSTE LA TAILLE DU TABLEAU POUR UN SOMMET DE PLUS
      IF ( MOXYZS .LT. L+3 ) THEN
         CALL TNMCAU( 'ENTIER', MOXYZS , L+3 , MOXYZS , MNXYZS)
         MOXYZS = L+3
         NBSOM  = MCN(MNXYZS+WNBSOM)
      ENDIF
C     ON RENTRE LE NOUVEAU POINT
      RMCN(MNXYZS+L  )=XYZPT(1)
      RMCN(MNXYZS+L+1)=XYZPT(2)
      RMCN(MNXYZS+L+2)=XYZPT(3)
      MCN(MNXYZS+WNBSOM)=NBSOM+1
C     LE DERNIER SOMMET DU MAILLAGE EST LE NOUVEAU SOMMET
      NOUSOM = MCN(MNXYZS+WNBSOM)
      CALL ECDATE( MCN(MNXYZS) )
      MCN(MNXYZS+MOTVAR(6))=NONMTD('~>>>XYZSOMMET')
C
C     MISE A JOUR DU TABLEAU NSEF
C     ===========================
      NBTQ1 = MCN(MNNSEF+WBEFOB)
      L     = WUSOEF + 4 * NBTQ1
      IF ( MONSEF .LT. L+12 ) THEN
         CALL TNMCAU( 'ENTIER', L, L+12, L, MNNSEF )
         MONSEF=L+12
      ENDIF
C
C     NUMERO DES SOMMETS DU TRIANGLE OU QUADRANGLE NUMTQ
      MN  = MNNSEF + WUSOEF - 4
      MNS = MN + 4 * NUMTQ
      NOSOEF(1)=MCN(MNS )
      NOSOEF(2)=MCN(MNS+1)
      NOSOEF(3)=MCN(MNS+2)
      NOSOEF(4)=MCN(MNS+3)
C
      IF( NOSOEF(4).EQ. 0 ) THEN
C
C        POINT CLIQUE DANS UN TRIANGLE NBTQ1 => 3 SOUS-TRIANGLES
         MNS = MN + 4 * NBTQ1 + 4
         MCN(MNS  )=NOSOEF(2)
         MCN(MNS+1)=NOSOEF(3)
         MCN(MNS+2)=NOUSOM
         MCN(MNS+3)=0
C
C        TROISIEME TRIANGLE NBTQ1+1 A MODIFIER
         MNS = MN + 4 * NBTQ1 + 8
         MCN(MNS  )=NOSOEF(3)
         MCN(MNS+1)=NOSOEF(1)
         MCN(MNS+2)=NOUSOM
         MCN(MNS+3)=0
C
C        3-EME SOMMET DU TRIANGLE INITIAL A MODIFIER
         MCN(MN+4*NUMTQ+2)=NOUSOM
C
C        NOUVEAU NOMBRE DE TRIANGLES QUADRANGLES
         MCN(MNNSEF+WBEFOB)=NBTQ1+2
C
      ELSE
C
C        POINT CLIQUE DANS UN QUADRANGLE => 4 SOUS-TRIANGLES
         MNS = MN + 4 * NBTQ1 + 4
         MCN(MNS  )=NOSOEF(1)
         MCN(MNS+1)=NOSOEF(2)
         MCN(MNS+2)=NOUSOM
         MCN(MNS+3)=0
C
         MNS = MN + 4 * NBTQ1 + 8
         MCN(MNS  )=NOSOEF(2)
         MCN(MNS+1)=NOSOEF(3)
         MCN(MNS+2)=NOUSOM
         MCN(MNS+3)=0
C
         MNS = MN + 4 * NBTQ1 + 12
         MCN(MNS  )=NOSOEF(3)
         MCN(MNS+1)=NOSOEF(4)
         MCN(MNS+2)=NOUSOM
         MCN(MNS+3)=0
C
         MNS = MN + 4 * NUMTQ
         MCN(MNS  )=NOSOEF(4)
         MCN(MNS+1)=NOSOEF(1)
         MCN(MNS+2)=NOUSOM
         MCN(MNS+3)=0
C
C        NOUVEAU NOMBRE DE TRIANGLES QUADRANGLES
         MCN(MNNSEF+WBEFOB)=NBTQ1+3
C
      ENDIF
C
C     MISE A JOUR FINALE DU TMS NSEF
      MCN(MNNSEF+WUTYOB)=3
      MCN(MNNSEF+WUTYMA)=0
      MCN(MNNSEF+WUTFMA)=0
      MCN(MNNSEF+WBSOEF)=4
      MCN(MNNSEF+WBTGEF)=0
      MCN(MNNSEF+WBEFTG)=0
      MCN(MNNSEF+WBEFAP)=0
C
      CALL ECDATE( MCN(MNNSEF) )
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )
C
      RETURN
      END
