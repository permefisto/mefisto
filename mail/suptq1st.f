       SUBROUTINE SUPTQ1ST( MNXYZ1, MNNSEF1,  NUMST, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER UN SOMMET CLIQUE DU MAILLAGE D'UNE SURFACE 3D en
C -----    SUPPRIMANT TOUS LES TRIANGLES QUADRANGLES DONT IL EST SOMMET

C ENTREES:
C --------
C MNXYZ1 : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF1: ADRESSE DU TABLEAU NSEF DE LA SURFACE

C SORTIES:
C --------
C NUMST  : >0 NUMERO DU SOMMET SUPPRIME
C          =0 SI PAS DE SOMMET CLIQUE
C IERR   : =0 SI PAS D'ERREUR
C           1 PAS DE POINT CLIQUE RETROUVE
C           2 PRESENCE DE QUADRANGLE
C           9 SATURATION DU SUPER-TABLEAU MCN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C MODIFS : PERRONNET ALAIN Saint PIERRE du PERRAY               Mai 2020
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))

      IERR   = 0
      MNNTAB = 0
      MNSUTQ = 0

C     NBSOM NOMBRE DE SOMMETS
      NBSOM  = MCN(MNXYZ1+WNBSOM)

C     NBEFOB NOMBRE DES EF
      NBEFOB = MCN(MNNSEF1+WBEFOB)

C     NBSOEF NOMBRE DE SOMMETS PAR EF
      NBSOEF = MCN(MNNSEF1+WBSOEF)

C     ADRESSE MCN DE DEBUT DES NO DE SOMMETS DES EF
      MNSOEL = MNNSEF1 + WUSOEF

C     SELECTIONNER A L'AIDE D'UN CLIC SOURIS UN SOMMET DU MAILLAGE
C     PRESENT DANS LA FENETRE GRAPHIQUE C-A-D
C     RETOURNER LE NUMERO NUMST0 DU SOMMET VISIBLE LE PLUS PROCHE
C     DU POINT CLIQUE EN NX, NY DANS UN MAILLAGE 3D
C     ------------------------------------------------------------
      CALL SESTCLIC( RMCN(MNXYZ1+WYZSOM), ITST, NUMST )
      IF( NUMST .LE. 0 ) GOTO 9990

C     RECHERCHE DES QUADRANGLES OU TRIANGLES CONTENANT LE SOMMET NUMST
C     ----------------------------------------------------------------
      CALL TNMCDC( 'ENTIER', NBSOM, MNNTAB )
      IF( MNNTAB .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'suptq1st: PAS ASSEZ de PLACE en MCN'
         ELSE
            KERR(1) = 'suptq1st: NOT ENOUGH MCN MEMORY'
         ENDIF
         CALL LEREUR
         IERR = 9
         RETURN
      ENDIF

C     NOMBRE D'EF TQ A SUPPRIMER
      NBTQAS=0
      MNEF = MNSOEL - 4
      DO NEF=1,NBEFOB
         MNEF = MNEF + 4
         IF( (MCN(MNEF)  .EQ.NUMST).OR.(MCN(MNEF+1).EQ.NUMST)
     %   .OR.(MCN(MNEF+2).EQ.NUMST).OR.(MCN(MNEF+3).EQ.NUMST)) THEN
C           NO NEF DU TQ CONTENANT LE SOMMET NUMST
            MCN( MNNTAB+NBTQAS ) = NEF
            NBTQAS = NBTQAS + 1
         ENDIF
      ENDDO
      IF( NBTQAS .LE. 0 ) GOTO 9900


C     LES NBTQAS QUADRANGLES ou TRIANGLES SONT SUPPRIMES
C     --------------------------------------------------
      CALL SUPEFMAR( NBTQAS, MCN(MNNTAB), NBSOEF, NBEFOB, MCN(MNSOEL) )

 9900 IF( MNNTAB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOM,  MNNTAB )

      IF(LANGAG .EQ. 0 ) THEN
         PRINT*,'suptq1st: le SOMMET',NUMST,' et ses EF sont SUPPRIMES'
      ELSE
         PRINT*,'suptq1st: the VERTEX',NUMST,' and its FE are DELETED'
      ENDIF
      GOTO 9999


C     PAS DE POINT CLIQUE et PAS DE SOMMET SUPPRIME
C     ---------------------------------------------
 9990 NUMST = 0
      IERR  = 1
      IF(LANGAG .EQ. 0 ) THEN
         PRINT*,'suptq1st: PAS de SOMMET SUPPRIME'
      ELSE
         PRINT*,'suptq1st: NO VERTEX DELETED'
      ENDIF

 9999 RETURN
      END
