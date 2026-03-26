       SUBROUTINE TUERTQ3D( MNXYZS, MNNSEF, NOSOTQ, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CLIQUER 1 POINT INTERNE A UN TRIANGLE OU UN QUADRANGLE 3D
C -----    VISIBLE ET DETRUIRE CE TRIANGLE OU QUADRANGLE
C ENTREES:
C --------
C MNXYZS : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF : ADRESSE DU TABLEAU NSEF      DE LA SURFACE
C
C SORTIES:
C --------
C NOSOTQ : NUMERO DES 4 SOMMETS DU TQ DETRUIT et NOSOTQ(1)=0 SINON
C IERR   : 0 SI PAS D'ERREUR et LE TQ DE NOSOTQ a ETE DETRUIT
C          1 SI POINT EN DEHORS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NOSOTQ(4)

      IERR = 0

C     DANS QUEL TRIANGLE ou QUADRANGLE?
C     =================================
      MNSOEL = MNNSEF + WUSOEF
      NBEFOB = MCN(MNNSEF+WBEFOB)
      MNXYZ  = MNXYZS + WYZSOM

C     SELECTIONNER A L'AIDE D'UN CLIC SOURIS UNE FACE DU MAILLAGE
C     PRESENTE DANS LA FENETRE GRAPHIQUE C-A-D
C     RETOURNER LE NUMERO NOFACL DE LA FACE VISIBLE DONT LE CLIC
C     EST INTERNE POUR UN MAILLAGE d'une SURFACE 2D ou 3D
C     -----------------------------------------------------------
      NBNOST = 0
      NOSOTQ(1) = 0
      CALL SEFACLIC( RMCN(MNXYZ), NCJAUN, NBNOST, NOSOTQ,
     %               MCN(MNSOEL), NOFACL )

      IF( NOFACL .EQ. 0 ) THEN
         NBLGRC(NRERR)=2
         IF( LANGAG .EQ. 0 ) THEN
          KERR(1)='POINT EXTERNE aux TRIANGLES ou QUADRANGLES VISIBLES'
          KERR(2)='CLIQUER un POINT a L''INTERIEUR d''un TRIANGLE ou QUA
     %DRANGLE'
         ELSE
            KERR(1)='POINT OUTSIDE TRIANGLES or QUADRANGLES'
            KERR(2)='CLICK a POINT INSIDE a TRIANGLE or QUADRANGLE'
         ENDIF
         CALL LEREUR
         IERR=1
         GOTO 9999
      ENDIF

C     SUPPRESSION DU TRIANGLE ou QUADRANGLE NOFACL
C     LE DERNIER TRIANGLE DEVIENT LE TRIANGLE ou QUADRANGLE NOFACL
C     ===========================================================
      MN  = MNSOEL + 4 * NBEFOB - 5
      MN1 = MNSOEL + 4 * NOFACL - 5
      PRINT*,'tuertq3d: SUPPRESSION du TQANGLE',NOFACL,' St:',
     %       (MCN(MN1+K),K=1,4)
      DO K=1,4
C        SAUVEGARDE DES 4 NUMEROS DES SOMMETS DU TQ SUPPRIME
         NOSOTQ(K)  = MCN(MN1+K)
C        LE TQ NOFACL PREND LES SOMMETS DU DERNIER TQ
         MCN(MN1+K) = MCN(MN+K)
C        LE TQ NBEFOB EST DESACTIVE
         MCN(MN +K) = 0
      ENDDO

C     UN EF DE MOINS
      NBEFOB = NBEFOB - 1
      MCN( MNNSEF+WBEFOB ) = NBEFOB

      CALL ECDATE( MCN(MNNSEF) )
      MCN( MNNSEF+MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

 9999 RETURN
      END
