       SUBROUTINE TUERTQ2D( MNXYZS, MNNSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CLIQUER 1 POINT INTERNE A UN TRIANGLE ET DETRUIRE CE TRIANGLE
C -----
C ENTREES:
C --------
C MNXYZS : ADRESSE DU TABLEAU XYZSOMMET DE LA SURFACE
C MNNSEF : ADRESSE DU TABLEAU NSEF      DE LA SURFACE
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SI POINT EN DEHORS DE LA TRIANGULATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Labo J-L. LIONS  UPMC  PARIS   SEPTEMBRE 2007
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
C
      IERR = 0
C
C     CLIC DU POINT EXTERIEUR AU MAILLAGE
      CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
C
C     DEJA DANS UN TRIANGLE ou QUADRANGLE?
C     ====================================
      MNSOEL = MNNSEF + WUSOEF
      NBTQ   = MCN(MNNSEF+WBEFOB)
      MNXYZ  = MNXYZS + WYZSOM
      CALL TQPTCLIC( NX, NY, 2, MCN(MNXYZ), NBTQ, MCN(MNSOEL),  NUMTQ )
      IF( NUMTQ .EQ. 0 ) THEN
         NBLGRC(NRERR)=2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='POINT EXTERNE AUX TRIANGLES VISIBLES'
            KERR(2)='CLIQUER UN POINT A L''INTERIEUR D''UN TRIANGLE'
         ELSE
            KERR(1)='POINT OUTSIDE THE TRIANGLES'
            KERR(2)='CLICK A POINT INSIDE ONE TRIANGLE'
         ENDIF
         CALL LEREUR
         IERR=1
         RETURN
      ENDIF
C
C     SUPPRESSION DU TRIANGLE ou QUADRANGLE NUMTQ
C     LE DERNIER EF DEVIENT L'EF NUMTQ
C     ===========================================
      MN   = MNSOEL + 4 * NBTQ  - 5
      MN1  = MNSOEL + 4 * NUMTQ - 5
      DO K=1,4
         MCN(MN1+K) = MCN(MN+K)
      ENDDO
C
C     UN EF DE MOINS
      NBTQ = NBTQ - 1
      MCN(MNNSEF+WBEFOB) = NBTQ
C
      CALL ECDATE( MCN(MNNSEF) )
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )
C
      RETURN
      END
