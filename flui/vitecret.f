      SUBROUTINE VITECRET( NDIM, NBNOVI, VX, VY, VZ, VitECR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECRETAGE DE LA VITESSE AUX NOEUDS 
C -----    DES TETRAEDRES de TAYLOR HOOD ou BREZZI-FORTIN
C
C ENTREES :
C ---------
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C         (POUR TAYLOR HOOD CE SONT LES SOMMETS ET MILIEUX DES ARETES)
C         (POUR BREZZI FORTIN CE SONT LES SOMMETS ET BARYCENTRES DES EF
C          ET IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBNOPR
C             LES AUTRES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF)
C VitECR : VITESSE D'ECRETAGE

C ENTREES et SORTIES:
C -------------------
C VX     : LA VITESSE EN X EN CHAQUE NOEUD VITESSE ET NBVECT FOIS
C VY     : LA VITESSE EN Y EN CHAQUE NOEUD VITESSE ET NBVECT FOIS
C VZ     : LA VITESSE EN Z EN CHAQUE NOEUD VITESSE ET NBVECT FOIS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Veulettes sur mer               Octobre 2020
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VX(NBNOVI), VY(NBNOVI), VZ(NBNOVI),
     %                  VitECR, VNORM, D, SQRT

      NBECRET = 0
      IF( VitECR .LE. 0D0 ) GOTO 9999

C     TAYLOR-HOOD ou BREZZI-FORTIN
C     ===========    =============
      DO I=1,NBNOVI

C        NORME DE LA VITESSE au TEMPS tn
         VNORM = VX(I)**2 + VY(I)**2
         IF( NDIM .EQ. 3 ) THEN
            VNORM = VNORM  + VZ(I)**2
         ENDIF
         VNORM = SQRT( VNORM )

         IF( VNORM .GT. VitECR ) THEN
C           ECRETAGE DE LA VITESSE AU NOEUD I
            NBECRET = NBECRET + 1
            D = VitECR / VNORM
            VX(I) = VX(I) * D
            VY(I) = VY(I) * D
            IF( NDIM .EQ. 3 ) VZ(I) = VZ(I) * D
         ENDIF

      ENDDO

 9999 PRINT*,'vitecret.f:',NBECRET,'/',NBNOVI,
     %       ' NOEUDS de VITESSE ECRETEE a',VitECR

      RETURN
      END
