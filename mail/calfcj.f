      SUBROUTINE CALFCJ(ID,I0,J0,AI,BI,NBS1,NBS2,FCJ)
C***********************************************************************
C BUT : CALCUL DES FONCTIONS DE CONTROLE SUR TOUS LES NOEUDS DU MAILLAGE
C***********************************************************************
C
C ENTREES:
C           I0     : ABSCISSE TOPO. DU NOEUD OU EST CENTREE LA F. CONT.
C           J0     : ORDONNEE TOPO. DU NOEUD OU EST CENTREE LA F. CONT.
C           AI     : AMPLITUDE DE LA FONCTION DE CONTROLE
C           BI     : VALEUR DE 'L ETENDUE' DE L ACTION DE LA F. CONT.
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C
C SORTIES:
C           FCJ    : VECTEUR DES FONCTIONS DE CONTROLE EN J
C
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT NOVEMBRE 1988
C****+7**************************************************************012
      DIMENSION FCJ(NBS1*NBS2)
      J1 = MAX((J0-10),2)
      J2 = MIN((J0+10),(NBS2-1))
      I1 = MAX((I0-10),2)
      I2 = MIN((I0+10),(NBS1-1))
      IF (ID.NE.0) THEN
        IF (ID.EQ.1) THEN
          J1 = MAX(J0,2)
        ELSE
          J2 = MIN(J0,(NBS2-1))
        ENDIF
      ENDIF
      DO 20 J=J1,J2
        DO 10 I=I1,I2
          NE = (J-1)*NBS1+I
          FCJ(NE)= FCJ(NE)+FCJPNT(I,J,I0,J0,AI,BI)
   10   CONTINUE
   20 CONTINUE
      END
