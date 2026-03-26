      FUNCTION FCJPNT(I,J,I0,J0,AI,BI)
C***********************************************************************
C BUT : CALCUL DE LA FONCTION DE CONTROLE ATTIRANT LES LIGNES ISO J
C       AUTOUR DU POINT (I0,J0)
C***********************************************************************
C
C ENTREES:
C           I      : ABSCISSE TOPO. DU NOEUD OU ON CALCULE LA F. CONT.
C           J      : ORDONNEE TOPO. DU NOEUD OU ON CALCULE LA F. CONT.
C           I0     : ABSCISSE TOPO. DU NOEUD OU EST CENTREE LA F. CONT.
C           J0     : ORDONNEE TOPO. DU NOEUD OU EST CENTREE LA F. CONT.
C           AI     : AMPLITUDE DE LA FONCTION DE CONTROLE
C           BI     : VALEUR DE 'L ETENDUE' DE L ACTION DE LA F. CONT.
C
C SORTIES:
C           FCJPNT : VALEUR DE LA FONCTION DE CONTROLE
C
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT NOVEMBRE 1988
C****+7**************************************************************012
      IF (J0.GT.J) THEN
        SGN = 1.
      ELSE
        IF (J0.EQ.J) THEN
          SGN = 0.
        ELSE
          SGN =-1.
        ENDIF
      ENDIF
      DIST = (I-I0)**2+(J-J0)**2
      FCJPNT = SGN*AI*EXP(-BI*SQRT(DIST))
      RETURN
      END
