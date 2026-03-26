      SUBROUTINE MEMSUI(SUITE1,SUITE2,LS1,LS2,IPOS)
C ======================================================================
C BUT : RETROUVER LA SUITE1 DANS LA SUITE2
C ======================================================================
C ENTREE :
C     SUITE1 : LA SUITE NUMERO 1
C     SUITE2 : LA SUITE NUMERO 2
C     LS1    : LONGUEUR DE LA SUITE 1
C     LS2    : LONGUEUR DE LA SUITE 2
C     IPOS   = 0 SEUL DANS LE SENS DIRECTE
C            = X  DANS LES DEUX SENS
C SORTIE :
C     IPOS   : POSITION DU DEBUT DE SUITE 1 DANS SUITE 2
C              0 SI PAS RETROUVE
C REMARQUE :
C     LA SIUTE 1 PEUT ETRE INTERROMPUE PAR LA FIN DE SUITE 2 ( BOUCLE )
C ======================================================================
C     A.GOLGOLAB  INRIA JUILLET 1988
C ======================================================================
C
      INTEGER SUITE1(LS1),SUITE2(LS2)
C
      DO 1 IS2 = 1 , LS2
         IPOSI = IS2
         DO 2 IS1 = 1 , LS1
            I = IPOSI + IS1 - 1
            ID = LS2 - I
            IF ( ID .LT. 0 ) THEN
               I = -ID
            ENDIF
            IF ( SUITE2(I) .NE. SUITE1(IS1) ) GOTO 1
    2    CONTINUE
         IPOS = IPOSI
         RETURN
    1 CONTINUE
C
      IF ( IPOS .EQ. 0 ) RETURN
C
      DO 3 IS2 = 1 , LS2
         IPOSI = IS2
         DO 4 IS1 =  LS1 , 1 , -1
            I = IPOSI + LS1 - IS1
            ID = LS2 - I
            IF ( ID .LT. 0 ) THEN
               I = -ID
            ENDIF
            IF ( SUITE2(I) .NE. SUITE1(IS1) ) GOTO 3
    4    CONTINUE
         IPOS = IPOSI
         RETURN
    3 CONTINUE
      IPOS = 0
      END
