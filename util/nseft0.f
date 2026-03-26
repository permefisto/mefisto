      SUBROUTINE NSEFT0( NTNSEF, MNNSEF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   SUPPRIMER LE NUMERO DES TANGENTES DES EF DU TABLEAU
C -----   'NSEF' POUR UN MAILLAGE STRUCTURE OU NON
C
C ENTREES:
C ---------
C NTNSEF :  NUMERO DE TMS DU TABLEAU 'NSEF'
C MNNSEF :  ADRESSE DU TABLEAU 'NSEF'
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS        MAI 1996
C2345X7..............................................................012
      include"./incl/a___nsef.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      IF( MNNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NSEFT0: ADRESSE NSEF INCORRECTE'
         CALL LEREUR
         RETURN
      ENDIF
C
C     EXISTE T IL DES TANGENTES ET UN TABLEAU DE POINTEURS SUR LES EF A TG? ?
      NBTGEF = MCN( MNNSEF + WBTGEF )
      IF( NBTGEF .EQ. 0 .AND. MCN(MNNSEF+WBEFAP) .EQ. 0 ) RETURN
C
C     LES PARAMETRES DU MAILLAGE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF,
     %             NBEFOB, NX    , NY    , NZ    ,
     %             IERR   )
C
C     IL FAUT SUPPRIMER LES TANGENTES EXISTANTES
      MCN( MNNSEF + WBTGEF ) = 0
      MCN( MNNSEF + WBEFTG ) = 0
      MCN( MNNSEF + WBEFAP ) = 0
C
C     REDUCTION DE LA TAILLE DU TMS NSEF
      CALL TAMSRA( NTNSEF, LDAPEF )
C
      RETURN
      END
