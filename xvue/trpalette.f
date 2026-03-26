      SUBROUTINE TRPALETTE( NOPALC )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA PALETTE NOPALC DES COULEURS (Cf xvue/palcde.f)
C -----
C ENTREE :
C --------
C NOPALC : NUMERO DE LA PALETTE DES COULEURS A GENERER
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Septembre 2012
C........................................................................
C     DECLARATIONS NECESSAIRES POUR LA MANIPULATION DES COULEURS
      include"./incl/trvari.inc"
      CHARACTER*3 KNMF
C
C     TRACE DE LA PALETTE NOPALC DE COULEURS
C     =====================================
      CALL EFFACEMEMPX
C
C     NOMBRE D'EPAISSEURS DES TRAITS POUR VISUALISER LES COULEURS
      NBEP = 3
      CALL XVEPAISSEUR( NBEP )
C
C     LE NUMERO DE LA PALETTE EST TRACE
      CALL XVCOULEUR( NCBLAN )
      WRITE( KNMF, '(I3)' ) NOPALC
      CALL XVTEXTE( 'PALETTE'//KNMF, 10, 0, 20 )
C
C     TRACE DE LA PALETTE NOPALC
      NXR = 0
      DO J=1,NDCOUL
         CALL XVCOULEUR( J )
         CALL XVTRAIT( NXR, 25+J*NBEP, NXR+300, 25+J*NBEP )
      ENDDO
C
      CALL MEMPXFENETRE
      CALL XVVOIR
C
C     ARRET AVANT LA FIN
C     ==================
      PRINT *
      PRINT *,'TRPALETTE: ENTRER UN CARACTERE POUR FINIR'
      CALL XVSOURIS( NOTYEV, NC, NPX, NPY )
C
      RETURN
      END
