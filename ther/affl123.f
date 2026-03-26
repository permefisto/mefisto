      SUBROUTINE AFFL123( NOMELE, NBELFX, NBPNFX, ndimco, NDIMFX,
     %                    NBCAS,  NCAS0, NCAS1,   COPNFX, FLUXNP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES VECTEURS FLUX DE CHALEUR AUX INTERFACES DES EF
C -----    D'UN TYPE D'EF DONNE POUR UN OBJET 1D ou 2D ou 3D
C
C ENTREES:
C --------
C NOMELE : NOM DU TYPE D'ELEMENT FINI
C NBELFX : NOMBRE D'ELEMENTS FINIS DU TYPE D'EF A TRAITER
C NBPNFX : NOMBRE DE POINTS DE CALCUL DES FLUX POUR 1 EF
C ndimco : NOMBRE DE COORDONNEES DU TABLEAU COPNFX
C NDIMFX : ESPACE DE TRAVAIL 1 OU 2 OU 3 DES EF
C NBCAS  : NOMBRE TOTAL DE CAS TRAITES
C NCAS0  : NUMERO DU PREMIER CAS A TRAITER
C NCAS1  : NUMERO DU DERNIER CAS A TRAITER
C COPNFX : LES NDIM COORDONNEES DES POINTS DE CALCUL DES FLUX NORMAUX
C FLUXNP : LES FLUX NORMAUX DE CHALEUR AUX POINTS DES FACES DES EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TIMS NTU TAIPEI TAIWAN           Octobre 2009
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL              COPNFX(1:ndimco,1:NBPNFX,1:NBELFX,1:2)
      DOUBLE PRECISION  FLUXNP(1:NBPNFX,1:NBELFX,1:NBCAS)
      CHARACTER*4       NOMELE(2)
      CHARACTER*2       EF
C
C     AFFICHAGE DES FLUX NORMAUX DES SOLUTIONS
C     ----------------------------------------
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         EF = 'EF'
         WRITE(IMPRIM,*) 'Les FLUX NORMAUX de CHALEUR aux POINTS des FAC
     %ES DES EF ', NOMELE
         WRITE(IMPRIM,*) '(XYZ COORDONNEES des POINTS  &  NX,NY,NZ COMPO
     %SANTES du VECTEUR NORMAL a la FACE)'
      ELSE
         EF = 'FE'
         WRITE(IMPRIM,*) 'The NORMAL HEAT FLUXES at FACE POINTS of FE ',
     %                    NOMELE
         WRITE(IMPRIM,*) '(XYZ POINT COORDINATES  &  NX,NY,NZ COMPONENTS
     %of FACE NORMAL VECTOR)'
      ENDIF
C
      IF( NDIMFX .EQ. 1 ) THEN
C
C        DIMENSION 1
         DO 15 K=1,NBELFX
            DO 10 L=1,NBPNFX
               WRITE(IMPRIM,10015) EF,K,(COPNFX(1,L,K,M),M=1,2),
     %                            (FLUXNP(L,K,NCAS),NCAS=NCAS0,NCAS1)
 10         CONTINUE
 15      CONTINUE
10015    FORMAT(A2,I6,': X=',G14.6,'  NX=',G14.6,(T46,' FLUX=',G14.6) )
C
      ELSE IF( NDIMFX .EQ. 2 ) THEN
C
C        DIMENSION 2
         DO 25 K=1,NBELFX
            DO 20 L=1,NBPNFX
               WRITE(IMPRIM,10025) EF,K,((COPNFX(N,K,L,M),N=1,2),M=1,2),
     %                            (FLUXNP(L,K,NCAS),NCAS=NCAS0,NCAS1)
 20         CONTINUE
 25      CONTINUE
10025  FORMAT(A2,I6,': X=',G14.6,' Y=',G14.6,'  NX=',G14.6,' NY=',G14.6,
     %        (T81,' FLUX=',G14.6) )
C
      ELSE IF( NDIMFX .EQ. 3 ) THEN
C
C        DIMENSION 3
         DO 35 K=1,NBELFX
            DO 30 L=1,NBPNFX
               WRITE(IMPRIM,10035) EF,K,((COPNFX(N,K,L,M),N=1,3),M=1,2),
     %                            (FLUXNP(L,K,NCAS),NCAS=NCAS0,NCAS1)
 30         CONTINUE
 35      CONTINUE
10035    FORMAT(A2,I6,': X=',G14.6,' Y=',G14.6,' Z=',G14.6,
     %         '  NX=',G14.6,' NY=',G14.6,' NZ=',G14.6,
     %         (T116,'  FLUX=',G14.6) )
      ENDIF
C
      RETURN
      END
