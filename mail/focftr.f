      SUBROUTINE FOCFTR( NBTRCF, NOTRCF, nbarpi, PXYD,   NOARST,
     %                   MOSOAR, MXSOAR, N1SOAR, NOSOAR,
     %                   MOARTR, N1ARTR, NOARTR,
     %                   NBARCF, N1ARCF, NOARCF, nbstpe, nostpe,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FORMER UN CONTOUR FERME (CF) AVEC LES ARETES SIMPLES DES
C -----    NBTRCF TRIANGLES DU TABLEAU NOTRCF
C          DESTRUCTION DES NBTRCF TRIANGLES DU TABLEAU NOARTR
C          DESTRUCTION DES ARETES DOUBLES   DU TABLEAU NOSOAR
C
C          ATTENTION: LE CHAINAGE LCHAIN DE NOSOAR DEVIENT CELUI DES CF
C
C ENTREES:
C --------
C NBTRCF : NOMBRE DE  TRIANGLES DU CF A FORMER
C NOTRCF : NUMERO DES TRIANGLES DANS LE TABLEAU NOARTR
c nbarpi : numero du dernier sommet frontalier ou interne impose
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE
C MXSOAR : NOMBRE MAXIMAL D'ARETES STOCKABLES DANS LE TABLEAU NOSOAR
C          ATTENTION: MXSOAR>3*MXSOMM OBLIGATOIRE!
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C
C ENTREES ET SORTIES :
C --------------------
C NOARST : NOARST(I) NUMERO D'UNE ARETE DE SOMMET I
C N1SOAR : NUMERO DE LA PREMIERE ARETE VIDE DANS LE TABLEAU NOSOAR
C          UNE ARETE I DE NOSOAR EST VIDE  <=>  NOSOAR(1,I)=0
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          HACHAGE DES ARETES = NOSOAR(1)+NOSOAR(2)*2
C N1ARTR : NUMERO DU PREMIER TRIANGLE VIDE DANS LE TABLEAU NOARTR
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOARTR(2,.)
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1 = 0 SI TRIANGLE VIDE => ARETE2 = TRIANGLE VIDE SUIVANT
C
C SORTIES:
C --------
C NBARCF : NOMBRE D'ARETES DU CF
C N1ARCF : NUMERO D'UNE ARETE DE CHAQUE CONTOUR
C NOARCF : NUMERO DES ARETES DE LA LIGNE DU CONTOUR FERME
C ATTENTION: CHAINAGE CIRCULAIRE DES ARETES
C            LES ARETES VIDES POINTES PAR N1ARCF(0) NE SONT PAS CHAINEES
c nbstpe : nombre de  sommets perdus dans la suppression des triangles
c nostpe : numero des sommets perdus dans la suppression des triangles
C IERR   :  0 SI PAS D'ERREUR
C          14 SI LES LIGNES FERMEES SE COUPENT => DONNEES A REVOIR
C          15 SI UNE SEULE ARETE SIMPLE FRONTALIERE
C          16 SI BOUCLE INFINIE CAR TOUTES LES ARETES SIMPLES
C                DE LA BOULE SONT FRONTALIERES!
C          17 SI BOUCLE INFINIE DANS CAETOI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE    UPMC PARIS  Mars    1997
C MODIFS : ALAIN PERRONNET Laboratoire JL LIONS UPMC PARIS  Octobre 2006
C....................................................................012
      PARAMETER        (LCHAIN=6, mxstpe=512)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  PXYD(3,*)
      INTEGER           NOTRCF(1:NBTRCF)
      INTEGER           NOSOAR(MOSOAR,MXSOAR),
     %                  NOARTR(MOARTR,*),
     %                  N1ARCF(0:*),
     %                  NOARCF(3,*),
     %                  NOARST(*),
     %                  nostpe(mxstpe),
     %                  nosotr(3)
C
C     FORMATION DES ARETES SIMPLES DU CF AUTOUR DE L'ARETE NS1-NS2
C     ATTENTION: LE CHAINAGE LCHAIN DU TABLEAU NOSOAR DEVIENT ACTIF
C     ============================================================
C     ICI TOUTES LES ARETES DU TABLEAU NOSOAR VERIFIENT NOSOAR(LCHAIN,I) = -1
C     CE QUI EQUIVAUT A DIRE QUE L'ETOILE DES ARETES SIMPLES EST VIDE
C     (INITIALISATION DANS LE SP INSOAR PUIS REMISE A -1 DANS LA SUITE!)
      N1AEOC = 0
      IERR   = 0
c
c     13/10/2006
c     nombre de sommets des triangles a supprimer sans repetition
      nbst = 0
c     13/10/2006
C
C     AJOUT A L'ETOILE DES ARETES SIMPLES DES 3 ARETES DES TRIANGLES A SUPPRIMER
C     SUPPRESSION DES TRIANGLES DE L'ETOILE POUR LES ARETES SIMPLES DE L'ETOILE
      DO 10 I=1,NBTRCF
C
C        AJOUT OU RETRAIT DES 3 ARETES DU TRIANGLE NOTRCF(I) DE L'ETOILE
         NT = NOTRCF( I )
c
c        13/10/2006  ...............................................
         call nusotr( nt, mosoar, nosoar, moartr, noartr, nosotr )
ccc         print *,'focftr: triangle supprime',nt,' st:',nosotr
c
c        ajout des numeros de sommets non encore vus dans l'etoile
         do 3 k=1,3
            do 2 j=1,nbst
               if( nosotr(k) .eq. nostpe(j) ) goto 3
 2          continue
c           ajout du sommet
            nbst = nbst + 1
            nostpe( nbst ) = nosotr(k)
 3       continue
c        13/10/2006 ................................................
C
         DO 5 J=1,3
C           L'ARETE DE NOSOAR A TRAITER
            NOAR = ABS( NOARTR(J,NT) )
            CALL CAETOI( NOAR,   MOSOAR, MXSOAR, N1SOAR, NOSOAR, NOARST,
     %                   N1AEOC, NBTRAR  )
            IF( NBTRAR .LE. 0 ) THEN
               IERR = 17
               RETURN
            ENDIF
C           SI ARETE SIMPLE ALORS SUPPRESSION DU NUMERO DE TRIANGLE
C           POUR CETTE ARETE
            IF( NBTRAR .EQ. 1 ) THEN
               IF( NOSOAR(4,NOAR) .EQ. NT ) THEN
                  NOSOAR(4,NOAR) = NOSOAR(5,NOAR)
               ENDIF
               NOSOAR(5,NOAR) = -1
C           ELSE
C              L'ARETE APPARTIENT A AUCUN TRIANGLE => ELLE EST VIDE
C              LES POSITIONS 4 ET 5 SERVENT MAINTENANT AUX CHAINAGES DES VIDES
            ENDIF
  5      CONTINUE
 10   CONTINUE
C
C     LES ARETES SIMPLES DE L'ETOILE SONT REORDONNEES POUR FORMER UNE
C     LIGNE FERMEE = UN CONTOUR FERME PERIPHERIQUE DE L'ETOILE ENCORE DIT 1 CF
C     ========================================================================
      N1AE00 = N1AEOC
 12   NA1    = N1AEOC
C     LA PREMIERE ARETE DU CONTOUR FERME
      NS0 = NOSOAR(1,NA1)
      NS1 = NOSOAR(2,NA1)
C
C     L'ARETE EST-ELLE DANS LE SENS DIRECT?
C     RECHERCHE DE L'ARETE DU TRIANGLE EXTERIEUR NT D'ARETE NA1
      NT = NOSOAR(4,NA1)
      IF( NT .LE. 0 ) NT = NOSOAR(5,NA1)
C
C     ATTENTION AU CAS DE L'ARETE INITIALE FRONTALIERE DE NO DE TRIANGLES 0 ET -
      IF( NT .LE. 0 ) THEN
C        PERMUTATION CIRCULAIRE DES ARETES SIMPLES CHAINEES
C        LA PREMIERE ARETE DOIT DEVENIR LA DERNIERE DU CHAINAGE,
C        LA 2=>1, LA 3=>2, ... , LA DERNIERE=>L'AVANT DERNIERE, 1=>DERNIERE
         N1AEOC = NOSOAR( LCHAIN, N1AEOC )
         IF( N1AEOC .EQ. N1AE00 ) THEN
C           ATTENTION: BOUCLE INFINIE SI TOUTES LES ARETES SIMPLES
C           DE LA BOULE SONT FRONTALIERES!... ARRETEE PAR CE TEST
            IERR = 16
            RETURN
         ENDIF
         NOAR = N1AEOC
         NA0  = 0
 14      IF( NOAR .GT. 0 ) THEN
C           LA SAUVEGARDE DE L'ARETE ET L'ARETE SUIVANTE
            NA0  = NOAR
            NOAR = NOSOAR(LCHAIN,NOAR)
            GOTO 14
         ENDIF
         IF( NA0 .LE. 0 ) THEN
C           UNE SEULE ARETE SIMPLE FRONTALIERE
            IERR = 15
            RETURN
         ENDIF
C        LE SUIVANT DE L'ANCIEN DERNIER EST L'ANCIEN PREMIER
         NOSOAR(LCHAIN,NA0) = NA1
C        LE NOUVEAU DERNIER EST L'ANCIEN PREMIER
         NOSOAR(LCHAIN,NA1) = 0
         GOTO 12
      ENDIF
C
C     ICI L'ARETE NA1 EST L'UNE DES ARETES DU TRIANGLE NT
      DO 15 I=1,3
         IF( ABS(NOARTR(I,NT)) .EQ. NA1 ) THEN
C           C'EST L'ARETE
            IF( NOARTR(I,NT) .GT. 0 ) THEN
C              ELLE EST PARCOURUE DANS LE SENS INDIRECT DE L'ETOILE
C             (CAR C'EST EN FAIT LE TRIANGLE EXTERIEUR A LA BOULE)
               NS0 = NOSOAR(2,NA1)
               NS1 = NOSOAR(1,NA1)
            ENDIF
            GOTO 17
         ENDIF
 15   CONTINUE
C
C     LE 1-ER SOMMET OU ARETE DU CONTOUR FERME
 17   N1ARCF( 1 ) = 1
C     LE NOMBRE DE SOMMETS DU CONTOUR FERME DE L'ETOILE
      NBARCF = 1
C     LE PREMIER SOMMET DE L'ETOILE
      NOARCF( 1, NBARCF ) = NS0
C     L'ARETE SUIVANTE DU CF
      NOARCF( 2, NBARCF ) = NBARCF + 1
C     LE NUMERO DE CETTE ARETE DANS LE TABLEAU NOSOAR
      NOARCF( 3, NBARCF ) = NA1
C     MISE A JOUR DU NUMERO D'ARETE DU SOMMET NS0
      NOARST(NS0) = NA1
C
C     TRACE DE L'ARETE
      CALL DVTRAR( PXYD, NS0, NS1, NCVERT, NCBLAN )
C
C     L'ARETE SUIVANTE A CHAINER
      N1AEOC = NOSOAR( LCHAIN, NA1 )
C     L'ARETE NA1 N'EST PLUS DANS L'ETOILE
      NOSOAR( LCHAIN, NA1 ) = -1
C
C     BOUCLE SUR LES ARETES SIMPLES DE L'ETOILE
 20   IF( N1AEOC .GT. 0 ) THEN
C
C        RECHERCHE DE L'ARETE DE 1-ER SOMMET NS1
         NA0 = -1
         NA1 = N1AEOC
 25      IF( NA1 .GT. 0 ) THEN
C
C           LE NUMERO DU DERNIER SOMMET DE L'ARETE PRECEDENTE
C           EST IL L'UN DES 2 SOMMETS DE L'ARETE NA1?
            IF ( NS1 .EQ. NOSOAR(1,NA1) ) THEN
C               L'AUTRE SOMMET DE L'ARETE NA1
                NS2 = NOSOAR(2,NA1)
            ELSE IF( NS1 .EQ. NOSOAR(2,NA1) ) THEN
C               L'AUTRE SOMMET DE L'ARETE NA1
                NS2 = NOSOAR(1,NA1)
            ELSE
C              NON: PASSAGE A L'ARETE SUIVANTE
               NA0 = NA1
               NA1 = NOSOAR( LCHAIN, NA1 )
               GOTO 25
            ENDIF
C
C           OUI: NA1 EST L'ARETE PERIPHERIQUE SUIVANTE
C                NA0 EST SA PRECEDENTE DANS LE CHAINAGE
C           UNE ARETE DE PLUS DANS LE CONTOUR FERME (CF)
            NBARCF = NBARCF + 1
C           LE PREMIER SOMMET DE L'ARETE NBARCF PERIPHERIQUE
            NOARCF( 1, NBARCF ) = NS1
C           L'ARETE SUIVANTE DU CF
            NOARCF( 2, NBARCF ) = NBARCF + 1
C           LE NUMERO DE CETTE ARETE DANS LE TABLEAU NOSOAR
            NOARCF( 3, NBARCF ) = NA1
C           MISE A JOUR DU NUMERO D'ARETE DU SOMMET NS1
            NOARST(NS1) = NA1
C
C           TRACE DE L'ARETE
            CALL DVTRAR( PXYD, NS1, NS2, NCVERT, NCBLAN )
C
C           SUPPRESSION DE L'ARETE DES ARETES SIMPLES DE L'ETOILE
            IF( N1AEOC .EQ. NA1 ) THEN
                N1AEOC = NOSOAR( LCHAIN, NA1 )
            ELSE
                NOSOAR( LCHAIN, NA0 ) = NOSOAR( LCHAIN, NA1 )
            ENDIF
C           L'ARETE N'EST PLUS UNE ARETE SIMPLE DE L'ETOILE
            NOSOAR( LCHAIN, NA1 ) = -1
C
C           LE SOMMET FINAL DE L'ARETE A RECHERCHER ENSUITE
            NS1 = NS2
            GOTO 20
         ENDIF
      ENDIF
C
C     VERIFICATION
      IF( NS1 .NE. NS0 ) THEN
C        ARETE NON RETROUVEE : L'ETOILE NE SE REFERME PAS
         NBLGRC(NRERR) = 3
         KERR(1) = 'FOCFTR: REVOYEZ VOS DONNEES'
         KERR(2) = 'LES LIGNES FERMEES DOIVENT ETRE DISJOINTES'
         KERR(3) = 'LA FONCTION TAILLE_IDEALE PAS TROP DISCONTINUE'
         CALL LEREUR
         IERR = 14
         RETURN
      ENDIF
C
C     L'ARETE SUIVANT LA DERNIERE ARETE DU CF EST LA PREMIERE DU CF
C     => REALISATION D'UN CHAINAGE CIRCULAIRE DES ARETES DU CF
      NOARCF( 2, NBARCF ) = 1
c
c     13/10/2006
c     existe t il des sommets perdus?
c     -------------------------------
      if( nbst .gt. mxstpe ) then
         write(imprim,*)'focftr: tableau nostfe(',mxstpe,') a augmenter'
         ierr = 15
         return
      endif
c     le nombre de sommets perdus
      nbstpe = nbst - nbarcf
      if( nbstpe .gt. 0 ) then
c        oui: stockage dans nostpe des sommets perdus
c        tout sommet des aretes de l'etoile est supprime
c        de la liste des sommets
         do 40 i=1,nbarcf
c           le numero du sommet de l'arete du cf
            ns1 = noarcf( 1, i )
            do 30 j=1,nbst
               if( ns1 .eq. nostpe(j) ) then
c                 le sommet peripherique est supprime
c                 de la liste des sommets perdus
                  nostpe(j) = 0
                  goto 40
               endif
 30         continue
 40      continue
c
c        compression
         n = 0
         do 45 i=1,nbst
            if( nostpe(i) .eq. 0 .or. nostpe(i) .gt. nbarpi ) then
c              un sommet de l'etoile ou perdu mais supprimable
c              ce qui apporte plus de qualites aux triangles a former
               n = n + 1
            else
c              un sommet perdu
               nostpe(i-n) = nostpe(i)
            endif
 45      continue
         nbstpe = nbst - n
ccc      print*,'focftr:',nbstpe,' sommets isoles:',(nostpe(k),k=1,nbstpe)
      endif
c     13/10/2006
C
C     DESTRUCTION DES TRIANGLES DE L'ETOILE DU TABLEAU NOARTR
C     -------------------------------------------------------
      DO 60 N=1,NBTRCF
C        LE NUMERO DU TRIANGLE DANS NOARTR
         NT0 = NOTRCF( N )
C        L'ARETE 1 DE NT0 DEVIENT NULLE
         NOARTR( 1, NT0 ) = 0
C        CHAINAGE DE NT0 EN TETE DU CHAINAGE DES TRIANGLES VIDES DE NOARTR
         NOARTR( 2, NT0 ) = N1ARTR
         N1ARTR = NT0
 60   CONTINUE
      END
