      SUBROUTINE T2DVCO( LECRIT, NBS,    NOPXYD,
     %                   PXYD,   N1TRVI, NOTRIA, NOTRSO,
     %                   LETREE, NOLESO,
     %                   MXETRI, NAETOI, NTETOI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES TRIANGLES DELAUNAY A PARTIR DES POINTS
C -----    1 A NBS DU TABLEAU NOPXYD A VALEUR DANS PXYD(3+1 A 3+NBS)
C
C ENTREES:
C --------
C LECRIT : NUMERO DU CRITERE DE TRIANGULATION ( 1,2,3,...)
C NBS    : NOMBRE DE SOMMETS A TRAITER LORS DU PREMIER JEU
C NOPXYD : NUMERO DANS PXYD DES POINTS A TRIANGULER
C          PERMET L'ORDRE SOMMETS DE TE D'ABORD, PUIS POINTS UTILISATEUR
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : NUMERO DU 1 PREMIER TRIANGLE VIDE DANS LE TABLEAU NOTRIA
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOTRIA(4,.)
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS I EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE I
C
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          3 SATURATION DES TRIANGLES
C          8 AUCUN CERCLE CIRCONSCRIT AUX TRIANGLES ACTUELS NE CONTIENT
C            LE POINT COURANT
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXETRI : NOMBRE MAXIMAL D'ARETES OU TRIANGLES DECLARABLES
C          DANS UNE ETOILE
C NAETOI : ENTIER (1:4,1:MXETRI)
C NTETOI : ENTIER (1:MXETRI)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JUILLET 1994
C....................................................................012
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,*),
     %                  NOTRSO(*),
     %                  NOPXYD(0:*),
     %                  NAETOI(4,MXETRI),
     %                  NTETOI(MXETRI),
     %                  LETREE(0:8,0:*),
     %                  NOLESO(*)
      DOUBLE PRECISION  PXYD(3,*), D, DD
      CHARACTER*6       KTXT
C
C     REINITIALISATION A VIDE DES ARETES DE L'ETOILE
      N1AEVI = 1
      N1AEOC = 0
      DO 10 I=1,MXETRI
C        NUMERO DU TRIANGLE DE L'AUTRE COTE DE L'ARETE
C         NAETOI(1,I) = 0
C        NUMERO LOCAL AU TRIANGLE DE L'ARETE
C         NAETOI(2,I) = 0
C        NUMERO DANS NAETOI DE L'ARETE SUIVANTE
         NAETOI(4,I) = I+1
  10  CONTINUE
      NAETOI(4,MXETRI) = 0
C
C     TRAITEMENT DES SOMMETS DE TE D'ABORD PUIS DES POINTS INTERNES ET
C     FRONTALIERS NUMEROTES DE 4 A 3+NBS
      DO 1000 NN=1,NBS
C
C        LE NUMERO DANS PXYD DU POINT NN
         N = NOPXYD( NN )
C
C        REINITIALISATION A VIDE DES ARETES ETOILANT LE POINT N-1
         CALL ETOVID( N1AEVI, N1AEOC, NAETOI )
C
C        LE TABLEAU DES TRIANGLES A DETRUIRE EST VIDE
         NBTRET = 0
C
C        RECHERCHE DU 1-ER TRIANGLE NOTRI1
C        CONTENANT LE POINT N DE COORDONNEES ( PXYD(1,N) , PXYD(2,N) )
C        =============================================================
         CALL  TRDUPT( N, PXYD, LETREE, NOLESO, NOTRIA, NOTRSO,
     %                 NOTRI1 )
C
C        LE TRIANGLE NOTRI1 CONTIENT LE SOMMET N OU BIEN
C        PRODUIT UNE MEILLEURE QUALITE DE L'ETOILE
C        N EST MIS DANS LE TABLEAU NTETOI POUR ETRE ENSUITE DETRUIT
 270     NBTRET = NBTRET + 1
         NTETOI( NBTRET ) = NOTRI1
C
C        AJOUT OU RETRAIT DES 3 ARETES DU TRIANGLE NOTRI1 A L'ETOILE
         CALL AJTRET( NOTRI1, NOTRIA, N1AEVI, N1AEOC, NAETOI )
C
C        PASSAGE AU TRIANGLE SUIVANT PAR LES ARETES RECENSEES
C        DES TRIANGLES ADJACENTS DES ARETES NON TRAITEES DE L'ETOILE
C        ===========================================================
         NA1 = N1AEOC
         GOTO 420
C        ARETE SUIVANTE DE L'ETOILE
 400     NA1 = NAETOI( 4, NA1 )
 420     IF( NA1 .GT. 0 ) THEN
C           LE NO LOCAL DE L'ARETE DE L'ETOILE
            NAOP = NAETOI(2,NA1)
            IF( NAOP .LT. 0 ) THEN
C              ARETE DEJA TRAITEE . PASSAGE A LA SUIVANTE
               GOTO 400
            ELSE
C              L'ARETE ET SON TRIANGLE SONT ANALYSES ENSUITE
C              LE TRIANGLE DE L'AUTRE COTE
               NOTRI1 = NOTRIA( NAOP+3 , NAETOI(1,NA1) )
               IF( NOTRI1 .LE. 0 ) THEN
C                 L'ARETE EST DONC TRAITEE
                  NAETOI(2,NA1) = -NAOP
C                 L'ARETE SUIVANTE
                  GOTO 400
               ENDIF
C              LE TEMOIN DE RECHERCHE EFFECTUEE
               NAETOI(2,NA1) = -NAOP
C
C              LA QUALITE SERA T ELLE MEILLEURE AVEC CE TRIANGLE ?
               CALL MAXQ2T( LECRIT, N, PXYD,
     &                      NAOP, NAETOI(1,NA1), NOTRI1, NOTRIA,
     &                      QUAL12, QUAL34 )
               IF( QUAL12 .LT. QUAL34 ) THEN
C                 MEILLEURE QUALITE DE L'ETOILE AVEC CE TRIANGLE EN PLUS
                  GOTO 270
               ELSE
C                 PASSAGE A L'ARETE SUIVANTE
                  GOTO 400
               ENDIF
            ENDIF
         ENDIF
C        L'ETOILE EST COMPLETE
C
C        RECHERCHE D'UNE EVENTUELLE ARETE DONT N SERAIT MILIEU ?
C        BUT : EVITER LA CREATION D'UN TRIANGLE PLAT
C        =======================================================
         NA1 = N1AEOC
C        BOUCLE SUR LES ARETES DE L'ETOILE
 650     IF( NA1 .GT. 0 ) THEN
C
C           LE NO DU TRIANGLE ET LOCAL DE L'ARETE
            NT = NAETOI(1,NA1)
            I  = ABS( NAETOI(2,NA1) )
C           LE NUMERO DU TRIANGLE AU DELA DE L'ARETE
            NT1 = NOTRIA(I+3,NT)
            IF( NT1 .LE. 0 ) GOTO 680
C           LES 2 SOMMETS DE L'ARETE I DU TRIANGLE NT
            NS1 = NOTRIA( I, NT )
            NS2 = NOTRIA( NOSUI3(I), NT )
C
C           LE POINT N EST IL MILIEU DE L'ARETE NS1-NS2 ?
C           SI N EST MILIEU DE L'ARETE NS1-NS2 => GENERATION 1 TRIANGLE PLAT
            D  = 0
            DD = 0
            DO 660 J=1,2
               D  = D +( (PXYD(J,NS1)+PXYD(J,NS2))*0.5 - PXYD(J,N) )**2
               DD = DD+( PXYD(J,NS2) - PXYD(J,NS1) ) ** 2
 660        CONTINUE
            IF( D .LE. 1D-4 * DD ) THEN
C              ICI LE POINT N EST MILIEU DE NS1-NS2
C              LE TRIANGLE OPPOSE EST AJOUTE A L'ETOILE
               NOTRI1 = NT1
               GOTO 270
C              AINSI PAS DE TRIANGLE PLAT
            ENDIF
C
C           PASSAGE A L'ARETE SUIVANTE DE L'ETOILE
 680        NA1 = NAETOI(4,NA1)
            GOTO 650
         ENDIF
C
C        LE NO DE TRIANGLE => NO TRIANGLE AU DELA DE L'ARETE
C        LE NO LOCAL DANS LE TRIANGLE => NO 1-ER SOMMET DE L'ARETE
C        LE NO INUTILISE              => NO 2-ME SOMMET DE L'ARETE
C        =========================================================
         NA1 = N1AEOC
C        BOUCLE SUR LES ARETES
 690     IF( NA1 .GT. 0 ) THEN
C           LE NO DU TRIANGLE ET LOCAL DE L'ARETE
            NOTRI1 = NAETOI(1,NA1)
            I      = ABS( NAETOI(2,NA1) )
C           LE NUMERO DU TRIANGLE AU DELA DE L'ARETE
            NT     = NOTRIA(I+3,NOTRI1)
            NAETOI(1,NA1) = NT
C           LE NUMERO DES 2 SOMMETS DANS LE SENS DIRECT
            IF( I .NE. 3 ) THEN
               NS2 = I + 1
            ELSE
               NS2 = 1
            ENDIF
C           NUMERO DU SOMMET 1 DE L'ARETE
            NAETOI(2,NA1) = NOTRIA(I  ,NOTRI1)
C           NUMERO DU SOMMET 2 DE L'ARETE
            NAETOI(3,NA1) = NOTRIA(NS2,NOTRI1)
C           PASSAGE A L'ARETE SUIVANTE
            NA1 = NAETOI(4,NA1)
            GOTO 690
         ENDIF
C
C        DESTRUCTION DES TRIANGLES DU TABLEAU NTETOI ETOILANT CE POINT N
C        ================================================================
         DO 700 I=1,NBTRET
            NT = NTETOI( I )
C           TRACE DU TRIANGLE DETRUIT
CCC            CALL  DEVO29( PXYD , NOTRIA , NT , 0 , 7 )
C           NT DEVIENT LE PREMIER TRIANGLE VIDE
            NOTRIA(1,NT) = 0
            NOTRIA(4,NT) = N1TRVI
            N1TRVI       = NT
 700     CONTINUE
CCCC        LE TRACE DU POINT
CCC         CALL VASAB2( PXYD(1,N) , PXYD(2,N) )
CCC         CALL ECRSYM
         NBTRET = 0
C
C        DECLARATION DES TRIANGLES ETOILANT LE POINT N
C        =============================================
         NA2  = N1AEOC
C        BOUCLE SUR LES ARETES DE L'ETOILE
 750     IF( NA2 .GT. 0 ) THEN
C           L'ARETE NA2 ET LE POINT N FORMENT UN NOUVEAU TRIANGLE A DECLARER
C           LES 2 SOMMETS DE L'ARETE
            NS1 = NAETOI(2,NA2)
            NS2 = NAETOI(3,NA2)
C           LE TRIANGLE EXTERIEUR A L'ARETE AU DELA DE L'ARETE
            NT1 = NAETOI(1,NA2)
            IF( NT1 .LE. 0 ) THEN
C
C              LE POINT N EST IL MILIEU DE L'ARETE FRONTALIERE NS1-NS2 ?
C              SI N EST MILIEU DE L'ARETE NS1-NS2 => GENERATION 1 TRIANGLE PLAT
               D  = 0
               DD = 0
               DO 755 J=1,2
                  D =D +( (PXYD(J,NS1)+PXYD(J,NS2))/2 - PXYD(J,N) )**2
                  DD=DD+( PXYD(J,NS2) - PXYD(J,NS1) ) ** 2
 755           CONTINUE
               IF( D .LE. 1D-4 * DD ) THEN
C                 ICI LE POINT N EST MILIEU DE NS1-NS2
C                 ARETE NS1-NS2 FRONTALIERE.
C                 LE TRIANGLE PLAT FRONTALIER N'EST PAS GENERE
                  GOTO 780
               ENDIF
C              AINSI PAS DE TRIANGLE PLAT DU AUX TE FRONTALIERS
            ENDIF
C
            IF( N1TRVI .LE. 0 ) THEN
C              SATURATION DES TRIANGLES
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES TRIANGLES'
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
C           MISE A JOUR DU 1-ER TRIANGLE VIDE
            NT     = N1TRVI
            N1TRVI = NOTRIA(4,N1TRVI)
C
C           REMPLISSAGE DU TRIANGLE NT DE SOMMETS CEUX DE L'ARETE NA2 ET N
C           NS1 NS2 SOMMETS DE L'ARETE DANS LE SENS DIRECT DU TRIANGLE NT
C           N EST TOUJOURS SUPERIEUR A NS1 NS2
            NOTRIA(1,NT) = NS1
            NOTRIA(2,NT) = NS2
            NOTRIA(3,NT) = N
C           LE CHAINAGE DES TRIANGLES PAR LES ARETES
            NOTRIA(4,NT) = NT1
            NOTRIA(5,NT) = 0
            NOTRIA(6,NT) = 0
C           RECHERCHE DU NO LOCAL A NT1 DE L'ARETE NS1 NS2
            IF( NT1 .GT. 0 ) THEN
               DO 760 J=1,3
                  IF( NOTRIA(J,NT1) .EQ. NS2 ) GOTO 770
 760           CONTINUE
 770           NOTRIA(J+3,NT1) = NT
            ENDIF
C
C           MISE A JOUR D'UN TRIANGLE CONTENANT CHAQUE SOMMET
C           TOUS LES TRIANGLES ONT PU DISPARAITRE
C           CAS DE L'ETOILE REDUITE AUX 2 TRIANGLES INITIAUX
            NOTRSO(NS1) = NT
            NOTRSO(NS2) = NT
C
C           LA PILE DES TRIANGLES CREES DANS L'ETOILE
            NBTRET = NBTRET + 1
            IF( NBTRET .GT. MXETRI ) THEN
C              SATURATION DES TRIANGLES DE L'ETOILE
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES TRIANGLES DE L''ETOILE'
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
            NTETOI(NBTRET) = NT
C           TRACE DU TRIANGLE CONSTRUIT
            CALL DVTRTR(PXYD,NOTRIA, NT, MOD(NT-N1COUL,NDCOUL)+N1COUL,7)
C           LE TRACE DU POINT
            WRITE( KTXT, '(I6)' ) N
            XX = REAL( PXYD(1,N) )
            YY = REAL( PXYD(2,N) )
            CALL SYMBOLE2D( NCMAGE, XX, YY, KTXT )
C
C           PASSAGE A L'ARETE SUIVANTE DE L'ETOILE
 780        NA2 = NAETOI(4,NA2)
            GOTO 750
         ENDIF
C
C        LE POINTEUR SOMMET => TRIANGLE
         NOTRSO(N) = NT
C
C        COMPLETION DES CHAINAGES DES TRIANGLES CREES DANS L'ETOILE
 800     IF( NBTRET .GT. 0 ) THEN
C            LE HAUT DE LA PILE
             NT = NTETOI(NBTRET)
C            LE TRIANGLE EST DEPILE
             NBTRET = NBTRET - 1
             IF( NT .LE. 0 ) GOTO 800
C
C            QUEL EST LE TRIANGLE OPPOSE A L'ARETE 2 DE SOMMETS NS1 NS2
             NS2  = NOTRIA(2,NT)
             DO 810 J=1,NBTRET
                NT1 = NTETOI(J)
                IF( NT1 .LE. 0 ) GOTO 810
                IF( NOTRIA(1,NT1) .NE. NS2 ) GOTO 810
C               L'ARETE EST RETROUVEE
                NOTRIA(5,NT ) = NT1
                NOTRIA(6,NT1) = NT
C               SI LE TRIANGLE NT1 EST TOTALEMENT CHAINE
C               IL EST RETIRE DE L'ETOILE
                IF( NOTRIA(5,NT1) .GT. 0 ) NTETOI(J)=0
                GOTO 820
 810         CONTINUE
C
C            SI LE TRIANGLE NT EST TOTALEMENT CHAINE SAUT DE L'ARETE 3
 820         IF( NOTRIA(6,NT) .GT. 0 ) GOTO 800
C
C            QUEL EST LE TRIANGLE OPPOSE A L'ARETE 3 DE SOMMETS N-NS1
             NS1  = NOTRIA(1,NT)
             DO 830 J=1,NBTRET
                NT1 = NTETOI(J)
                IF( NT1 .EQ. 0 ) GOTO 830
                IF( NOTRIA(2,NT1) .NE. NS1 ) GOTO 830
C               L'ARETE EST RETROUVEE
                NOTRIA(6,NT ) = NT1
                NOTRIA(5,NT1) = NT
C               SI LE TRIANGLE NT1 EST TOTALEMENT CHAINE
C               IL EST RETIRE DE L'ETOILE
                IF( NOTRIA(6,NT1) .GT. 0 ) NTETOI(J)=0
                GOTO 800
 830         CONTINUE
C            RETOUR EN HAUT DE PILE
             GOTO 800
         ENDIF
C
C        L'ETOILE EST TRAITEE . PASSAGE AU POINT N SUIVANT
 1000 CONTINUE
C
C     FIN DU JEU DE POINTS 3+1 A 3+NBS
 9999 RETURN
      END
