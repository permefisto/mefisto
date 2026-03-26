      SUBROUTINE EDGE45 ( MCSEG, SIZESG, MCNXY, MCEDGE, MCXYED, NBEDGE,
     &                    MCCPCX, NBCPCX )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ALLOUER ET ETABLIR DES TABLEAUX MCN(MCEDGE) ET MCN(MCXYED) CONTENANT
C -----    LES POINTS DU BORD D'UNE SURFACE CARACTERISEE PAR LA HASH TABLE DE
C          SES SEGMENTS (cf. hash45.f) DANS MCN(MCSEG) ET PAR LES COORDONNEES
C          DE SES POINTS DANS MCN(MCNXY)
C
C ENTREES:
C --------
C MCNSEG : ADRESSE DE LA HASH TABLE DANS MCN (cf hash45.f)
C SIZESG : NOMBRE DE LIGNES DE LA HASH TABLE
C MCNXY  : ADRESSE DE LA TABLE DES COORDONNEES XYZ DANS RMCN
C          cf td/d/a___xyzsommet
C
C SORTIES:
C --------
C MCEDGE : ADRESSE DU TABLEAU CONTENANT LES POINTS DU BORD DE STRUCTURE :
C
C     MCN(MCEDGE + (i-1) * 4) = No DU ieme POINT TROUVE SUR LE BORD DANS
C                               LE TABLEAU MCN(MCNXY).
C
C     MCN(MCEDGE + (i-1) * 4 + 1) = POINTEUR SUR LA POSITION OCCUPEE DANS
C                                   LE TABLEAU PAR LE POINT SUIVANT DU BORD
C                                   PARCOURU DANS LE SENS DIRECT (INTERIEUR
C                                   VU SUR LA GAUCHE DU PARCOURS).
C
C     MCN(MCEDGE + (i-1) * 4 + 2) = POINTEUR SUR LA POSITION OCCUPEE DANS
C                                   LE TABLEAU PAR LE POINT SUIVANT DU BORD
C                                   PARCOURU DANS LE SENS INDIRECT (INTERIEUR
C                                   VU SUR LA DROITE DU PARCOURS).
C
C     MCN(MCEDGE + (i-1) * 4 + 3) = SI LE POINT EN COURS EST UN POINT
C                                   D'INTERSECTION AVEC UN CONTOUR D'UNE AUTRE
C                                   SURFACE, POINTEUR INDIQUANT LA LIGNE DU
C                                   TABLEAU MCN(MCEDGE) DE L'AUTRE SURFACE
C                                   CORRESPONDANT AU MEME POINT.
C
C                                   Rem : CETTE CASE EST LAISSEE A 0 DANS
C                                   CE SOUS-PROGRAMME. ELLE N'EST REMPLIE QUE
C                                   DANS cut45.f QUI CALCULE LES INTERSECTIONS
C                                   ENTRE DEUX CONTOURS.
C
C MCXYED : TABLEAU CONTENANT LES COORDONNEES DES POINTS DU BORD :
C
C     RMCN(MCXYED + (i-1) * 2) = X DU POINT LATERAL CORRESPONDANT A LA LIGNE i
C                               DU TABLEAU MCN(MCEDGE).
C
C     RMCN(MCXYED + (i-1) * 2 + 1) = Y DU POINT LATERAL CORRESPONDANT A LA
C                                   LIGNE i DU TABLEAU MCN(MCEDGE).
C
C NBEDGE : NOMBRE DE LIGNES ALLOUEES EN MEMOIRE POUR STOCKER LE BORD
C          => TAILLE DU TABLEAU MCN(MCEDGE) = 4 * NBEDGE
C          => TAILLE DU TABLEAU RMCN(MCEDGE) = 2 * NBEDGE
C
C MCCPCX : TABLEAU CONTENANT LES NUMEROS D'UN POINT DE CHAQUE COMPOSANTE CONNEXE
C          DU CONTOUR TRANSMIS.
C
C NBCPCX : NOMBRE DE COMPOSANTES CONNEXES DU CONTOUR.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : BALDENSPERGER-MATHIEU DEA ANALYSE NUMERIQUE UPMC JANVIER 1999
C2345X7..............................................................012
      INCLUDE "./incl/a___xyzsommet.inc"
C
      include"./incl/pp.inc"
      COMMON       MCN(MOTMCN)
      REAL        RMCN(1)
      EQUIVALENCE (MCN(1), RMCN(1))
C
      INTEGER MCSEG,  SIZESG, MCNXY
      INTEGER MCEDGE, MCXYED, NBEDGE
      INTEGER MCCPCX, NBCPCX
C
      INTEGER MCTEMP
      INTEGER NBPTS, I, IEDGE, PNTARR, POSPNT, LIGNE1, LIGNE2, POS
      INTEGER DEPART, CURRNT, ICPCX
C
C     Nombre de points du maillage
      NBPTS = MCN(MCNXY + WNBSOM)
C
C     Le nombre de points sur le bord est exactement egal au nombre de
C     segments sur le bord. Pour determiner la taille a allouer pour le
C     tableau MCN(MCEDGE), on parcourt la hash table des segments et on
C     compte le nombre de segments qui sont au bord.
      NBEDGE = 0
      DO I = 0, SIZESG - 1
        IF ( MCN(MCSEG + I * 3 + 2) .NE. 0 ) THEN
          NBEDGE = NBEDGE + 1
        ENDIF
      ENDDO
C
C     On alloue les tableaux MCN(MCEDGE) et RMCN(MCXYED) avec le nombre de
C     point lateraux NBEDGE qu'on vient de determiner et un tableau qui
C     va servir de crible pour la suite de taille le nombre de
      CALL TNMCDC ('ENTIER', NBEDGE * 4, MCEDGE)
      CALL TNMCDC ('REEL', NBEDGE * 2, MCXYED)
      CALL TNMCDC ('ENTIER', SIZESG, MCTEMP)
      CALL AZEROI (NBEDGE * 4, MCN(MCEDGE))
      CALL AZEROI (SIZESG, MCN(MCTEMP))
C
C     On parcourt les NBPTS premieres lignes de la hash table
      IEDGE = 0
      DO I = 0, NBPTS - 1
         POS = MCSEG + I * 3
C        Si il y a bien un segment qui part du point numero I vers un
C        point de numero superieur.
         IF ( (MCN(POS) .NE. 0) ) THEN
C
C           On boucle sur tous les segments partant du point I vers un point de
C           numero superieur (on suit la chaine de la hash table).
 100        IF ( MCN(POS + 2) .NE. 0 ) THEN
C              Si on n'a pas deja rencontre le point de depart I, on l'ajoute da
C              les tableaux des points du bord et on le biffe dans le crible.
C              Sinon, on recupere juste son numero de ligne LIGNE1 dans le
C                     tableau MCEDGE.
               IF ( MCN(MCTEMP + I) .EQ. 0 ) THEN
                  MCN(MCEDGE + IEDGE * 4) = I + 1
                  POSPNT = MCNXY + WYZSOM + 3 * I
                  RMCN(MCXYED + IEDGE * 2) = RMCN(POSPNT)
                  RMCN(MCXYED + IEDGE * 2 + 1) = RMCN(POSPNT + 1)
                  MCN(MCTEMP + I) = IEDGE + 1
                  LIGNE1 = IEDGE + 1
                  IEDGE = IEDGE + 1
               ELSE
                  LIGNE1 = MCN(MCTEMP + I)
               ENDIF
C           Pour le point d'arrivee, on effectue la meme operation que pour le
C           point de depart : si on l'a deja rencontre (presence dans le crible
C           MCTEMP), on ne note que la ligne a laquelle il se trouve dans MCEDGE
C           dans LIGNE2 sinon on l'ajoute dans MCEDGE et MCXYED et on le biffe
C           dans MCTEMP.
               PNTARR = MCN(POS)
               IF ( MCN(MCTEMP + PNTARR - 1) .EQ. 0 ) THEN
                  MCN(MCEDGE + IEDGE * 4) = PNTARR
                  POSPNT = MCNXY + WYZSOM + 3 * (PNTARR - 1)
                  RMCN(MCXYED + IEDGE * 2) = RMCN(POSPNT)
                  RMCN(MCXYED + IEDGE * 2 + 1) = RMCN(POSPNT + 1)
                  MCN(MCTEMP + PNTARR - 1) = IEDGE + 1
                  LIGNE2 = IEDGE + 1
                  IEDGE = IEDGE + 1
               ELSE
                  LIGNE2 = MCN(MCTEMP + PNTARR - 1)
               ENDIF
C           Si le segment parcouru dans l'ordre croissant des numeros de ses
C           extremites dans MCN(MCNXY) laisse la surface a sa gauche on fait
C           pointer le point I sur PNTARR dans le sens direct et PNTARR sur
C           le point I dans le sens indirect.
               IF ( MCN(POS + 2) .LT. 0 ) THEN
                  MCN(MCEDGE + (LIGNE1 - 1) * 4 + 1) = LIGNE2
                  MCN(MCEDGE + (LIGNE2 - 1) * 4 + 2) = LIGNE1
C           Si le segment laisse la surface a sa droite, on fait l'inverse
               ELSE
                  MCN(MCEDGE + (LIGNE1 - 1) * 4 + 2) = LIGNE2
                  MCN(MCEDGE + (LIGNE2 - 1) * 4 + 1) = LIGNE1
               ENDIF
C           Si la hash table pointe sur un nouveau segment issu du point I,
C           on boucle vers 100, sinon on passe au point I+1.
            ENDIF
            IF ( MCN(POS + 1) .NE. 0 ) THEN
               POS = MCSEG + (MCN(POS + 1) - 1) * 3
               GOTO 100
            ENDIF
         ENDIF
C     On est arrive au bout de la chaine de segments partant du point I
C     vers un point de numero plus eleve. On passe au point I+1.
      ENDDO
C
C     Recherche du nombre de composantes connexes: on parcourt le tableau
C     MCEDGE en suivant le bord dans le sens qui laisse la surface sur la
C     gauche du parcours. Lorsqu'on retombe sur le point de depart, on
C     ajoute 1 au nombre de composantes connexes. On biffe au fur et a
C     mesure dans le crible les points sur lesquels on est deja passe.
      CALL AZEROI (SIZESG, MCN(MCTEMP))
      NBCPCX = 0
      DEPART = 1
 200  CURRNT = DEPART
 300  MCN(MCTEMP + MCN(MCEDGE + (CURRNT - 1) * 4) - 1) = 1
      CURRNT = MCN(MCEDGE + (CURRNT - 1) * 4 + 1)
      IF ( CURRNT .EQ. DEPART ) THEN
        NBCPCX = NBCPCX + 1
      ELSE
        GOTO 300
      ENDIF
C
C     On cherche un nouveau point de depart sur lequel on n'est pas encore
C     passe. On s'arrete si on arrive au nombre de points du bord sans en
C     avoir trouve un.
 400  DEPART = DEPART + 1
      IF ( DEPART .LE. NBEDGE ) THEN
        IF ( MCN(MCTEMP+MCN(MCEDGE+(DEPART-1)*4)-1) .EQ. 0 ) THEN
          GOTO 200
        ELSE
          GOTO 400
        ENDIF
      ENDIF
C
C     On connait maintenant le nombre de composantes connexes avec CMPCNX.
C     On alloue le tableau MCN(MCCPCX) pour stocker le point de depart de
C     chaque composante connexe. Comment le trouve-t-on ? Mais betement
C     en refaisant exactement la meme chose que pour trouver NBCPCX sans
C     oublier de mettre chaque depart dans MCCPCX.
C     Avis aux amateurs pour trouver une methode plus elegante...
      CALL TNMCDC ('ENTIER', NBCPCX, MCCPCX)
      CALL AZEROI (SIZESG, MCN(MCTEMP))
      ICPCX = 0
      DEPART = 1
 500  CURRNT = DEPART
 600  MCN(MCTEMP + MCN(MCEDGE + (CURRNT - 1) * 4) - 1) = 1
      CURRNT = MCN(MCEDGE + (CURRNT - 1) * 4 + 1)
      IF ( CURRNT .EQ. DEPART ) THEN
         MCN(MCCPCX + ICPCX) = MCN(MCEDGE + (DEPART - 1) * 4)
         ICPCX = ICPCX + 1
      ELSE
         GOTO 600
      ENDIF
C
 700  DEPART = DEPART + 1
      IF ( DEPART .LE. NBEDGE ) THEN
         IF ( MCN(MCTEMP+MCN(MCEDGE+(DEPART-1)*4)-1) .EQ. 0 ) THEN
            GOTO 500
         ELSE
            GOTO 700
         ENDIF
      ENDIF
C
C     Desallocation du tableau temporaire MCN(MCTEMP)
      CALL TNMCDS ('ENTIER', SIZESG, MCTEMP)
C
      RETURN
      END
