      SUBROUTINE HASH45 ( MNNSEF, MNXYZS, MNHSEG, TAILLE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ETABLIR UNE TABLE DE HACHAGE DES SEGMENTS D'UN MAILLAGE D'UNE
C -----    SURFACE DEFINIE PAR LES TMS XYZSOMMET ET NSEF
C
C ENTREES:
C --------
C MNNSEF : ADRESSE DU TMS NSEF      DANS MCN (cf td/d/a___nsef)
C MNXYZS : ADRESSE DU TMS XYZSOMMET DANS MCN (cf td/d/a___xyzsommet)
C
C SORTIES:
C --------
C MNHSEG : ADRESSE DE LA HASHING TABLE DANS MCN CONTENANT :
C
C  MCN(MNHSEG + (i-1) * 3) = No DU 2eme SOMMET D'UN SEGMENT PARTANT
C                            DU POINT DONT LE NUMERO EST CELUI DE LA LIGNE
C                            A L'ORIGINE DE LA CHAINE DE POINTEURS VERS i.
C                            L'ORIGINE EST LE POINT i SI i<NBSOM, SINON, IL
C                            FAUT REMONTER LA CHAINE.
C
C  MCN(MNHSEG + (i-1) * 3 + 1) = POINTEUR SUR UNE AUTRE LIGNE DECRIVANT
C                                UN SEGMENT ISSU DU MEME POINT OU 0 SI IL
C                                N'Y EN A PLUS APRES.
C
C  MCN(MNHSEG + (i-1) * 3 + 2) = VAUT 0 SI LE SEGMENT EST INTERNE AU
C                                MAILLAGE, 1 SI LE SEGMENT EST AU BORD ET
C                                VOIT LA SURFACE A DROITE, -1 SI IL VOIT
C                                LA SURFACE A GAUCHE
C
C TAILLE : NOMBRE DE LIGNES RESERVEES EN MEMOIRE POUR LA HASHING TABLE
C          => TAILLE DE MCN(MNHSEG) = 3 * TAILLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : BALDENSPERGER-MATHIEU DEA ANALYSE NUMERIQUE UPMC JANVIER 1999
C2345X7..............................................................012
      include "./incl/langue.inc"
      include "./incl/gsmenu.inc"
      include "./incl/a___nsef.inc"
      include "./incl/a___xyzsommet.inc"
C
      include"./incl/pp.inc"
      COMMON MCN(MOTMCN)
      REAL RMCN(1)
      EQUIVALENCE (MCN(1), RMCN(1))
C
      INTEGER MNNSEF, MNXYZS
      INTEGER MNHSEG, TAILLE
C
      INTEGER NBTRI, NBSOM, NEXT
      INTEGER I, J, NCOGEL, POSPT1, POSPT2, POSPT3
      INTEGER NOSOEL(64)
      INTEGER ORIENT, INVERS, SOMMIN, SOMMAX, SUIVNT, POINT2
C
      REAL XPNT1, YPNT1, XPNT2, YPNT2, XPNT3, YPNT3
      REAL XVECT1, YVECT1, XVECT2, YVECT2, DETERM
C
C     LES CARACTERISTIQUES DES NSEF DE CETTE SURFACE CF LE TMS ~/td/d/a___nsef
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBTRI,
     %             NX    , NY    , NZ    ,
     %             IERR   )
C     NUTYMA : NUMERO DE TYPE DU MAILLAGE
C              0 : 'NON STRUCTURE'      , 2 : 'SEGMENT    STRUCTURE',
C              3 : 'TRIANGLE  STRUCTURE', 4 : 'QUADRANGLE STRUCTURE',
C              5 : 'TETRAEDRE STRUCTURE', 6 : 'PENTAEDRE  STRUCTURE',
C              7 : 'HEXAEDRE  STRUCTURE'
C     NBSOEL : NOMBRE DE SOMMETS DES NSEF
C              0 SI MAILLAGE NON STRUCTURE
C     NBSOEF : NOMBRE DE SOMMETS DE STOCKAGE DES NSEF
C              ( TRIANGLE NBSOEL=3  NBSOEF=4 )
C     NBTRI  : NOMBRE DE NSEF (TRIANGLES ET QUADRANGLES) DU MAILLAGE
C     NX, NY, NZ : LE NOMBRE D'ARETES DANS LES DIRECTIONS X Y Z
C                  SI LE MAILLAGE EST STRUTURE
C
C     NOMBRE DE SOMMETS DU MAILLAGE
      NBSOM = MCN( MNXYZS + WNBSOM )
C
C     POINTEUR AU DELA DU DERNIER SOMMET DU MAILLAGE
      NEXT  = NBSOM + 1
C
C     MAJORATION DU NOMBRE DE SEGMENTS DU CONTOUR
      TAILLE = NBTRI + 2 * NBSOM - 1
C
C     Reservation place memoire (majoree) pour la table et mise a zero
      CALL TNMCDC ('ENTIER', TAILLE * 3, MNHSEG)
      CALL AZEROI (TAILLE * 3, MCN(MNHSEG))
C
C     Boucle sur les elements finis
      DO 500 I = 0, NBTRI - 1
C
C        LE NUMERO DES SOMMETS DE L'EF I+1
         CALL NSEFNS( I+1,    NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C        NCOGEL: 3 POUR UN TRIANGLE ou 4 POUR UN QUADRANGLE
C        NOSOEL: LE NUMERO DES SOMMETS ET TANGENTES DE L'EF
         IF( IERR .NE. 0 ) RETURN
         IF( NCOGEL .NE. 3 .AND. NCOGEL .NE. 4 ) THEN
C           ERREUR
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'SURFACE AVEC AU MOINS UN EF DIFFERENT'
               KERR(2) = 'DE TRIANGLE OU QUADRANGLE'
            ELSE
               KERR(1) = 'SURFACE WITH AT LEAST ONE FINITE ELEMENT'
               KERR(2) = 'WHICH IS NOT a TRIANGLE or a QUADRANGLE'
            ENDIF
            CALL LEREUR
            IERR = 5
            RETURN
         ENDIF
C
C        boucle sur les sommets de l'ef
         NOSOEL(NCOGEL + 1) = NOSOEL(1)
C
C        Calcul de l'orientation du triangle ou quadrangle
         POSPT1 = MNXYZS + WYZSOM + 3 * (NOSOEL(1) - 1)
         XPNT1  = RMCN(POSPT1)
         YPNT1  = RMCN(POSPT1 + 1)
C
         POSPT2 = MNXYZS + WYZSOM + 3 * (NOSOEL(2) - 1)
         XPNT2  = RMCN(POSPT2)
         YPNT2  = RMCN(POSPT2 + 1)
C
         POSPT3 = MNXYZS + WYZSOM + 3 * (NOSOEL(NCOGEL) - 1)
         XPNT3  = RMCN(POSPT3)
         YPNT3  = RMCN(POSPT3 + 1)
C
         XVECT1 = XPNT2 - XPNT1
         YVECT1 = YPNT2 - YPNT1
C
         XVECT2 = XPNT3 - XPNT1
         YVECT2 = YPNT3 - YPNT1
C
C        LA SURFACE ORIENTEE
         DETERM = XVECT1 * YVECT2 - YVECT1 * XVECT2
         ORIENT = -INT(SIGN(1.0E0,DETERM))
C
C        Boucle sur les segments de l'element
         DO 400 J = 1, NCOGEL
C
C           Orientation du segment dans l'ordre des numeros des points
            IF ( NOSOEL(J) .LT. NOSOEL(J + 1) ) THEN
               SOMMIN = NOSOEL(J)
               SOMMAX = NOSOEL(J + 1)
               INVERS = 1
            ELSE
               SOMMIN = NOSOEL(J + 1)
               SOMMAX = NOSOEL(J)
               INVERS = -1
            ENDIF
C
C           On cherche l'emplacement a remplir dans la table
 100        SUIVNT = MCN(MNHSEG + (SOMMIN - 1) * 3 + 1)
            POINT2 = MCN(MNHSEG + (SOMMIN - 1) * 3)
            IF ( (SUIVNT .NE. 0) .AND. (POINT2 .NE. SOMMAX) ) THEN
              SOMMIN = SUIVNT
              GOTO 100
            ENDIF
C
C           On est soit a la fin de la chaine, soit sur le meme segment
            IF ( POINT2 .EQ. SOMMAX ) THEN
C              On est sur le meme segment
               MCN(MNHSEG + (SOMMIN - 1) * 3 + 2) = 0
            ELSE
C              On a un nouveau segment qu'on place sur la case en cours si elle
C              est vide et sinon sur la premiere case vide. Orientation.
               IF ( MCN(MNHSEG + (SOMMIN - 1) * 3) .EQ. 0 ) THEN
                  MCN(MNHSEG + (SOMMIN - 1) * 3) = SOMMAX
                  MCN(MNHSEG + (SOMMIN - 1) * 3 + 2) = ORIENT * INVERS
               ELSE
                  MCN(MNHSEG + (SOMMIN - 1) * 3 + 1) = NEXT
                  MCN(MNHSEG + (NEXT - 1) * 3) = SOMMAX
                  MCN(MNHSEG + (NEXT - 1) * 3 + 2) = ORIENT * INVERS
                  NEXT = NEXT + 1
               ENDIF
            ENDIF
C
C       Fin de la boucle sur les segments de l'element fini
 400    CONTINUE
C
C     Fin de la boucle sur les elements finis
 500  CONTINUE
C
      RETURN
      END
