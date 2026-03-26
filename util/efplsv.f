      SUBROUTINE EFPLSV( MNNPEF, NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RETROUVER LE NUMERO DE VOLUME DE L'ELEMENT FINI
C -----             LE NUMERO DE SURFACE DES FACES   DE L'EF
C                   LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C                   LE NUMERO DE POINT   DES SOMMETS DE L'EF
C
C     ATTENTION :   LE PARCOURS DES ELEMENTS FINIS DOIT ETRE SEQUENTIEL
C     ===========   DU PREMIER AU DERNIER !  ( SINON ERREUR )
C
C ENTREES :
C ---------
C MNNPEF : ADRESSE MCN DU TABLEAU 'NPEF"NOM_ELEMENT_FINI'
C NUELEM : NUMERO DE L'ELEMENT FINI DANS CE TYPE
C
C MODIFIES :
C ----------
C NOVCEL : POINTE SUR LE 1-ER VOLUME
C NOSFEL : POINTE SUR LA 1-ERE SURFACE D'UNE FACE   DES EF
C NOLAEL : POINTE SUR LA 1-ERE LIGNE   D'UNE ARETE  DES EF
C NOPSEL : POINTE SUR LE 1-ER  POINT   D'UN  SOMMET DES EF
C
C       SI NUELEM=1 ALORS NOVCEL=NOSFEL=NOLAEL=NOPSEL SONT IMPOSES A 0
C       ENSUITE LA VALEUR SORTIE EST CELLE D'ENTREE DE L'EF NUELEM+1
C
C SORTIES :
C ---------
C NOOBVC : NUMERO DU VOLUME  DE L'EF
C NOOBSF : NUMERO DE SURFACE DES NFACE  FACES   DE L'EF
C NOOBLA : NUMERO DE LIGNE   DES NARET  ARETES  DE L'EF
C NOOBPS : NUMERO DE POINT   DES NBNSOM SOMMETS DE L'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/a___npef.inc"
      include"./incl/ponoel.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOOBSF(1:*),NOOBLA(1:*),NOOBPS(1:*)
C
C     LE NOMBRE D'ELEMENTS FINIS
      NBELEM = MCN( MNNPEF + WBELEM )
C
C     LE NOMBRE DE NOEUDS DE L'ELEMENT FINI
      NBNDEL = MCN( MNNPEF + WBNDEL )
C
C     LE NOMBRE DE POINTS GEOMETRIQUES DE L'ELEMENT FINI
      NBPGEL = MCN( MNNPEF + WBPGEL )
C
C     NBPSEL   NOMBRE DE POINTS UTILISATEUR SOMMET D UN ELEMENT
C     MOPSEL   NOMBRE DE VARIABLES DES TABLEAUX NLPSEL ET NEPSEL
C     NBLAEL   NOMBRE DE LIGNES UTILISATEUR ARETE D UN ELEMENT
C     MOLAEL   NOMBRE DE VARIABLES DES TABLEAUX NLLAEL ET NELAEL
C     NBSFEL   NOMBRE DE SURFACES UTILISATEUR FACE D UN ELEMENT
C     MOSFEL   NOMBRE DE VARIABLES DES TABLEAUX NLSFEL ET NESFEL
C     NBVCEL   NOMBRE DE VOLUMES UTILISATEUR CUBE D UN ELEMENT
C     MOVCEL   NOMBRE DE VARIABLES DES TABLEAUX NLVCEL ET NEVCEL
C
      NBPSEL = MCN( MNNPEF + WBPSEL )
      MOPSEL = MCN( MNNPEF + WOPSEL )
      NBLAEL = MCN( MNNPEF + WBLAEL )
      MOLAEL = MCN( MNNPEF + WOLAEL )
      NBSFEL = MCN( MNNPEF + WBSFEL )
      MOSFEL = MCN( MNNPEF + WOSFEL )
      NBVCEL = MCN( MNNPEF + WBVCEL )
      MOVCEL = MCN( MNNPEF + WOVCEL )
C
C     LES ADRESSES MCN DU DEBUT DES TABLEAUX
      MNNDEL = MNNPEF  + WUNDEL
      IF( NBPGEL .GT. 0 ) THEN
C        NOEUDS#POINTS
         MNPGEL = MNNDEL + NBELEM * NBNDEL
         MNPSEL = MNPGEL + NBELEM * NBPGEL
      ELSE
C        NOEUDS=POINTS
         MNPGEL = MNNDEL
         MNPSEL = MNPGEL + NBELEM * NBNDEL
      ENDIF
C
      MNLAEL = MNPSEL + NBPSEL + MOPSEL + MOPSEL
      MNSFEL = MNLAEL + NBLAEL + MOLAEL + MOLAEL
      MNVCEL = MNSFEL + NBSFEL + MOSFEL + MOSFEL
C
      IF( NUELEM .LE. 1 ) THEN
C        LE NOMBRE D'ITEMS PARCOURUS DANS LES TABLEAUX EST INITIALISE
C        POUR LE PREMIER ELEMENT FINI DE CE TYPE
         NOVCEL = 0
         NOSFEL = 0
         NOLAEL = 0
         NOPSEL = 0
      ENDIF
C
C     LE NUMERO DE VOLUME UTILISATEUR DU VOLUME DE L'ELEMENT
C     ------------------------------------------------------
C     ATTENTION : CE PROGRAMME SOUS ENTEND QUE
C                 LES ELEMENTS FINIS SONT RANGES PAR ORDRE CROISSANT
      IF( NBVCEL .GT. 0 ) THEN
         IF( MOVCEL .GT. 0  ) THEN
C           LE NUMERO DE L'ELEMENT EST IL NUELEM ?
            MN1 = MNVCEL + NOVCEL
            MN  = MN1    + NBVCEL + MOVCEL
            IF( MCN(MN) .EQ. NUELEM ) THEN
C              LE NUMERO DU VOLUME UTILISATEUR
               NOOBVC = MCN( MN1 )
C              L'ITEM EST INCREMENTE
               NOVCEL = NOVCEL + 1
            ELSE
               NOOBVC = 0
            ENDIF
         ELSE
C           PAS DE NUMERO LOCAL => TOUS SONT DONNES DE 1 A NBELEM
            NOOBVC = MCN( MNVCEL + NUELEM - 1 )
         ENDIF
      ELSE
         NOOBVC = 0
      ENDIF
C
C     LE NUMERO DE SURFACE UTILISATEUR DES FACES DE L'ELEMENT
C     -------------------------------------------------------
C     ATTENTION : CE PROGRAMME SOUS ENTEND QUE LES ITEMS
C                 DES ELEMENTS SONT RANGES PAR ORDRE CROISSANT
C                 DES ELEMENTS, DES FACES DES ELEMENTS
C
      IF( NBSFEL .GT. 0 ) THEN
         IF( MOSFEL .GT. 0 ) THEN
C           RESTE T IL DES FACES SUR UNE SURFACE ET
C           LE NUMERO DE L'ELEMENT EST IL NUELEM ?
            MN = MNSFEL + NBSFEL + NOSFEL
            IF(NOSFEL .LT. NBSFEL .AND. MCN(MN+MOSFEL) .EQ. NUELEM) THEN
C              BOUCLE SUR LES FACES
               DO 350 K=1,NFACE
C                 LA FACE APPARTIENT ELLE A UNE SURFACE ?
                  IF( MCN(MN+MOSFEL) .EQ. NUELEM .AND.
     %                      MCN(MN)  .EQ. K ) THEN
C                    OUI
                     NOOBSF(K) = MCN( MNSFEL + NOSFEL )
C                    PASSAGE A L'ITEM SUIVANT
                     NOSFEL    = NOSFEL + 1
                     MN        = MN + 1
                  ELSE
                     NOOBSF(K) = 0
                  ENDIF
 350           CONTINUE
            ELSE
               DO 360 K=1,NFACE
                  NOOBSF(K) = 0
 360           CONTINUE
            ENDIF
         ELSE
C           PAS DE NUMERO LOCAL => TOUS SONT DONNES DE 1 A NBELEM
            NOOBSF(1) = MCN( MNSFEL + NUELEM - 1 )
         ENDIF
      ELSE
         DO 361 K=1,NFACE
            NOOBSF(K) = 0
 361     CONTINUE
      ENDIF
C
C     LE NUMERO DE LIGNE UTILISATEUR DES ARETES DE L'ELEMENT
C     ------------------------------------------------------
C     ATTENTION : CE PROGRAMME SOUS ENTEND QUE LES ITEMS
C                 DES ELEMENTS SONT RANGES PAR ORDRE CROISSANT
C                 DES ELEMENTS,DES ARETES DES ELEMENTS,DES SOMMETS DES EF
C
      IF( NBLAEL .GT. 0 ) THEN
         IF( MOLAEL .GT. 0 ) THEN
C           RESTE T IL UNE ARETE SUR UNE LIGNE?
C           LE NUMERO DE L'ELEMENT EST IL NUELEM ?
            MN = MNLAEL + NBLAEL + NOLAEL
            IF(NOLAEL .LT. NBLAEL .AND. MCN(MN+MOLAEL) .EQ. NUELEM) THEN
C              BOUCLE SUR LES ARETES
               DO 370 K=1,NARET
C                 L'ARETE APPARTIENT ELLE A UNE LIGNE ?
                  IF( MCN(MN+MOLAEL) .EQ. NUELEM .AND.
     %                      MCN(MN)  .EQ. K ) THEN
C                    OUI
                     NOOBLA(K) = MCN( MNLAEL + NOLAEL )
C                    PASSAGE A L'ITEM SUIVANT
                     NOLAEL    = NOLAEL + 1
                     MN        = MN + 1
                  ELSE
                     NOOBLA(K) = 0
                  ENDIF
 370           CONTINUE
            ELSE
               DO 380 K=1,NARET
                  NOOBLA(K) = 0
 380           CONTINUE
            ENDIF
         ELSE
C           PAS DE NUMERO LOCAL => TOUS SONT DONNES DE 1 A NBELEM
            NOOBLA(1) = MCN( MNLAEL + NUELEM - 1 )
         ENDIF
      ELSE
         DO 381 K=1,NARET
            NOOBLA(K) = 0
 381     CONTINUE
      ENDIF
C
C     LE NUMERO DE POINT UTILISATEUR DES SOMMETS DE L'ELEMENT
C     -------------------------------------------------------
C     ATTENTION : CE PROGRAMME SOUS ENTEND QUE LES ITEMS
C                 DES ELEMENTS SONT RANGES PAR ORDRE CROISSANT
C                 DES ELEMENTS, DES SOMMETS DES ELEMENTS
C
      IF( NBPSEL .GT. 0 ) THEN
         IF( MOPSEL .GT. 0 ) THEN
C           RESTE T IL UN SOMMET QUI SOIT UN POINT?
C           LE NUMERO DE L'ELEMENT EST IL NUELEM ?
            MN = MNPSEL + NBPSEL + NOPSEL
            IF(NOPSEL .LT. NBPSEL .AND. MCN(MN+MOPSEL) .EQ. NUELEM) THEN
C              BOUCLE SUR LES SOMMETS
               DO 390 K=1,NBNSOM
C                 LE SOMMET APPARTIENT IL A UN POINT ?
                  IF( MCN(MN+MOPSEL) .EQ. NUELEM .AND.
     %                      MCN(MN)  .EQ. K ) THEN
C                    OUI
                     NOOBPS(K) = MCN( MNPSEL + NOPSEL )
C                    PASSAGE A L'ITEM SUIVANT
                     NOPSEL    = NOPSEL + 1
                     MN        = MN + 1
                  ELSE
                     NOOBPS(K) = 0
                  ENDIF
 390           CONTINUE
            ELSE
               DO 395 K=1,NBNSOM
                  NOOBPS(K) = 0
 395           CONTINUE
            ENDIF
         ELSE
C           PAS DE NUMERO LOCAL => TOUS SONT DONNES DE 1 A NBELEM
            NOOBPS(1) = MCN( MNPSEL + NUELEM - 1 )
         ENDIF
      ELSE
         DO 396 K=1,NBNSOM
            NOOBPS(K) = 0
 396     CONTINUE
      ENDIF
      END
