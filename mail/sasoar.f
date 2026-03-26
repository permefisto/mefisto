      SUBROUTINE SASOAR( NOAR, MOSOAR, MXSOAR, N1SOAR, NOSOAR, NOARST )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER L'ARETE NOAR DU TABLEAU NOSOAR
C -----    SI CELLE CI N'EST PAS UNE ARETE DES LIGNES DE LA FONTIERE
C
C          LA METHODE EMPLOYEE ICI EST CELLE DU HACHAGE
C          AVEC POUR FONCTION D'ADRESSAGE H = MIN( NU2SAR(1), NU2SAR(2) )
C
C          ATTENTION: IL FAUT METTRE A JOUR LE NO D'ARETE DES 2 SOMMETS
C                     DE L'ARETE SUPPRIMEE DANS LE TABLEAU NOARST!
C
C ENTREES:
C --------
C NOAR   : NUMERO DE L'ARETE DE NOSOAR A SUPPRIMER
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE H
C MXSOAR : NOMBRE MAXIMAL D'ARETES STOCKABLES DANS LE TABLEAU NOSOAR
C          ATTENTION: MXSOAR>3*MXSOMM OBLIGATOIRE!
C
C MODIFIES:
C ---------
C N1SOAR : NO DE L'EVENTUELLE PREMIERE ARETE LIBRE DANS LE TABLEAU NOSOAR
C          CHAINAGE DES VIDES SUIVANT EN 3 ET PRECEDANT EN 2 DE NOSOAR
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          UNE ARETE I DE NOSOAR EST VIDE <=> NOSOAR(1,I)=0 ET
C          NOSOAR(4,ARETE VIDE)=L'ARETE VIDE QUI PRECEDE
C          NOSOAR(5,ARETE VIDE)=L'ARETE VIDE QUI SUIT
C NOARST : NUMERO D'UNE ARETE DE NOSOAR POUR CHAQUE SOMMET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE    UPMC PARIS  MARS    1997
C MODIFS : ALAIN PERRONNET LABORATOIRE JL LIONS UPMC PARIS  OCTOBRE 2006
C ...................................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NOSOAR(MOSOAR,MXSOAR), NOARST(*), NS(2)
c
c     13/10/2006
c     mise a jour de noarst pour les 2 sommets de l'arete a supprimer
c     necessaire uniquement pour les sommets frontaliers et internes imposes
c     le numero des 2 sommets de l'arete noar a supprimer
      ns(1) = nosoar(1,noar)
      ns(2) = nosoar(2,noar)
      do 8 k=1,2
         if( noarst(ns(k)) .eq. noar ) then
c           il faut remettre a jour le pointeur sur une arete
            if(nosoar(1,ns(k)).eq.ns(k) .and. nosoar(2,ns(k)).gt.0
     %         .and. nosoar(4,ns(k)) .gt. 0 ) then
c              arete active de sommet ns(k)
               noarst( ns(k) ) = ns(k)
            else
               do 5 i=1,mxsoar
                  if( nosoar(1,i).gt.0 .and. nosoar(4,i).gt.0 ) then
c                    arete non vide
                     if( nosoar(2,i).eq.ns(k) .or.
     %                  (nosoar(1,i).eq.ns(k).and.nosoar(2,i).gt.0))then
c                       arete active de sommet ns(k)
                        noarst( ns(k) ) = i
                        goto 8
                     endif
                  endif
 5             continue
ccc               print *,'sasoar: PB pas d arete de sommet',ns(k)
            endif
         endif
 8    continue
c     13/10/2006
C
      IF( NOSOAR(3,NOAR) .LE. 0 ) THEN
C
C        L'ARETE N'EST PAS FRONTALIERE => ELLE DEVIENT UNE ARETE VIDE
C
C        RECHERCHE DE L'ARETE QUI PRECEDE DANS LE CHAINAGE DU HACHAGE
         NOAR0 = 0
         NOAR1 = NOSOAR(1,NOAR)
C
C        PARCOURS DU CHAINAGE DU HACHAGE JUSQU'A RETROUVER L'ARETE NOAR
 10      IF( NOAR1 .NE. NOAR ) THEN
C
C           L'ARETE SUIVANTE PARMI CELLES AYANT MEME FONCTION D'ADRESSAGE
            NOAR0 = NOAR1
            NOAR1 = NOSOAR( MOSOAR, NOAR1 )
            IF( NOAR1 .GT. 0 ) GOTO 10
C
C           L'ARETE NOAR N'A PAS ETE RETROUVEE DANS LE CHAINAGE => ERREUR
            WRITE(IMPRIM,*) 'ERREUR SASOAR:ARETE NON DANS LE CHAINAGE '
     %                      ,NOAR
            WRITE(IMPRIM,*) 'ARETE DE ST1=',NOSOAR(1,NOAR),
     %      ' ST2=',NOSOAR(2,NOAR),' LIGNE=',NOSOAR(3,NOAR),
     %      ' TR1=',NOSOAR(4,NOAR),' TR2=',NOSOAR(5,NOAR)
            WRITE(IMPRIM,*) 'CHAINAGES=',(NOSOAR(I,NOAR),I=6,MOSOAR)
            CALL XVPAUSE
C           L'ARETE N'EST PAS DETRUITE
            RETURN
C
         ENDIF
C
         IF( NOAR .NE. NOSOAR(1,NOAR) ) THEN
C
C           SAUT DE L'ARETE NOAR DANS LE CHAINAGE DU HACHAGE
C           NOAR0 INITIALISEE EST ICI L'ARETE QUI PRECEDE NOAR DANS CE CHAINAGE
            NOSOAR( MOSOAR, NOAR0 ) = NOSOAR( MOSOAR, NOAR )
C
C           LE CHAINAGE DU HACHAGE N'EXISTE PLUS POUR NOAR
C           PAS UTILE CAR MISE A ZERO FAITE DANS LE SP HASOAR
CCC         NOSOAR( MOSOAR, NOAR ) = 0
C
C           NOAR DEVIENT LA NOUVELLE PREMIERE ARETE DU CHAINAGE DES VIDES
            NOSOAR( 4, NOAR ) = 0
            NOSOAR( 5, NOAR ) = N1SOAR
C           LA NOUVELLE PRECEDE L'ANCIENNE PREMIERE
            NOSOAR( 4, N1SOAR ) = NOAR
            N1SOAR = NOAR
C
CCC      ELSE
C
C           NOAR EST LA PREMIERE ARETE DU CHAINAGE DU HACHAGE H
C           CETTE ARETE NE PEUT ETRE CONSIDEREE DANS LE CHAINAGE DES VIDES
C           CAR LE CHAINAGE DU HACHAGE DOIT ETRE CONSERVE (SINON PERTE...)
C
         ENDIF
C
C        LE TEMOIN D'ARETE VIDE
         NOSOAR( 1, NOAR ) = 0
      ENDIF
      END
