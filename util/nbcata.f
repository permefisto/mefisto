      INTEGER FUNCTION NBCATA( NBCACH , NBCHAI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE NOMBRE NBCATA DE CARACTERES REELS D UN TABLEAU DE
C ----- NBCHAI CHAINES DE NBCACH CARACTERES
C
C ENTREES :
C ---------
C NBCACH : NOMBRE DE CARACTERES PAR CHAINE
C NBCHAI : NOMBRE DE CHAINES DE NBCACH CARACTERES
C SORTIE  :
C ---------
C NBCATA : NOMBRE REEL DE CARACTERES DU TABLEAU
C***********************************************************************
C ATTENTION: SOUS-PROGRAMME DEPENDANT MACHINE ICI VERSION IBM
C***********************************************************************
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  DECEMBRE 1983
C.......................................................................
      include"./incl/msvaau.inc"
C
C     ATTENTION :
C     =========
C
C     SI L ORDINATEUR EST A ADRESSAGE PAR MOT
C     LA SEQUENCE SUIVANTE EST CORRECTE
C     ========================================================
C
C     LE NOMBRE DE MOTS NECESSAIRES POUR UNE CHAINE DE NBCACH CARACTERES
C     NBCATA = ( NBCACH - 1 ) / NBCHMO + 1
C
C     LE NOMBRE DE CARACTERES REELS DU TABLEAU
C     NBCATA = NBCATA * NBCHMO * NBCHAI
C
C     SI L ORDINATEUR EST A ADRESSAGE PAR OCTET( VERSION IBM )
C     LA SEQUENCE SUIVANTE EST CORRECTE
C     ========================================================
      NBCATA = NBCACH * NBCHAI
      RETURN
      END
