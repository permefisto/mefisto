      FUNCTION NBELFI( MNNSEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE NOMBRE D'ELEMENTS FINIS DU MAILLAGE 'NSEF'
C -----    D'UN PLSV
C
C ENTREES:
C --------
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF' A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1994
C....................................................................012
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON    MCN(MOTMCN)
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBELFI,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      END
