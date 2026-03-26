      SUBROUTINE MEFIFOAM( NF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : ECRIRE UN ENTETE D'UN FICHIER OpenFoam GENERE A PARTIR D'UN
C ----- MAILLAGE DE MEFISTO
C
C ENTREE:
C -------
C NF    : NUMERO D'UNITE DU FICHIER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     Mars 2011
C23456---------------------------------------------------------------012
C
10000 FORMAT(
     %'/*--------------------------------*- C++ -*----------------------
     %------------*\'/
     %'| =========  MEFISTO        | From a MEFISTO MESH to an OpenFOAM
     %MESH         |'/
     %'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolb
     %ox           |'/
     %'|  \\    /   O peration     | Version:  2.0
     %             |'/
     %'|   \\  /    A nd           | Web:      http://www.OpenFOAM.org
     %             |'/
     %'|    \\/     M anipulation  |
     %             |'/
     %'\*---------------------------------------------------------------
     %------------*/')
C
      WRITE(NF,10000)
C
      RETURN
      END
