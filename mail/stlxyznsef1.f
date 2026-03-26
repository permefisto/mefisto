      SUBROUTINE STLXYZNSEF1( KNMSFFR, MNXYZS, MNNSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   EXPORTER SUR UN FICHIER .stl en CARACTERES ASCII
C -----   la TRIANGULATION d'une SURFACE 3D FERMEE
C         c-a-d Le nom de la surface fermee
C               le nombre de triangles
C               Pour chaque triangle
C                   XYZN du VECTEUR NORMAL dirige vers l'exterieur
C                   XYZ1 du SOMMET 1 du triangle
C                   XYZ1 du SOMMET 2 du triangle
C                   XYZ1 du SOMMET 3 du triangle dans le sens direct
C
C         obtenus a partir des TMS XYZSOMMET et NSEF Mefisto
C         de la triangulation de la SURFACE 3D FERMEE
C
C ENTREES:
C --------
C KNMSFFR : NOM DU PLSV DANS LE LEXIQUE DES PLSV
C MNXYZS  : >0 ADRESSE MCN du TMS 'XYZSOMMET'
C MNNSEF  : =0 ADRESSE MCN du TMS 'NSEF' pour un POINT
C           >0 ADRESSE MCN du TMS 'NSEF' pour une LIGNE ou SURFACE ou VOLUME

C SORTIE :
C --------
C IERR   : =0  le FICHIER KNMSFFR.stl est ECRIT dans
C              le REPERTOIRE du PROJET
C          >0  ERREUR RENCONTREE dans le TMS 'NSEF'. PAS de FICHIER CREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY          Decembre 2024
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*24       KNMSFFR
      CHARACTER*28       KNMFIC
      INTEGER            NOSOEF(4)
      REAL               XYZVN(3), XYZST(3,3)
      LOGICAL            LEXIST, LOPEN

      NUTYOB = 3
      IERR   = 0
      IF( MNXYZS .LE. 0 .OR. MNNSEF .LE.0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *, 'Le fichier xyznsef ou nsef N''EXISTE PAS'
         ELSE
            PRINT *, 'The file xyznsef or nsef IS UNKNOWN'
         ENDIF
         IERR = 1
         RETURN
      ENDIF

C     DECLARATION OUVERTURE DU FICHIER KNMSFFR.stl (en lettres minuscules)
C     ============================================
C     LE NOM DU FICHIER a ECRIRE dans le REPERTOIRE du PROJET
      CALL MINUSC( KNMSFFR )
      I      = NUDCNB( KNMSFFR )
      KNMFIC = KNMSFFR(1:I) // '.stl'

C     SI LE FICHIER KNMFIC EXISTE ALORS IL EST DETRUIT PUIS RECONSTRUIT
      INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( LEXIST ) THEN
C        LE FICHIER KNMFIC EXISTE
         IF( .NOT. LOPEN ) THEN
C           OUVERTURE DU FICHIER KNMFIC
            CALL TRUNIT( NF )
            OPEN( FILE=KNMFIC, UNIT=NF, STATUS='OLD' )
         ENDIF
C        DESTRUCTION DU FICHIER
         CLOSE( NF, STATUS='DELETE' )
      ENDIF

C     CREATION DU FICHIER KNMFIC = KNMSFFR.stl
C     ----------------------------------------
      CALL TRUNIT( NF )
      OPEN( UNIT=NF, ERR=9900, STATUS='NEW',
     %      FILE=KNMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )

C     ECRITURE DU NOM de la SURFACE en PREMIERE LIGNE du FICHIER.stl
C     --------------------------------------------------------------
      WRITE(NF,10001) KNMSFFR
10001 FORMAT('solid ', A24)

C     LE NOMBRE DE SOMMETS DU MAILLAGE
      NBSOM = MCN( MNXYZS + WNBSOM )

      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,' NOM de la SURFACE FERMEE ', KNMSFFR
      ELSE
         PRINT *,' CLOSED SURFACE NAME ', KNMSFFR
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
C        ECRITURE DU NOMBRE DE SOMMETS 
         PRINT *,' NOMBRE de SOMMETS de la SURFACE=',NBSOM
      ELSE
         PRINT *,' SURFACE VERTICES NUMBER=',NBSOM
      ENDIF

C     LA TOPOLOGIE : RECUPERATION du TMS NSEF NON STRUCTURE
C     ============== --------------------------------------
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN

      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,' NOMBRE de TRIANGLES de la SURFACE=',NBEFOB
      ELSE
         PRINT *,' SURFACE TRIANGLES NUMBER=',NBEFOB
      ENDIF

C     LA BOUCLE SUR LES NBEFOB TRIANGLES DU MAILLAGE
C     ==============================================
C     ADRESSE-4 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS  = MNXYZS + WYZSOM - 4
      IERR = 0

      DO 100 NEF = 1, NBEFOB

C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( NEF   , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEF, IERR )

C        NBSOEF NOMBRE de SOMMETS DE L'EF
C        WRITE(NF,11108) (NOSOEF(K),K=1,NBSOEF)
         IF( NOSOEF(4) .NE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT *,' EF ',NEF,' de SOMMETS ', NOSOEF,
     %                 ' NON un TRIANGLE'
            ELSE
               PRINT *,' FE ',NEF,' of VERTICES',NOSOEF,
     %                 ' is NOT a TRIANGLE'
            ENDIF
            IERR = IERR + 1
            GOTO 100
         ENDIF

C        RECUPERATION DES XYZ DES 3 SOMMETS DU TRIANGLE
         DO NS = 1, 3
            MNST = MNS + 3 * NOSOEF( NS )
            DO J=1,3
               XYZST( J, NS ) = RMCN( MNST + J )
            ENDDO
         ENDDO

C        CALCUL du VECTEUR NORMAL UNITAIRE XYZVN au TRIANGLE NEF
         CALL NORFAC( 3, XYZST, XYZVN, IER )
C        IER : =0 SI VECTEUR NORMAL CORRECTEMENT CALCULE
C              >0 SINON   et XYZVN=(0,0,0)

C        ECRITURE du VECTEUR NORMAL UNITAIRE SUR LE FICHIER
10010    FORMAT('  facet normal', 3(1X,E12.6) )
         WRITE(NF,10010) XYZVN

10020    FORMAT('    outer loop')
         WRITE(NF,10020)

10030    FORMAT('      vertex',  3(1X,E12.6) )
         DO J=1,3
            WRITE(NF,10030) (XYZST(I,J),I=1,3)
         ENDDO

10040    FORMAT('    endloop')
         WRITE(NF,10040)

10050    FORMAT('  endfacet')
         WRITE(NF,10050)

 100  ENDDO
      IF( IERR .NE. 0 ) GOTO 9910


C     ECRITURE DU NOM de la SURFACE en DERNIERE LIGNE du FICHIER
C     ----------------------------------------------------------
      WRITE(NF,10100) KNMSFFR
10100 FORMAT('endsolid ', A24)

      CLOSE( NF )

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'Le fichier ',KNMFIC,' est ECRIT avec ',
     %             NBSOM,' sommets et ',NBEFOB,' facettes'
      ELSE
         PRINT *, 'The file ',KNMFIC,' is WRITTEN with ',
     %             NBSOM,' vertices and ',NBEFOB,' facets'
      ENDIF
      PRINT*
      RETURN

C     PB A L'OUVERTURE DU FICHIER xyznsef
 9900 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'PB a l''OUVERTURE DU FICHIER xyznsef'
      ELSE
         KERR(1) = 'FILE xyznsef CAN NOT BE OPENED'
      ENDIF
      CALL LEREUR
      RETURN

C     DESTRUCTION du FICHIER
 9910 CLOSE( NF, STATUS='DELETE' )

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *, IERR,' EF NON des TRIANGLES'
         PRINT *, 'Le fichier ',KNMFIC,' N"" est PAS CREE'
      ELSE
         PRINT *, IERR,' FE NOT TRIANGLES'
         PRINT *, 'The file ',KNMFIC,' is NOT CREATED'
      ENDIF
      PRINT*

      RETURN
      END
