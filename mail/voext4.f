        SUBROUTINE VOEXT4( NBSA, NPSOM, NSENS, XYZA, XYZ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECOPIER LES COORDONNEES DES SOMMETS DE LA FACE
C -----    D'UN TETRAEDRE DANS L'ORDRE DE NUMEROTATION NATUREL
C
C ENTREES:
C --------
C NBSA   : NOMBRE DE SOMMETS SUR CHAQUE ARETE
C NPSOM  : NUMERO DU SOMMET COINCIDANT AVEC LE SOMMET 1 DE LA FACE
C NSENS  : SENS DE PARCOURS DE LA FACE
C XYZA   : COORDONNEES DES SOMMETS DE LA FACE AVANT RENUMEROTATION
C
C SORTIES:
C --------
C XYZ    : COORDONNEES DES SOMMETS DE LA FACE APRES RENUMEROTATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS            MARS 1989
C.......................................................................
      REAL XYZA(3,*),XYZ(3,*)
C
      IF ( NSENS .EQ. 1 ) THEN
         IF ( NPSOM .EQ. 1 ) THEN
            DO N1=1,NBSA
               DO N2=1,N1
                  N3 = N1
                  N4 = N2
                  NUM  = (N1*N1-N1)/2 + N2
                  NUMA = (N3*N3-N3)/2 + N4
                  DO J=1,3
                     XYZ(J,NUM) = XYZA(J,NUMA)
                  ENDDO
               ENDDO
            ENDDO
         ELSE IF ( NPSOM .EQ. 2 ) THEN
            DO N1=1,NBSA
               DO N2=1,N1
                  N3 = NBSA + 1 - N2
                  N4 = N1 - N2 + 1
                  NUM  = (N1*N1-N1)/2 + N2
                  NUMA = (N3*N3-N3)/2 + N4
                  DO J=1,3
                     XYZ(J,NUM) = XYZA(J,NUMA)
                  ENDDO
               ENDDO
            ENDDO
         ELSE IF ( NPSOM .EQ. 3 ) THEN
            DO N1=1,NBSA
               DO N2=1,N1
                  N3 = NBSA - N1 + N2
                  N4 = NBSA + 1 - N1
                  NUM  = (N1*N1-N1)/2 + N2
                  NUMA = (N3*N3-N3)/2 + N4
                  DO J=1,3
                     XYZ(J,NUM) = XYZA(J,NUMA)
                  ENDDO
               ENDDO
            ENDDO
         END IF
      ELSE
C        SENS RETROGRADE
         IF ( NPSOM .EQ. 1 ) THEN
            DO N1=1,NBSA
               DO N2=1,N1
                  N3 = N1 + 1 - N2
                  N4 = N1
                  NUM  = (N1*N1-N1)/2 + N2
                  NUMA = (N4*N4-N4)/2 + N3
                  DO J=1,3
                     XYZ(J,NUM) = XYZA(J,NUMA)
                  ENDDO
               ENDDO
            ENDDO
         ELSE IF ( NPSOM .EQ. 2 ) THEN
            DO N1=1,NBSA
               DO N2=1,N1
                  N3 = N2
                  N4 = NBSA +N2 - N1
                  NUM  = (N1*N1-N1)/2 + N2
                  NUMA = (N4*N4-N4)/2 + N3
                  DO J=1,3
                     XYZ(J,NUM) = XYZA(J,NUMA)
                  ENDDO
               ENDDO
            ENDDO
         ELSE IF ( NPSOM .EQ. 3 ) THEN
            DO N1=1,NBSA
               DO N2=1,N1
                  N3 = NBSA + 1 - N1
                  N4 = NBSA + 1 - N2
                  NUM  = (N1*N1-N1)/2 + N2
                  NUMA = (N4*N4-N4)/2 + N3
                  DO J=1,3
                     XYZ(J,NUM) = XYZA(J,NUMA)
                  ENDDO
               ENDDO
            ENDDO
         END IF
      END IF
C
      RETURN
      END
