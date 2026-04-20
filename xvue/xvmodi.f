      SUBROUTINE XVUE_MODULE_INIT( NAME, NAME_LEN )
!GCC$ ATTRIBUTES WEAK :: XVUE_MODULE_INIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     HOOK D'INITIALISATION PAR MODULE (mail, elas, flui, ther,
C -----     nlse) POUR LA COUCHE GRAPHIQUE X11. NO-OP STUB.
C
C           The Qt backend (xvue/qt/src/xvue_qt_api.cpp) provides the
C           real (strong) body; the X11 backend has no QAction menu
C           system so this stub simply returns. The CALL site in
C           prpr/pp*.f is identical for both backends -- only the
C           linked library differs.
C
C           IMPORTANT: the GCC ATTRIBUTES WEAK directive marks this
C           symbol as weak so the Qt-side strong definition in
C           xvue/qt/build/libxvueqt.a wins at link time when both
C           xvue/lib AND libxvueqt.a are linked together (Qt build).
C           The X11 build links only xvue/lib so this weak symbol is
C           the one that resolves and runs as a no-op.
C
C ENTREES :
C ---------
C NAME     : module name (ignored in X11 build)
C NAME_LEN : length of NAME (ignored in X11 build)
C
C SORTIES :  aucune
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  AVRIL 2026
C2345X7..............................................................012
      CHARACTER*(*)  NAME
      INTEGER        NAME_LEN
C
C     Reference the args once so -Wall does not complain under -Wunused.
      IF( .FALSE. ) THEN
         WRITE(*,*) NAME, NAME_LEN
      ENDIF
C
      RETURN
      END
