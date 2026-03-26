       subroutine PetscIni()
       implicit none       
#include "finclude/petsc.h"

       PetscErrorCode ierr
       call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
       end
