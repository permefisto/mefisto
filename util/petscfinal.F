       subroutine PetscFinal()
       implicit none       
#include "finclude/petsc.h"

       PetscErrorCode ierr
       call PetscFinalize(ierr)
       end
