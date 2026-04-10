# Technology Stack

**Analysis Date:** 2026-04-10

## Languages

**Primary:**
- Fortran 77 (fixed-form) - Core solver and mesh generation modules (`mail/`, `elas/`, `flui/`, `ther/`, `reso/`, `util/`, `xvue/`)
- Fortran 95 with OpenMP - Parallel variants of core modules (files with `.f95` extension for threading support)

**Secondary:**
- C (single file) - X11 graphics interface layer at `xvue/xvuelc.c`

## Runtime

**Environment:**
- Linux 64-bit (primary target: `-m64` compilation flag)
- Fortran runtime: `libgfortran`

**Build System:**
- Shell scripts (bash) - Primary build mechanism
- Manual compilation scripts in `bin/` directory
- Static library archiving with `ar` (GNU ar)

## Compilers

**Fortran:**
- `gfortran` (GNU Fortran) - Primary compiler for Fortran 77 and Fortran 95 code
- Compilation flags: `-O` (optimization), `-fopenmp` (OpenMP threading), `-m64` (64-bit), `-mcmodel=large` (large memory model), `-fPIC` (position-independent code), `-Wall` (all warnings)

**C:**
- `gcc` (GNU C Compiler) - For X11 graphics interface
- Compilation flags: `-O` (optimization), `-m64` (64-bit), `-openmp` (OpenMP), `-fPIC` (position-independent code)

## Build System Workflow

**Entry Point:** `bin/cbl_tout` (compile everything)

**Sequence:**
1. Generate `incl/homdir.inc` - Encodes `$MEFISTO` path as Fortran PARAMETER (dynamically created per installation)
2. Compile C graphics layer: `bin/ccxvue` → `xvue/xvuelc.o`
3. Compile all module libraries using `bin/cbl_tous_f`:
   - `elas/lib` - Elasticity solver
   - `flui/lib` - Fluid dynamics solver
   - `ther/lib` - Thermal solver
   - `mail/lib` - Mesh generation
   - `reso/lib` - Linear solver
   - `util/lib` - Shared utilities
   - `xvue/lib` - X11 graphics
4. Link individual module executables: `bin/cbinit`, `bin/cbmail`, `bin/cbelas`, `bin/cbflui`, `bin/cbther`, `bin/cbnlse`, `bin/cbpoba`
5. Output: executables in `pp/` directory

**Module-specific build scripts:**
- `bin/cbinit` - Compiles `prpr/ppinit.f` → `pp/ppinit`
- `bin/cbmail` - Compiles `prpr/ppmail.f` → `pp/ppmail`
- `bin/cbelas` - Compiles `prpr/ppelas.f` → `pp/ppelas`
- `bin/cbflui` - Compiles `prpr/ppflui.f` → `pp/ppflui`
- `bin/cbther` - Compiles `prpr/ppther.f` → `pp/ppther`
- `bin/cbnlse` - Compiles nonlinear solver → `pp/ppnlse`
- `bin/cbpoba` - Compiles post-processing tool → `pp/pppoba`

## External Libraries

**X11 Graphics:**
- System library: `/usr/X11R6/lib64/libX11.a` (or `/usr/X11R6/lib/libX11.a` on 32-bit)
- Required header `<X11/Xlib.h>` (bundled in `incl/Xlib.h` for portability)
- Other X11 headers in `incl/`: `Xatom.h`, `Xutil.h`
- Link flag: `-lX11`

**OpenMP Threading:**
- GNU OpenMP (`libgomp`) - Implicit linking via `gfortran -fopenmp`
- Thread parameters defined in `incl/threads.inc`
- Minimum threshold for parallelization: 512 iterations minimum per thread
- Environment variable: `OMP_NUM_THREADS` (user-configurable at runtime)

**C Standard Libraries:**
- `libc` - Standard C runtime
- Includes: `limits.h`, `unistd.h`, `stdio.h`, `string.h`, `stdlib.h`, `math.h`, `time.h`, `setjmp.h`, `ctype.h`, `sys/types.h`, `sys/time.h`

## Configuration Files

**Build Configuration:**
- `incl/pp.inc` - Memory allocation configuration (MOTMCN parameter sets maximum memory pool)
  - Configurable from 28,000,000 to 2,130,000,000 words (64-bit limit due to 32-bit integer addressing)
  - Current value: 1,532,000,000 words (approximately 6 GB when accounting for 4-byte words)
- `incl/threads.inc` - OpenMP thread configuration (MININDSTHR, MININD1THR, NBTHREADS)
- `incl/homdir.inc` - **Generated at build time** - Contains MEFISTO source directory path

**Runtime Configuration:**
- Location: `incl/` directory contains shared include files for all modules
- Key includes define:
  - Common block structures (MCN - shared memory pool)
  - Physical constants and problem parameters
  - Language settings (FRANCAISE/ANGLAISE) via `incl/langue.inc`
  - Menu system via `incl/gsmenu.inc`
  - Type definitions via `incl/typnoobj.inc`

## Environment Variables (Required)

**Build Time:**
- `MEFISTO` - Source directory root (used by all `bin/cb*` scripts)
- `$MEFISTO` is embedded in generated `incl/homdir.inc` as Fortran PARAMETER

**Runtime:**
- `MEFISTOX` - User project working directory (where simulations are run)
- `PATH` - Must include `$MEFISTO/bin` for launcher scripts (INITIER, MAILLER, ELASTICER, etc.)
- `CDPATH` - Optional but recommended for convenience (includes `$MEFISTO` and `$MEFISTOX`)
- `OMP_NUM_THREADS` - OpenMP thread count (exported by launcher scripts like `bin/ELASTICER`)

## Platform Requirements

**Development (Build):**
- Linux 64-bit OS
- `gfortran` (Fortran 77/95 compiler with OpenMP)
- `gcc` (C compiler)
- `libX11-dev` (X11 development headers and libraries)
- `ar` (GNU archive tool)
- `bash` shell
- System paths: `/usr/X11R6/lib64/` or `/usr/X11R6/lib/` for X11 libraries

**Production (Runtime):**
- Linux 64-bit OS
- X11 server (for interactive graphical display)
- `libgfortran` runtime library (from gfortran)
- `libX11` runtime library
- User project directory structure under `$MEFISTOX/ProjectName/`

## Memory Model

**Virtual Memory:**
- Large memory model (`-mcmodel=large`) - Supports up to 2 TB address space
- 64-bit pointer support via `-m64` flag
- Maximum workable memory: ~6 GB per process (limited by 32-bit integer indexing in MCN array)

**Memory Allocation:**
- Static global memory pool (MCN common block) sized at build time via `incl/pp.inc`
- No dynamic allocation via Fortran 90 allocatable arrays in main solvers
- Traditional Fortran 77 approach: all arrays allocated from static MCN pool at startup

## Legacy Notes

**Previous Build System References (not used):**
- `Makefile`, `MakefileIBM`, `MakefileMefisto` in `bin/` - Historical, superseded by bash scripts
- `bin.lnx64/` directory - Prebuilt Linux 64-bit binaries and legacy scripts
- Comments reference obsolete compilers: `g77`, `f77`, `HP`, `SUN`, `SGI` machines

---

*Stack analysis: 2026-04-10*
