# External Integrations

**Analysis Date:** 2026-04-10

## Graphics & Display

**X11 Display Server:**
- Purpose: Interactive graphical visualization of 2D/3D meshes and simulation results
- Implementation: C-Fortran bridge in `xvue/xvuelc.c`
- Compiled object: `xvue/xvuelc.o` (linked into all solver executables)
- Protocol: Native X Window System (X11R6 or later)
- X Server requirement: User must have X11 server running (typically `Xvfb` for headless, native display for interactive)
- X11 Libraries:
  - `-lX11` (core X11 library from `/usr/X11R6/lib64/` or equivalent)
  - Headers bundled in `incl/`: `Xlib.h`, `Xatom.h`, `Xutil.h` (for compatibility across Linux distributions)
- Features implemented in `xvue/`:
  - Window creation and management
  - Color palette support (256 colors, 8-bit planes)
  - PostScript printing (via X11 graphics context)
  - Interactive event handling (mouse, keyboard)
  - Font loading and text rendering

**Related modules:**
- `xvue/xvouvrir.f` - Opens X11 display and window
- `xvue/xvinit.f` - Initializes X11 graphics context
- `xvue/palcde.f` - Color palette setup
- `xvue/videofin.f` - Animation/video output to GIF via external tool
- `xvue/video1.f` - Frame capture for animation

## Image & Animation Output

**ImageMagick (convert utility):**
- Purpose: Convert X11 window dumps (.xwd files) to standard image formats (PNG, GIF)
- Used in: `xvue/videofin.f`, `xvue/video1.f`, `util/trtable.f`
- Example workflow:
  ```fortran
  CALL SYSTEM('convert ' // NMFVIDEO(1:L) // '.xwd ' // NMFVIDEO(1:L) // '.png')
  ```
- Animation creation (commented in code): `convert -rotate 90 zfxy*.eps -delay 10 -extent 980x550 zfxy.gif`
- System call format: `SYSTEM()` (Fortran intrinsic executes shell commands)
- Install requirement: `ImageMagick` package for `convert` command

## Parallel Computing

**OpenMP (Open Multi-Processing):**
- Purpose: Shared-memory parallelization of compute-intensive loops
- Implementation: Fortran 95 variants (`.f95` files) with OpenMP directives
- Compiler flag: `-fopenmp` (gfortran)
- OpenMP library: `libgomp` (linked implicitly via gfortran)
- Directives used:
  - `!$OMP PARALLEL` / `!$OMP END PARALLEL` - Parallel regions
  - `!$OMP MASTER` / `!$OMP END MASTER` - Master thread only
  - `!$OMP BARRIER` - Synchronization point
  - Thread-local constructs for reduction operations

**Configuration:**
- Environment variable: `OMP_NUM_THREADS` (set by launcher scripts in `bin/`)
- Detection: `OMP_GET_NUM_THREADS()` at runtime (from `OMP_LIB` module)
- Parallelization thresholds in `incl/threads.inc`:
  - `MININDSTHR = 512` - Minimum loop iterations to enable parallelization
  - `MININD1THR = 128` - Minimum iterations per thread
- Example from `prpr/ppelas.f`:
  ```fortran
  !     USE OMP_LIB
        include"./incl/threads.inc"
        ...
        NBTHREADS = 1
  !$OMP PARALLEL
  !$OMP MASTER
  !$    NBTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
  ```

## Project Workflow & File I/O

**Project Directory Structure:**
- Location: `$MEFISTOX/ProjectName/` (user working directory)
- Initialized by: `bin/INITIER` launcher script
- Contains: Project-specific mesh files, input data, simulation results

**Data File Formats:**
- Mesh format: MEFISTO native (internal binary/text format)
- Input files: Fortran-based configuration (read via `OPEN`, `READ` statements)
- Output files: Results stored in project directory after each solver run
- Interactive project name: Passed via command line to launcher scripts

**Launcher Scripts (bin/):**
- `INITIER` - Initialize new project (creates `$MEFISTOX/$ProjectName/`)
- `MAILLER` - Interactive mesh generation (opens X11 window)
- `ELASTICER` - Run elasticity solver (links `pp/ppelas` executable)
- `FLUIDER` - Run fluid dynamics solver (links `pp/ppflui` executable)
- `THERMICER` - Run thermal solver (links `pp/ppther` executable)
- `NLSER` - Run nonlinear solver (links `pp/ppnlse` executable)
- `WAVER` - Run wave equation solver (if available)

**Environment Setup (via CLAUDE.md):**
```bash
export MEFISTO=/path/to/mefistosource   # Source root
export MEFISTOX=$HOME/mefistox          # User working directory
export PATH=.:$PATH:$MEFISTO/bin        # Adds launcher scripts to PATH
export CDPATH=.:$HOME:$MEFISTO:$MEFISTOX
```

## Fortran Runtime Features

**System Calls:**
- `CALL SYSTEM('command')` - Execute shell commands (used for `convert`, `pwd`, file operations)
- `CALL GETARG(N, ARGUMENT)` - Retrieve command-line arguments
- `CALL IARGC()` - Get count of command-line arguments

**I/O Subsystem:**
- Unit numbering: Managed via `incl/lu.inc` (logical unit definitions)
- COMMON block: `/UNITES/` defines LECTEU (reader), IMPRIM (printer), INTERA (interactive), NUNITE array (29 units)
- Format statements: Fortran 77 style with labeled FORMAT statements
- Interactive menu system: Via `incl/gsmenu.inc`

**Language Support:**
- Bilingual interface: French (default) and English
- Detection: File at `td/m/anglais` determines menu language
- Variables: `LANGAG` (language flag), `LANGUE` (language name string)

## Hardware Acceleration

**None detected:**
- No GPU acceleration (CUDA, OpenCL)
- No MPI (Message Passing Interface) for distributed computing
- Parallel support limited to shared-memory OpenMP threading on single machine

## Monitoring & Debugging

**Error Handling:**
- Custom error logging system via `incl/ppmck.inc`
- KERR array for error messages
- NBLGRC array for error code tracking
- Error function: `CALL LEREUR` (in various modules)

**Runtime Diagnostics:**
- Print statements for progress indication
- Common block MCN traces for memory allocation (`incl/a___trace.inc`)
- No external logging framework (built-in Fortran mechanisms only)

## Documentation & Help

**Included Documentation:**
- `doc/normes.ps` - Programming standards (PostScript format)
- `doca/` - English documentation (symlinked to `doc/`)
- `docf/` - French documentation (symlinked to `doc/`)
- `td/` - Tutorials and demo data:
  - `td/da/`, `td/df/` - Demo cases (English, French)
  - `td/ia/`, `td/if/` - Initialization cases (English, French)
  - `td/ma/`, `td/mf/` - Menu examples (English, French)

**Test Suites:**
- `testa/` - English test cases
- `testf/` - French test cases
- Used for validation after compilation

## External Dependencies Summary

| Component | Type | Required For | Location | Status |
|-----------|------|--------------|----------|--------|
| gfortran | Compiler | Building all modules | System (gcc-gfortran) | Essential |
| gcc | Compiler | Building X11 interface | System (gcc) | Essential |
| libX11 | Library | Graphics rendering | `/usr/X11R6/lib64/` | Essential |
| libX11-dev | Headers | Compiling X11 interface | `/usr/X11R6/include/` (bundled) | Essential |
| ImageMagick convert | Utility | Animated GIF output | System (imagemagick) | Optional* |
| OpenMP (libgomp) | Library | Parallelization | System (via gfortran) | Optional** |
| X11 Server | Runtime | Interactive display | User's display | Optional*** |

*Optional: GIF output only; all core functionality works without it
**Optional: Sequential execution possible with OMP_NUM_THREADS=1
***Optional: Batch mode via redirected output possible without X11 Server

---

*Integration audit: 2026-04-10*
