# LEXICON-AUDIT-flui.md

**Phase:** 06.3-fluid-flui-menu-wiring
**Generated:** 2026-04-27
**Status:** FROZEN (ratified by user 2026-04-27)
**Requirement:** UX-05 (flui slice) Success Criterion #1

## Scope

Full recursive LIMTCL tree walk from `td/m/debuflui` via `flui/` + `util/` +
`reso/` + `prpr/ppflui.f` call-sites. **42 sub-menus** (re-verified
2026-04-27 via `grep -rEoh "LIMTCL\( *'[^']+'" flui/ util/ prpr/ppflui.f reso/`
— 0 drift from planner's interfaces list). De-duplication rule from 6.1
Plan 01 / 6.2 Plan 01:

- **Flui-specific sub-menus fully expanded** — `cl_flui`, `methrecr`,
  `nbvitpre`, `noxyzmod`, `particle`, `tracmovi`, `tracoura`, `tracpr2d`,
  `tracpr3d`, `tractem2`, `tractem3`, `traerrpr`, `traerrvi`, `trarot2d`,
  `trarot3d`, `trflevit`, `vitepr2d`, `vitepr3d`, `vitetype`, `vp0flui`
  (20 sub-menus, full leaf expansion, 123 leaves total).
- **Shared util sub-menus listed ONCE** with notes cell citing 6.1 mail
  / 6.2 elas audits for leaf-level detail — `affiche`, `couleur0`,
  `couleurs`, `entrsort`, `fichier`, `lecteur`, `managfen`, `managmef`,
  `methvvpr`, `mode_es`, `opt_lign`, `opt_surf`, `sectvopl`, `suivi_ms`,
  `suivitms`, `tracaxes`, `tuer`, `typ_objt`, `typtrait`, `vuesplan`,
  `zeros` (21 sub-menus, compressed to 1 row each).

**Total row count: 159** — within validator bound `[80, 250]` (same as
6.1/6.2). 15 debuflui top-level + 123 flui-unique leaves + 21 shared
util compressed.

Frequency bucketing from `testa/*/*.{stoke56cr, ns7cr, ns8cr, stokes6cr,
ns7cg, ns8cg, t3060, t030}` grep evidence (13 fixtures across cavity2d,
cavity3d, marche, nstg3d, pan2d). Top-level evidence (comment markers):
`1;`=2, `2;`=4, `3;`=2, `5;`=3, `6;`=2, `7;`=3, `8;`=3, `9;`=2, `10;`=7,
`20;`=6 — supplemented by domain review for sparse signals (small flui
testa corpus). HIGH/MED/LOW rubric: HIGH ≥5, MED 2-4, LOW 0-1; HIGH set
domain-promoted to cover canonical Navier-Stokes workflow (1; object,
2; physical data, 5; or 7; solver, 10; visual, 99; save). Per CONTEXT
carryover from 6.1/6.2: rows with `qaction=yes` have `frequency` in
{high, med}; exactly 5 rows have `toolbar=yes`.

Schema: 9 columns, enforced by `tools/validate_audit_md.py` (shipped in
6.1 Plan 01; module-aware ICONS resolver since 6.2 Plan 02 — reused
verbatim, no tooling changes). Rule 9 runs in WARN mode for new flui
SVGs (Plan 02 ships them); `mesh-draw.svg` resolves cleanly into
`icons/mail/` (no WARN — already shipped by 6.1).

### Menu taxonomy (proposed — user ratifies at Task 2)

The flui module's Qt menu bar uses **`{File, Fluid, View, Help}`** per
ROADMAP Phase 6.3 line 201 + 06.0-UI-SPEC §Per-Module Conformance
Contract (flui `<Module>` menu title is `Fluid / Fluide`). This differs
from 6.1's `{File, Mesh, View, Help}` and 6.2's `{File, Solve, View,
Help}` — "Fluid" replaces "Solve" / "Mesh".

Content distribution (proposed):

- **File**: codes `70;` `99;` + shared 6.0 File actions
- **Fluid**: codes `1;` `2;` `3;` `5;` `6;` `7;` `8;` `9;` + a `Parameters`
  sub-menu with codes `20;` `60;`
- **View**: codes `10;` `19;` + shared 6.0 View actions
- **Help**: code `97;` + shared 6.0 Help/About

## Known data quirks

- `prpr/ppflui.f` GOTO table at line 241 dispatches codes 1..10 plus
  hidden codes 19, 20, 35, 36, 60, 70, 97, 99. Codes `35;` (NS2DTGV
  Taylor-Green 2D) and `36;` (NS3DTGV Taylor-Green 3D) ARE in the GOTO
  table at lines 235-236 but NOT in `td/m/debuflui` — treated as hidden
  lexicons with frequency=low and qaction=no unless user review at
  Task 2 promotes them. NS3DTGV is invoked by `nstg3d` testa fixtures
  through chained scripts (not as a top-level interactive command).
- Code `4` is absent from both `td/m/debuflui` AND the GOTO table
  (slot 4 in `(100, 200, 300, 30, 500, 600, 700, 800, 900, 1000)` is
  `30` — a no-op input-loop reset). NO `4;` row.
- `prpr/ppflui.f` lines 379-394 contain commented-out `ccc` blocks for
  PISO+Charac (NAVSTO=4) and PISO+(V.D)V (NAVSTO=5) — alternate
  algorithms for codes 7;,8;. Dead code; NO audit rows added.
- `flui/fluidens.f:1460` `CALL LIMTCL( 'vp0flui', ...)` is also
  commented `ccc`. The `td/m/vp0flui` file exists with 3 leaves but
  `vp0flui` is NEVER reached at runtime (dead lexicon). Audit lists
  vp0flui rows with frequency=low and notes cell flagging the dead
  call site — Plan 02 may exclude vp0flui from QAction set entirely;
  user can ratify exclusion at Task 2.
- `cl_flui` is reached from `flui/dfflui.f:269` under code `2;` (DFFLUI
  prompts physical-data boundary conditions). It is NOT also reached
  from DFTHER (code 3;) — DFTHER does NOT call `LIMTCL('cl_flui')`
  directly. Canonical parent prefix = `2;` (no co-parent ambiguity,
  unlike 6.2 elas's contdepl which had two parents).
- `methrecr` is reached from `flui/chmereso.f:42` which is called by
  `STOKESTA` (code 5;) and `FLUIDENS` (codes 6;, 7;, 8;, 9;). Used
  canonical parent `5;` (steady Stokes — primary first-solver path);
  cited co-reach via 6;/7;/8;/9; in `notes`. Synthetic prefix `5;81;`
  used to avoid collision with `cl_flui` (which sits under 2; not 5;)
  and to follow 6.2 D-04 high-number-offset convention.
- `nbvitpre`, `vitepr2d`, `vitepr3d`, `tracpr2d`, `tracpr3d`, `tractem2`,
  `tractem3`, `traerrpr`, `traerrvi`, `trarot2d`, `trarot3d`, `trflevit`,
  `tracmovi`, `tracoura`, `vitetype`, `noxyzmod`, `particle` are all
  reached from `flui/tfluide.f` which is called by ppflui code `10;`
  (TFLUIDE — DRAWING of VELOCITIES PRESSURES). Canonical parent prefix
  for all visualization sub-menus = `10;`. Synthetic high-number suffix
  used per leaf to keep path uniqueness (e.g. `10;1;` for nbvitpre,
  `10;20;` for vitepr2d, `10;30;` for vitepr3d, etc.). Hierarchy
  documented in `notes`.
- `td/m/vp0flui` line 2 has unbalanced quote `' DEPART avec le VECTEUR
  VITESSE+PRESSION deja calcule au TEMPS INITIAL}` (closing brace `}`
  instead of closing quote `'`). The LIMTCL parser tolerates this in
  practice (no runtime crash); audit honours the FR text as written
  and flags in `notes`. Apply 6.1 Pitfall 5 log-and-fallback discipline.
- `td/m/vitepr2d` and `td/m/vitepr3d` use FR-only text `Mi Mx`
  (Min/Max abbrev) which has no exact EN cognate; td/ma/* preserves
  the same shorthand. Treated as bilingual carry-over.

## Legend

- `lexicon_path` — single-line semicolon-separated, no spaces
- `frequency` — `high` / `med` / `low` (from testa/ counts + domain review)
- `qaction` — `yes` iff frequency is high or med
- `toolbar` — `yes` iff in the top-5 toolbar slice (exactly 5 rows)
- `icon_source` — Qt `SP_*` enum name, custom `.svg` filename, or `—`
- `shortcut` — Qt accelerator or `—`
- `notes` — clarifications, deferred items, source flags

## Audit Table

| lexicon_path | description_fr | description_en | frequency | qaction | toolbar | icon_source | shortcut | notes |
|--------------|----------------|----------------|-----------|---------|---------|-------------|----------|-------|
| 1; | NOM de l'OBJET a traiter | OBJECT NAME to treat | high | yes | no | flui-object.svg | — | debuflui leaf; testa count=2 (cavity2d/cavity3d set OBJECT first); domain-promoted HIGH (mandatory first command per ppflui.f:243-318); Fluid menu root item |
| 2; | ENTRER les DONNEES PHYSIQUES du FLUIDE | GIVE the INPUT PHYSICAL DATA of FLUID | high | yes | yes | flui-input-physics.svg | — | debuflui leaf; testa count=4 (every fluid project sets density+viscosity); dispatches to DFFLUI which prompts cl_flui sub-menu; toolbar top-5 (Input data setup) |
| 3; | ENTRER les DONNEES THERMIQUES du FLUIDE | GIVE the INPUT HEAT DATA of FLUID | med | yes | no | flui-input-heat.svg | — | debuflui leaf; testa count=2 (only Boussinesq fixtures); dispatches to DFTHER (heat boundary conditions) |
| 5; | RESOUDRE le PROBLEME de STOKES STATIONNAIRE | INCOMPRESSIBLE STEADY STOKES SOLVER | high | yes | yes | solve-stokes.svg | — | debuflui leaf; testa count=3 (cavity2d/cavity3d steady Stokes baseline); dispatches to STOKESTA (which prompts methrecr); domain-promoted HIGH (canonical first solver run); toolbar top-5 |
| 6; | RESOUDRE le PROBLEME de STOKES INSTATIONNAIRE | INCOMPRESSIBLE UNSTEADY STOKES SOLVER | med | yes | no | solve-stokes.svg | — | debuflui leaf; testa count=2 (cavity2d/nstg3d unsteady Stokes); dispatches to FLUIDENS(NAVSTO=0); icon REUSED — same flow type as 5; |
| 7; | RESOUDRE NAVIER-STOKES INSTATIONNAIRE IMPLICITE | IMPLICIT UNSTEADY NAVIER-STOKES SOLVER | high | yes | yes | solve-navier-stokes.svg | — | debuflui leaf; testa count=3 (marche/nstg3d implicit NS production); dispatches to FLUIDENS(NAVSTO=1); domain-promoted HIGH (primary production solver); toolbar top-5 |
| 8; | RESOUDRE NAVIER-STOKES INSTATIONNAIRE FRACTIONNAIRE | FRACTIONAL UNSTEADY NAVIER-STOKES SOLVER | med | yes | no | solve-navier-stokes.svg | — | debuflui leaf; testa count=3 (marche/nstg3d fractional NS variant); dispatches to FLUIDENS(NAVSTO=2); icon REUSED — same NS family as 7; |
| 9; | RESOUDRE NAVIER-STOKES THERMIQUE BOUSSINESQ INSTATIONNAIRE | UNSTEADY BOUSSINESQ HEAT NAVIER-STOKES SOLVER | med | yes | no | solve-navier-stokes.svg | — | debuflui leaf; testa count=2 (Boussinesq HEAT NS variant); dispatches to FLUIDENS(NAVSTO=3); icon REUSED — NS family |
| 10; | DESSINER des VITESSES PRESSIONS TEMPERATURES | DRAWING of VELOCITIES PRESSURES TEMPERATURES | high | yes | yes | draw-velocity.svg | — | debuflui leaf; testa count=7 (every result-producing project draws results); dispatches to TFLUIDE (which prompts the entire visualization tree: nbvitpre, vitepr2d/3d, tracpr2d/3d, tractem2/3, traerrpr/vi, trarot2d/3d, trflevit, tracmovi, tracoura, vitetype, noxyzmod, particle); toolbar top-5 (primary result visual) |
| 19; | TRACER des MAILLAGES des PLSVO | DRAWING of PLSVO MESHES | med | yes | no | mesh-draw.svg | — | debuflui leaf; testa count=0 (no comment-marker hits but ppflui.f:1900 dispatches TRMAIL); domain-promoted MED (mesh inspection routine); REUSES 6.1 mail/mesh-draw.svg (no copy needed) |
| 20; | PRECISION pour INVERSER A x = b | PRECISION to SOLVE A x = b | med | yes | no | — | — | debuflui leaf; testa count=6 (every CG-based fixture sets precision); dispatches to ZEROGC shared util (zeros compressed row below); Parameters submenu of Fluid |
| 35; | TABLE des RESULTATS 2D du TOURBILLON de TAYLOR-GREEN | 2D TAYLOR-GREEN VORTEX RESULTS TABLE | low | no | no | — | — | HIDDEN code — NOT in td/m/debuflui menu but live in ppflui.f:235 GOTO table; dispatches to NS2DTGV; reachable via typed `35;` only |
| 36; | TABLE des RESULTATS 3D du TOURBILLON de TAYLOR-GREEN | 3D TAYLOR-GREEN VORTEX RESULTS TABLE | low | no | no | — | — | HIDDEN code — NOT in td/m/debuflui menu but live in ppflui.f:236 GOTO table; dispatches to NS3DTGV; reachable via typed `36;` only; nstg3d testa fixtures use this through chained scripts |
| 60; | GERER la FENETRE GRAPHIQUE de Mefisto | MANAGE the Mefisto GRAPHIC WINDOW | low | no | no | SP_ComputerIcon | — | debuflui leaf; dispatches to MANAGFEN shared util (managfen compressed row); Parameters submenu of Fluid |
| 70; | GERER les TMS Fichiers Unites de Mefisto | MANAGE the Mefisto TMS Files Units | low | no | no | SP_DirIcon | — | debuflui leaf; dispatches to MANAGMEF shared util (managmef compressed row); File menu |
| 97; | AFFICHER la DATE de la version de Mefisto | Mefisto VERSION NAME | low | no | no | SP_MessageBoxInformation | — | debuflui leaf; dispatches to VRSION; Help menu item |
| 99; | SAUVEGARDER les DONNEES et FIN du TRAITEMENT | SAVE DATA and QUIT | high | yes | yes | SP_DialogCloseButton | Ctrl+Q | debuflui leaf; domain-promoted HIGH (every session closes here); toolbar top-5 (shared convention with 6.1/6.2); File menu |
| 2;1; | FORCE Force externe sur la frontiere | BOUNDARY EXTERIOR FORCE | med | yes | no | — | — | cl_flui leaf; reached from DFFLUI under debuflui code 2;; cavity2d uses `1;` (Mass Density via WTBLVI) — note: cavity2d FORCE-IMPOSED line uses code 1 in the typed lexicon for mass density during DFFLUI's outer prompt loop, distinct from cl_flui inner loop |
| 2;2; | BLocage de la VITESSE imposee | BLOCKING-UP the IMPOSED VELOCITY | high | yes | no | — | — | cl_flui leaf; testa count=many (cavity2d/cavity3d/marche/nstg3d all fix velocity with `2;`); canonical boundary condition — domain HIGH despite indented status |
| 2;3; | BLocage de la PRESSION imposee | BLOCKING-UP the IMPOSED PRESSURE | med | yes | no | — | — | cl_flui leaf; pressure boundary condition — used in marche outflow |
| 2;4; | Vitesses initiales | INITIAL VELOCITIES | low | no | no | — | — | cl_flui leaf; initial condition for unsteady NS solve |
| 2;5; | Pressions initiales | INITIAL PRESSURES | low | no | no | — | — | cl_flui leaf; initial condition for unsteady NS solve |
| 5;81;1; | FACTORISATION CROUT avec stockage PROFIL | FACTORIZATION CROUT with SKYLINE STORAGE | med | yes | no | — | — | methrecr leaf — synthetic prefix 5;81; to avoid cl_flui collision; reached from CHMERESO under codes 5; (STOKESTA) AND 6;/7;/8;/9; (FLUIDENS); canonical Crout factorization |
| 5;81;2; | GRADIENT CONJUGUE avec stockege CONDENSE | CONJUGATE GRADIENT with CONDENSED STORAGE | med | yes | no | — | — | methrecr leaf; CG with incomplete Crout preconditioner |
| 5;81;3; | DOUBLE GRADIENT ACCELERE CGS | CONJUGATE GRADIENT SQUARED CGS | low | no | no | — | — | methrecr leaf; CGS variant |
| 5;81;4; | DOUBLE GRADIENT STABILISE BiCGStab | Bi-CONJUGATE GRADIENT Stabilised BiCGStab | low | no | no | — | — | methrecr leaf; BiCGStab variant |
| 10;1;1; | Numeros INITIAL et FINAL des VECTEURS SOLUTIONS a TRACER | INITIAL and FINAL NUMBER of SOLUTIONS VECTORS to DRAW | med | yes | no | — | — | nbvitpre leaf — synthetic prefix 10;1; (nbvitpre is the FIRST prompt under TFLUIDE); reached from flui/tfluide.f:725 under debuflui code 10;; canonical "by index" branch |
| 10;1;2; | TEMPS INITIAL et FINAL des VECTEURS SOLUTIONS a TRACER | INITIAL and FINAL TIME of SOLUTIONS VECTORS to DRAW | med | yes | no | — | — | nbvitpre leaf; "by time" branch |
| 10;20;1; | Traces de la PRESSION | The PRESSURE drawings | med | yes | no | — | — | vitepr2d leaf — synthetic prefix 10;20; for 2D variant; reached from flui/tfluide.f:1507; dispatches to tracpr2d sub-menu |
| 10;20;2; | Traces de la VITESSE sous FORME de FLECHES | The VELOCITY is drawn with ARROWS | high | yes | no | — | — | vitepr2d leaf; testa count=high (every visualization fixture draws velocity arrows via 90; under trflevit); dispatches to trflevit sub-menu |
| 10;20;3; | Traces de la COMPOSANTE X Y Z ou MODULE de la VITESSE | One VELOCITY COMPONENT X Y Z or the VELOCITY MAGNITUDE | med | yes | no | — | — | vitepr2d leaf; dispatches to noxyzmod -> tracmovi sub-menu chain |
| 10;20;4; | VITESSE Moy Mx PRESSION Mi Mx TEMPERATURE Moy Mi Mx Temps | VELOCITY Mean Max PRESSURE Max-Min TEMPERATURE Mean Time | low | no | no | — | — | vitepr2d leaf; printing summary statistics |
| 10;20;5; | Traces de la TEMPERATURE | The TEMPERATURE drawings | med | yes | no | — | — | vitepr2d leaf; dispatches to tractem2 sub-menu (only meaningful for code 9; Boussinesq) |
| 10;20;6; | Traces des FLUX de TEMPERATURE | The TEMPERATURE FLUX drawings | low | no | no | — | — | vitepr2d leaf; thermal flux variant |
| 10;20;7; | Calcul et Trace du PARCOURS de PARTICULES dans le FLUIDE | COMPUTE and DRAW the RUN of PARTICLES in the FLUID | low | no | no | — | — | vitepr2d leaf; dispatches to particle sub-menu |
| 10;20;9; | FONCTION COURANT Psi ou Vx Vy de dPsi sur dy | The STREAM FUNCTION Psi | med | yes | no | — | — | vitepr2d leaf; 2D-only stream function; dispatches to tracoura sub-menu |
| 10;20;10; | VORTICITE TOURBILLON ROTATIONNEL dVy sur dx moins dVx sur dy | The VORTICITY Rotational Velocity 2D | med | yes | no | — | — | vitepr2d leaf; 2D vorticity; dispatches to trarot2d sub-menu |
| 10;20;90; | Affichage VITESSE et PRESSION et ERREUR aux Noeuds | The PRINTING VELOCITY PRESSURE ERRORS at Nodes | med | yes | no | — | — | vitepr2d leaf; canonical print-summary trigger |
| 10;20;91; | Trace ERREURS VITESSE si Vitesse exacte connue | The VELOCITY ERRORS from Velocity Exact GIVEN | low | no | no | — | — | vitepr2d leaf; dispatches to traerrvi sub-menu |
| 10;20;92; | Trace ERREURS PRESSION si Pression Exacte connue | The PRESSURE ERRORS from Pressure Exact is GIVEN | low | no | no | — | — | vitepr2d leaf; dispatches to traerrpr sub-menu |
| 10;30;1; | Traces de la PRESSION | The PRESSURE drawings | med | yes | no | — | — | vitepr3d leaf — synthetic prefix 10;30; for 3D variant; reached from flui/tfluide.f:1509; dispatches to tracpr3d sub-menu |
| 10;30;2; | Traces de la VITESSE sous FORME de FLECHES | The VELOCITY is drawn with ARROWS | high | yes | no | — | — | vitepr3d leaf; dispatches to trflevit sub-menu |
| 10;30;3; | Traces de la COMPOSANTE X Y Z ou MODULE de la VITESSE | One VELOCITY COMPONENT X Y Z or VELOCITY MAGNITUDE | med | yes | no | — | — | vitepr3d leaf; dispatches to noxyzmod -> tracmovi |
| 10;30;4; | VITESSE Moy Mx PRESSION Mi Mx TEMPERATURE Moy Mi Mx Temps | VELOCITY Mean Max PRESSURE Min Max TEMPERATURE Mean Time | low | no | no | — | — | vitepr3d leaf; statistics summary |
| 10;30;5; | Traces de la TEMPERATURE | The TEMPERATURE drawings | med | yes | no | — | — | vitepr3d leaf; dispatches to tractem3 |
| 10;30;6; | Traces des FLUX de TEMPERATURE | The TEMPERATURE FLUX drawings | low | no | no | — | — | vitepr3d leaf; thermal flux 3D |
| 10;30;7; | Calcul et Trace du PARCOURS de PARTICULES dans le FLUIDE | COMPUTE and DRAW the RUN of PARTICLES in the FLUID | low | no | no | — | — | vitepr3d leaf; dispatches to particle |
| 10;30;10; | VORTICITE TOURBILLON ROTATIONNEL de la VITESSE | The VORTICITY Rotational Velocity 3D | med | yes | no | — | — | vitepr3d leaf; dispatches to trarot3d |
| 10;30;12; | INTEGRALE de la PRESSION sur les SURFACES de l OBJET 3D | The INTEGRAL of PRESSURE on OBJECT SURFACES | low | no | no | — | — | vitepr3d leaf; 3D-only pressure integral |
| 10;30;13; | FLUX NORMAL de la VITESSE sur les SURFACES de l OBJET 3D | The VELOCITY NORMAL FLUX on OBJECT SURFACES | low | no | no | — | — | vitepr3d leaf; 3D-only velocity flux |
| 10;30;90; | Affichage VITESSE et PRESSION et ERREUR aux Noeuds | The PRINTING VELOCITY PRESSURE ERRORS at Nodes | med | yes | no | — | — | vitepr3d leaf; canonical print-summary trigger |
| 10;30;91; | Trace ERREURS VITESSE si Vitesse exacte connue | The VELOCITY ERRORS from Velocity Exact GIVEN | low | no | no | — | — | vitepr3d leaf; dispatches to traerrvi |
| 10;30;92; | Trace ERREURS PRESSION si Pression Exacte connue | The PRESSURE ERRORS from Pressure Exact GIVEN | low | no | no | — | — | vitepr3d leaf; dispatches to traerrpr |
| 10;40;1; | COURBES ISO-PRESSIONS | ISO-PRESSURE CURVES | med | yes | no | — | — | tracpr2d leaf — synthetic prefix 10;40; for 2D pressure trace; reached from flui/tfluide.f:1541; |
| 10;40;2; | SURFACES COULEURS ISO-PRESSIONS | ISO-PRESSURES COLOR ZONES | med | yes | no | — | — | tracpr2d leaf; canonical 2D color rendering |
| 10;40;9; | En 2D SURFACE X Y PRESSION X Y | In 2D SURFACE X Y PRESSURE X Y | low | no | no | — | — | tracpr2d leaf; 3D-perspective lift |
| 10;50;1; | SURFACES ISO-PRESSIONS en COULEURS | ISO-PRESSURE SURFACES | med | yes | no | — | — | tracpr3d leaf — synthetic prefix 10;50; for 3D pressure trace; reached from flui/tfluide.f:1543; |
| 10;50;2; | ZONES COULEURS ISO-PRESSIONS | ISO-PRESSURES COLOR ZONES | med | yes | no | — | — | tracpr3d leaf |
| 10;50;3; | ZONES par SECTIONS X ou Y ou Z=CTE | CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | tracpr3d leaf; dispatches to sectvopl |
| 10;50;4; | PROFILS dans PLANS X ou Y ou Z=CTE | PROFILE on CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | tracpr3d leaf |
| 10;50;5; | LE LONG d'une DROITE de 2 POINTS | Pressure ALONG a LINE DEFINED by 2 POINTS | low | no | no | — | — | tracpr3d leaf |
| 10;50;6; | PRESSION sur des SURFACES de l'OBJET | Pressure on a SURFACE or more of the OBJECT | low | no | no | — | — | tracpr3d leaf |
| 10;50;7; | PRESSION sur les SURFACES FRONTALIERES 3D | Pressure on 3D BOUNDARY SURFACES | low | no | no | — | — | tracpr3d leaf |
| 10;50;8; | Calcule INTEGRALE sur SURFACES de l'OBJET | Compute Pressure INTEGRAL on OBJECT SURFACES | low | no | no | — | — | tracpr3d leaf |
| 10;55;1; | ISOTHERMES | ISOTHERMAL LINES | low | no | no | — | — | tractem2 leaf — synthetic prefix 10;55; for 2D temperature trace; reached from flui/tfluide.f:2105 (only when code 9; Boussinesq runs) |
| 10;55;2; | ZONES COULEURS ISOTHERMES | ZONES of ISOTHERMAL COLORS | low | no | no | — | — | tractem2 leaf |
| 10;55;7; | SURFACE X Y TEMPERATURE X Y | SURFACE X Y TEMPERATURE X Y | low | no | no | — | — | tractem2 leaf; 3D-perspective lift |
| 10;55;8; | ERREUR si TEMPERATURE EXACTE | ERROR from EXACT TEMPERATURE | low | no | no | — | — | tractem2 leaf; error mode |
| 10;57;1; | ISO-SOLUTIONS | ISO-SOLUTION SURFACES | low | no | no | — | — | tractem3 leaf — synthetic prefix 10;57; for 3D temperature trace; reached from flui/tfluide.f:2107 (Boussinesq 3D) |
| 10;57;2; | ZONES COULEURS ISO-SOLUTIONS | ZONES of ISO-SOLUTION COLORS | low | no | no | — | — | tractem3 leaf |
| 10;57;3; | ZONES par SECTIONS X ou Y ou Z=CTE | ZONES by CUTTING PLANES X or Y or Z=Const | low | no | no | — | — | tractem3 leaf; dispatches to sectvopl |
| 10;57;4; | PROFILS dans PLANS X ou Y ou Z=CTE | PROFILES by CUTTING PLANES X or Y or Z=Const | low | no | no | — | — | tractem3 leaf |
| 10;57;5; | LE LONG d'une DROITE def 2 POINTS | ALONG A LINE DEFINED BY 2 POINTS | low | no | no | — | — | tractem3 leaf |
| 10;57;8; | ERREUR ABSOLUE avec SOLUTION EXACTE | ABSOLUTE ERROR from EXACT SOLUTION GIVEN | low | no | no | — | — | tractem3 leaf |
| 10;60;1; | ISO-ERREURS sur la PRESSION | ISO-PRESSURE ERROR SURFACES or CURVES | low | no | no | — | — | traerrpr leaf — synthetic prefix 10;60; for pressure-error trace; reached from flui/tfluide.f:3209 |
| 10;60;2; | ZONES de COULEURS des ERREURS sur la PRESSION | COLOR ISO-PRESSURE ERROR on BOUNDARY | low | no | no | — | — | traerrpr leaf |
| 10;60;3; | ZONES des ERREURS par SECTIONS X ou Y ou Z=CTE | ERROR CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | traerrpr leaf |
| 10;60;4; | PROFILS des ERREURS dans PLANS X ou Y ou Z=CTE | ERROR PROFILE on CROSS SECTION | low | no | no | — | — | traerrpr leaf |
| 10;60;5; | ERREURS le LONG d'une DROITE de 2 POINTS | ERROR ALONG a LINE DEFINED by 2 POINTS | low | no | no | — | — | traerrpr leaf |
| 10;60;6; | ERREURS PRESSION sur des SURFACES de l'OBJET | ERROR on ONE SURFACE or more of the OBJECT | low | no | no | — | — | traerrpr leaf |
| 10;60;7; | ERREURS PRESSION sur les SURFACES FRONTALIERES | ERROR on the 3D BOUNDARY SURFACES | low | no | no | — | — | traerrpr leaf |
| 10;60;9; | En 2D SURFACE X Y ErreurPression X Y | In 2D SURFACE X Y PressureError X Y | low | no | no | — | — | traerrpr leaf; 2D-only |
| 10;65;1; | ISO-ERREURS du module VITESSE | ISO-VELOCITY ERROR SURFACES or CURVES | low | no | no | — | — | traerrvi leaf — synthetic prefix 10;65; for velocity-error trace; reached from flui/tfluide.f:3042 |
| 10;65;2; | ZONES de COULEURS de l'ERREUR VITESSE | COLOR ZONES of ISO-VELOCITY ERRORS on BOUNDARY | low | no | no | — | — | traerrvi leaf |
| 10;65;3; | ZONES ERREUR VITESSE par SECTIONS X ou Y ou Z=CTE | VELOCITY ERROR CROSS SECTIONS | low | no | no | — | — | traerrvi leaf |
| 10;65;4; | PROFILS ERREUR VITESSE dans PLANS X ou Y ou Z=CTE | VELOCITY ERROR PROFILES on CROSS SECTION | low | no | no | — | — | traerrvi leaf |
| 10;65;5; | ERREURS VITESSE le LONG d'une DROITE de 2 POINTS | VELOCITY ERROR ALONG A LINE DEFINED BY 2 POINTS | low | no | no | — | — | traerrvi leaf |
| 10;65;6; | ERREURS VITESSE sur UNE SURFACE ou PLUS de l'OBJET | VELOCITY ERROR on ONE SURFACE or more of the OBJECT | low | no | no | — | — | traerrvi leaf |
| 10;65;9; | En 2D SURFACE X Y Erreur VitesseExact moins Calc X Y | In 2D SURFACE X Y Error VelocityExact-Comp X Y | low | no | no | — | — | traerrvi leaf; 2D-only |
| 10;70;1; | ISOVALEURS de Rot VITESSE | Rot VELOCITY ISO-VALUES | med | yes | no | — | — | trarot2d leaf — synthetic prefix 10;70; for 2D rotational trace; reached from flui/tfluide.f:2787 |
| 10;70;2; | ISOVALEURS Rot VITESSE en ZONES de COULEURS | Rot VELOCITY ISO-COLOR ZONES | med | yes | no | — | — | trarot2d leaf; canonical 2D vorticity color zones |
| 10;70;9; | TRACE la SURFACE X Y Rot VITESSE X Y | DRAWING of SURFACE X Y Rot VELOCITY X Y | low | no | no | — | — | trarot2d leaf; 2D 3D-perspective lift |
| 10;75;0; | No de la COMPOSANTE du ROTATIONNEL X Y Z ou MODULE | COMPONENT NUMBER X Y Z or MAGNITUDE of Rotational | med | yes | no | — | — | trarot3d leaf — synthetic prefix 10;75; for 3D rotational trace; reached from flui/tfluide.f:2789; canonical entry numeric prompt |
| 10;75;1; | ISO-MODULE de Rot VITESSE | ISO-Rot VELOCITY MAGNITUDES | med | yes | no | — | — | trarot3d leaf |
| 10;75;2; | ISO-MODULE Rot VITESSE en ZONES de COULEURS | Rot VELOCITY COLOR ZONES | med | yes | no | — | — | trarot3d leaf |
| 10;75;3; | ZONES Rot VITESSE par SECTIONS X ou Y ou Z=CTE | Rot VELOCITY on CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | trarot3d leaf |
| 10;75;4; | PROFILS Rot VITESSE dans PLANS X ou Y ou Z=CTE | Rot VELOCITY PROFILE on CROSS SECTION | low | no | no | — | — | trarot3d leaf |
| 10;75;5; | Rot VITESSE sur une DROITE definie par 2 POINTS 3D | Rot VELOCITY ALONG A 3D LINE | low | no | no | — | — | trarot3d leaf |
| 10;75;6; | Rot VITESSE sur UNE SURFACE 3D de l'OBJET 3D | Rot VELOCITY On ONE 3D SURFACE | low | no | no | — | — | trarot3d leaf |
| 10;75;11; | TRACE de la COMPOSANTE X de Rot VITESSE | Drawing of X-COMPONENT of Rot Velocity | low | no | no | — | — | trarot3d leaf |
| 10;75;12; | TRACE de la COMPOSANTE Y de Rot VITESSE | Drawing of Y-COMPONENT of Rot Velocity | low | no | no | — | — | trarot3d leaf |
| 10;75;13; | TRACE de la COMPOSANTE Z de Rot VITESSE | Drawing of Z-COMPONENT of Rot Velocity | low | no | no | — | — | trarot3d leaf |
| 10;75;14; | TRACE de la NORME Rot VITESSE | Drawing of MAGNITUDE of Rot Velocity | low | no | no | — | — | trarot3d leaf |
| 10;80;2; | FLECHE MAX EN CM | CM of MAX VELOCITY ARROW | med | yes | no | — | — | trflevit leaf — synthetic prefix 10;80; for velocity-arrow trace; reached from flui/tfluide.f:1698; testa: cavity2d/nstg3d use 90 trigger after configuring 2-9 |
| 10;80;3; | Un CM VAUT en NORME de la VITESSE | VELOCITY NORM of 1CM of ARROW | med | yes | no | — | — | trflevit leaf; scale prompt |
| 10;80;4; | ECHELLE PRECEDENTE de la VITESSE | PREVIOUS SCALE of VELOCITY | low | no | no | — | — | trflevit leaf |
| 10;80;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | trflevit leaf; prompts couleur0 |
| 10;80;6; | TYPE TRAIT des ARETES MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | trflevit leaf; prompts typtrait |
| 10;80;7; | COULEUR des FLECHES de la VITESSE | VELOCITY ARROW COLOR | low | no | no | — | — | trflevit leaf; prompts couleurs |
| 10;80;8; | EPAISSEUR des FLECHES de la VITESSE | THICKNESS of VELOCITY ARROWS | low | no | no | — | — | trflevit leaf |
| 10;80;9; | PAS DE TRACE des FLECHES-VITESSE | Drawing STEP of VELOCITY ARROWS | low | no | no | — | — | trflevit leaf |
| 10;80;50; | EFFACER le trace actuel | ERASE the WINDOW | low | no | no | — | — | trflevit leaf |
| 10;80;90; | TRACE des FLECHES de la VITESSE | DRAWING of VELOCITY ARROWS | high | yes | no | — | — | trflevit leaf; testa count=high (cavity2d/nstg3d 90; canonical exec trigger); domain HIGH |
| 10;85;0; | CHOIX de la COMPOSANTE X Y Z ou MODULE de la VITESSE | CHOICE of X Y Z COMPONENT or VELOCITY MAGNITUDE | med | yes | no | — | — | tracmovi leaf — synthetic prefix 10;85; for component-velocity trace; reached from flui/tfluide.f:1981; numeric-prompt entry |
| 10;85;1; | ISO-LIGNES ou ISO-SURFACES de la VITESSE dans l'OBJET | ISO-VELOCITY LINES or SURFACES in the OBJECT | med | yes | no | — | — | tracmovi leaf |
| 10;85;2; | ISO-COULEURS ZONES de la VITESSE | ISO-VELOCITY COLORS ZONES | med | yes | no | — | — | tracmovi leaf |
| 10;85;3; | ISO-COULEURS de la VITESSE par SECTIONS X ou Y ou Z=CTE | ISO-VELOCITY COLOR on X or Y or Z=Const CROSS SECTIONS | low | no | no | — | — | tracmovi leaf |
| 10;85;4; | ISO-PROFILS de la VITESSE sur des PLANS X ou Y ou Z=CTE | ISO-VELOCITY PROFILE on X or Y or Z=Const CROSS SECTIONS | low | no | no | — | — | tracmovi leaf |
| 10;85;5; | ISO-VITESSE le LONG d'une DROITE definie par 2 POINTS | ISO-VELOCITY ALONG A LINE DEFINED BY 2 POINTS | low | no | no | — | — | tracmovi leaf |
| 10;85;6; | ISO-VITESSE sur UNE SURFACE ou PLUS de l'OBJET | ISO-VELOCITY on ONE SURFACE or more of the OBJECT | low | no | no | — | — | tracmovi leaf |
| 10;85;9; | En 2D SURFACE X Y Vitesse X Y | In 2D SURFACE X Y Velocity X Y | low | no | no | — | — | tracmovi leaf; 2D-only |
| 10;87;1; | ISOVALEUR-COURBES de la FONCTION COURANT | ISOVALUE CURVES of STREAM FUNCTION | med | yes | no | — | — | tracoura leaf — synthetic prefix 10;87; for stream-function trace; reached from flui/tfluide.f:2647 (2D-only) |
| 10;87;2; | ZONES de COULEURS ISO-VALEURS | COLOR ZONES of ISO-STREAM FUNCTION | med | yes | no | — | — | tracoura leaf |
| 10;87;9; | SURFACE 3d X Y FonctionCourant X Y | 3d-SURFACE X Y StreamFunction X Y | low | no | no | — | — | tracoura leaf; 3D-perspective lift |
| 10;88;1; | Trace des VITESSES des EF SECTIONNES par un PLAN | Drawing of VELOCITY ARROWS in FE CUT by a PLANE | low | no | no | — | — | vitetype leaf — synthetic prefix 10;88; for arrow-type filter; reached from flui/tfluide.f:1828 |
| 10;88;2; | Trace de toutes les FLECHES des VITESSES | Drawing of all VELOCITY ARROWS | low | no | no | — | — | vitetype leaf |
| 10;88;3; | Intervalle de la NORME d'une VITESSE a tracer | VELOCITY MAGNITUDE Min-Max INTERVAL to draw | low | no | no | — | — | vitetype leaf |
| 10;88;4; | Trace des VITESSES de NORME dans l'INTERVALLE | Drawing of VELOCITY ARROWS of MAGNITUDE in the INTERVAL | low | no | no | — | — | vitetype leaf |
| 10;89;1; | CHOIX de la COMPOSANTE X de la Vitesse | Choice of Velocity X-Component | low | no | no | — | — | noxyzmod leaf — synthetic prefix 10;89; for component selector; reached from flui/tfluide.f:1988 (and :2803 under trarot3d) |
| 10;89;2; | CHOIX de la COMPOSANTE Y de la Vitesse | Choice of Velocity Y-Component | low | no | no | — | — | noxyzmod leaf |
| 10;89;3; | CHOIX de la COMPOSANTE Z de la Vitesse | Choice of Velocity Z-Component | low | no | no | — | — | noxyzmod leaf |
| 10;89;4; | CHOIX du MODULE de la Vitesse | Choice of Velocity Magnitude | low | no | no | — | — | noxyzmod leaf |
| 10;95;1; | Donnees Longitudes Latitudes Rayon Temps0 des PARTICULES | Give LONGITUDES LATITUDES RADIUS TIME0 of PARTICLES to RUN | low | no | no | — | — | particle leaf — synthetic prefix 10;95; for particle-tracker; reached from flui/tfluide.f:2315 via PARTICULE subroutine |
| 10;95;2; | Donnees XYZ VitesseXYZ Rayon Temps0 des PARTICULES | Give XYZ VelocityXYZ RADIUS TIME0 of PARTICLES to RUN | low | no | no | — | — | particle leaf |
| 10;95;3; | CALCUL du PARCOURS des PARTICULES donnees | Compute the RUN of GIVEN PARTICLES | low | no | no | — | — | particle leaf |
| 10;95;4; | CALCUL du PARCOURS des PARTICULES-COVID dans l'AIR | Compute the RUN of GIVEN PARTICLES-DROPLETS COVID in AIR | low | no | no | — | — | particle leaf; COVID droplet variant |
| 10;95;5; | TRACE du PARCOURS deja CALCULE de CHAQUE PARTICULE | Draw the already COMPUTED RUN of GIVEN PARTICLES | low | no | no | — | — | particle leaf |
| 99;81;1; | DEPART avec le VECTEUR VITESSE+PRESSION deja calcule au TEMPS INITIAL | START from the V+P VECTOR COMPUTED at INITIAL TIME | low | no | no | — | — | vp0flui leaf — synthetic prefix 99;81; (vp0flui is DEAD CODE: flui/fluidens.f:1460 commented `ccc`); listed for completeness; Plan 02 may exclude vp0flui from QAction set |
| 99;81;2; | DEPART avec la DONNEE sur l'OBJET du VECTEUR VITESSE PRESSION a T0 | START from OBJECT INITIAL DATA to CONSTRUCT INITIAL V+P VECTOR | low | no | no | — | — | vp0flui leaf — DEAD CODE; see flui/fluidens.f:1460 |
| 99;81;3; | REPRISE avec VECTEUR VITESSE PRESSION de TEMPS le PLUS PROCHE de T0 | RESTART from V+P VECTOR of TIME NEAREST to INITIAL TIME | low | no | no | — | — | vp0flui leaf — DEAD CODE; see flui/fluidens.f:1460 |
| 10;80;5;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleur0 sub-menu — shared utility, 16 colour names; reached from flui/tfluide.f:1759 (under trflevit 5;); see 6.1 LEXICON-AUDIT-mail.md rows for full enumeration |
| 10;80;7;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleurs sub-menu — shared utility, 16 colour names; reached from flui/tfluide.f:1782 (under trflevit 7;); see 6.1 LEXICON-AUDIT-mail.md rows for full enumeration |
| 10;80;6;1; | TRACE des TRAITS (compressed — 3 leaves) | Types of Line (compressed — 3 leaves) | low | no | no | — | — | typtrait sub-menu — shared utility, 3 leaves; reached from flui/tfluide.f:1775 (under trflevit 6;); see 6.1 LEXICON-AUDIT-mail.md row for full enumeration |
| 70;1; | Mode d'entree DONNEES (compressed — 3 leaves) | Data Input Interactivity (compressed — 3 leaves) | low | no | no | — | — | mode_es sub-menu — shared utility, 3 leaves; reached from util/suives.f via debuflui 70;; see 6.1/6.2 audits for full enumeration |
| 70;2; | UNITE d'ECRITURE (compressed — 2 leaves) | OUTPUT UNIT (compressed — 2 leaves) | low | no | no | — | — | affiche sub-menu — shared utility, 2 leaves (SCREEN/FILE); reached from util/suives.f via debuflui 70;; see 6.1/6.2 audits |
| 70;3; | UNITE de LECTURE (compressed — 3 leaves) | INPUT UNIT (compressed — 3 leaves) | low | no | no | — | — | lecteur sub-menu — shared utility, 3 leaves; reached from util/suives.f via debuflui 70;; see 6.1/6.2 audits |
| 70;4; | Gestion des UNITES I/O (compressed — 6 leaves) | Input Output units management (compressed — 6 leaves) | low | no | no | — | — | entrsort sub-menu — shared utility, 6 leaves; reached from util/suives.f via debuflui 70;; see 6.1/6.2 audits |
| 70;5; | descriptif FICHIER MS (compressed — 3 leaves) | Characteristics MS file (compressed — 3 leaves) | low | no | no | — | — | fichier sub-menu — shared utility, 3 leaves; reached from util/ajfich.f via debuflui 70;; see 6.1/6.2 audits |
| 60;1; | LARGEUR HAUTEUR PIXELS (compressed — 2 leaves) | WINDOW WIDTH & HEIGHT (compressed — 2 leaves) | low | no | no | — | — | managfen sub-menu — shared utility, 2 leaves; reached from util/managfen.f via debuflui 60;; see 6.1/6.2 audits |
| 70;6; | TUER TMS PLSVO (compressed — 4 leaves) | DELETE TMS PLSVO (compressed — 4 leaves) | low | no | no | — | — | tuer sub-menu — shared utility, 4 leaves; reached from util/tuer via debuflui 70;; see 6.1/6.2 audits |
| 70;7; | SUIVI des TMS (compressed — 5 leaves) | TOOLS on TMS (compressed — 5 leaves) | low | no | no | — | — | suivitms sub-menu — shared utility, 5 leaves; reached from util/suitms via debuflui 70;; see 6.1/6.2 audits |
| 70;8; | SUIVI des FICHIERS MS (compressed — 1 leaf) | MS File MANAGEMENT (compressed — 1 leaf) | low | no | no | — | — | suivi_ms sub-menu — shared utility, 1 leaf; reached from util/suifms via debuflui 70;; see 6.1/6.2 audits |
| 70;9; | GESTION de la MEMOIRE (compressed — 4 leaves) | MS MANAGEMENT (compressed — 4 leaves) | low | no | no | — | — | managmef sub-menu — shared utility, 4 leaves; reached from util/managmef.f via debuflui 70; dispatch; see 6.1/6.2 audits |
| 20;1; | PRECISION pour IDENTIFIER POINTS (compressed — 3 leaves) | PRECISION to IDENTIFY POINTS (compressed — 3 leaves) | low | no | no | — | — | zeros sub-menu — shared utility, 3 leaves; reached from util/zeros.f via debuflui 20;; see 6.1/6.2 audits |
| 19;1; | TYPES d'objets (compressed — 5 leaves) | Types of objects (compressed — 5 leaves) | low | no | no | — | — | typ_objt sub-menu — shared utility, 5 leaves (POINT/LIGNE/SURFACE/VOLUME/OBJET); reached via TRMAIL path under debuflui 19;; see 6.1 LEXICON-AUDIT-mail.md |
| 19;2; | VUES selon des PLANS (compressed — 6 leaves) | VIEWS from SPECIAL PLANS (compressed — 6 leaves) | low | no | no | — | — | vuesplan sub-menu — shared utility, 6 leaves; reached via TRMAIL path under debuflui 19;; see 6.1 audit |
| 19;3; | TRACE des LIGNES options (compressed — 33 leaves) | Line drawing types (compressed — 33 leaves) | low | no | no | — | — | opt_lign sub-menu — shared utility, 33 leaves; reached via TRMAIL path under debuflui 19;; see 6.1 audit |
| 19;4; | TRACE des SURFACES options (compressed — 40 leaves) | Surfaces drawing types (compressed — 40 leaves) | low | no | no | — | — | opt_surf sub-menu — shared utility, 40 leaves; reached via TRMAIL path under debuflui 19;; see 6.1 audit |
| 19;5; | TRACE des AXES (compressed — 4 leaves) | Drawing of AXES (compressed — 4 leaves) | low | no | no | — | — | tracaxes sub-menu — shared utility, 4 leaves; reached via TRMAIL path under debuflui 19;; see 6.1 audit |
| 19;6; | PLAN SECTION VOLUMES (compressed — 4 leaves) | Section Planes for VOLUMES (compressed — 4 leaves) | low | no | no | — | — | sectvopl sub-menu — shared utility, 4 leaves; reached via TRMAIL or TFLUIDE 3D paths; see 6.1 audit |
| 6;1; | SOUS ESPACES (compressed — 2 leaves) | SUB SPACES (compressed — 2 leaves) | low | no | no | — | — | methvvpr sub-menu — shared utility, 2 leaves (eigen-method); reached from reso/calvvp.f:229 — flui never reaches eigenmodes solver but methvvpr is in the cross-module LIMTCL grep set; see 6.2 LEXICON-AUDIT-elas.md rows 6;1; / 6;2; for full enumeration |

<!-- End of audit table — validator will count rows above this line. -->

## Summary Statistics

- **Total rows:** 159
- **By frequency:**
  - `high`: 6 — `1;`, `2;`, `5;`, `7;`, `10;`, `99;` (top-level domain-promoted Navier-Stokes workflow) plus `2;2;` (cl_flui velocity fixation), `10;20;2;` (vitepr2d velocity arrows), `10;30;2;` (vitepr3d velocity arrows), `10;80;90;` (trflevit exec trigger) — total **10** counting indented HIGH leaves
  - `med`: 32 — debuflui codes 3;, 6;, 8;, 9;, 19;, 20; + sub-menu workflow leaves: cl_flui 2;1;, 2;3;; methrecr 5;81;1;, 5;81;2;; nbvitpre 10;1;1;, 10;1;2;; vitepr2d 10;20;1;, 10;20;3;, 10;20;5;, 10;20;9;, 10;20;10;, 10;20;90;; vitepr3d 10;30;1;, 10;30;3;, 10;30;5;, 10;30;10;, 10;30;90;; tracpr2d 10;40;1;, 10;40;2;; tracpr3d 10;50;1;, 10;50;2;; trarot2d 10;70;1;, 10;70;2;; trarot3d 10;75;0;, 10;75;1;, 10;75;2;; trflevit 10;80;2;, 10;80;3;; tracmovi 10;85;0;, 10;85;1;, 10;85;2;; tracoura 10;87;1;, 10;87;2;
  - `low`: 117 — remainder (long-tail visualization leaves + 21 shared util compressed rows + hidden codes 35;/36; + dead-code vp0flui rows)
- **By qaction:** yes count == high+med == **42** (5 toolbar=yes rows are a subset of qaction=yes; QAction set for Plan 02 = 42 rows — larger than 6.2 elas's 20-row set because the visualization tree under code 10; has more workflow-relevant leaves)
- **By toolbar:** exactly 5 `yes` — `2;`, `5;`, `7;`, `10;`, `99;`

## Top-5 Toolbar (Draft — to be ratified at Task 2)

1. `2;` — GIVE PHYSICAL DATA / ENTRER les DONNEES PHYSIQUES du FLUIDE — icon `flui-input-physics.svg`
2. `5;` — INCOMPRESSIBLE STEADY STOKES SOLVER / RESOUDRE le PROBLEME de STOKES STATIONNAIRE — icon `solve-stokes.svg`
3. `7;` — IMPLICIT UNSTEADY NAVIER-STOKES SOLVER / RESOUDRE NAVIER-STOKES INSTATIONNAIRE IMPLICITE — icon `solve-navier-stokes.svg`
4. `10;` — DRAWING of VELOCITIES PRESSURES TEMPERATURES / DESSINER des VITESSES PRESSIONS TEMPERATURES — icon `draw-velocity.svg`
5. `99;` — SAVE DATA and QUIT / SAUVEGARDER les DONNEES et FIN du TRAITEMENT — icon `SP_DialogCloseButton`

Rationale: `2;` is the mandatory data-input gate before any solver runs (canonical first action after `1;` Object); `5;` is the canonical first solver run (steady Stokes baseline); `7;` is the production unsteady Navier-Stokes solver (testa marche/nstg3d evidence); `10;` is the universal result-visualization entry (every result-producing fixture draws results); `99;` is the shared 6.1/6.2 close-session convention. The user may at Task 2 checkpoint:
- Swap `5;` for `1;` if "Object first" workflow is more common (testa shows OBJECT NAME always FIRST but `1;` is implicit-mandatory rather than user-toggleable)
- Swap `7;` for `8;` if FRACTIONAL Navier-Stokes is the production solver (testa marche/nstg3d use BOTH 7; and 8; equally)
- Swap `10;` for `19;` if mesh inspection dominates (testa evidence: result-drawing dominates mesh-inspection in flui workflows)

Any change must keep the count at exactly 5.

## Flui-unique SVG icon set (Plan 02 ships)

Custom .svg filenames introduced by this audit that Plan 02 must ship under
`xvue/qt/resources/icons/flui/`:

- `flui-object.svg` — OBJECT NAME (code 1;)
- `flui-input-physics.svg` — PHYSICAL DATA (code 2;)
- `flui-input-heat.svg` — HEAT DATA (code 3;)
- `solve-stokes.svg` — STEADY/UNSTEADY STOKES SOLVER (codes 5;, 6;)
- `solve-navier-stokes.svg` — NAVIER-STOKES SOLVER family (codes 7;, 8;, 9;)
- `draw-velocity.svg` — DRAWING of VELOCITIES PRESSURES TEMPERATURES (code 10;)

Total: **6 new flui-specific SVGs** (one fewer than 6.2 elas's 7 because the
Stokes/Navier-Stokes solver families share icons across their unsteady
variants). User may consolidate further at Task 2 — e.g., a single
`solve-flow.svg` covering both Stokes and Navier-Stokes drops the count
to 5 new SVGs.

## Shared SVG reuses (from 6.1 mail/)

- `mesh-draw.svg` — reused from `xvue/qt/resources/icons/mail/` for flui
  code `19;` (DRAWING of PLSVO MESHES); the qrc prefix `/xvue/qt/icons`
  resolves the mail path regardless of module context. No file copy
  needed. Same reuse pattern as 6.2 elas code `10;`.

## Cross-References

- 06.1 `LEXICON-AUDIT-mail.md` — 9-column schema template (this audit
  mirrors the structure verbatim) + leaf-level enumeration of shared
  util sub-menus (compressed here)
- 06.2 `LEXICON-AUDIT-elas.md` — same template, second-iteration
  example; demonstrated the synthetic high-number-offset prefix
  pattern (3;81;, 3;91;, 4;50;) replicated here as 5;81;, 10;20;,
  10;30;, etc.
- 06.1 `06.1-01-SUMMARY.md` — de-duplication rule (module-specific
  fully expanded; shared util compressed to 1 row each)
- 06.2 `06.2-01-SUMMARY.md` — domain-review promotion pattern for
  sparse testa corpora (replicated here for the small flui fixture
  set)
- 06.0-UI-SPEC §Per-Module Conformance Contract — flui module-title
  Fluid/Fluide
- ROADMAP Phase 6.3 line 201 — `{File, Fluid, View, Help}` taxonomy
- `tools/validate_audit_md.py` — 9-rule schema validator; module-aware
  ICONS resolver since 6.2 Plan 02 (no validator changes needed in 6.3)
- Plan 02 consumes: rows with `qaction=yes` as the
  `registerFluiActions_stub_` QAction set; rows with `toolbar=yes` as
  the toolbar append set; `icon_source` ending `.svg` (filename not
  already in `icons/mail/` or `icons/elas/`) as the `xvue_icons.qrc`
  append set for `icons/flui/`
