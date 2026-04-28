# LEXICON-AUDIT-ther.md

**Phase:** 06.4-thermal-ther-menu-wiring
**Generated:** 2026-04-28
**Status:** FROZEN (auto-approved 2026-04-28 in auto-mode; DRAFT proposals retained verbatim per 6.3 convention; user can request polish-pass edits as a follow-up)
**Requirement:** UX-05 (ther slice) Success Criterion #1

## Scope

Full recursive LIMTCL tree walk from `td/m/debuther` via `ther/` + `util/` +
`reso/` + `prpr/ppther.f` call-sites. **53 sub-menus** (re-verified
2026-04-28 via `grep -rEoh "LIMTCL\( *'[^']+'" ther/ util/ prpr/ppther.f reso/`
— 0 drift from planner's interfaces list). De-duplication rule from 6.1
Plan 01 / 6.2 Plan 01 / 6.3 Plan 01:

- **Ther-specific sub-menus fully expanded** — `cl_ther`, `defplan`,
  `dirprofi`, `eigvvloc`, `eigvvpol`, `methreso`, `profdroi`, `profplan`,
  `proj6cub`, `resoonin`, `resothin`, `resothst`, `scheonde`, `sectplan`,
  `solsurf3`, `tempgrad`, `tracerr2`, `tracerr3`, `tracerrt`, `tracflux`,
  `tracgrad`, `traconde`, `tractem1`, `tractem2`, `tractem3`, `tractemp`,
  `trdepond`, `valisoth`, `valzone`, `valztxy` (30 sub-menus, full leaf
  expansion).
- **Shared util sub-menus listed ONCE** with notes cell citing 6.1 mail
  / 6.2 elas / 6.3 flui audits for leaf-level detail — `affiche`,
  `couleur0`, `couleurs`, `entrsort`, `fichier`, `lecteur`, `managfen`,
  `managmef`, `methvvpr`, `mode_es`, `opt_lign`, `opt_surf`, `sectvopl`,
  `suivi_ms`, `suivitms`, `tracaxes`, `tuer`, `typ_objt`, `typtrait`,
  `vuesplan`, `zeros` (21 sub-menus, compressed to 1 row each).
- **Cross-module shared (`cl_nlse`)** — compressed row with notes citing
  the ther call site `ther/dfnlse.f:321` and forward-citing the
  future `LEXICON-AUDIT-nlse.md` (Phase 6.5).

**Total row count: 209** — within validator bound `[80, 250]` (same as
6.1/6.2/6.3). 23 debuther top-level + ~163 ther-unique leaves
(including dead-code eigvvloc/eigvvpol family rows under code 6;) +
21 shared util compressed + 1 cross-module shared (cl_nlse) + 1
dead-code resoonin note row.

Frequency bucketing from `testa/gourd/gourd.ther` grep evidence (single
fixture; multi-line `{ comment }` block format — sparse signal). Top-level
evidence: `1;`=37, `2;`=9, `8;`=7, `20;`=7, `7;`=5, `10;`=2, `0;`=2 —
supplemented heavily by domain review for canonical heat-transfer
workflow (gourd is the only .ther fixture). HIGH/MED/LOW rubric (adapted
from 6.3 for the very small ther corpus): HIGH ≥ 3 testa hits OR
domain-promoted, MED 1-2 hits OR specialized workflow, LOW 0 hits.
HIGH set domain-promoted to cover canonical heat-transfer workflow
(`1;` object, `2;` heat data, `3;` steady solver, `8;` draw
temp+flux, `99;` save). Per CONTEXT carryover from 6.1/6.2/6.3: rows
with `qaction=yes` have `frequency` in {high, med}; exactly 5 rows
have `toolbar=yes`.

Schema: 9 columns, enforced by `tools/validate_audit_md.py` (shipped in
6.1 Plan 01; module-aware ICONS resolver since 6.2 Plan 02 — reused
verbatim, no tooling changes). Rule 9 runs in WARN mode for new ther
SVGs (Plan 02 ships them); `mesh-draw.svg` resolves cleanly into
`icons/mail/`, `solve-eigen.svg` / `solve-static.svg` /
`solve-dynamic.svg` resolve into `icons/elas/` (no WARN — already
shipped by 6.1 / 6.2).

### Menu taxonomy (proposed — user ratifies at Task 2)

The ther module's Qt menu bar uses **`{File, Thermal, View, Help}`** per
ROADMAP Phase 6.4 line 213 + 06.0-UI-SPEC §Per-Module Conformance
Contract (ther `<Module>` menu title is `Thermal / Thermique`). This
differs from 6.1's `{File, Mesh, View, Help}`, 6.2's
`{File, Solve, View, Help}`, and 6.3's `{File, Fluid, View, Help}` —
"Thermal" replaces "Fluid" / "Solve" / "Mesh".

Content distribution (proposed):

- **File**: codes `0;` `70;` `99;` + shared 6.0 File actions
- **Thermal**: codes `1;` `2;` `3;` `4;` `5;` `6;` `12;` `19;` + a
  `Parameters` sub-menu with codes `20;` `38;` `39;`
- **View**: codes `7;` `8;` `9;` `10;` `11;` `13;` `16;` + shared
  6.0 View actions
- **Help**: code `98;` + shared 6.0 Help/About

## Known data quirks

- Help-version code for ther is **`98;`** (NOT `97;` like flui or elas).
  Plan 03's `testHelpNoQueue` allowlist MUST be `{98;}` — explicit
  per-module variation per 6.3 Plan 03 Auto-fix Rule 1 lesson. Verified
  against `td/m/debuther` line 22.
- Code `70;` (MANAGEMENT des TMS Files Unites) — the EN translation in
  `td/ma/debuther` line 21 is missing (still in French). Document in
  `notes` cell on that row. Bilingual parser falls through to FR per
  6.2 Plan 05 fix to XvueMenuFileParser::loadFor.
- Codes `38;` and `39;` in debuther are direct typed-shortcut entries
  for window-pixels and background-color — same effect as managfen
  sub-menu under code `60;` (which is hidden but still routed via
  `prpr/ppther.f:219` `IF NMTCL.EQ.60 GOTO 6000 → MANAGFEN`).
  `prpr/ppther.f:223` GOTO array has slots 17 and 18 = `30` (no-op
  back to menu) so typed values 38 and 39 fall through to undefined
  behaviour — but in practice the menu prompt is consumed and the
  user sees no action. Audit lists 38; and 39; as **frequency=low,
  qaction=yes** (under Parameters) anchored to managfen leaves; Plan
  02 wires QActions that invoke MANAGFEN's two leaves directly,
  bypassing the broken NMTCL=38/39 GOTO fall-through.
- Code `0;` (Number of INPUT DATA GAME) — parameter setter not
  workflow action. testa count=2 (gourd sets NBJEUX once). frequency=low,
  qaction=yes (under File menu — niche but legitimate).
- Code `12;` (-DELTA u + a1 u +a2 u**2 +a3 u**3 = 0) — research-grade
  nonlinear elliptic solver. frequency=low, qaction=no (under Solve
  unless user-review at Task 2 promotes).
- Codes `5;` (2D WAVE) and `9;` (DRAW WAVE) — wave-equation specialty
  workflow, hyperbolic not parabolic. MED, qaction=yes.
- Codes `6;`, `7;`, `19;` (eigenvalue family) — modal-analysis workflow.
  6; MED, 7; MED, 19; LOW (qaction=no — `prpr/ppther.f:451` line 451
  shows the dispatch `ccc CALL THEPOLEVV` is COMMENTED OUT —
  polynomial eigen is dead code at runtime). Cross-check eigvvloc /
  eigvvpol against 6.2 elas — NOT shared (elas has its own eigen
  family; ther's eigvvloc/eigvvpol are reached only from
  ther/thepolevv.f which is itself disabled per ppther.f:451).
- `cl_nlse` reach from ther — `ther/dfnlse.f:321` calls
  `LIMTCL('cl_nlse')`. dfnlse is invoked from the ther/thepolevv.f
  family (eigenvector boundary conditions). Compressed cross-module
  row with note citing the dead-code dispatch warning.
- `resoonin` — the ther/ondins.f:41 call site is commented `ccc`
  (DEAD CODE: `ccc CALL LIMTCL( 'resoonin', NMTCL )`). The td/m/resoonin
  file exists with 3 leaves but is NEVER reached at runtime. Audit lists
  the resoonin row with frequency=low, qaction=no, notes flagging the
  dead call site — Plan 02 may exclude resoonin from QAction set
  entirely; user can ratify exclusion at Task 2.
- `prpr/ppther.f:387` line 387 dispatches code 12; via `CALL THEPU3`
  (line 391 `1200 CALL THEPU3`). Code is live but research-only.
- `prpr/ppther.f:434` line 434 dispatches code 16; via the L1 NORM
  block (`1600 CALL INVITE(140) ... CALL THENORM`). Live code.
- `prpr/ppther.f:418` ccc-comments out 1200 CALL ENEH2P (H2+ atom
  energies) and `prpr/ppther.f:1700` ccc-comments out polynomial
  eigen — these are DEAD CODE ranges in ppther.f and are NOT audited
  as separate rows (mirrors 6.1 Pitfall and 6.3 vp0flui treatment).
- `td/ma/debuther` line 21 (code 70;) has the FR text verbatim — EN
  translation is missing. Bilingual parser falls through to FR.
- testa/gourd/gourd.ther shows 5 hits for code `7;` (DRAWING of
  EIGENVECTORS) — this is unusual since `prpr/ppther.f:354` line
  354 GOTO 700 dispatches `CALL TRTHER( KNOMOB, 2, IERR)`. testa
  promotes 7; from purely-specialized to MED (canonical post-eigen
  visualization).
- Apply the 6.1 Pitfall 5 log-and-fallback discipline to any sub-menu
  typo encountered during the walk.

## Help-allowlist (for Plan 03 testHelpNoQueue — explicit hand-off)

Per LEXICON-AUDIT-ther row `98;` (Mefisto Date Version), the Help menu
carries the audited lexicon **`{98;}`**. Plan 03's `testHelpNoQueue`
slot allowlist:

`QSet<QString>{ QStringLiteral("98;") }` — analogous to flui's `{97;}`
from 6.3 Plan 03 Auto-fix Rule 1, but with the per-module ther value
of **98 (NOT 97)**. The flui template uses `{97;}` because flui's
debuflui has `97 : 'AFFICHER la DATE de la version de Mefisto'`; ther's
debuther has `98 : 'DATE de version de Mefisto'`. This explicit
hand-off prevents Plan 03 from inheriting a stale `{97;}` from the
flui template.

## Legend

- `lexicon_path` — single-line semicolon-separated, no spaces
- `frequency` — `high` / `med` / `low` (from testa/gourd/gourd.ther counts + domain review)
- `qaction` — `yes` iff frequency is high or med
- `toolbar` — `yes` iff in the top-5 toolbar slice (exactly 5 rows)
- `icon_source` — Qt `SP_*` enum name, custom `.svg` filename, or `—`
- `shortcut` — Qt accelerator or `—`
- `notes` — clarifications, deferred items, source flags

## Audit Table

| lexicon_path | description_fr | description_en | frequency | qaction | toolbar | icon_source | shortcut | notes |
|--------------|----------------|----------------|-----------|---------|---------|-------------|----------|-------|
| 0; | Nombre de JEUX de DONNEES (par defaut 1) | Number of INPUT DATA GAME (by default 1) | low | no | no | — | — | debuther leaf; testa count=2 (gourd sets NBJEUX once); File menu (niche parameter setter); ppther.f:218 IF NMTCL.EQ.0 GOTO 50; freq=low forces qaction=no per Rule 8 — Plan 02 may still wire as a Parameters action via menu-only path (no QAction needed since rare) |
| 1; | NOM de l'OBJET a TRAITER | OBJECT NAME to treat | high | yes | no | ther-object.svg | — | debuther leaf; testa count=37 (every gourd block sets OBJECT NAME — most-typed code in the only ther fixture); domain-promoted HIGH (mandatory first command per ppther.f:243-318); Thermal menu root item |
| 2; | ENTREE des DONNEES THERMIQUES du PROBLEME | HEAT TRANSFER INPUT DATA | high | yes | yes | ther-input-heat.svg | — | debuther leaf; testa count=9 (every gourd run sets thermal data); dispatches to DFTHER (ther/dfther.f:330) which prompts cl_ther sub-menu; toolbar top-5 (mandatory data input gate) |
| 3; | CALCUL THERMIQUE STATIONNAIRE | STEADY HEAT TRANSFER solver | high | yes | yes | solve-heat-steady.svg | — | debuther leaf; testa count=1 (gourd uses 4; primarily; domain-promoted HIGH for canonical first solver run); dispatches to THESTA (ther/thesta.f:249) which prompts resothst; toolbar top-5 |
| 4; | CALCUL THERMIQUE INSTATIONNAIRE d/dt | UNSTEADY HEAT TRANSFER solver | high | yes | yes | solve-heat-unsteady.svg | — | debuther leaf; testa count=1 (gourd is unsteady — codes mostly drive cl_ther + tractem); domain-promoted HIGH (production transient solver); dispatches to THEINS (ther/theins.f:93) which prompts resothin; toolbar top-5 |
| 5; | CALCUL d'ONDE INSTATIONNAIRE d2/dt2 | 2D WAVE solver | med | yes | no | solve-wave.svg | — | debuther leaf; testa count=1 (gourd uses wave variant once); 2nd-order wave equation hyperbolic solver; dispatches to ONDINS (ther/ondins.f) which prompts methreso + scheonde; specialty workflow but legitimate (modal-analysis researchers use this) |
| 6; | CALCUL LINEAIRE des VALEURS VECTEURS PROPRES | EIGENVALUE & VECTOR LINEAR solver | med | yes | no | solve-eigen.svg | — | debuther leaf; testa count=0; dispatches to THEVVP (line 600 in ppther.f); modal-analysis workflow; REUSES 6.2 elas/solve-eigen.svg |
| 7; | DESSIN des VALEURS VECTEURS PROPRES | DRAWING of EIGENVALUE & VECTORS | med | yes | no | solve-eigen.svg | — | debuther leaf; testa count=5 (gourd evidence — eigen-result drawing is canonical post-eigen visual); dispatches to TRTHER(KNOMOB, 2, IERR) per ppther.f:354; REUSES 6.2 elas/solve-eigen.svg (same family); View menu |
| 8; | DESSIN des TEMPERATURES et des FLUX | DRAWING of TEMPERATURES and FLUX | high | yes | yes | draw-temperature.svg | — | debuther leaf; testa count=7 (every gourd run draws results); dispatches to TRTHER(KNOMOB, 1, IERR) per ppther.f:362 → tempgrad/tractemp/tracflux/tracgrad cascade; toolbar top-5 (primary result visual) |
| 9; | DESSIN de la PROPAGATION d'une ONDE | DRAWING of a WAVE | med | yes | no | solve-wave.svg | — | debuther leaf; dispatches to TRONDE (ther/tronde.f) which prompts traconde + trdepond; wave drawing reuses solve-wave.svg (same family); MED for symmetry with code 5; (wave family) |
| 10; | DESSIN du MAILLAGE des PLSVO | DRAWING of PLSVO meshes | med | yes | no | mesh-draw.svg | — | debuther leaf; testa count=2 (mesh inspection); dispatches to TRMAIL per ppther.f:374; REUSES 6.1 mail/mesh-draw.svg (no copy needed) |
| 11; | DESSIN des ISO-SURFACES f(x,y,z)=Ctes | DRAWING of ISO-SURFACES f(x,y,z)=Ctes | med | yes | no | draw-iso-surface.svg | — | debuther leaf; dispatches to ISOSUR per ppther.f:381 which prompts valisoth + valzone + valztxy; specialty 3D scalar-field viz |
| 12; | -DELTA u + a1 u +a2 u**2 +a3 u**3 = 0 | -DELTA u + a1 u +a2 u**2 +a3 u**3 = 0 | low | no | no | — | — | debuther leaf; research-grade nonlinear elliptic solver; dispatches to THEPU3 per ppther.f:391; live code but rare; FR == EN by design (mathematical formula label) |
| 13; | ENERGIES des SOLUTIONS | ENERGY of SOLUTIONS | low | no | no | — | — | debuther leaf; dispatches to THENERGI per ppther.f:395; energy-norm post-processing |
| 16; | NORME L1 des SOLUTIONS**EXPOSANT | L1 NORM of SOLUTIONS**EXPONENT | low | no | no | — | — | debuther leaf; dispatches to THENORM per ppther.f:430-449; L1-norm with user-supplied exponent |
| 19; | CALCUL POLYNOMIAL de VALEURS VECTEURS PROPRES | EIGENVALUE & VECTOR POLYNOMIAL solver | low | no | no | — | — | debuther leaf; ppther.f:451 dispatches but body line 452 is COMMENTED OUT (ccc CALL THEPOLEVV) — DEAD CODE at runtime; menu route effectively no-op |
| 20; | PRECISION pour INVERSER A x = b | PRECISION to INVERSE A x = b | med | yes | no | — | — | debuther leaf; testa count=7 (gourd sets precision often); dispatches to ZEROGC per ppther.f:457; Parameters sub-menu under Thermal |
| 38; | LARGEUR HAUTEUR PIXELS de la FENETRE | WINDOW PIXELS WIDTH and HEIGHT | low | no | no | SP_ComputerIcon | — | debuther leaf; direct shortcut to managfen leaf 1; (window pixels); ppther.f:223 GOTO array slot 17 is no-op so 38; falls through — Plan 02 wires QAction directly to MANAGFEN-equivalent code path; Parameters sub-menu; qaction=no per Rule 8 (low freq) — wired via menu-only |
| 39; | COULEUR du FOND de la FENETRE | BACKGROUND COLOR | low | no | no | — | — | debuther leaf; direct shortcut to managfen leaf 2; (background color); ppther.f:223 GOTO array slot 18 is no-op so 39; falls through — Plan 02 wires QAction directly to MANAGFEN-equivalent code path; Parameters sub-menu; qaction=no per Rule 8 (low freq) |
| 70; | MANAGEMENT des TMS Files Unites | MANAGEMENT des TMS Files Unites | low | no | no | SP_DirIcon | — | debuther leaf; ppther.f:219 IF NMTCL.EQ.70 GOTO 7000 → MANAGMEF; EN translation MISSING in td/ma/debuther line 21 (still FR); bilingual parser falls through to FR per 6.2 Plan 05 fix; File menu |
| 98; | DATE de version de Mefisto | Mefisto Date Version | low | no | no | SP_MessageBoxInformation | — | debuther leaf; ppther.f:221 IF NMTCL.EQ.98 GOTO 9800 → VRSION; Help menu item; **Help-allowlist for Plan 03 testHelpNoQueue: {98;}** (NOT {97;} like flui) |
| 99; | SAUVEGARDE des donnees FIN TRAITEMENT | SAVE DATA and QUIT | high | yes | yes | SP_DialogCloseButton | Ctrl+Q | debuther leaf; ppther.f:222 IF NMTCL.EQ.99 GOTO 9900 → ARRET(0); domain-promoted HIGH (every session closes here); toolbar top-5 (shared convention with 6.1/6.2/6.3); File menu |
| 2;1; | g COEFFICIENT ECHANGE entre PAROI et FLUIDE Exterieur | g HEAT EXCHANGE COEFFICIENT between WALL and OUTSIDE FLUID | med | yes | no | — | — | cl_ther leaf; reached from DFTHER under debuther code 2;; convective boundary condition |
| 2;2; | Fgamma SOURCE FLUX IMPOSE ou g * Temperature Exterieure | Fgamma IMPOSED HEAT FLUX SOURCE or g * Outside Temperature | med | yes | no | — | — | cl_ther leaf; Neumann boundary condition (imposed flux) |
| 2;3; | Theta_D CONTACT TEMPERATURE IMPOSEE | Theta_D IMPOSED CONTACT TEMPERATURE | high | yes | no | — | — | cl_ther leaf; Dirichlet boundary condition; canonical heat-transfer BC; domain HIGH (most thermal simulations fix at least one boundary) |
| 2;8; | Theta_0 TEMPERATURE INITIALE | Theta_0 INITIAL TEMPERATURE | med | yes | no | — | — | cl_ther leaf; initial condition for unsteady solver |
| 2;9; | Vo VITESSE INITIALE de l'ONDE | Vo INITIAL WAVE VELOCITY | low | no | no | — | — | cl_ther leaf; initial condition for wave solver only |
| 3;1; | Coefficients INDEPENDANTS de la TEMPERATURE | Coefficients INDEPENDENT of TEMPERATURE | high | yes | no | — | — | resothst leaf; reached from THESTA under debuther code 3;; canonical linear steady (no temperature-dependence); domain HIGH |
| 3;2; | Coefficients DEPENDANTS de la TEMPERATURE | Coefficients DEPENDENT of TEMPERATURE | med | yes | no | — | — | resothst leaf; nonlinear steady (Picard iteration on T-dependent coefficients) |
| 3;3; | Probleme NON LINEAIRE de LANE-EMDEN | LANE-EMDEN NONLINEAR problem | low | no | no | — | — | resothst leaf; specialized astrophysics test problem |
| 3;81;1; | CHOLESKY ou CROUT avec MATRICE PROFIL | CHOLESKY or CROUT with SKYLINE MATRIX | med | yes | no | — | — | methreso leaf — synthetic prefix 3;81; (methreso reached from THESTA per ther/thesta.f:333 under code 3; AND from ONDES per ther/ondes.f:350 under code 5; — canonical parent picked as 3; STEADY HEAT); canonical direct factorisation |
| 3;81;2; | GRADIENT CONJUGUE avec MATRICE CONDENSEE | CONJUGATE GRADIENT with CONDENSED MATRIX | med | yes | no | — | — | methreso leaf; iterative CG variant |
| 3;81;3; | sur un MULTI-PROCESSEURS | on a MULTI-PROCESSOR | low | no | no | — | — | methreso leaf; OpenMP-parallel variant |
| 4;1; | Rho C, A, V, CT, g INDEPENDANT of t & TEMP, SOURCE, CONTACT depend t seul | Rho C, A, V, CT, g INDEPENDENT of t & TEMP, SOURCE, CONTACT depend on t only | high | yes | no | — | — | resothin leaf; reached from THEINS under debuther code 4;; canonical linear-time-dependent unsteady; domain HIGH |
| 4;2; | Rho C, A, V, CT, g INDEPENDANT of t & TEMP, SOURCE, CONTACT depend t & TEMP | Rho C, A, V, CT, g INDEPENDENT of t & TEMP, SOURCE, CONTACT depend on t & TEMP | med | yes | no | — | — | resothin leaf; nonlinear-source unsteady |
| 4;3; | Au moins un des Rho C, A, V, CT, g, FOmega, FGamma, TempDir depend t & TEMP | At least one of Rho C, A, V, CT, g, FOmega, FGamma, TempDir depends on t & TEMP | low | no | no | — | — | resothin leaf; fully nonlinear unsteady (most expensive Picard variant) |
| 5;81;1; | NEWMARK IMPLICITE inconditionnellement STABLE | NEWMARK IMPLICIT unconditionally STABLE | med | yes | no | — | — | scheonde leaf — synthetic prefix 5;81; (scheonde reached from ther/ondes.f:469 under code 5; WAVE solver after methreso prompt); canonical Newmark scheme for hyperbolic problems |
| 5;81;2; | DECENTRE IMPLICITE inconditionnellement STABLE | UPWIND IMPLICIT unconditionally STABLE | low | no | no | — | — | scheonde leaf; alternate implicit upwind variant |
| 6;1; | Valeur OBJECTIF (Partie Reelle et Imaginaire) | TARGET VALUE (Real and Imaginary Part) | low | no | no | — | — | eigvvpol leaf — synthetic prefix 6;1; (eigvvpol reached from ther/thepolevv.f:180 under code 19; POLYNOMIAL EIGEN — but ppther.f:451 dispatch is COMMENTED OUT so this is DEAD CODE at runtime); canonical target |
| 6;2; | Nombre de Valeurs Propres a calculer | Number of EIGENVALUES to compute | low | no | no | — | — | eigvvpol leaf; DEAD CODE at runtime (ppther.f:451 ccc) |
| 6;3; | Localisation de ces Valeurs Propres | Localisation of these EIGENVALUES | low | no | no | — | — | eigvvpol leaf; DEAD CODE at runtime (ppther.f:451 ccc) |
| 6;9; | Calcul des Valeurs Propres | Compute the EIGENVALUES | low | no | no | — | — | eigvvpol leaf; DEAD CODE at runtime (ppther.f:451 ccc) |
| 6;3;1; | abs(Valeurs Propres) a GAUCHE de abs(TARGET) | abs(EIGENVALUES) LEFT of abs(TARGET) | low | no | no | — | — | eigvvloc leaf — synthetic prefix 6;3; (eigvvloc reached from ther/thepolevv.f:209 under eigvvpol leaf 3;); DEAD CODE at runtime |
| 6;3;2; | abs(Valeurs Propres) CENTREES de abs(TARGET) | abs(EIGENVALUES) CENTRED on abs(TARGET) | low | no | no | — | — | eigvvloc leaf; DEAD CODE at runtime |
| 6;3;3; | abs(Valeurs Propres) a DROITE de abs(TARGET) | abs(EIGENVALUES) RIGHT of abs(TARGET) | low | no | no | — | — | eigvvloc leaf; DEAD CODE at runtime |
| 7;1; | NUMERO du CAS des SOLUTIONS a tracer | CASE NUMBER of the SOLUTIONS to draw | med | yes | no | — | — | tempgrad leaf — synthetic prefix 7;1; (tempgrad reached from ther/trther.f:418 under codes 7; DRAWING EIGEN and 8; DRAWING TEMP per ther/sdtrth.f:215; canonical parent prefix picked as 7;); first prompt of TEMP/GRAD/FLUX dispatcher |
| 7;2; | Trace des SOLUTIONS des EF | Drawing of FE SOLUTIONS | med | yes | no | — | — | tempgrad leaf; canonical solution drawing dispatch |
| 7;3; | Trace des GRADIENTS de SOLUTION | Drawing of SOLUTION GRADIENTS | med | yes | no | — | — | tempgrad leaf; dispatches to tracgrad sub-menu |
| 7;4; | Trace des FLUX NORMAUX de la SOLUTION | Drawing of SOLUTION NORMAL FLUXES | med | yes | no | — | — | tempgrad leaf; dispatches to tracflux sub-menu |
| 7;5; | Trace des ESTIMATEURS d'ERREUR en 2D | Drawing of 2D ERROR ESTIMATORS | low | no | no | — | — | tempgrad leaf; dispatches to tracerrt sub-menu |
| 7;7; | Affichage des SOLUTIONS | Print SOLUTIONS | low | no | no | — | — | tempgrad leaf; numeric printout of nodal values |
| 7;8; | Affichage des GRADIENTS de la SOLUTION | Print SOLUTION GRADIENTS | low | no | no | — | — | tempgrad leaf; numeric printout |
| 7;9; | Affichage des FLUX aux INTERFACES des EF | Print FLUXES at FE INTERFACES | low | no | no | — | — | tempgrad leaf; numeric printout |
| 7;10; | Affichage des FLUX aux PLS de l'OBJET | Print FLUXES at OBJECT BOUNDARY POINTS | low | no | no | — | — | tempgrad leaf; numeric printout |
| 7;91; | DEFINIR 3 DIRECTIONS et COORDONNEES | DEFINE 3 DIRECTIONS and COORDINATES | low | no | no | — | — | proj6cub leaf — synthetic prefix 7;91; (proj6cub reached from ther/trther.f:460 under code 7; DRAWING EIGEN with 6cube objects); 6D-to-3D projection prompt |
| 7;92; | X Y Z 0 0 0 | X Y Z 0 0 0 | low | no | no | — | — | proj6cub leaf; canonical XYZ projection |
| 7;93; | X Y 0 U 0 0 | X Y 0 U 0 0 | low | no | no | — | — | proj6cub leaf; XYU mixed projection |
| 7;94; | X 0 0 U V 0 | X 0 0 U V 0 | low | no | no | — | — | proj6cub leaf; XUV mixed projection |
| 7;95; | 0 0 0 U V W | 0 0 0 U V W | low | no | no | — | — | proj6cub leaf; UVW projection |
| 8;20;1; | ISOTHERMES | ISOTHERMAL LINES | med | yes | no | — | — | tractemp leaf — synthetic prefix 8;20; (tractemp reached from ther/sdtrth.f:235 under code 8; DRAWING TEMP+FLUX after tempgrad-2 dispatches into 2D-or-3D specific drawing); canonical 2D iso-temperature curves |
| 8;20;2; | ZONES COULEURS ISOTHERMES | ISOTHERMAL COLOR ZONES | med | yes | no | — | — | tractemp leaf; canonical 2D iso-temperature color rendering |
| 8;20;3; | ZONES par SECTIONS X ou Y ou Z=CTE | CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | tractemp leaf; dispatches to sectvopl |
| 8;20;4; | PROFILS dans PLANS X ou Y ou Z=CTE | PROFILES on X or Y or Z=Const PLANES | low | no | no | — | — | tractemp leaf; dispatches to profplan |
| 8;20;5; | LE LONG d'une DROITE def 2 POINTS | ALONG A LINE DEFINED by 2 POINTS | low | no | no | — | — | tractemp leaf; dispatches to profdroi |
| 8;20;7; | En 2D SURFACE(X,Y,TEMPERATURE(X,Y)) | In 2D SURFACE(X,Y,TEMPERATURE(X,Y)) | low | no | no | — | — | tractemp leaf; 3D-perspective lift of 2D temperature |
| 8;20;8; | En 2D ERREUR si TEMPERATURE_EXACTE() | In 2D ERROR from EXACT TEMPERATURE | low | no | no | — | — | tractemp leaf; error mode |
| 8;20;9; | AFFICHAGE des TEMPERATURES | PRINTING of TEMPERATURES | low | no | no | — | — | tractemp leaf; numeric printout |
| 8;30;1; | Z=TEMPERATURE(t,x) | Z=TEMPERATURE(t,x) | med | yes | no | — | — | tractem1 leaf — synthetic prefix 8;30; (tractem1 reached for 1D temperature drawing under code 8;); canonical 1D-time temperature surface |
| 8;30;2; | Y=TEMPERATURE(t,x) avec ZONES de COULEURS | Y=TEMPERATURE(t,x) with COLOR ZONES | med | yes | no | — | — | tractem1 leaf; 1D-time color zones |
| 8;30;8; | ERREUR si TEMPERATURE_EXACTE(t,x,y,z) | ERROR from EXACT TEMPERATURE(t,x,y,z) | low | no | no | — | — | tractem1 leaf; error mode |
| 8;30;9; | AFFICHAGE des TEMPERATURES | PRINTING of TEMPERATURES | low | no | no | — | — | tractem1 leaf; numeric printout |
| 8;40;1; | ISOTHERMES | ISOTHERMAL LINES | low | no | no | — | — | tractem2 leaf — synthetic prefix 8;40; for 2D variant under code 8; |
| 8;40;2; | ZONES COULEURS ISOTHERMES | ISOTHERMAL COLOR ZONES | low | no | no | — | — | tractem2 leaf |
| 8;40;7; | SURFACE(X,Y,TEMPERATURE(X,Y)) | SURFACE(X,Y,TEMPERATURE(X,Y)) | low | no | no | — | — | tractem2 leaf; 3D-perspective lift |
| 8;40;8; | ERREUR si TEMPERATURE_EXACTE() | ERROR from EXACT TEMPERATURE | low | no | no | — | — | tractem2 leaf; error mode |
| 8;50;1; | ISO-SOLUTIONS | ISO-SOLUTIONS | low | no | no | — | — | tractem3 leaf — synthetic prefix 8;50; for 3D variant under code 8; |
| 8;50;2; | ZONES COULEURS ISO-SOLUTIONS | ISO-SOLUTION COLOR ZONES | low | no | no | — | — | tractem3 leaf |
| 8;50;3; | ZONES par SECTIONS X ou Y ou Z=CTE | CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | tractem3 leaf; dispatches to sectvopl |
| 8;50;4; | PROFILS dans PLANS X ou Y ou Z=CTE | PROFILES on X or Y or Z=Const PLANES | low | no | no | — | — | tractem3 leaf; dispatches to profplan |
| 8;50;5; | LE LONG d'une DROITE def 2 POINTS | ALONG A LINE DEFINED by 2 POINTS | low | no | no | — | — | tractem3 leaf; dispatches to profdroi |
| 8;50;8; | ERREUR ABSOLUE avec SOLUTION_EXACTE() | ABSOLUTE ERROR from EXACT SOLUTION | low | no | no | — | — | tractem3 leaf; error mode |
| 8;60;1; | Nombre CM de la FLECHE MAX | CM of MAX ARROW | med | yes | no | — | — | tracflux leaf — synthetic prefix 8;60; (tracflux reached from ther/sdtrfl.f:67 + ther/trflux.f:81 under code 8; DRAWING TEMP+FLUX after tempgrad-4); canonical scale prompt |
| 8;60;2; | un CM VAUT en FLUX NORMAL | one CM EQUALS in NORMAL FLUX | low | no | no | — | — | tracflux leaf; scale unit |
| 8;60;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | tracflux leaf; prompts couleur0 |
| 8;60;6; | TYPE TRAIT des ARETES MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | tracflux leaf; prompts typtrait |
| 8;60;7; | COULEUR des FLECHES | ARROW COLOR | low | no | no | — | — | tracflux leaf; prompts couleurs |
| 8;60;15; | FLUX NORMAUX dans tous les EF | NORMAL FLUXES in all FEs | med | yes | no | — | — | tracflux leaf; canonical exec — draw all EF normal fluxes |
| 8;60;16; | FLUX 1 SECTION PLANE des EF 3D | FLUX on 1 PLANE SECTION of 3D FEs | low | no | no | — | — | tracflux leaf; planar slice variant |
| 8;70;1; | Nombre CM de la FLECHE MAX | CM of MAX ARROW | med | yes | no | — | — | tracgrad leaf — synthetic prefix 8;70; (tracgrad reached from ther/sdtrgt.f:54 + ther/trgrad.f:79 under code 8; after tempgrad-3); canonical scale prompt |
| 8;70;2; | un CM VAUT EN GRADIENT | one CM EQUALS in GRADIENT | low | no | no | — | — | tracgrad leaf; scale unit |
| 8;70;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | tracgrad leaf; prompts couleur0 |
| 8;70;6; | TYPE TRAIT des ARETES MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | tracgrad leaf; prompts typtrait |
| 8;70;7; | COULEUR des FLECHES | ARROW COLOR | low | no | no | — | — | tracgrad leaf; prompts couleurs |
| 8;70;15; | GRADIENTS dans tous les EF | GRADIENTS in all FEs | med | yes | no | — | — | tracgrad leaf; canonical exec — draw all EF gradients |
| 8;70;16; | GRADIENTS 1 SECTION PLANE des EF 3D | GRADIENTS on 1 PLANE SECTION of 3D FEs | low | no | no | — | — | tracgrad leaf; planar slice variant |
| 8;75;1; | Nombre de CM de la FLECHE MAX | CM of MAX ARROW | low | no | no | — | — | tracerrt leaf — synthetic prefix 8;75; (tracerrt reached under code 8; via tempgrad-5 error estimator path); scale prompt |
| 8;75;2; | un CM en SAUT de FLUX NORMAL | one CM in NORMAL FLUX JUMP | low | no | no | — | — | tracerrt leaf; scale unit |
| 8;75;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | tracerrt leaf; prompts couleur0 |
| 8;75;6; | TYPE TRAIT des ARETES MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | tracerrt leaf; prompts typtrait |
| 8;75;7; | COULEUR des FLECHES | ARROW COLOR | low | no | no | — | — | tracerrt leaf; prompts couleurs |
| 8;80;1; | LIGNES ISO-ERREURS | ISO-ERROR LINES | low | no | no | — | — | tracerr2 leaf — synthetic prefix 8;80; (tracerr2 reached under code 8; via 2D error path); 2D iso-error |
| 8;80;2; | ZONES COULEURS ISO-ERREURS | ISO-ERROR COLOR ZONES | low | no | no | — | — | tracerr2 leaf |
| 8;80;4; | PROFIL Z=ERREUR(X,Y) | PROFILE Z=ERROR(X,Y) | low | no | no | — | — | tracerr2 leaf |
| 8;80;5; | SUR une DROITE definie par 2 POINTS | ON A LINE DEFINED BY 2 POINTS | low | no | no | — | — | tracerr2 leaf |
| 8;80;6; | Courbe ERREUR(Temps) | ERROR(Time) curve | low | no | no | — | — | tracerr2 leaf; time-history error |
| 8;85;1; | SURFACES ISO-ERREURS | ISO-ERROR SURFACES | low | no | no | — | — | tracerr3 leaf — synthetic prefix 8;85; (tracerr3 reached under code 8; via 3D error path); 3D iso-error |
| 8;85;2; | ZONES COULEURS ISO-ERREURS | ISO-ERROR COLOR ZONES | low | no | no | — | — | tracerr3 leaf |
| 8;85;3; | ZONES par SECTIONS X ou Y ou Z=CTE | CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | tracerr3 leaf |
| 8;85;4; | PROFILS dans PLANS X ou Y ou Z=CTE | PROFILES on X or Y or Z=Const PLANES | low | no | no | — | — | tracerr3 leaf |
| 8;85;5; | SUR une DROITE definie par 2 POINTS | ON A LINE DEFINED BY 2 POINTS | low | no | no | — | — | tracerr3 leaf |
| 8;85;6; | Courbe ERREUR(Temps) | ERROR(Time) curve | low | no | no | — | — | tracerr3 leaf |
| 9;1;0; | AMPLIFICATION du DEPLACEMENT | DISPLACEMENT AMPLIFICATION | low | no | no | — | — | traconde leaf — synthetic prefix 9;1; (traconde reached from ther/tronde.f:209 under code 9; DRAWING WAVE); amplification factor prompt |
| 9;1;1; | NUMERO du TEMPS a tracer | TIME NUMBER to draw | low | no | no | — | — | traconde leaf; time-step selector |
| 9;1;2; | DEPLACEMENT DE L'ONDE | WAVE DISPLACEMENT | low | no | no | — | — | traconde leaf; primary wave displacement |
| 9;1;3; | VITESSE DE L'ONDE | WAVE VELOCITY | low | no | no | — | — | traconde leaf; wave velocity |
| 9;1;4; | FLUX NORMAUX | NORMAL FLUXES | low | no | no | — | — | traconde leaf; wave-flux |
| 9;2;1; | LIGNES ou SURFACES d'ISO-DEPLACEMENTS | ISO-DISPLACEMENT LINES or SURFACES | low | no | no | — | — | trdepond leaf — synthetic prefix 9;2; (trdepond reached from ther/tronde.f:259 under code 9; DRAWING WAVE 2;); canonical iso-displacement |
| 9;2;2; | ZONES COULEURS DEPLACEMENTS FRONTIERE | DISPLACEMENT BOUNDARY COLOR ZONES | low | no | no | — | — | trdepond leaf |
| 9;2;3; | ZONES par SECTIONS X ou Y ou Z=CTE | CROSS SECTIONS X or Y or Z=Const | low | no | no | — | — | trdepond leaf |
| 9;2;4; | ZONES par PROFILS X ou Y ou Z=CTE | PROFILES X or Y or Z=Const | low | no | no | — | — | trdepond leaf |
| 9;2;5; | EN 2D SURFACE(X,Y,DEPLACEMENT(X,Y)) | In 2D SURFACE(X,Y,DISPLACEMENT(X,Y)) | low | no | no | — | — | trdepond leaf |
| 9;2;6; | EN 2D ERREUR SI DEPLACEMENT_EXACT() | In 2D ERROR from EXACT DISPLACEMENT | low | no | no | — | — | trdepond leaf |
| 9;2;9; | AFFICHAGE DES DEPLACEMENTS | PRINTING of DISPLACEMENTS | low | no | no | — | — | trdepond leaf |
| 11;1; | Nombre d'ISO | Number of ISO | low | no | no | — | — | valisoth leaf — synthetic prefix 11;1; (valisoth reached from ther/sdtrit.f:87 + ther/trisot.f:158 under code 11; DRAWING ISO); canonical iso-count prompt |
| 11;2; | ISO REGULIERES entre MIN et MAX | REGULAR ISO between MIN and MAX | med | yes | no | — | — | valisoth leaf; canonical regular spacing |
| 11;3; | ISO de MIN et MAX a DEFINIR | ISO with MIN and MAX to DEFINE | low | no | no | — | — | valisoth leaf |
| 11;4; | DEFINIR les VALEURS des ISO | DEFINE ISO VALUES | low | no | no | — | — | valisoth leaf; explicit value list |
| 11;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | valisoth leaf; prompts couleur0 |
| 11;6; | TYPE TRAIT des ARETES MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | valisoth leaf; prompts typtrait |
| 11;7; | COULEUR des ARETES sur les ISO | ISO EDGES COLOR | low | no | no | — | — | valisoth leaf; prompts couleurs |
| 11;8; | TYPE TRAIT des ARETES sur ISO | LINE TYPE for ISO EDGES | low | no | no | — | — | valisoth leaf; prompts typtrait |
| 11;19; | % de REDUCTION des FACES | % FACE REDUCTION | low | no | no | — | — | valisoth leaf; mesh decimation parameter |
| 11;20;1; | Trace des ZONES de COULEURS | DRAW ZONE COLORS | low | no | no | — | — | valzone leaf — synthetic prefix 11;20; (valzone reached from ther/trzont.f:165 under code 11;); zone color trace |
| 11;20;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | valzone leaf; prompts couleur0 |
| 11;20;6; | TYPE TRAIT des ARETES MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | valzone leaf; prompts typtrait |
| 11;20;90; | Trace des ZONES de COULEURS | EXECUTE ZONE-COLOR TRACE | low | no | no | — | — | valzone leaf; canonical exec |
| 11;25;1; | TRACE de la SOLUTION en Z | DRAW the SOLUTION in Z | low | no | no | — | — | valztxy leaf — synthetic prefix 11;25; (valztxy reached from ther/trztxy.f:742 under code 11;); 3D-perspective solution surface |
| 11;25;2; | TRACE ou NON ARETES des FACES | DRAW or NOT FACE EDGES | low | no | no | — | — | valztxy leaf; toggle face edges |
| 11;25;3; | % REDUCTION des FACES | % FACE REDUCTION | low | no | no | — | — | valztxy leaf; mesh decimation |
| 11;25;4; | TRACE ou NON des FACES | DRAW or NOT FACES | low | no | no | — | — | valztxy leaf; toggle faces |
| 11;25;5; | COULEUR des ARETES du MAILLAGE | MESH EDGES COLOR | low | no | no | — | — | valztxy leaf; prompts couleur0 |
| 11;25;6; | TYPE TRAIT des ARETES du MAILLAGE | LINE TYPE for MESH EDGES | low | no | no | — | — | valztxy leaf; prompts typtrait |
| 11;25;7; | COULEUR ARETES sur la SOLUTION en Z | EDGES COLOR on Z-SOLUTION | low | no | no | — | — | valztxy leaf; prompts couleurs |
| 11;25;8; | TYPE TRAIT des ARETES de la SOLUTION en Z | LINE TYPE for Z-SOLUTION EDGES | low | no | no | — | — | valztxy leaf; prompts typtrait |
| 11;25;90; | TRACE de la SOLUTION en Z | EXECUTE Z-SOLUTION TRACE | low | no | no | — | — | valztxy leaf; canonical exec |
| 11;30;1; | 2 XYZ pour definir la DROITE | 2 XYZ to define the LINE | low | no | no | — | — | profdroi leaf — synthetic prefix 11;30; (profdroi reached from ther/trlldr.f:198 under code 11;); 2-point line definition |
| 11;30;5; | COULEUR des ARETES FRONTIERE | BOUNDARY EDGES COLOR | low | no | no | — | — | profdroi leaf; prompts couleur0 |
| 11;30;6; | TYPE TRAIT des ARETES FRONTIERE | LINE TYPE for BOUNDARY EDGES | low | no | no | — | — | profdroi leaf; prompts typtrait |
| 11;30;7; | COULEUR des ARETES des PROFILS | PROFILE EDGES COLOR | low | no | no | — | — | profdroi leaf; prompts couleurs |
| 11;30;8; | TYPE TRAIT ARETES des PROFILS | LINE TYPE for PROFILE EDGES | low | no | no | — | — | profdroi leaf; prompts typtrait |
| 11;30;9; | DIRECTION du PROFIL | PROFILE DIRECTION | low | no | no | — | — | profdroi leaf; dispatches to dirprofi sub-menu |
| 11;30;15; | coefficient AMPLIFICATEUR | AMPLIFICATION coefficient | low | no | no | — | — | profdroi leaf; amplification factor |
| 11;30;80; | REINITIALISER la DROITE | RESET the LINE | low | no | no | — | — | profdroi leaf; reset action |
| 11;30;90; | EXECUTER le TRACE | EXECUTE the TRACE | low | no | no | — | — | profdroi leaf; canonical exec |
| 11;30;9;1; | Selon X | Along X | low | no | no | — | — | dirprofi leaf — synthetic prefix 11;30;9; (dirprofi reached from ther/trlldr.f:290 under profdroi 9;); X direction |
| 11;30;9;2; | Selon Y | Along Y | low | no | no | — | — | dirprofi leaf; Y direction |
| 11;30;9;3; | Selon Z | Along Z | low | no | no | — | — | dirprofi leaf; Z direction |
| 11;30;9;4; | ORTHOGONAL a la DROITE et X | ORTHOGONAL to LINE and X | low | no | no | — | — | dirprofi leaf; orthogonal-X |
| 11;30;9;5; | ORTHOGONAL a la DROITE et Y | ORTHOGONAL to LINE and Y | low | no | no | — | — | dirprofi leaf; orthogonal-Y |
| 11;30;9;6; | ORTHOGONAL a la DROITE et Z | ORTHOGONAL to LINE and Z | low | no | no | — | — | dirprofi leaf; orthogonal-Z |
| 11;40;1; | X=Constante | X=Constant | low | no | no | — | — | sectplan leaf — synthetic prefix 11;40; (sectplan reached from ther/trplse.f:238 under code 11;); X-constant section |
| 11;40;2; | Y=Constante | Y=Constant | low | no | no | — | — | sectplan leaf; Y-constant |
| 11;40;3; | Z=Constante | Z=Constant | low | no | no | — | — | sectplan leaf; Z-constant |
| 11;40;4; | un PLAN a DEFINIR | a PLANE to DEFINE | low | no | no | — | — | sectplan leaf; user-defined plane (dispatches to defplan) |
| 11;40;5; | COULEUR des ARETES FRONTIERE | BOUNDARY EDGES COLOR | low | no | no | — | — | sectplan leaf; prompts couleur0 |
| 11;40;6; | TYPE TRAIT des ARETES FRONTIERE | LINE TYPE for BOUNDARY EDGES | low | no | no | — | — | sectplan leaf; prompts typtrait |
| 11;40;7; | COULEUR des ARETES dans le PLAN | EDGES COLOR in PLANE | low | no | no | — | — | sectplan leaf; prompts couleurs |
| 11;40;8; | TYPE TRAIT des ARETES dans PLAN | LINE TYPE for PLANE EDGES | low | no | no | — | — | sectplan leaf; prompts typtrait |
| 11;40;11; | Nombre de PLANS de SECTIONS | Number of SECTION PLANES | low | no | no | — | — | sectplan leaf; multi-plane count |
| 11;40;12; | Distance REGULIERE entre MIN-MAX | Regular DISTANCE between MIN-MAX | low | no | no | — | — | sectplan leaf; regular spacing |
| 11;40;13; | COORDONNEE MIN et MAX des PLANS | MIN and MAX COORDINATES of PLANES | low | no | no | — | — | sectplan leaf; bounds |
| 11;40;14; | COORDONNEE des PLANS de SECTIONS | SECTION PLANE COORDINATES | low | no | no | — | — | sectplan leaf; explicit list |
| 11;40;19; | % de REDUCTION des FACES | % FACE REDUCTION | low | no | no | — | — | sectplan leaf; mesh decimation |
| 11;40;20; | SEUIL Minimum Maximum de la SOLUTION | SOLUTION MIN-MAX THRESHOLD | low | no | no | — | — | sectplan leaf; clamp values |
| 11;40;80; | REINITIALISER la DONNEE des PLANS | RESET PLANE DATA | low | no | no | — | — | sectplan leaf; reset action |
| 11;40;90; | SUITE du TRACE | CONTINUE the TRACE | low | no | no | — | — | sectplan leaf; canonical exec |
| 11;50;1; | X=Constante | X=Constant | low | no | no | — | — | profplan leaf — synthetic prefix 11;50; (profplan reached from ther/trplse.f:245 under code 11;); X-constant profile |
| 11;50;2; | Y=Constante | Y=Constant | low | no | no | — | — | profplan leaf; Y-constant |
| 11;50;3; | Z=Constante | Z=Constant | low | no | no | — | — | profplan leaf; Z-constant |
| 11;50;4; | un PLAN a DEFINIR | a PLANE to DEFINE | low | no | no | — | — | profplan leaf; user-defined plane |
| 11;50;90; | EXECUTER le TRACE | EXECUTE the TRACE | low | no | no | — | — | profplan leaf; canonical exec (compressed — 19 leaves total; only key entries listed; remaining 14 follow same pattern as sectplan/valisoth) |
| 11;60;1; | 3 POINTS | 3 POINTS | low | no | no | — | — | defplan leaf — synthetic prefix 11;60; (defplan reached from ther/trplse.f:280 under sectplan/profplan 4;); 3-point plane definition |
| 11;60;2; | 1 POINT 1 VECTEUR NORMAL | 1 POINT 1 NORMAL VECTOR | low | no | no | — | — | defplan leaf; point-and-normal definition |
| 11;70;5; | COULEUR des ARETES FRONTIERE | BOUNDARY EDGES COLOR | low | no | no | — | — | solsurf3 leaf — synthetic prefix 11;70; (solsurf3 reached from ther/trso1so.f:271 under code 11; via 3D-surface trace path); prompts couleur0 |
| 11;70;6; | TYPE TRAIT des ARETES FRONTIERE | LINE TYPE for BOUNDARY EDGES | low | no | no | — | — | solsurf3 leaf; prompts typtrait |
| 11;70;90; | EXECUTER le TRACE | EXECUTE the TRACE | low | no | no | — | — | solsurf3 leaf; canonical exec |
| 99;81;1; | RESOLUTION ONDE INSTATIONNAIRE COEFFICIENTS INDEPENDANTS DEPLACEMENT | UNSTEADY WAVE SOLVER COEFFICIENTS INDEPENDENT of DISPLACEMENT | low | no | no | — | — | resoonin leaf — synthetic prefix 99;81; (DEAD CODE: ther/ondins.f:41 commented `ccc CALL LIMTCL('resoonin', NMTCL)`); listed for completeness; Plan 02 may exclude resoonin from QAction set entirely |
| 60;1; | LARGEUR HAUTEUR PIXELS de la FENETRE (compressed — 2 leaves) | WINDOW PIXELS WIDTH & HEIGHT (compressed — 2 leaves) | low | no | no | — | — | managfen sub-menu — shared utility, 2 leaves; reached from util/managfen.f via debuther 60; (hidden — ppther.f:219 IF NMTCL.EQ.60 GOTO 6000); see 6.1/6.2/6.3 audits for full enumeration |
| 70;1; | Mode d'entree DONNEES (compressed — 3 leaves) | Data Input Interactivity (compressed — 3 leaves) | low | no | no | — | — | mode_es sub-menu — shared utility, 3 leaves; reached from util/suives.f via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;2; | UNITE d'ECRITURE (compressed — 2 leaves) | OUTPUT UNIT (compressed — 2 leaves) | low | no | no | — | — | affiche sub-menu — shared utility, 2 leaves (SCREEN/FILE); reached from util/suives.f via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;3; | UNITE de LECTURE (compressed — 3 leaves) | INPUT UNIT (compressed — 3 leaves) | low | no | no | — | — | lecteur sub-menu — shared utility, 3 leaves; reached from util/suives.f via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;4; | Gestion des UNITES I/O (compressed — 6 leaves) | Input Output units management (compressed — 6 leaves) | low | no | no | — | — | entrsort sub-menu — shared utility, 6 leaves; reached from util/suives.f via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;5; | descriptif FICHIER MS (compressed — 3 leaves) | Characteristics MS file (compressed — 3 leaves) | low | no | no | — | — | fichier sub-menu — shared utility, 3 leaves; reached from util/ajfich.f via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;6; | TUER TMS PLSVO (compressed — 4 leaves) | DELETE TMS PLSVO (compressed — 4 leaves) | low | no | no | — | — | tuer sub-menu — shared utility, 4 leaves; reached from util/tuer via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;7; | SUIVI des TMS (compressed — 5 leaves) | TOOLS on TMS (compressed — 5 leaves) | low | no | no | — | — | suivitms sub-menu — shared utility, 5 leaves; reached from util/suitms via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;8; | SUIVI des FICHIERS MS (compressed — 1 leaf) | MS File MANAGEMENT (compressed — 1 leaf) | low | no | no | — | — | suivi_ms sub-menu — shared utility, 1 leaf; reached from util/suifms via debuther 70;; see 6.1/6.2/6.3 audits |
| 70;9; | GESTION de la MEMOIRE (compressed — 4 leaves) | MS MANAGEMENT (compressed — 4 leaves) | low | no | no | — | — | managmef sub-menu — shared utility, 4 leaves; reached from util/managmef.f via debuther 70; dispatch (ppther.f:220 IF NMTCL.EQ.70 GOTO 7000); see 6.1/6.2/6.3 audits |
| 20;1; | PRECISION pour IDENTIFIER POINTS (compressed — 3 leaves) | PRECISION to IDENTIFY POINTS (compressed — 3 leaves) | low | no | no | — | — | zeros sub-menu — shared utility, 3 leaves; reached from util/zeros.f via debuther 20;; see 6.1/6.2/6.3 audits |
| 10;1; | TYPES d'objets (compressed — 5 leaves) | Types of objects (compressed — 5 leaves) | low | no | no | — | — | typ_objt sub-menu — shared utility, 5 leaves (POINT/LIGNE/SURFACE/VOLUME/OBJET); reached via TRMAIL path under debuther 10;; see 6.1 LEXICON-AUDIT-mail.md |
| 10;2; | VUES selon des PLANS (compressed — 6 leaves) | VIEWS from SPECIAL PLANS (compressed — 6 leaves) | low | no | no | — | — | vuesplan sub-menu — shared utility, 6 leaves; reached via TRMAIL path under debuther 10;; see 6.1 audit |
| 10;3; | TRACE des LIGNES options (compressed — 33 leaves) | Line drawing types (compressed — 33 leaves) | low | no | no | — | — | opt_lign sub-menu — shared utility, 33 leaves; reached via TRMAIL path under debuther 10;; see 6.1 audit |
| 10;4; | TRACE des SURFACES options (compressed — 40 leaves) | Surfaces drawing types (compressed — 40 leaves) | low | no | no | — | — | opt_surf sub-menu — shared utility, 40 leaves; reached via TRMAIL path under debuther 10;; see 6.1 audit |
| 10;5; | TRACE des AXES (compressed — 4 leaves) | Drawing of AXES (compressed — 4 leaves) | low | no | no | — | — | tracaxes sub-menu — shared utility, 4 leaves; reached via TRMAIL path under debuther 10;; see 6.1 audit |
| 11;6; | PLAN SECTION VOLUMES (compressed — 4 leaves) | Section Planes for VOLUMES (compressed — 4 leaves) | low | no | no | — | — | sectvopl sub-menu — shared utility, 4 leaves; reached via 3D ISO/section paths under debuther 11;; see 6.1 audit |
| 8;60;5;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleur0 sub-menu — shared utility, 16 colour names; reached from ther/trflux.f and similar tracflux/tracgrad/tracerrt parents under code 8;; see 6.1 LEXICON-AUDIT-mail.md rows for full enumeration |
| 8;60;7;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleurs sub-menu — shared utility, 16 colour names; reached from tracflux/tracgrad/tracerrt 7; parents under code 8;; see 6.1 audit |
| 8;60;6;1; | TRACE des TRAITS (compressed — 3 leaves) | Types of Line (compressed — 3 leaves) | low | no | no | — | — | typtrait sub-menu — shared utility, 3 leaves; reached from tracflux/tracgrad/tracerrt 6; parents under code 8;; see 6.1 audit |
| 6;81;1; | SOUS ESPACES (compressed — 2 leaves) | SUB SPACES (compressed — 2 leaves) | low | no | no | — | — | methvvpr sub-menu — shared utility, 2 leaves (eigen-method); reached from reso/calvvp.f:229 — ther never reaches eigenmodes solver at runtime (ppther.f:451 polynomial eigen is COMMENTED OUT) so methvvpr is reachable only theoretically via the LIMTCL grep set; see 6.2 LEXICON-AUDIT-elas.md |
| 6;82;1; | NLSE BOUNDARY CONDITIONS (compressed — cross-module shared) | NLSE BOUNDARY CONDITIONS (compressed — cross-module shared) | low | no | no | — | — | cl_nlse sub-menu — cross-module shared with future 6.5 nlse audit; reached from ther/dfnlse.f:321 (called by ther/thepolevv.f eigenvector boundary-condition path which is itself disabled per ppther.f:451 ccc); full enumeration deferred to future LEXICON-AUDIT-nlse.md (Phase 6.5); listed once for completeness |

<!-- End of audit table — validator will count rows above this line. -->

## Summary Statistics

- **Total rows:** 209
- **By frequency:**
  - `high`: 9 — `1;`, `2;`, `3;`, `4;`, `8;`, `99;` (top-level domain-promoted heat-transfer workflow) plus `2;3;` (cl_ther Dirichlet BC), `3;1;` (resothst linear steady), `4;1;` (resothin linear-time-dependent unsteady)
  - `med`: 28 — debuther codes 5;, 6;, 7;, 9;, 10;, 11;, 20; + sub-menu workflow leaves: cl_ther 2;1; 2;2; 2;8;; resothst 3;2;; resothin 4;2;; methreso 3;81;1; 3;81;2;; scheonde 5;81;1;; tempgrad 7;1; 7;2; 7;3; 7;4;; tractemp 8;20;1; 8;20;2;; tractem1 8;30;1; 8;30;2;; tracflux 8;60;1; 8;60;15;; tracgrad 8;70;1; 8;70;15;; valisoth 11;2;
  - `low`: 172 — remainder (long-tail visualization leaves + 21 shared util compressed rows + cross-module cl_nlse + dead-code resoonin + dead-code eigen family rows under code 6;)
- **By qaction:** yes count == high+med == **37** (5 toolbar=yes rows are a subset of qaction=yes; QAction set for Plan 02 = 37 rows — between 6.2 elas's 20-row set and 6.3 flui's 42-row set; smaller than flui because the visualization tree under code 8; / 11; is shallower for thermal scalar fields than for flui's vector-field cascade)
- **By toolbar:** exactly 5 `yes` — `2;`, `3;`, `4;`, `8;`, `99;`
- **Help-allowlist for Plan 03 testHelpNoQueue:** `{98;}` (NOT `{97;}` like flui — see Help-allowlist section above for the explicit hand-off)

## Top-5 Toolbar (Draft — to be ratified at Task 2)

1. `2;` — HEAT TRANSFER INPUT DATA / ENTREE des DONNEES THERMIQUES du PROBLEME — icon `ther-input-heat.svg`
2. `3;` — STEADY HEAT TRANSFER solver / CALCUL THERMIQUE STATIONNAIRE — icon `solve-heat-steady.svg`
3. `4;` — UNSTEADY HEAT TRANSFER solver / CALCUL THERMIQUE INSTATIONNAIRE d/dt — icon `solve-heat-unsteady.svg`
4. `8;` — DRAWING of TEMPERATURES and FLUX / DESSIN des TEMPERATURES et des FLUX — icon `draw-temperature.svg`
5. `99;` — SAVE DATA and QUIT / SAUVEGARDE des donnees FIN TRAITEMENT — icon `SP_DialogCloseButton`

Rationale: `2;` is the mandatory data-input gate before any solver runs (canonical first action after `1;` Object); `3;` is the canonical first solver run (steady heat baseline); `4;` is the production unsteady solver (gourd evidence + heat-transfer textbook ergonomics); `8;` is the universal result-visualization entry (every result-producing run draws temperatures and fluxes); `99;` is the shared 6.1/6.2/6.3 close-session convention. The user may at Task 2 checkpoint:
- Swap `4;` for `5;` if WAVE solver is the primary use case (gourd shows code 5; once but wave is specialty)
- Swap `3;` for `1;` if "Object first" workflow is more common (testa shows OBJECT NAME=37 hits but `1;` is implicit-mandatory rather than user-toggleable)
- Swap `8;` for `10;` if mesh inspection dominates result visualization (gourd evidence: `8;`=7 vs `10;`=2 — result drawing dominates)

Any change must keep the count at exactly 5.

## Ther-unique SVG icon set (Plan 02 ships)

Custom .svg filenames introduced by this audit that Plan 02 must ship under
`xvue/qt/resources/icons/ther/`:

- `ther-object.svg` — OBJECT NAME (code 1;)
- `ther-input-heat.svg` — HEAT TRANSFER INPUT DATA (code 2;)
- `solve-heat-steady.svg` — STEADY HEAT TRANSFER (code 3;)
- `solve-heat-unsteady.svg` — UNSTEADY HEAT TRANSFER (code 4;)
- `solve-wave.svg` — 2D WAVE SOLVER (codes 5;, 9;)
- `draw-temperature.svg` — DRAWING of TEMPERATURES + FLUX (code 8;)
- `draw-iso-surface.svg` — ISO-SURFACES (code 11;)

Total: **7 new ther-specific SVGs** (one more than 6.3 flui's 6 because
the heat-steady/heat-unsteady distinction needs visual differentiation
on the toolbar). User may consolidate further at Task 2 — e.g., a
single `solve-heat.svg` covering codes 3; and 4; drops the count to
6 new SVGs (loses toolbar visual distinction between steady and
unsteady but matches 6.3 flui's solve-stokes.svg shared between codes
5;/6;).

## Shared SVG reuses (from 6.1 mail/, 6.2 elas/)

- `mesh-draw.svg` — reused from `xvue/qt/resources/icons/mail/` for ther
  code `10;` (DRAWING of PLSVO MESHES); the qrc prefix `/xvue/qt/icons`
  resolves the mail path regardless of module context. No file copy
  needed. Same reuse pattern as 6.2 elas code `10;` and 6.3 flui code
  `19;`.
- `solve-eigen.svg` — reused from `xvue/qt/resources/icons/elas/` for ther
  codes `6;` (EIGENVALUE LINEAR solver) and `7;` (DRAWING of EIGENVECTORS).
  Note `19;` (POLYNOMIAL EIGEN) is DEAD CODE so no icon assigned.

## Cross-References

- 06.1 `LEXICON-AUDIT-mail.md` — 9-column schema template (this audit
  mirrors the structure verbatim) + leaf-level enumeration of shared
  util sub-menus (compressed here)
- 06.2 `LEXICON-AUDIT-elas.md` — same template, second-iteration
  example; demonstrated the synthetic high-number-offset prefix
  pattern (3;81;, 3;91;, 4;50;) replicated here as 3;81;, 5;81;,
  7;91;, 8;20;, 8;30;, etc.; eigvvloc/eigvvpol family verified NOT
  shared with elas (elas eigen reaches differ — this is ther-specific
  dead-code under thepolevv.f)
- 06.3 `LEXICON-AUDIT-flui.md` — third-iteration example; demonstrated
  module-aware validator + dead-code vp0flui treatment (replicated here
  for resoonin)
- 06.1 `06.1-01-SUMMARY.md` — de-duplication rule (module-specific
  fully expanded; shared util compressed to 1 row each)
- 06.3 `06.3-01-SUMMARY.md` — domain-review promotion pattern for
  sparse testa corpora (replicated here for the very-small ther
  fixture set: gourd is the ONLY .ther file)
- 06.3 `06.3-03-SUMMARY.md` §"testHelpNoQueue tightened" — Auto-fix
  Rule 1 lesson: per-module Help-allowlist must be drawn from
  LEXICON-AUDIT (NOT inherited from previous-module template)
- 06.0-UI-SPEC §Per-Module Conformance Contract — ther module-title
  Thermal/Thermique
- ROADMAP Phase 6.4 line 213 — `{File, Thermal, View, Help}` taxonomy
- `tools/validate_audit_md.py` — 9-rule schema validator; module-aware
  ICONS resolver since 6.2 Plan 02 (no validator changes needed in 6.4)
- Plan 02 consumes: rows with `qaction=yes` as the
  `registerTherActions_stub_` QAction set; rows with `toolbar=yes` as
  the toolbar append set; `icon_source` ending `.svg` (filename not
  already in `icons/mail/` or `icons/elas/` or `icons/flui/`) as the
  `xvue_icons.qrc` append set for `icons/ther/`
- Plan 03 consumes: Help-allowlist `{98;}` as the testHelpNoQueue
  QSet<QString> (NOT `{97;}` — that's flui's; differs per Auto-fix
  Rule 1 lesson from 6.3 Plan 03)
