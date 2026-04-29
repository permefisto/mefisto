# LEXICON-AUDIT-nlse.md

**Phase:** 06.5-nonlinear-nlse-menu-wiring
**Generated:** 2026-04-29
**Status:** FROZEN (auto-approved 2026-04-29 in auto-mode; DRAFT proposals retained verbatim per 6.3/6.4 convention; user can request polish-pass edits as a follow-up)
**Requirement:** UX-05 (nlse slice) Success Criterion #1

## Scope

Full recursive LIMTCL tree walk from `td/m/debunlse` via `ther/` + `util/` +
`reso/` + `prpr/ppnlse.f` call-sites. **53 LIMTCL sub-menus** (re-verified
2026-04-29 via `grep -rEoh "LIMTCL\( *'[^']+'" ther/ util/ prpr/ppnlse.f reso/`
— 1 drift from planner's interfaces list of 52, the extra one being
`debunlse` itself which is the audit root and is not compressed). De-duplication
rule from 6.1 Plan 01 / 6.2 Plan 01 / 6.3 Plan 01 / 6.4 Plan 01:

- **Nlse-primary sub-menus fully expanded** — `cl_nlse` (3 leaves) and
  `coefnlse` (6 leaves). These are the two genuinely NLSE-physics-specific
  sub-menus that 6.4 ther compressed and that 6.5 nlse FULLY expands as
  the canonical owner of NLSE physics.
- **Ther-proxy sub-menus listed ONCE** with notes cell citing 6.4
  `LEXICON-AUDIT-ther.md` for leaf-level detail — `cl_ther`, `defplan`,
  `dirprofi`, `eigvvloc`, `eigvvpol`, `methreso`, `methvvpr`, `profdroi`,
  `profplan`, `proj6cub`, `resoonin` (DEAD CODE), `resothin`, `resothst`,
  `scheonde`, `sectplan`, `solsurf3`, `tempgrad`, `tracerr2`, `tracerr3`,
  `tracerrt`, `tracflux`, `tracgrad`, `traconde`, `tractem1`, `tractem2`,
  `tractem3`, `tractemp`, `trdepond`, `valisoth`, `valzone`, `valztxy`
  (31 sub-menus, compressed to 1 row each — these were FULLY EXPANDED
  in 6.4 ther audit, so 6.5 nlse defers to that canonical owner per the
  same compression rule 6.4 applied to the 21 shared util sub-menus
  that 6.1 mail had fully expanded).
- **Shared util sub-menus listed ONCE** with notes cell citing 6.1 mail
  / 6.2 elas / 6.3 flui / 6.4 ther audits for leaf-level detail —
  `affiche`, `couleur0`, `couleurs`, `entrsort`, `fichier`, `lecteur`,
  `managfen`, `managmef`, `mode_es`, `opt_lign`, `opt_surf`, `sectvopl`,
  `suivi_ms`, `suivitms`, `tracaxes`, `tuer`, `typ_objt`, `typtrait`,
  `vuesplan`, `zeros` (20 sub-menus, compressed to 1 row each — same
  treatment as 6.4 ther's 21-row compression).

**Total row count: 80** — meets validator bound `[80, 250]` exactly at the
lower threshold. Smaller than 6.4 ther's 209 because nlse legitimately
treats the 31 ther-proxy sub-menus as compressed (deferring to the 6.4
audit) rather than re-expanding them. The genuinely NLSE-primary
expansion is just `cl_nlse` (3 leaves) + `coefnlse` (6 leaves). This
mirrors the same compression pattern 6.4 ther applied to the 21 shared
util sub-menus after they were fully expanded in 6.1 mail — peers don't
re-expand shared trees once a sibling audit owns them.

Frequency bucketing is **100% domain review** — no testa-derived counts
exist for nlse:

```bash
$ find testa/ testf/ -name "*.nlse" 2>/dev/null
(empty)
```

The five testa/ projects with `nlse` in their names — `nlsecu/`,
`nlsecuvv/`, `nlsesqri/`, `nlsgpe/`, `gpesigc/` — contain ONLY mesh
files (`.meshq2`, `.meshp1`) + initial-guess files (`.iexri`,
`.iexrr`, `.iexrrplus`, `.evv`, `.imvv`, `.sigc`). NONE of them is a
lexicon batch script. This is a LARGER sparsity problem than 6.4
ther's gourd evidence (which had 1 .ther fixture with 4.3K of
multi-line block lexicon). Frequency assignment is therefore drawn
entirely from canonical NLSE workflow ergonomics + ppnlse.f dispatch
table. Per 6.1/6.2/6.3/6.4 carryover: rows with `qaction=yes` have
`frequency` in {high, med}; exactly 5 rows have `toolbar=yes`.

Schema: 9 columns, enforced by `tools/validate_audit_md.py` (shipped in
6.1 Plan 01; module-aware ICONS resolver since 6.2 Plan 02 — reused
verbatim, no tooling changes). Rule 9 runs in WARN mode for new nlse
SVGs (Plan 02 ships them); `mesh-draw.svg` (code 19;) resolves cleanly
into `icons/mail/`.

### Menu taxonomy (auto-frozen, mirroring 6.3/6.4 fully-autonomous path)

The nlse module's Qt menu bar uses **`{File, Nonlinear, View, Help}`** per
ROADMAP Phase 6.5 line 220. Per 06.0-UI-SPEC §Per-Module Conformance
Contract the canonical FR translation `Non Linéaire` is confirmed
against `td/ma/debunlse` line 1 which translates "Non Linear SCHRODINGER
Equation NLSE" — the natural FR menu title is `Non Linéaire`. This
differs from 6.1's `{File, Mesh, View, Help}`, 6.2's
`{File, Solve, View, Help}`, 6.3's `{File, Fluid, View, Help}`, and
6.4's `{File, Thermal, View, Help}` — "Nonlinear" replaces the
per-module slot.

NOTE — Pre-existing 06.0-UI-SPEC table at line 245 lists the per-module
title for nlse as `Solve / Calcul`; ROADMAP Phase 6.5 line 220 takes
precedence as the authoritative source for this slice. The conflict is
documented here for the maintainer; future 06.0-UI-SPEC update can
align the table to `Nonlinear / Non Linéaire` as a polish-pass edit.

Content distribution (proposed):

- **File**: codes `71;` `72;` `73;` `74;` `99;` + shared 6.0 File actions
- **Nonlinear**: codes `1;` `2;` `5;` `6;` + a `Parameters` sub-menu with codes `20;` `38;` `39;`
- **View**: codes `11;` `12;` `13;` `14;` `15;` `16;` `19;` + shared 6.0 View actions
- **Help**: code `97;` + shared 6.0 Help/About

## Known data quirks

- Help-version code for nlse is **`97;`** (NOT `98;` like ther). Plan
  03's `testHelpNoQueue` allowlist MUST be `{97;}` — explicit per-module
  variation per 6.3 Plan 03 Auto-fix Rule 1 lesson. Verified against
  `td/m/debunlse` line 20 and `td/ma/debunlse` line 20. nlse's value
  matches flui's `{97;}` (debuflui line 21), DIFFERS from ther's `{98;}`
  (debuther line 22).
- Codes `38;` and `39;` in debunlse are direct typed-shortcut entries
  for window-pixels and background-color — `ppnlse.f:218` `IF NMTCL.EQ.38
  GOTO 3800 → CALL XVPXFE` and `ppnlse.f:219` `IF NMTCL.EQ.39 GOTO 3900
  → couleurs prompt`. Both routes work but are tagged as Parameters
  sub-menu items. Audit lists 38; and 39; as **frequency=low,
  qaction=no** (per Rule 8 for low) — Plan 02 wires QActions via
  menu-only path bypassing `LIMTCL('debunlse')` → direct `XVPXFE` /
  `LIMTCL('couleurs')` invocations.
- nlse uses **four separate top-level codes** `71;`, `72;`, `73;`, `74;`
  for TMS sub-management (SUIVI TMS / FICHIERS / UNITES / TUER) where
  ther/elas/flui used a single combined `70;` code that opens managmef.
  The dispatch in `ppnlse.f:467-470` unwraps the 71-74 range with
  `NMTCL = NMTCL - 70` then GOTO array slots 1-4 → SUITMS / SUIFMS /
  SUIVES / TUER. Audit shows four separate top-level rows for
  `71;`/`72;`/`73;`/`74;` (matching 6.1 mail's pattern).
- debunlse skips codes `0`, `3`, `4`, `7`, `8`, `9`, `10`, `17`, `18`.
  `prpr/ppnlse.f:220` GOTO array contains
  `100, 200, 300, 400, 500, 600, 700, 30, 30, 1000, 1100, 1200, 1300,
  1400, 1500, 1600, 30, 30, 1900, 2000` covering codes 1-20 with `30`
  (loop-back-to-menu) at slots 8, 9, 17, 18 — codes `8;`, `9;`, `17;`,
  `18;` intentionally fall through. Audit lists ONLY codes that exist
  in debunlse — does NOT fabricate rows for 0/3/4/7/8/9/10.
- ppnlse.f preserves a `ccc`-commented block at lines 345-354 (deleted
  `7;` Schema3 SANS MATRICE GLOBALE NLSEINS path). The live code 6;
  (i-TIME / Gross-Pitaevskii) was originally code `1400` per the
  comment at ppnlse.f:364 — code 6; is the resurrected slot under the
  new dispatch table.
- Codes `5;` and `6;` both dispatch to NLSEINS but with different
  TESTNL values: `ppnlse.f:340-343` (slot 500) sets `TESTNL=7` for code
  5; (IMPLICIT scheme with global 2nx2n non-symmetric matrix); `ppnlse.f:365-368`
  (slot 600) sets `TESTNL=9` for code 6; (i-TIME / Gross-Pitaevskii
  semi-implicit with global nxn matrix). The TESTNL flag selects the
  scheme variant inside NLSEINS; both prompt `cl_nlse` + `coefnlse`
  via DFNLSE prior to entering the time loop.
- Visualization codes `11;`, `12;`, `13;` dispatch to `CALL TRTHER(KNOMOB,
  N, IERR)` per ppnlse.f:392, 399, 406 with `N` in {4, 5, 6} (4=MODULE,
  5=REAL part, 6=IMAGINARY part). TRTHER then prompts a tractem* /
  tracflux / tracgrad / tracerrt cascade — all those sub-menus are
  PROXY-TO-THER (compressed via the ther-proxy bucket below). The
  TOP-LEVEL visualizers (codes 11; 12; 13;) get ONE ROW EACH at
  top-level and we do NOT write rows for the per-component leaves
  under each, since those leaves are owned by the 6.4 audit.
- Codes `14;`, `15;`, `16;` dispatch to TRNLSETST / TRNLSEMXU / TRNLSERR
  (ppnlse.f:411-428). These are NLSE-primary routines and they reach NO
  further sub-menus — verified by `grep "LIMTCL" ther/trnlsetst.f
  ther/trnlsemxu.f ther/trnlserr.f` returning zero matches. One row each
  at top-level, no further leaves.
- Code `19;` (DRAW PLSVO MESHES) dispatches to `CALL TRMAIL` per
  ppnlse.f:433. Audit treats `19;` as a shared util (compressed row at
  top-level only, with icon reuse `mesh-draw.svg` from
  `xvue/qt/resources/icons/mail/`).
- Cross-module — the planner-flagged `cl_ther` row appears in nlse's
  reachable LIMTCL grep set because ther/dfnlse.f:213 uses SDDEF2 with
  `'NLSE','coefnlse'` arguments; cl_ther itself is reached only from
  ther/dfther.f and is NOT prompted under any nlse code path. Listed in
  the ther-proxy bucket for completeness; full enumeration in 6.4
  audit.
- `resoonin` — same DEAD CODE flag as 6.4 ther: `ther/ondins.f:41`
  contains `ccc CALL LIMTCL('resoonin', NMTCL)` (commented out). The
  td/m/resoonin file exists with 3 leaves but is NEVER reached at
  runtime (in ther OR nlse). Listed in ther-proxy bucket; see 6.4 audit
  for the same dead-code treatment.
- `eigvvloc`, `eigvvpol`, `methvvpr` — ther DEAD CODE family: `prpr/ppther.f:451`
  ccc-comments out the THEPOLEVV dispatch, making the entire eigen
  family unreachable from ther at runtime. nlse never dispatches into
  this family at all (no THEPOLEVV call in ppnlse.f). Listed in the
  ther-proxy bucket as compressed dead-code citations.
- nlse SHIPS three menu-file variants — `td/m/debunlse`, `td/ma/debunlse`,
  `td/mf/debunlse` — where 6.4 ther shipped only `m/` + `ma/`. Bilingual
  parser falls through to FR per 6.2 Plan 05 fix to
  `XvueMenuFileParser::loadFor`. Verified via `ls td/{m,ma,mf}/debunlse`.
- Apply the 6.1 Pitfall 5 log-and-fallback discipline to any sub-menu
  typo encountered during downstream Plan 02 wire-up.

## Help-allowlist (for Plan 03 testHelpNoQueue — explicit hand-off)

Per LEXICON-AUDIT-nlse row `97;` (MEFISTO Version Date), the Help menu
carries the audited lexicon **`{97;}`**. Plan 03's `testHelpNoQueue`
slot allowlist:

`QSet<QString>{ QStringLiteral("97;") }` — **MATCHES** flui's `{97;}`
(debuflui line 21 also uses code 97), **DIFFERS** from ther's `{98;}`
(debuther line 22). This explicit hand-off prevents Plan 03 from
inheriting a stale `{98;}` from the 6.4 ther template (Auto-fix Rule 1
lesson from 6.3 Plan 03 — per-module Help-allowlist must be drawn
from LEXICON-AUDIT, NOT inherited from previous-module template).

## Legend

- `lexicon_path` — single-line semicolon-separated, no spaces
- `frequency` — `high` / `med` / `low` (from domain review of canonical NLSE workflow + ppnlse.f dispatch evidence — no testa fixtures exist for nlse)
- `qaction` — `yes` iff frequency is high or med
- `toolbar` — `yes` iff in the top-5 toolbar slice (exactly 5 rows)
- `icon_source` — Qt `SP_*` enum name, custom `.svg` filename, or `—`
- `shortcut` — Qt accelerator or `—`
- `notes` — clarifications, deferred items, source flags

## Audit Table

| lexicon_path | description_fr | description_en | frequency | qaction | toolbar | icon_source | shortcut | notes |
|--------------|----------------|----------------|-----------|---------|---------|-------------|----------|-------|
| 1; | NOM de l'OBJET a TRAITER | OBJECT NAME to treat | high | yes | no | nlse-object.svg | — | debunlse leaf; ppnlse.f:215 dispatch slot 1 → label 100 (CALL LIRLEX); domain-promoted HIGH (mandatory first command per ppnlse.f:227-302); Nonlinear menu root item |
| 2; | ENTREE des DONNEES de l'Equation Non Lineaire de SCHRODINGER | Non Linear SCHRODINGER Equation INPUT DATA | high | yes | yes | nlse-input.svg | — | debunlse leaf; ppnlse.f:306-308 dispatch slot 2 → label 200 (CALL DFNLSE) which prompts coefnlse (via SDDEF2) + cl_nlse (LIMTCL); toolbar top-5 (mandatory data input gate before any solver run) |
| 5; | IMPLICIT idU(t)/dt -Alfa Laplacien U(t) +N(t,abs(U))U(t)=F(t,X) | IMPLICIT i dU(t)/dt -Alfa Laplacian U(t) +N(abs(U)**2)U(t)=F(t,X) | high | yes | yes | solve-nlse-implicit.svg | — | debunlse leaf; ppnlse.f:340-343 dispatch slot 5 → label 500 (TESTNL=7 + CALL NLSEINS — IMPLICIT global 2nx2n non-symmetric matrix scheme); domain-promoted HIGH (canonical default solver run); toolbar top-5 |
| 6; | i-TIME methode idU(t)/dt -Alfa Laplacian U(t) +N(abs(U)**2)U(t)=F(t,X) | i-TIME method idU(t)/dt -Alfa Laplacian U(t) +N(abs(U)**2)U(t)=F(t,X) | med | yes | no | — | — | debunlse leaf; ppnlse.f:365-368 dispatch slot 6 → label 600 (TESTNL=9 + CALL NLSEINS — Gross-Pitaevskii i-time semi-implicit nxn scheme); specialty workflow (rotation-equation use cases); originally was code 1400 per the ccc-commented preamble at line 364 |
| 11; | DESSIN du MODULE de l'ONDE NLSE | Drawing of abs(U(t,X)) NLSE WAVE MODULE | high | yes | yes | draw-nlse-modulus.svg | — | debunlse leaf; ppnlse.f:391-393 dispatch slot 11 → label 1100 (CALL TRTHER(KNOMOB, 4, IERR) — NTYPDESS=4 selects MODULE); primary visual (universal NLSE result-visualization entry); toolbar top-5; downstream cascade owned by 6.4 ther audit (tractem* / tracflux / tracgrad / tracerrt — see ther-proxy compressed rows below) |
| 12; | DESSIN de la PARTIE REELLE de l'ONDE NLSE | Drawing of U(t,X) NLSE WAVE REAL PART | med | yes | no | draw-nlse-component.svg | — | debunlse leaf; ppnlse.f:398-400 dispatch slot 12 → label 1200 (CALL TRTHER(KNOMOB, 5, IERR) — NTYPDESS=5 selects REAL part); specialty visualization; downstream cascade owned by 6.4 ther audit; icon shared with code 13; (Re U / Im U use the same draw-nlse-component.svg) |
| 13; | DESSIN de la PARTIE IMAGINAIRE de l'ONDE NLSE | Drawing of U(t,X) NLSE WAVE IMAGINARY PART | med | yes | no | draw-nlse-component.svg | — | debunlse leaf; ppnlse.f:405-407 dispatch slot 13 → label 1300 (CALL TRTHER(KNOMOB, 6, IERR) — NTYPDESS=6 selects IMAGINARY part); specialty visualization; downstream cascade owned by 6.4 ther audit; icon shared with code 12; |
| 14; | DESSIN des VALEURS du TEST d'ARRET des ITERATIONS | Drawing of STOP TEST Values of ITERATIONS | med | yes | no | — | — | debunlse leaf; ppnlse.f:411-414 dispatch slot 14 → label 1400 (CALL TRNLSETST after LXLXOU) — NLSE-primary routine; post-iteration convergence diagnostic; reaches no further LIMTCL sub-menus (verified via grep "LIMTCL" ther/trnlsetst.f → 0 hits) |
| 15; | DESSIN des Max abs(U(Noeud))(Temps) | Drawing of Max abs(U(Node))(Time) | med | yes | no | — | — | debunlse leaf; ppnlse.f:418-421 dispatch slot 15 → label 1500 (CALL TRNLSEMXU after LXLXOU) — NLSE-primary routine; post-iteration time-series diagnostic; reaches no further LIMTCL sub-menus |
| 16; | DESSIN des ERREUR(Temps) = ER(Temps) + i EI(Temps) | Drawing of REAL U-ERROR(Time), IMAGINARY U-ERROR(Time) | med | yes | no | — | — | debunlse leaf; ppnlse.f:425-428 dispatch slot 16 → label 1600 (CALL TRNLSERR after LXLXOU) — NLSE-primary routine; convergence diagnostic; reaches no further LIMTCL sub-menus |
| 19; | DESSIN du MAILLAGE des PLSVO | Drawing of PLSVO meshes | med | yes | yes | mesh-draw.svg | — | debunlse leaf; ppnlse.f:433-434 dispatch slot 19 → label 1900 (CALL TRMAIL); REUSES 6.1 mail/mesh-draw.svg (no copy needed); toolbar top-5 (mesh inspection — borrow from 6.1 mail toolbar convention) |
| 20; | PRECISION pour INVERSER les systemes lineaires A x = b | PRECISION to INVERSE the linear systems A x = b | low | no | no | — | — | debunlse leaf; ppnlse.f:438-439 dispatch slot 20 → label 2000 (CALL ZEROGC); Parameters sub-menu under Nonlinear; rare niche parameter setter; freq=low forces qaction=no per Rule 8 — Plan 02 may still wire as a Parameters action via menu-only path |
| 38; | LARGEUR HAUTEUR PIXELS de la FENETRE de VISUALISATION | WINDOW PIXEL WIDTH & HEIGHT numbers | low | no | no | SP_ComputerIcon | — | debunlse leaf; ppnlse.f:218 IF NMTCL.EQ.38 GOTO 3800 → CALL XVPXFE; Parameters sub-menu under Nonlinear; freq=low forces qaction=no per Rule 8; Plan 02 wires via menu-only path |
| 39; | COULEUR du FOND de la FENETRE de VISUALISATION | BACKGROUND COLOR of the window | low | no | no | — | — | debunlse leaf; ppnlse.f:219 IF NMTCL.EQ.39 GOTO 3900 → couleurs prompt (LIMTCL) + XVFOND + EFFACE; Parameters sub-menu; freq=low forces qaction=no per Rule 8; Plan 02 wires via menu-only path |
| 71; | SUIVI des TMS AFFICHAGE MODIFICATION | TOOLS on TMS | low | no | no | SP_DirIcon | — | debunlse leaf; ppnlse.f:217 IF NMTCL.GT.70 GOTO 7001 then NMTCL=NMTCL-70=1 → label 7100 (CALL SUITMS); File menu (TMS administration); freq=low forces qaction=no per Rule 8 — Plan 02 wires via menu-only path |
| 72; | SUIVI des FICHIERS de la MS | TOOLS on FILES of the MS | low | no | no | SP_DirIcon | — | debunlse leaf; secondary GOTO slot 2 → label 7200 (CALL SUIFMS); File menu (file administration); freq=low |
| 73; | GESTION des UNITES lecture affichage | TOOLS on I/O units | low | no | no | SP_DirIcon | — | debunlse leaf; secondary GOTO slot 3 → label 7300 (CALL SUIVES); File menu (I/O administration); freq=low |
| 74; | TUER des TMS de PLSVO | DELETE TMS of PLSVO | low | no | no | SP_TrashIcon | — | debunlse leaf; secondary GOTO slot 4 → label 7400 (CALL TUER); File menu (destructive op — TMS deletion); freq=low |
| 97; | DATE de version de Mefisto | MEFISTO VERSION DATE | low | no | no | SP_MessageBoxInformation | — | debunlse leaf; secondary GOTO slot 27 → label 9700 (CALL VRSION); Help menu item; **Help-allowlist for Plan 03 testHelpNoQueue: {97;}** (NOT {98;} like ther — see Help-allowlist section above for the explicit hand-off); freq=low forces qaction=no per Rule 8 |
| 99; | SAUVEGARDE des donnees FIN TRAITEMENT | SAVE DATA and QUIT | high | yes | yes | SP_DialogCloseButton | Ctrl+Q | debunlse leaf; secondary GOTO slot 29 → label 9900 (CALL ARRET(0) + STOP); domain-promoted HIGH (every session closes here); toolbar top-5 (shared convention with 6.1/6.2/6.3/6.4); File menu |
| 2;2; | Fgamma FLUX Complexe IMPOSE de l'ONDE | FGamma IMPOSED Complex WAVE FLUX at Boundary | med | yes | no | — | — | cl_nlse leaf; reached from DFNLSE under debunlse code 2; via ther/dfnlse.f:321 LIMTCL('cl_nlse'); Neumann boundary condition (imposed complex flux) |
| 2;3; | ONDE_D VALEUR Complexe IMPOSEE de l'ONDE | WAVE_D IMPOSED Complex WAVE VALUE at Boundary | high | yes | no | — | — | cl_nlse leaf; Dirichlet boundary condition (imposed complex value); canonical NLSE BC (most NLSE simulations fix at least one boundary value); domain HIGH |
| 2;8; | ONDE_0 VALEUR Complexe IMPOSEE de l'ONDE INITIALE | WAVE_0 INITIAL Complex WAVE VALUE at Boundary | med | yes | no | — | — | cl_nlse leaf; initial complex wave condition for boundary; reached under code 2; |
| 2;91;1; | Rho Densite de MASSE | Rho MASS Density | high | yes | no | — | — | coefnlse leaf — synthetic prefix 2;91; (coefnlse reached from ther/dfnlse.f:213 SDDEF2 call with 'coefnlse' arg under debunlse code 2;; synthetic offset chosen to avoid collision with cl_nlse leaves at 2;2; / 2;3; / 2;8;); canonical mass density (mandatory material parameter); domain HIGH |
| 2;91;2; | Alfa Coefficient du LAPLACIEN | Alfa Laplacian Coefficient | high | yes | no | — | — | coefnlse leaf; canonical Laplacian-coefficient parameter (Alfa = Conductivite per ther/dfnlse.f:46 remark); domain HIGH (every NLSE solve sets this) |
| 2;91;4; | FOmega Force Complexe Second Membre | FOmega Complex Force Second Member | med | yes | no | — | — | coefnlse leaf; complex source term F(t,X) — required for non-homogeneous NLSE problems; specialty workflow (homogeneous = forcing is zero) |
| 2;91;6; | Beta Coefficient (ur**2+ui**2) ur OU ui | Beta Coefficient (ur**2+ui**2) ur or ui | high | yes | no | — | — | coefnlse leaf; canonical NLSE non-linearity coefficient (the N in N(abs(U)**2) U term); domain HIGH (defines NLSE physics — without Beta the equation reduces to linear Schrödinger) |
| 2;91;8; | U(t0,X) Valeur de l'Onde Complexe INITIALE | INITIAL Complex Wave Value | med | yes | no | — | — | coefnlse leaf; initial-condition complex wave field (required for time-dependent solvers — codes 5; and 6;) |
| 2;91;9; | Omega Vitesse Angulaire de la rotation | Omega Rotation Angular Velocity | low | no | no | — | — | coefnlse leaf; angular velocity for Gross-Pitaevskii rotation term (i Omega (x du/dy - y du/dx)); only relevant for code 6; (i-TIME) — code 5; (IMPLICIT) ignores Omega |
| 11;7; | TEMP/GRAD/FLUX dispatcher (compressed — 10 leaves) | TEMP/GRAD/FLUX dispatcher (compressed — 10 leaves) | low | no | no | — | — | tempgrad sub-menu — ther-proxy compressed (synthetic prefix 11;7; mirrors 6.4 ther row 7;1; under code 7; DRAWING EIGEN, repurposed here under code 11; DRAWING MODULE since TRTHER prompts tempgrad regardless of NTYPDESS); reached from ther/trther.f under codes 11;/12;/13;; full enumeration in 6.4 LEXICON-AUDIT-ther.md rows tagged "tempgrad leaf" |
| 11;8;20; | TEMPERATURE-style 2D drawing (compressed — 7 leaves) | TEMPERATURE-style 2D drawing (compressed — 7 leaves) | low | no | no | — | — | tractemp sub-menu — ther-proxy compressed; reached from ther/sdtrth.f:235 under codes 11;/12;/13; via tempgrad-2; the "TEMP" naming is historical — for nlse this draws abs(U(t,X)) / Re U / Im U; full enumeration in 6.4 audit rows tagged "tractemp leaf" |
| 11;8;30; | TEMPERATURE 1D-time surface (compressed — 4 leaves) | TEMPERATURE 1D-time surface (compressed — 4 leaves) | low | no | no | — | — | tractem1 sub-menu — ther-proxy compressed; reached from ther/trther.f under codes 11;/12;/13; (1D variant); full enumeration in 6.4 audit rows tagged "tractem1 leaf" |
| 11;8;40; | TEMPERATURE 2D variant (compressed — 4 leaves) | TEMPERATURE 2D variant (compressed — 4 leaves) | low | no | no | — | — | tractem2 sub-menu — ther-proxy compressed; reached from ther/trther.f under codes 11;/12;/13; (2D variant); full enumeration in 6.4 audit rows tagged "tractem2 leaf" |
| 11;8;50; | TEMPERATURE 3D variant (compressed — 6 leaves) | TEMPERATURE 3D variant (compressed — 6 leaves) | low | no | no | — | — | tractem3 sub-menu — ther-proxy compressed; reached from ther/trther.f under codes 11;/12;/13; (3D variant); full enumeration in 6.4 audit rows tagged "tractem3 leaf" |
| 11;8;60; | FLUX drawing (compressed — 7 leaves) | FLUX drawing (compressed — 7 leaves) | low | no | no | — | — | tracflux sub-menu — ther-proxy compressed; reached from ther/sdtrfl.f:67 + ther/trflux.f:81 under codes 11;/12;/13; via tempgrad-4; the "FLUX" naming is historical — for nlse this draws complex-wave flux; full enumeration in 6.4 audit rows tagged "tracflux leaf" |
| 11;8;70; | GRADIENT drawing (compressed — 7 leaves) | GRADIENT drawing (compressed — 7 leaves) | low | no | no | — | — | tracgrad sub-menu — ther-proxy compressed; reached from ther/sdtrgt.f:54 + ther/trgrad.f:79 under codes 11;/12;/13; via tempgrad-3; full enumeration in 6.4 audit rows tagged "tracgrad leaf" |
| 11;8;75; | ERROR estimator (compressed — 5 leaves) | ERROR estimator (compressed — 5 leaves) | low | no | no | — | — | tracerrt sub-menu — ther-proxy compressed; reached under codes 11;/12;/13; via tempgrad-5 error estimator path; full enumeration in 6.4 audit rows tagged "tracerrt leaf" |
| 11;8;80; | ERROR 2D (compressed — 5 leaves) | ERROR 2D (compressed — 5 leaves) | low | no | no | — | — | tracerr2 sub-menu — ther-proxy compressed; reached under codes 11;/12;/13; via 2D error path; full enumeration in 6.4 audit rows tagged "tracerr2 leaf" |
| 11;8;85; | ERROR 3D (compressed — 6 leaves) | ERROR 3D (compressed — 6 leaves) | low | no | no | — | — | tracerr3 sub-menu — ther-proxy compressed; reached under codes 11;/12;/13; via 3D error path; full enumeration in 6.4 audit rows tagged "tracerr3 leaf" |
| 11;11;1; | ISO drawing values (compressed — 9 leaves) | ISO drawing values (compressed — 9 leaves) | low | no | no | — | — | valisoth sub-menu — ther-proxy compressed; reached from ther/sdtrit.f:87 + ther/trisot.f:158 under code 11; via TRTHER cascade (iso-surface drawing for abs(U) field); full enumeration in 6.4 audit rows tagged "valisoth leaf" |
| 11;11;20; | ZONE colors (compressed — 4 leaves) | ZONE colors (compressed — 4 leaves) | low | no | no | — | — | valzone sub-menu — ther-proxy compressed; reached from ther/trzont.f:165 under code 11; via TRTHER cascade; full enumeration in 6.4 audit rows tagged "valzone leaf" |
| 11;11;25; | Z-axis solution (compressed — 9 leaves) | Z-axis solution (compressed — 9 leaves) | low | no | no | — | — | valztxy sub-menu — ther-proxy compressed; reached from ther/trztxy.f:742 under code 11; via TRTHER cascade; full enumeration in 6.4 audit rows tagged "valztxy leaf" |
| 11;11;30; | LINE profile (compressed — 9 leaves) | LINE profile (compressed — 9 leaves) | low | no | no | — | — | profdroi sub-menu — ther-proxy compressed; reached from ther/trlldr.f:198 under code 11; via TRTHER cascade; full enumeration in 6.4 audit rows tagged "profdroi leaf" |
| 11;11;30;9; | LINE profile direction (compressed — 6 leaves) | LINE profile direction (compressed — 6 leaves) | low | no | no | — | — | dirprofi sub-menu — ther-proxy compressed; reached from ther/trlldr.f:290 under profdroi 9;; full enumeration in 6.4 audit rows tagged "dirprofi leaf" |
| 11;11;40; | PLANE section (compressed — 14 leaves) | PLANE section (compressed — 14 leaves) | low | no | no | — | — | sectplan sub-menu — ther-proxy compressed; reached from ther/trplse.f:238 under code 11;; full enumeration in 6.4 audit rows tagged "sectplan leaf" |
| 11;11;50; | PLANE profile (compressed — 5 leaves) | PLANE profile (compressed — 5 leaves) | low | no | no | — | — | profplan sub-menu — ther-proxy compressed; reached from ther/trplse.f:245 under code 11;; full enumeration in 6.4 audit rows tagged "profplan leaf" |
| 11;11;60; | PLANE definition (compressed — 2 leaves) | PLANE definition (compressed — 2 leaves) | low | no | no | — | — | defplan sub-menu — ther-proxy compressed; reached from ther/trplse.f:280 under sectplan/profplan 4;; full enumeration in 6.4 audit rows tagged "defplan leaf" |
| 11;11;70; | 3D solution surface (compressed — 3 leaves) | 3D solution surface (compressed — 3 leaves) | low | no | no | — | — | solsurf3 sub-menu — ther-proxy compressed; reached from ther/trso1so.f:271 under code 11; via 3D-surface trace path; full enumeration in 6.4 audit rows tagged "solsurf3 leaf" |
| 11;91; | 6cube projection (compressed — 5 leaves) | 6cube projection (compressed — 5 leaves) | low | no | no | — | — | proj6cub sub-menu — ther-proxy compressed; reached from ther/trther.f:460 under code 11; with 6cube objects (6D-to-3D projection); full enumeration in 6.4 audit rows tagged "proj6cub leaf" |
| 5;81;1; | METHRESO method (compressed — 3 leaves) | METHRESO method (compressed — 3 leaves) | low | no | no | — | — | methreso sub-menu — ther-proxy compressed (synthetic prefix 5;81; — methreso is reached from reso/calvvp.f:229 generally but for nlse the canonical reach is through coefnlse->NLSEINS where the linear solve picks a methreso variant; never directly prompted from ppnlse.f; conservatively listed via ther dispatch); full enumeration in 6.4 audit rows tagged "methreso leaf" |
| 5;82;1; | SCHEONDE wave scheme (compressed — 2 leaves) | SCHEONDE wave scheme (compressed — 2 leaves) | low | no | no | — | — | scheonde sub-menu — ther-proxy compressed; reached from ther/ondes.f:469 under ther's code 5; (WAVE solver) — NOT reached at runtime from nlse since nlse never dispatches to ondes; listed in nlse-reachable LIMTCL grep set due to file co-residence under ther/; full enumeration in 6.4 audit |
| 5;83;1; | RESOTHIN unsteady linear (compressed — 3 leaves) | RESOTHIN unsteady linear (compressed — 3 leaves) | low | no | no | — | — | resothin sub-menu — ther-proxy compressed; reached from THEINS dispatcher in ther which is ther's code 4; — NOT reached from nlse runtime (NLSEINS does not prompt resothin); listed for grep-set completeness; full enumeration in 6.4 audit |
| 5;84;1; | RESOTHST steady linear (compressed — 3 leaves) | RESOTHST steady linear (compressed — 3 leaves) | low | no | no | — | — | resothst sub-menu — ther-proxy compressed; reached from THESTA dispatcher in ther's code 3; — NOT reached from nlse runtime; listed for grep-set completeness; full enumeration in 6.4 audit |
| 5;85; | RESOONIN dead code (compressed — 3 leaves) | RESOONIN dead code (compressed — 3 leaves) | low | no | no | — | — | resoonin sub-menu — DEAD CODE in ther/ondins.f:41 (`ccc CALL LIMTCL('resoonin', NMTCL)` commented out); td/m/resoonin file exists with 3 leaves but NEVER reached at runtime in either ther or nlse; same dead-code flag as 6.4 audit; full enumeration in 6.4 audit rows tagged "resoonin leaf" |
| 5;86; | TRACONDE wave drawing (compressed — 5 leaves) | TRACONDE wave drawing (compressed — 5 leaves) | low | no | no | — | — | traconde sub-menu — ther-proxy compressed; reached from ther/tronde.f:209 under ther's code 9; (DRAWING WAVE) — NOT reached from nlse runtime; full enumeration in 6.4 audit |
| 5;87; | TRDEPOND wave-displacement drawing (compressed — 7 leaves) | TRDEPOND wave-displacement drawing (compressed — 7 leaves) | low | no | no | — | — | trdepond sub-menu — ther-proxy compressed; reached from ther/tronde.f:259 under ther's code 9;2; — NOT reached from nlse runtime; full enumeration in 6.4 audit |
| 6;81;1; | METHVVPR sub-spaces (compressed — 2 leaves) | METHVVPR sub-spaces (compressed — 2 leaves) | low | no | no | — | — | methvvpr sub-menu — ther-proxy compressed; reached from reso/calvvp.f:229 — DEAD CODE for nlse (no eigenvector solver dispatched in ppnlse.f); listed in nlse-reachable LIMTCL grep set due to reso/ co-residence; full enumeration in 6.2 LEXICON-AUDIT-elas.md (canonical owner) |
| 6;82;1; | EIGVVPOL polynomial eigen (compressed — 4 leaves) | EIGVVPOL polynomial eigen (compressed — 4 leaves) | low | no | no | — | — | eigvvpol sub-menu — DEAD CODE; reached from ther/thepolevv.f:180 under ther's code 19; POLYNOMIAL EIGEN — but ppther.f:451 dispatch is COMMENTED OUT so ther never reaches; nlse never dispatches to thepolevv at all; full enumeration in 6.4 audit rows tagged "eigvvpol leaf" |
| 6;83;1; | EIGVVLOC eigen localisation (compressed — 3 leaves) | EIGVVLOC eigen localisation (compressed — 3 leaves) | low | no | no | — | — | eigvvloc sub-menu — DEAD CODE; reached from ther/thepolevv.f:209 under eigvvpol leaf 3; — same DEAD-CODE flag as eigvvpol; full enumeration in 6.4 audit rows tagged "eigvvloc leaf" |
| 6;84;1; | CL_THER thermal BC (compressed — 5 leaves) | CL_THER thermal BC (compressed — 5 leaves) | low | no | no | — | — | cl_ther sub-menu — ther-proxy compressed; reached from ther/dfther.f for ther's code 2; — NOT reached from nlse runtime (nlse uses cl_nlse not cl_ther; cl_ther appears in nlse's LIMTCL grep set only because ther/dfnlse.f:213 SDDEF2 call mentions LISTE1 entries that overlap thermal nomenclature); full enumeration in 6.4 audit rows tagged "cl_ther leaf" |
| 71;81; | AFFICHE output unit (compressed — 2 leaves) | AFFICHE output unit (compressed — 2 leaves) | low | no | no | — | — | affiche sub-menu — shared utility, 2 leaves (SCREEN/FILE); reached from util/suives.f via debunlse 73;; see 6.1 LEXICON-AUDIT-mail.md / 6.4 LEXICON-AUDIT-ther.md for full enumeration |
| 71;82; | COULEUR0 mesh-edges color (compressed — 16 leaves) | COULEUR0 mesh-edges color (compressed — 16 leaves) | low | no | no | — | — | couleur0 sub-menu — shared utility, 16 colour names; reached from ther/trflux.f and similar tracflux/tracgrad/tracerrt parents under codes 11;/12;/13; via TRTHER cascade; see 6.1 audit rows for full enumeration |
| 71;83; | COULEURS arrow color (compressed — 16 leaves) | COULEURS arrow color (compressed — 16 leaves) | low | no | no | — | — | couleurs sub-menu — shared utility, 16 colour names; reached from tracflux/tracgrad/tracerrt 7; parents under codes 11;/12;/13;; ALSO reached from ppnlse.f:449 directly under code 39; (BACKGROUND COLOR); see 6.1 audit |
| 71;84; | ENTRSORT I/O units (compressed — 6 leaves) | ENTRSORT I/O units (compressed — 6 leaves) | low | no | no | — | — | entrsort sub-menu — shared utility, 6 leaves; reached from util/suives.f via debunlse 73;; see 6.1 audit |
| 71;85; | FICHIER MS file (compressed — 3 leaves) | FICHIER MS file (compressed — 3 leaves) | low | no | no | — | — | fichier sub-menu — shared utility, 3 leaves; reached from util/ajfich.f via debunlse 72;; see 6.1 audit |
| 71;86; | LECTEUR input unit (compressed — 3 leaves) | LECTEUR input unit (compressed — 3 leaves) | low | no | no | — | — | lecteur sub-menu — shared utility, 3 leaves; reached from util/suives.f via debunlse 73;; see 6.1 audit |
| 71;87; | MANAGFEN window-state (compressed — 2 leaves) | MANAGFEN window-state (compressed — 2 leaves) | low | no | no | — | — | managfen sub-menu — shared utility, 2 leaves (window pixels + background color); for nlse reached only indirectly — codes 38; and 39; bypass managfen and call XVPXFE / couleurs prompt directly; listed for grep-set completeness; see 6.1 audit |
| 71;88; | MANAGMEF resource management (compressed — 4 leaves) | MANAGMEF resource management (compressed — 4 leaves) | low | no | no | — | — | managmef sub-menu — shared utility, 4 leaves; for nlse reached via util/suitms (debunlse 71;), util/suifms (debunlse 72;), util/suives (debunlse 73;), util/tuer (debunlse 74;) — these dispatch into managmef internally; see 6.1 audit |
| 71;89; | MODE_ES interactivity (compressed — 3 leaves) | MODE_ES interactivity (compressed — 3 leaves) | low | no | no | — | — | mode_es sub-menu — shared utility, 3 leaves; reached from util/suives.f via debunlse 73;; see 6.1 audit |
| 71;90; | OPT_LIGN line drawing (compressed — 33 leaves) | OPT_LIGN line drawing (compressed — 33 leaves) | low | no | no | — | — | opt_lign sub-menu — shared utility, 33 leaves; reached via TRMAIL path under debunlse 19; (DRAW PLSVO MESHES); see 6.1 audit |
| 71;91; | OPT_SURF surface drawing (compressed — 40 leaves) | OPT_SURF surface drawing (compressed — 40 leaves) | low | no | no | — | — | opt_surf sub-menu — shared utility, 40 leaves; reached via TRMAIL path under debunlse 19;; see 6.1 audit |
| 71;92; | SECTVOPL volume sections (compressed — 4 leaves) | SECTVOPL volume sections (compressed — 4 leaves) | low | no | no | — | — | sectvopl sub-menu — shared utility, 4 leaves; reached via 3D ISO/section paths under debunlse codes 11;/12;/13; (TRTHER cascade); see 6.1 audit |
| 71;93; | SUIVI_MS MS file mgmt (compressed — 1 leaf) | SUIVI_MS MS file mgmt (compressed — 1 leaf) | low | no | no | — | — | suivi_ms sub-menu — shared utility, 1 leaf; reached from util/suifms via debunlse 72;; see 6.1 audit |
| 71;94; | SUIVITMS TMS tools (compressed — 5 leaves) | SUIVITMS TMS tools (compressed — 5 leaves) | low | no | no | — | — | suivitms sub-menu — shared utility, 5 leaves; reached from util/suitms via debunlse 71;; see 6.1 audit |
| 71;95; | TRACAXES axes drawing (compressed — 4 leaves) | TRACAXES axes drawing (compressed — 4 leaves) | low | no | no | — | — | tracaxes sub-menu — shared utility, 4 leaves; reached via TRMAIL path under debunlse 19;; see 6.1 audit |
| 71;96; | TUER PLSVO deletion (compressed — 4 leaves) | TUER PLSVO deletion (compressed — 4 leaves) | low | no | no | — | — | tuer sub-menu — shared utility, 4 leaves; reached from util/tuer via debunlse 74; (TUER des TMS de PLSVO); see 6.1 audit |
| 71;97; | TYP_OBJT object types (compressed — 5 leaves) | TYP_OBJT object types (compressed — 5 leaves) | low | no | no | — | — | typ_objt sub-menu — shared utility, 5 leaves (POINT/LIGNE/SURFACE/VOLUME/OBJET); reached via TRMAIL path under debunlse 19;; see 6.1 audit |
| 71;98; | TYPTRAIT line type (compressed — 3 leaves) | TYPTRAIT line type (compressed — 3 leaves) | low | no | no | — | — | typtrait sub-menu — shared utility, 3 leaves; reached from tracflux/tracgrad/tracerrt 6; parents under codes 11;/12;/13;; see 6.1 audit |
| 71;99; | VUESPLAN special-plane views (compressed — 6 leaves) | VUESPLAN special-plane views (compressed — 6 leaves) | low | no | no | — | — | vuesplan sub-menu — shared utility, 6 leaves; reached via TRMAIL path under debunlse 19;; see 6.1 audit |
| 20;1; | ZEROS precision (compressed — 3 leaves) | ZEROS precision (compressed — 3 leaves) | low | no | no | — | — | zeros sub-menu — shared utility, 3 leaves; reached from util/zeros.f via debunlse 20; (PRECISION pour INVERSER A x = b dispatcher ZEROGC); see 6.1 audit |

<!-- End of audit table — validator will count rows above this line. -->

## Summary Statistics

- **Total rows:** 80
- **By frequency:**
  - `high`: 9 — top-level: `1;`, `2;`, `5;`, `11;`, `99;` (5 domain-promoted canonical NLSE workflow codes) plus sub-leaves: `2;3;` (cl_nlse Dirichlet BC) plus `2;91;1;`, `2;91;2;`, `2;91;6;` (coefnlse mandatory material parameters Rho / Alfa / Beta) = 4 HIGH leaves; total **9 HIGH rows**
  - `med`: 11 — top-level: `6;`, `12;`, `13;`, `14;`, `15;`, `16;`, `19;` (7 top-level MED) plus sub-leaves: `2;2;` (cl_nlse FGamma) plus `2;8;` (cl_nlse WAVE_0) plus `2;91;4;` (coefnlse FOmega) plus `2;91;8;` (coefnlse U(t0,X)) = 4 MED leaves; total **11 MED rows**
  - `low`: 60 — remainder: 8 top-level LOW (codes 20;, 38;, 39;, 71;, 72;, 73;, 74;, 97;) + 31 ther-proxy compressed + 20 shared-util compressed + 1 zeros compressed (under code 20;) = **60 LOW rows**
- **By qaction:** `yes` count == high+med == **20** (Plan 02's `registerNlseActions_stub_` QAction set; smaller than 6.4 ther's 37 because the visualization tree under codes 11;/12;/13; is fully delegated to the 6.4 ther audit via compressed proxy rows; nlse-primary expansion is just cl_nlse 3 + coefnlse 6 = 9 sub-leaves)
- **By toolbar:** exactly 5 `yes` — `2;`, `5;`, `11;`, `19;`, `99;`
- **Help-allowlist for Plan 03 testHelpNoQueue:** `{97;}` (NOT `{98;}` like ther — see Help-allowlist section above for the explicit hand-off; matches flui's `{97;}`)

## Top-5 Toolbar (Auto-frozen — mirrors 6.3/6.4 fully-autonomous path)

1. `2;` — NLSE INPUT DATA / ENTREE des DONNEES de l'Equation Non Lineaire de SCHRODINGER — icon `nlse-input.svg`
2. `5;` — IMPLICIT NLSE solver / IMPLICIT idU(t)/dt -Alfa Laplacien U(t) +N(t,abs(U))U(t)=F(t,X) — icon `solve-nlse-implicit.svg`
3. `11;` — DRAW NLSE WAVE MODULE / DESSIN du MODULE de l'ONDE NLSE — icon `draw-nlse-modulus.svg`
4. `19;` — DRAW PLSVO MESHES / DESSIN du MAILLAGE des PLSVO — icon `mesh-draw.svg` (REUSE from 6.1 mail/)
5. `99;` — SAVE DATA and QUIT / SAUVEGARDE des donnees FIN TRAITEMENT — icon `SP_DialogCloseButton`

Rationale: `2;` mandatory data-input gate (every NLSE simulation prompts
DFNLSE → cl_nlse + coefnlse cascade before solver dispatch); `5;`
canonical first solver run (IMPLICIT scheme is the default — Gross-Pitaevskii
i-time `6;` is specialty); `11;` primary result visual (the abs(U(t,X))
modulus is the universal NLSE result-visualization entry — Re/Im parts
`12;`/`13;` are specialty); `19;` mesh inspection (borrow from 6.1 mail
toolbar convention — preserves 6.1/6.4 reuse of `mesh-draw.svg`); `99;`
shared 6.1/6.2/6.3/6.4 close-session convention. The user may at a
polish-pass:
- Swap `19;` for `1;` if "Object first" workflow is more critical (testa
  has no nlse fixtures so the choice is domain-driven)
- Swap `11;` for `15;` if convergence diagnostic dominates result
  visualization (Max abs(U(Node))(Time) is the time-history diagnostic)

Any change must keep the count at exactly 5.

## Nlse-unique SVG icon set (Plan 02 ships)

Custom .svg filenames introduced by this audit that Plan 02 must ship under
`xvue/qt/resources/icons/nlse/`:

- `nlse-object.svg` — OBJECT NAME (code 1;) — abstract wave glyph
- `nlse-input.svg` — NLSE INPUT DATA (code 2;) — wave-equation glyph
- `solve-nlse-implicit.svg` — IMPLICIT NLSE solver (code 5;) — wave + integral sign
- `draw-nlse-modulus.svg` — DRAW abs(U(t,X)) (code 11;) — wave packet glyph
- `draw-nlse-component.svg` — DRAW Re U / Im U (codes 12; AND 13; share — single
  icon between the two leaves to keep new-SVG count modest, mirroring 6.4
  consolidation pattern)

Total: **5 new nlse-specific SVGs.** User may add `solve-nlse-itime.svg`
for code 6; in a polish-pass if Gross-Pitaevskii branch becomes a
documented use case (currently code 6; has `icon_source=—` since it is
not on the toolbar; menu-only entry is sufficient).

## Shared SVG reuses (from 6.1 mail/)

- `mesh-draw.svg` — reused from `xvue/qt/resources/icons/mail/` for nlse
  code `19;` (DRAWING of PLSVO MESHES); the qrc prefix `/xvue/qt/icons`
  resolves the mail path regardless of module context. No file copy
  needed. Same reuse pattern as 6.2 elas code `10;`, 6.3 flui code
  `19;`, and 6.4 ther code `10;`.

## Cross-References

- 06.1 `LEXICON-AUDIT-mail.md` — 9-column schema template (this audit
  mirrors the structure verbatim) + leaf-level enumeration of shared
  util sub-menus (compressed here)
- 06.2 `LEXICON-AUDIT-elas.md` — second-iteration example; demonstrated
  the synthetic high-number-offset prefix pattern (3;81;, 3;91;, 4;50;)
  replicated here as 2;91; (coefnlse), 5;81;..5;87; (ther-proxy
  solver/wave family), 6;81;..6;84; (ther-proxy eigen/cl_ther dead-code
  family), 71;81;..71;99; (shared-util compressed bucket)
- 06.3 `LEXICON-AUDIT-flui.md` — third-iteration example; demonstrated
  module-aware validator + Help-allowlist `{97;}` shape (matches
  nlse's `{97;}`)
- 06.4 `LEXICON-AUDIT-ther.md` — fourth-iteration example AND CANONICAL
  OWNER of the 31 ther-proxy compressed rows here; Plan 02 of 6.5 must
  read 6.4's audit to get leaf-level detail for every row tagged
  "compressed — N leaves" in the ther-proxy bucket (rows 11;7; through
  6;84;1; in this audit)
- 06.1 `06.1-01-SUMMARY.md` — de-duplication rule (module-specific
  fully expanded; shared util compressed to 1 row each)
- 06.4 `06.4-01-SUMMARY.md` — domain-review-only frequency bucketing for
  sparse testa corpora (replicated here for the EVEN-SPARSER nlse
  fixture set: ZERO .nlse files exist, 5 testa/ projects with nlse-named
  mesh-only data)
- 06.3 `06.3-03-SUMMARY.md` §"testHelpNoQueue tightened" — Auto-fix
  Rule 1 lesson: per-module Help-allowlist must be drawn from
  LEXICON-AUDIT (NOT inherited from previous-module template)
- 06.0-UI-SPEC §Per-Module Conformance Contract — nlse module-title
  `Nonlinear / Non Linéaire` per ROADMAP Phase 6.5 line 220 (overrides
  the 06.0-UI-SPEC table at line 245 which lists `Solve / Calcul`; see
  Menu taxonomy section above for the conflict resolution note)
- ROADMAP Phase 6.5 line 220 — `{File, Nonlinear, View, Help}` taxonomy
- `tools/validate_audit_md.py` — 9-rule schema validator; module-aware
  ICONS resolver since 6.2 Plan 02 (no validator changes needed in 6.5)
- Plan 02 consumes: rows with `qaction=yes` as the
  `registerNlseActions_stub_` QAction set (20 rows); rows with
  `toolbar=yes` as the toolbar append set (5 rows); `icon_source` ending
  `.svg` (filename not already in `icons/mail/` or `icons/elas/` or
  `icons/flui/` or `icons/ther/`) as the `xvue_icons.qrc` append set
  for `icons/nlse/` (5 new SVGs to ship)
- Plan 03 consumes: Help-allowlist `{97;}` as the testHelpNoQueue
  QSet<QString> (NOT `{98;}` — that's ther's; matches flui's `{97;}`;
  differs per Auto-fix Rule 1 lesson from 6.3 Plan 03)
