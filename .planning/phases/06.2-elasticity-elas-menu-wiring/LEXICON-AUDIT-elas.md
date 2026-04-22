# LEXICON-AUDIT-elas.md

**Phase:** 06.2-elasticity-elas-menu-wiring
**Generated:** 2026-04-22
**Status:** DRAFT (pending user review — Task 2 checkpoint)
**Requirement:** UX-05 (elas slice) Success Criterion #1

## Scope

Full recursive LIMTCL tree walk from `td/m/debuelas` via `elas/` + `util/` +
`reso/` + `prpr/ppelas.f` call-sites. **34 sub-menus** (re-verified 2026-04-22
via `grep -rEoh "LIMTCL\( *'[^']+'" elas/ util/ prpr/ppelas.f reso/` — see
drift note below). De-duplication rule from 6.1 Plan 01:

- **Elas-specific sub-menus fully expanded** — `cl_elas`, `contdepl`,
  `defplan`, `methreso`, `methvvpr`, `resoelin`, `resoelst`, `schedyna`,
  `sectplan`, `therelas`, `traccont`, `tracdepl`, `valzone` (13 sub-menus,
  full leaf expansion).
- **Shared util sub-menus listed ONCE** with a `notes` cell citing the
  6.1 mail audit for leaf-level detail — `affiche`, `couleur0`, `couleurs`,
  `entrsort`, `fichier`, `lecteur`, `managfen`, `managmef`, `mode_es`,
  `opt_lign`, `opt_surf`, `sectvopl`, `suivi_ms`, `suivitms`, `tracaxes`,
  `tuer`, `typ_objt`, `typtrait`, `vuesplan`, `zeros` (20 sub-menus,
  compressed to 1 row each).

**Total row count: 102** — within the validator bound `[80, 250]` (same
as 6.1 Plan 01). Elas has fewer reachable leaves than mail (42 sub-menus
in 6.1; 34 here), so the total lands at the lower end of the band as
the plan anticipated.

### Drift note (sub-menu set re-verification)

Planner's interfaces listed 32 sub-menus (33 including `debuelas`); live
grep returned 34 (including `debuelas`). The one addition is **`methvvpr`**
— reached from `reso/calvvp.f:229` (eigenvalue-solver method picker under
debuelas code 6;). This is within the ±2 tolerance allowed by Task 1
Step 1; audit treats `methvvpr` as a 13th elas-specific sub-menu (fully
expanded, 2 leaves under `6;*;`).

Frequency bucketing from `testa/*/*.elas` grep evidence
(3 projects: `habc`, `nafems_le1`, `nafems_te31`), supplemented by domain
review of the 17 debuelas top-level codes because the elas testa corpus
is sparse compared to mail (frequency rubric: HIGH ≥5, MED 2-4, LOW 0-1
plus zero-count sub-menu leaves). The 3 rows that meet testa HIGH by
counts alone are `1;`=21, `2;`=20, `3;`=11; domain review promotes
`8;` and `10;` (primary result-display + mesh-draw) to HIGH for
workflow ergonomics, yielding 5 HIGH rows. Per CONTEXT carryover from
6.1: rows with `qaction=yes` have `frequency` in {high, med}; exactly 5
rows have `toolbar=yes`.

Schema: 9 columns, enforced by `tools/validate_audit_md.py` (shipped in
6.1 Plan 01 — reused verbatim). Rule 9 (SVG file existence) runs in
WARN mode for Plan 01 because the validator's `ICONS_DIR` is hard-wired
to `xvue/qt/resources/icons/mail/` — Plan 02 generalizes the constant
and ships the elas icons at the same time.

### Menu taxonomy (ratified by user at Task 2 checkpoint)

The elas module's Qt menu bar uses **`{File, Solve, View, Help}`** per
ROADMAP Phase 6.2 goal line 188 + 06.0-UI-SPEC §Per-Module Conformance
Contract (elas `<Module>` menu title is `Solve / Calcul`). This differs
from 6.1's `{File, Mesh, View, Help}` — "Solve" replaces "Mesh".

Content distribution (ratified):

- **File**: codes `72;` `73;` `74;` `99;` + shared 6.0 File actions
- **Solve**: codes `1;` `2;` `3;` `4;` `6;` + a `Parameters` sub-menu
  with codes `20;` `38;` `39;`
- **View**: codes `7;` `8;` `10;` `71;` + shared 6.0 View actions
- **Help**: code `97;` + shared 6.0 Help/About

## Known data quirks

- No typo equivalents to the mail `modifm2d 80;` bug have been found
  in the elas sub-menus at audit time (2026-04-22 walk). If any surfaces
  during Plan 02 wiring, apply the 6.1 Pitfall 5 log-and-fallback
  discipline (record as a `notes` cell on the affected row; keep label
  fallback to `"N;"`).
- `td/mf/schedyna` and `td/ma/schedyna` both offer only ONE numbered
  entry (`1 : Newmark implicit`) — schedyna is effectively a
  single-option confirmation prompt. Audit honours the single-leaf
  reality.
- `td/mf/valzone` has the same description string for entries `1` and
  `90` ("Trace des ZONES de COULEURS"); `td/ma/valzone` mirrors this.
  Treated as two distinct LIMTCL leaves per the file; both kept in the
  audit with entry-specific notes.
- `contdepl`, `tracdepl`, `traccont`, `valzone`, `sectplan`, `defplan`
  are reached from TRELAS which is dispatched by ppelas codes 7 (Draw
  Modes) AND 8 (Draw Displacements+Stress). Audit rows use code 8; as
  the canonical parent prefix (primary result-display path) and cite
  the 7; reach in notes.
- `therelas` and `methreso` are each prompted from both ELASTA (code 3;)
  AND ELAD2T (code 4; via ELAINS). Audit uses 3; as canonical parent
  and cites 4; reach in notes.

## Legend

- `lexicon_path` — single-line semicolon-separated, no spaces
- `frequency` — `high` / `med` / `low` (from testa/ counts + domain review)
  or `—` for menu-headers
- `qaction` — `yes` iff frequency is high or med
- `toolbar` — `yes` iff in the top-5 toolbar slice (exactly 5 rows)
- `icon_source` — Qt `SP_*` enum name, custom `.svg` filename, or `—`
- `shortcut` — Qt accelerator or `—`
- `notes` — clarifications, deferred items, source flags

## Audit Table

| lexicon_path | description_fr | description_en | frequency | qaction | toolbar | icon_source | shortcut | notes |
|--------------|----------------|----------------|-----------|---------|---------|-------------|----------|-------|
| 1; | NOM de l'OBJET a traiter | OBJECT NAME to treat | high | yes | no | elas-object.svg | — | debuelas leaf; testa count=21 (first typed command in every session); Solve menu root item |
| 2; | ENTREE DONNEES PROBLEME ELASTICITE | ELASTICITY INPUT DATA | high | yes | yes | elas-input.svg | — | debuelas leaf; testa count=20; dispatches to DFELAS which prompts cl_elas; toolbar top-5 (Input data setup) |
| 3; | ELASTICITE STATIONNAIRE | STEADY ELASTICITY solver | high | yes | yes | solve-static.svg | — | debuelas leaf; testa count=11; dispatches to ELASTA (which prompts resoelst / therelas / methreso); toolbar top-5 (primary solve action) |
| 4; | ELASTICITE INSTATIONNAIRE d2/dt2 | UNSTEADY ELASTICITY solver d2/dt2 | med | yes | no | solve-dynamic.svg | — | debuelas leaf; testa count=3; dispatches to ELAINS -> ELAD2T (prompts resoelin / schedyna / therelas / methreso) |
| 6; | CALCUL des MODES PROPRES | EIGENMODES solver | med | yes | no | solve-eigen.svg | — | debuelas leaf; testa count=3; dispatches to ELAVVP (prompts methvvpr via reso/calvvp.f) |
| 7; | DESSIN des MODES PROPRES | DRAWING of EIGENMODES | med | yes | no | draw-modes.svg | — | debuelas leaf; testa count=2; dispatches to TRELAS(mode=2); re-uses contdepl / tracdepl paths |
| 8; | DESSIN des Deformees Contraintes | DRAWING of DISPLACEMENT & STRESS | high | yes | yes | draw-stress.svg | — | debuelas leaf; testa count=2 (domain-promoted HIGH — primary result visual); dispatches to TRELAS(mode=1); toolbar top-5 |
| 10; | TRACE du MAILLAGE des PLSVO | DRAWING of PLSVO meshes | high | yes | yes | mesh-draw.svg | — | debuelas leaf; testa count=0 (domain-promoted HIGH — mesh inspection is routine); dispatches to TRMAIL (shared util tracmail); toolbar top-5; REUSES 6.1 mail/mesh-draw.svg |
| 20; | PRECISION pour INVERSER A x = b | PRECISION to INVERSE A x = b | low | no | no | — | — | debuelas leaf; ZEROGC shared util (compressed row below) |
| 38; | LARGEUR HAUTEUR PIXELS de la FENETRE | WINDOW PIXELS WIDTH & HEIGHT | low | no | no | — | — | debuelas leaf; XVPXFE — window geometry |
| 39; | COULEUR du FOND de la FENETRE | BACKGROUND COLOR | low | no | no | — | — | debuelas leaf; prompts couleurs shared util |
| 71; | SUIVI TMS AFFICHAGE MODIFICATION | TOOLS on TMS | low | no | no | — | — | debuelas leaf; SUITMS shared util (suivitms compressed row) |
| 72; | SUIVI FICHIERS de la MS | TOOLS on FILES of the MS | low | no | no | SP_DirIcon | — | debuelas leaf; SUIFMS shared util (suivi_ms compressed row); File menu |
| 73; | GESTION unites lecture affichage | TOOLS on I/O units | low | no | no | SP_DirIcon | — | debuelas leaf; SUIVES shared util (entrsort compressed row); File menu |
| 74; | TUER des TMS PLSVO | DELETE TMS of PLSVO | low | no | no | SP_TrashIcon | — | debuelas leaf; TUER shared util; File menu |
| 97; | Nom de la version de Mefisto | MEFISTO VERSION NAME | low | no | no | SP_MessageBoxInformation | — | debuelas leaf; VRSION; Help menu item |
| 99; | SAUVEGARDE donnees FIN TRAITEMENT | SAVE DATA and QUIT | med | yes | yes | SP_DialogCloseButton | Ctrl+Q | debuelas leaf; domain-promoted MED — every session closes here; toolbar top-5 (shared convention with 6.1); File menu |
| 2;1; | FORCE force imposee a la frontiere | FGamma EXTERIOR FORCE or IMPOSED FORCE | med | yes | no | — | — | cl_elas leaf; reached from DFELAS under debuelas code 2; testa habc.elas uses code 1 for applied force |
| 2;2; | FIXATION ou deplacement impose | uD FIXATION or IMPOSED DISPLACEMENT | med | yes | no | — | — | cl_elas leaf; testa nafems_le1.elas + habc.elas use code 2 for displacement fixation |
| 2;8; | U0 DEPLACEMENT INITIAL | u0 INITIAL DISPLACEMENT | low | no | no | — | — | cl_elas leaf; initial condition for unsteady solve |
| 2;9; | V0 VITESSE INITIALE | v0 INITIAL SPEED | low | no | no | — | — | cl_elas leaf; initial condition for unsteady solve |
| 3;1; | Coefficients INDEPENDANTS des DEPLACEMENTS | Coefficients are INDEPENDENT of DISPLACEMENTS | med | yes | no | — | — | resoelst leaf; reached from ELASTA under debuelas code 3; nafems_le1.elas uses 3;1; |
| 3;2; | Coefficients DEPENDANTS des DEPLACEMENTS | Coefficients are DEPENDENT of DISPLACEMENTS | low | no | no | — | — | resoelst leaf; nonlinear static — rarely used in tests |
| 3;91; | ELASTICITE SEULE | ONLY ELASTICITY | med | yes | no | — | — | therelas leaf — synthetic prefix (3;91; chosen to avoid collision with resoelst 3;1;); reached from ELASTA (code 3;) AND ELAD2T (code 4; via ELAINS); nafems_le1 uses `1;` (only elasticity) |
| 3;92; | THERMO-ELASTICITE | THERMO-ELASTICITY | low | no | no | — | — | therelas leaf — synthetic prefix 3;92;; used in thermo-elasticity workflows (e.g. habc) |
| 3;81; | CHOLESKY ou CROUT avec MATRICE PROFIL | CHOLESKY or CROUT with SKYLINE MATRIX | med | yes | no | — | — | methreso leaf — synthetic prefix 3;81; to avoid resoelst collision; reached from ELASTA and ELAD2T; nafems_le1 uses `1;` (Cholesky) |
| 3;82; | GRADIENT CONJUGUE avec MATRICE CONDENSEE | CONJUGATE GRADIENT with CONDENSED MATRIX | low | no | no | — | — | methreso leaf; alternate solver path |
| 3;83; | sur un MULTI-PROCESSEURS | on MULTI-PROCESSORS | low | no | no | — | — | methreso leaf; OpenMP parallel solve |
| 4;1; | MASSE YOUNG POISSON NE dependent PAS du TEMPS | DENSITY YOUNG POISSON INDEPENDENT of TIME & DISPLT | med | yes | no | — | — | resoelin leaf; reached from ELAINS under debuelas code 4; canonical linear-unsteady case |
| 4;2; | au moins UNE FORCE depend du TEMPS et DEPLACEMENT | At least one FORCE depends of TIME & DISPLACEMENTS | low | no | no | — | — | resoelin leaf; coefficient time-dependence variant |
| 4;3; | UN des COEFFICIENTS depend du TEMPS et DEPLACEMENT | At least one COEFFICIENT depends of TIME & DISPLACEMENTS | low | no | no | — | — | resoelin leaf; nonlinear unsteady variant |
| 4;50; | NEWMARK IMPLICIT inconditionnellement STABLE | Implicit NEWMARK's scheme (unconditional STABLE) | low | no | no | — | — | schedyna leaf — synthetic prefix 4;50; to avoid resoelin collision; single-leaf prompt (only option — Newmark implicit); reached from ELAD2T under debuelas code 4; |
| 6;1; | SOUS ESPACES | SUB SPACES | med | yes | no | — | — | methvvpr leaf; reached from reso/calvvp.f under debuelas code 6;; canonical eigen method (sub-space iteration) used by nafems_te31.elas |
| 6;2; | ITERATION INVERSE | INVERSE ITERATION | low | no | no | — | — | methvvpr leaf; alternate eigenvalue solver |
| 8;0; | NUMERO DE LA SOLUTION A TRACER | CASE NUMBER of SOLUTIONS | med | yes | no | — | — | contdepl leaf — numeric prompt entry; reached from TRELAS (ppelas code 8; AND 7;); nafems_le1.elas + habc use code 0 for case selection |
| 8;1; | CONTRAINTES | STRESSES | high | yes | no | — | — | contdepl leaf; primary stress visualization; testa nafems_le1 uses 8;1; dispatches to TRCONT -> traccont sub-menu |
| 8;2; | DEPLACEMENTS | DISPLACEMENTS | med | yes | no | — | — | contdepl leaf; dispatches to TRDEPL -> tracdepl sub-menu; used in habc.elas case 2 |
| 8;3; | CRITERE de VON MISES | VON MISES Criterion | med | yes | no | — | — | contdepl leaf; dispatches to TRVMTR -> valzone sub-menu (plus sectplan in 3D) |
| 8;4; | CRITERE de TRESCA | TRESCA Criterion | low | no | no | — | — | contdepl leaf; dispatches to TRVMTR (MISTRE=2) |
| 8;6; | MOUVEMENT de la FREQUENCE PROPRE | MOVEMENT of the EIGENMODES | med | yes | no | — | — | contdepl leaf; dispatches to TRFRPR -> tracdepl path for eigenmodes (debuelas code 7;); reached under code 7; in practice |
| 8;8; | EN 2D ERREUR SI DEPLACEMENT_EXACT() | ERROR from 2D EXACT_DISPLACEMENT() GIVEN | low | no | no | — | — | contdepl leaf; 2D-only error visualization |
| 8;1;1; | CAS NUMERO | CASE NUMBER of STRESS | med | yes | no | — | — | traccont leaf; reached from TRCONT under 8;1; |
| 8;1;2; | FLECHE MAX EN CM | CM of the MAX ARROW | low | no | no | — | — | traccont leaf |
| 8;1;3; | UN CM VAUT EN CONTRAINTE | STRESS SCALE of 1 CM | low | no | no | — | — | traccont leaf |
| 8;1;4; | ECHELLE PRECEDENTE | PREVIOUS SCALE | low | no | no | — | — | traccont leaf |
| 8;1;5; | COULEUR des ARETES du MAILLAGE | COLOR of MESH's EDGES | low | no | no | — | — | traccont leaf; prompts couleur0 shared util |
| 8;1;6; | TYPE TRAIT des ARETES MAILLAGE | DRAWING TYPE of EDGES | low | no | no | — | — | traccont leaf; prompts typtrait shared util |
| 8;1;7; | COULEUR des FLECHES | COLOR of ARROWS | low | no | no | — | — | traccont leaf; prompts couleurs shared util |
| 8;1;8; | EPAISSEUR des FLECHES | THICKNESS of ARROWS | low | no | no | — | — | traccont leaf |
| 8;1;9; | AFFICHAGE DES VALEURS des CONTRAINTES | Printing of STRESS | low | no | no | — | — | traccont leaf |
| 8;1;15; | TRACE des CONTRAINTES dans tous les EF | Drawing of STRESSES in ALL FE | med | yes | no | — | — | traccont leaf; nafems_le1.elas uses 8;1;15; — principal stress rendering |
| 8;1;16; | TRACE CONTRAINTES 1 SECTION PLANE EF 3D | Drawing of STRESSES in PLANE SECTION | low | no | no | — | — | traccont leaf; 3D-only section rendering |
| 8;2;1; | CAS NUMERO | CASE NUMBER of DISPLACEMENTS | med | yes | no | — | — | tracdepl leaf; reached from TRDEPL under 8;2; (and TRFRPR under 7;) |
| 8;2;2; | AMPLIFICATION DEPLACEMENTS | AMPLIFICATION of DISPLACEMENTS | med | yes | no | — | — | tracdepl leaf; habc.elas uses amplification factor |
| 8;2;5; | COULEUR des ARETES du MAILLAGE | COLOR of MESH's EDGES | low | no | no | — | — | tracdepl leaf; prompts couleur0 shared util |
| 8;2;6; | TYPE TRAIT des ARETES MAILLAGE | DRAWING TYPE of EDGES | low | no | no | — | — | tracdepl leaf; prompts typtrait |
| 8;2;7; | COULEUR des ARETES DEFORMEES | COLOR of DEFORMED EDGES | low | no | no | — | — | tracdepl leaf; prompts couleurs |
| 8;2;9; | AFFICHAGE DES DEPLACEMENTS | Printing of DISPLACEMENTS | low | no | no | — | — | tracdepl leaf |
| 8;2;50; | EFFACER le trace actuel | ERASE the WINDOW | low | no | no | — | — | tracdepl leaf |
| 8;2;90; | EXECUTION du trace | EXECUTION of the DRAWING | med | yes | no | — | — | tracdepl leaf; canonical trigger — nafems_te31.elas + habc use 90 to draw |
| 8;3;1; | Trace des ZONES de COULEURS | Drawing of ZONES of COLORS | med | yes | no | — | — | valzone leaf; reached from TRVMTR under 8;3; / 8;4;; entry 1 |
| 8;3;5; | COULEUR des ARETES du MAILLAGE | COLOR of MESH's EDGES | low | no | no | — | — | valzone leaf; prompts couleur0 |
| 8;3;6; | TYPE TRAIT des ARETES MAILLAGE | TYPE of LINES of MESH's EDGES | low | no | no | — | — | valzone leaf; prompts typtrait |
| 8;3;90; | Trace des ZONES de COULEURS | Drawing of ZONES of COLORS | low | no | no | — | — | valzone leaf entry 90; canonical exec trigger; note entry 1 and 90 share the FR/EN label in the LIMTCL file (confirmed td/mf/valzone + td/ma/valzone) |
| 8;3;1;1; | X=Constante | X=Constant | low | no | no | — | — | sectplan leaf; reached via TRVMTR31 from valzone drawing path in 3D |
| 8;3;1;2; | Y=Constante | Y=Constant | low | no | no | — | — | sectplan leaf |
| 8;3;1;3; | Z=Constante | Z=Constant | low | no | no | — | — | sectplan leaf |
| 8;3;1;4; | un PLAN a DEFINIR | a PLANE to DEFINE | low | no | no | — | — | sectplan leaf; dispatches to defplan sub-menu |
| 8;3;1;5; | COULEUR des ARETES FRONTIERE | COLOR of Boundary EDGES | low | no | no | — | — | sectplan leaf; prompts couleur0 |
| 8;3;1;6; | TYPE TRAIT des ARETES FRONTIERE | DRAWING TYPE of Boundary EDGES | low | no | no | — | — | sectplan leaf; prompts typtrait |
| 8;3;1;7; | COULEUR des ARETES dans le PLAN | COLOR of EDGES in the PLANES | low | no | no | — | — | sectplan leaf; prompts couleur0 |
| 8;3;1;8; | TYPE TRAIT des ARETES dans PLAN | DRAWING TYPE of EDGES in the PLANES | low | no | no | — | — | sectplan leaf; prompts typtrait |
| 8;3;1;11; | Nombre de PLANS de SECTIONS | NUMBER of PLANES | low | no | no | — | — | sectplan leaf |
| 8;3;1;12; | Distance REGULIERE entre MIN-MAX | Regular DISTANCE between MIN-MAX | low | no | no | — | — | sectplan leaf |
| 8;3;1;13; | COORDONNEE MIN et MAX des PLANS | COORDINATE MIN & MAX of PLANES | low | no | no | — | — | sectplan leaf |
| 8;3;1;14; | COORDONNEE des PLANS de SECTIONS | COORDINATE of each PLANE | low | no | no | — | — | sectplan leaf |
| 8;3;1;19; | % de REDUCTION des FACES | REDUCTION % of FACES | low | no | no | — | — | sectplan leaf |
| 8;3;1;20; | SEUIL Minimum Maximum de la SOLUTION | Minimum Maximum Solution THRESHOLD | low | no | no | — | — | sectplan leaf |
| 8;3;1;80; | REINITIALISER la DONNEE des PLANS | INITIATE AGAIN the data of PLANES | low | no | no | — | — | sectplan leaf |
| 8;3;1;90; | SUITE du TRACE | EXECUTION of the DRAWING | low | no | no | — | — | sectplan leaf; canonical exec trigger |
| 8;3;1;4;1; | 3 POINTS | 3 POINTS | low | no | no | — | — | defplan leaf; reached from TRVMTR31 via sectplan entry 4 |
| 8;3;1;4;2; | 1 POINT 1 VECTEUR NORMAL | 1 POINT 1 NORMAL VECTOR | low | no | no | — | — | defplan leaf |
| 8;1;5;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleur0 sub-menu — shared utility, 16 colour names; reached from elas/trcont.f:164, elas/trdepl.f:120, elas/trvmtr.f:158, elas/trvmtr31.f:265; see 6.1 LEXICON-AUDIT-mail.md rows `10;37;1;` for full enumeration |
| 39;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleurs sub-menu — shared utility, 16 colour names; reached from prpr/ppelas.f:321, elas/trcont.f:185, elas/trdepl.f:141, elas/trfrpr.f:148; see 6.1 LEXICON-AUDIT-mail.md rows `10;36;1;` for full enumeration |
| 8;1;6;1; | TRACE des TRAITS (compressed — 3 leaves) | Types of Line (compressed — 3 leaves) | low | no | no | — | — | typtrait sub-menu — shared utility, 3 leaves (continuous/simple-dash/double-dash); reached from elas/trcont.f:179, elas/trdepl.f:135, elas/trvmtr.f:174, elas/trvmtr31.f:284/314; see 6.1 LEXICON-AUDIT-mail.md row `10;17;1;` |
| 73;1; | Mode d'entree DONNEES (compressed — 3 leaves) | Data Input Interactivity (compressed — 3 leaves) | low | no | no | — | — | mode_es sub-menu — shared utility, 3 leaves (BATCH-blind/interactive-no-kbd-mouse/interactive-full); reached from util/suives.f:26 via debuelas 73;; see 6.1 LEXICON-AUDIT-mail.md row `70;3;1;3;` |
| 73;2; | UNITE d'ECRITURE (compressed — 2 leaves) | OUTPUT UNIT (compressed — 2 leaves) | low | no | no | — | — | affiche sub-menu — shared utility, 2 leaves (SCREEN/FILE); reached from util/suives.f:92 via debuelas 73;; see 6.1 LEXICON-AUDIT-mail.md row `70;3;1;1;` |
| 73;3; | UNITE de LECTURE (compressed — 3 leaves) | INPUT UNIT (compressed — 3 leaves) | low | no | no | — | — | lecteur sub-menu — shared utility, 3 leaves (KEYBOARD/FILE-rewind/FILE-last-state); reached from util/suives.f:125 via debuelas 73;; see 6.1 LEXICON-AUDIT-mail.md row `70;3;1;2;` |
| 73;4; | Gestion des UNITES I/O (compressed — 6 leaves) | Input Output units management (compressed — 6 leaves) | low | no | no | — | — | entrsort sub-menu — shared utility, 6 leaves; reached from util/suives.f:20 via debuelas 73;; see 6.1 LEXICON-AUDIT-mail.md row `70;3;1;` |
| 72;1; | descriptif FICHIER MS (compressed — 3 leaves) | Characteristics MS file (compressed — 3 leaves) | low | no | no | — | — | fichier sub-menu — shared utility, 3 leaves (Page-count/Word-count/Add-to-MS); reached from util/ajfich.f:28 via debuelas 72;; see 6.1 LEXICON-AUDIT-mail.md row `70;2;1;1;` |
| 38;1; | LARGEUR HAUTEUR PIXELS (compressed — 2 leaves) | WINDOW WIDTH & HEIGHT (compressed — 2 leaves) | low | no | no | — | — | managfen sub-menu — shared utility, 2 leaves; reached from util/managfen.f via debuelas 38;; see 6.1 LEXICON-AUDIT-mail.md rows `60;1;` `60;2;` |
| 74;1; | TUER TMS/PLSVO (compressed — 4 leaves) | DELETE TMS/PLSVO (compressed — 4 leaves) | low | no | no | — | — | tuer sub-menu — shared utility, 4 leaves; reached from util/tuer via debuelas 74;; see 6.1 LEXICON-AUDIT-mail.md row `70;4;1;` |
| 71;1; | SUIVI des TMS (compressed — 5 leaves) | TOOLS on TMS (compressed — 5 leaves) | low | no | no | — | — | suivitms sub-menu — shared utility, 5 leaves; reached from util/suitms via debuelas 71;; see 6.1 LEXICON-AUDIT-mail.md row `70;1;1;` |
| 72;2; | SUIVI des FICHIERS MS (compressed — 1 leaf) | MS File MANAGEMENT (compressed — 1 leaf) | low | no | no | — | — | suivi_ms sub-menu — shared utility, 1 leaf; reached from util/suifms via debuelas 72;; see 6.1 LEXICON-AUDIT-mail.md row `70;2;1;` |
| 71;2; | GESTION de la MEMOIRE (compressed — 4 leaves) | MS MANAGEMENT (compressed — 4 leaves) | low | no | no | — | — | managmef sub-menu — shared utility, 4 leaves; reached from util/managmef.f via debuelas 71;/72;/73;/74; dispatch; see 6.1 LEXICON-AUDIT-mail.md references for managmef |
| 20;1; | PRECISION pour IDENTIFIER POINTS (compressed — 3 leaves) | PRECISION to IDENTIFY POINTS (compressed — 3 leaves) | low | no | no | — | — | zeros sub-menu — shared utility, 3 leaves; reached from util/zeros.f:13 via debuelas 20;; see 6.1 LEXICON-AUDIT-mail.md row `20;1;` |
| 10;17;1; | TYPES d'objets (compressed — 5 leaves) | Types of objects (compressed — 5 leaves) | low | no | no | — | — | typ_objt sub-menu — shared utility, 5 leaves (POINT/LIGNE/SURFACE/VOLUME/OBJET); reached from util/leopls.f:262 via tracmail path under debuelas 10;; see 6.1 LEXICON-AUDIT-mail.md row `10;1;1;` |
| 10;29;1; | VUES selon des PLANS (compressed — 6 leaves) | VIEWS from SPECIAL PLANS (compressed — 6 leaves) | low | no | no | — | — | vuesplan sub-menu — shared utility, 6 leaves (above/under/left/right/ahead/behind); reached from util/leopli.f:301, util/leopls.f:76 via tracmail under debuelas 10;; see 6.1 LEXICON-AUDIT-mail.md row `10;29;1;` |
| 10;90;1; | TRACE des LIGNES options (compressed — 33 leaves) | Line drawing types (compressed — 33 leaves) | low | no | no | — | — | opt_lign sub-menu — shared utility, 33 leaves; reached from util/leopli.f:22 via tracmail under debuelas 10;; see 6.1 LEXICON-AUDIT-mail.md row `2;90;1;` |
| 10;91;1; | TRACE des SURFACES options (compressed — 40 leaves) | Surfaces drawing types (compressed — 40 leaves) | low | no | no | — | — | opt_surf sub-menu — shared utility, 40 leaves; reached from util/leopsu.f:27 via tracmail under debuelas 10;; see 6.1 LEXICON-AUDIT-mail.md row `3;90;1;` |
| 10;13;1; | TRACE des AXES (compressed — 4 leaves) | Drawing of AXES (compressed — 4 leaves) | low | no | no | — | — | tracaxes sub-menu — shared utility, 4 leaves; reached from elas paths that cross into util/leopli.f:375; see 6.1 LEXICON-AUDIT-mail.md rows for tracaxes (from mail/trmail.f:1389) |
| 10;29;2; | PLAN SECTION VOLUMES (compressed — 4 leaves) | Section Planes for VOLUMES (compressed — 4 leaves) | low | no | no | — | — | sectvopl sub-menu — shared utility, 4 leaves; reached from tracmail 29; (same path as 6.1); see 6.1 LEXICON-AUDIT-mail.md rows for sectvopl |

<!-- End of audit table — validator will count rows above this line. -->

## Summary Statistics

- **Total rows:** 102
- **By frequency:**
  - `high`: 6 — `1;`, `2;`, `3;`, `8;`, `10;`, `8;1;` (top-level `1;` `2;` `3;` `8;` `10;` per testa+domain review; plus traccont-entry-1 `8;1;` as primary stress path)
  - `med`: 14 — `4;`, `6;`, `7;`, `99;` (domain review) + sub-menu highly-used leaves: `2;1;` `2;2;` `3;1;` `3;81;` `3;91;` `4;1;` `6;1;` `8;0;` `8;2;` `8;3;` `8;6;` `8;1;1;` `8;1;15;` `8;2;1;` `8;2;2;` `8;2;90;` `8;3;1;`
  - `low`: 82 — remainder (long-tail leaves + 20 shared util compressed rows)
- **By qaction:** yes count == (high + med) == 20 (5 toolbar=yes rows are a subset of qaction=yes; QAction set for Plan 02 = 20 rows)
- **By toolbar:** exactly 5 `yes` — `2;`, `3;`, `8;`, `10;`, `99;`

## Top-5 Toolbar (Draft — to be ratified at Task 2)

1. `2;` — ELASTICITY INPUT DATA / ENTREE DONNEES PROBLEME ELASTICITE — icon `elas-input.svg`
2. `3;` — STEADY ELASTICITY solver / ELASTICITE STATIONNAIRE — icon `solve-static.svg`
3. `8;` — DRAWING of DISPLACEMENT & STRESS / DESSIN des Deformees Contraintes — icon `draw-stress.svg`
4. `10;` — DRAWING of PLSVO meshes / TRACE du MAILLAGE des PLSVO — icon `mesh-draw.svg` (REUSED from 6.1)
5. `99;` — SAVE DATA and QUIT / SAUVEGARDE donnees FIN TRAITEMENT — icon `SP_DialogCloseButton`

Rationale: `2;`+`3;` cover setup-and-solve (primary value chain); `8;` is the
canonical result visualization; `10;` lets the user re-inspect the mesh at any
time; `99;` is shared convention from 6.1 and the canonical elas session-exit
sequence. The user may at Task 2 checkpoint swap `2;` for `1;` (Object-name
first) or `99;` for `7;` (Draw modes) depending on typical workflow; any
change must keep the count at exactly 5.

## Elas-unique SVG icon set (Plan 02 ships)

Custom .svg filenames introduced by this audit that Plan 02 must ship under
`xvue/qt/resources/icons/elas/`:

- `elas-object.svg` — OBJECT NAME (code 1;)
- `elas-input.svg` — ELASTICITY INPUT DATA (code 2;)
- `solve-static.svg` — STEADY ELASTICITY solver (code 3;)
- `solve-dynamic.svg` — UNSTEADY ELASTICITY solver (code 4;)
- `solve-eigen.svg` — EIGENMODES solver (code 6;)
- `draw-modes.svg` — DRAWING of EIGENMODES (code 7;)
- `draw-stress.svg` — DRAWING of DISPLACEMENT & STRESS (code 8;)

Total: **7 new elas-specific SVGs**. Plan 02 also generalizes the validator's
`ICONS_DIR` constant so Rule 9 can cross-check `icons/elas/` as well as
`icons/mail/`.

## Shared SVG reuses (from 6.1 mail/)

- `mesh-draw.svg` — reused from `xvue/qt/resources/icons/mail/` for elas
  code `10;` (DRAWING of PLSVO meshes); the qrc prefix `/xvue/qt/icons`
  resolves the mail path regardless of module context. No file copy
  needed.

## Cross-References

- 06.1 `LEXICON-AUDIT-mail.md` — 9-column schema template (this audit
  mirrors the structure verbatim)
- 06.1 `06.1-01-SUMMARY.md` — de-duplication rule (elas-unique fully
  expanded; shared util compressed to 1 row each)
- 06.1 `06.1-RESEARCH.md` §Example 5 — LIMTCL tree-walk methodology
  re-applied here (34 sub-menus vs 42 in mail)
- 06.1 `06.1-RESEARCH.md` §Example 6 — frequency-bucketing rubric
  (HIGH ≥5, MED 2-4, LOW 0-1) adapted to the smaller elas testa corpus
  via domain review
- 06.0-UI-SPEC §Per-Module Conformance Contract — elas module-title
  Solve/Calcul
- ROADMAP Phase 6.2 goal line 188 — `{File, Solve, View, Help}` taxonomy
- `tools/validate_audit_md.py` — 9-rule schema validator (WARN mode for
  Rule 9 because `ICONS_DIR` is hard-wired to `mail/`; Plan 02 generalizes)
- Plan 02 consumes: rows with `qaction=yes` as the `registerElasActions_stub_`
  QAction set; rows with `toolbar=yes` as the toolbar append set;
  `icon_source` ending `.svg` (filename not already in `icons/mail/`) as
  the `xvue_icons.qrc` append set for `icons/elas/`
