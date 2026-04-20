# LEXICON-AUDIT-mail.md

**Phase:** 06.1-mesher-mail-menu-wiring
**Generated:** 2026-04-20
**Status:** DRAFT (pending user review — Task 3 checkpoint)
**Requirement:** UX-05 (mail slice) Success Criterion #1

## Scope

Full recursive LIMTCL tree walk from `td/m/debut` via mail/ + util/
call-sites. 42 sub-menus (RESEARCH §Example 5, re-verified 2026-04-20
at 42 matches — no drift). De-duplication rule:

- **Mail-specific sub-menus fully expanded** — every leaf gets its own
  row with the full parent-prefixed `lexicon_path` (e.g. `1;1;`,
  `5;2;3;`, `10;37;`). Applies to: `saisi_pt`, `saisipt2`, `saisipt3`,
  `modifm2d`, `modifm3d`, `modifobj`, `tracmail`, `tracpoin`, `traclign`,
  `tracsurf`, `tracvolu`, `tracobje`, `tracaxes`, `typ3face`, `typ_plsv`,
  `typ_sv`, `tgoupas`, `sectvopl`, `pasftrcl`, `fidap`, `importer`,
  `exporter`, `interpol`.
- **Shared util sub-menus listed ONCE** with a `notes` cell noting
  "shared utility; reached from many parents" — compression form
  `<parent>;N;*` where N is the parent-menu entry and `*` marks the
  compressed body. Applies to: `couleur0`, `couleurs`, `non_oui`,
  `zeros`, `affiche`, `entrsort`, `lecteur`, `managfen`, `managmef`
  (expanded — 4 leaves only, all mail-reachable), `mode_es`, `opt_lign`,
  `opt_surf`, `suivi_ms`, `suivitms`, `tuer`, `typ_objt`, `typtrait`,
  `vuesplan`, `fichier`.

**Total row count: 211** — within the validator bound `[80, 250]`. This
widens CONTEXT D-11's original `100-150` target to accommodate the
~230 reachable leaves documented in RESEARCH §Example 5 while the
de-duplication rule trims ~40 rows off the full expansion. The widened
bound is ratified at Task 3 user-review checkpoint.

Frequency bucketing from `testa/*/*.mesh` grep evidence (RESEARCH
§Example 6, re-verified 2026-04-20 — the 4 mesh files scanned produced
the same frequency ranking). Per CONTEXT D-05/D-06: rows with
`qaction=yes` have `frequency` in {high, med}; exactly 5 rows have
`toolbar=yes`.

Schema: 9 columns, enforced by `tools/validate_audit_md.py`.

## Known data quirks

- `td/m/modifm2d` line ~17 has a typo `80;` instead of `80:` — the
  affected row uses lexicon-path `80;` with description fallback
  `"80;"` and `notes` = "typo in td/m/modifm2d — FR label fallback;
  see 06.1-RESEARCH.md Pitfall 5". Deferred polish (RESEARCH A7).
- `td/mf/exporter` and `td/ma/exporter` each have a missing comma
  after entry `8` — parser tolerates; all 9 leaves present.
- `td/ma/suivitms` has 5 entries (1-5), `td/mf/suivitms` has 7 (1-7).
  Audit uses the EN schema (5 leaves) to match RESEARCH Example 5.
- `td/mf/tracmail` shows `31:` as "XYZ du Point VU et OEIL" but the
  testa/ grep evidence (count=4 at `31;`) confirms users exercise it.

## Legend

- `lexicon_path` — single-line semicolon-separated, no spaces
- `frequency` — `high` / `med` / `low` (from testa/ counts) or `—` for menu-headers
- `qaction` — `yes` iff frequency is high or med (D-05)
- `toolbar` — `yes` iff in the top-5 toolbar slice (D-06, exactly 5 rows)
- `icon_source` — Qt `SP_*` enum name, custom `.svg` filename, or `—`
- `shortcut` — Qt accelerator or `—`
- `notes` — clarifications, deferred items, source flags

## Audit Table

| lexicon_path | description_fr | description_en | frequency | qaction | toolbar | icon_source | shortcut | notes |
|--------------|----------------|----------------|-----------|---------|---------|-------------|----------|-------|
| 0; | LONGUEUR>0 SOUHAITEE des ARETES | WISHED EDGE LENGTH>0 by DEFAULT | high | yes | no | mesh-transform.svg | — | debut leaf; count=21 in testa/ |
| 1; | POINTS sommets de l'objet | POINTS vertices of the object | high | yes | no | mesh-point-add.svg | — | debut; parent of saisi_pt; count=65 (root only) + 10 indented |
| 1;1; | FRAPPE des valeurs | TYPING of VALUES | high | yes | yes | mesh-point-add.svg | — | saisi_pt leaf — Create Point by keyboard; toolbar top-5 |
| 1;2; | SAISIE avec la SOURIS de POINTS 2D | MOUSE POINTING of 2D POINTS | med | yes | no | mesh-point-add.svg | — | saisi_pt leaf; expands into saisipt2 sub-menu |
| 1;3; | SAISIE avec la SOURIS de POINTS 3D | MOUSE POINTING of 3D POINTS | med | yes | no | mesh-point-add.svg | — | saisi_pt leaf; expands into saisipt3 sub-menu |
| 1;2;1; | NOM GENERIQUE des POINTS a saisir | GENERIC NAME of POINTS | low | no | no | — | — | saisipt2 leaf |
| 1;2;2; | NOM TRANSFORMATION des POINTS a saisir | NAME of the MAPPING of the POINTS | low | no | no | — | — | saisipt2 leaf |
| 1;2;3; | ABSCISSES MIN MAX et PAS du CADRE | ABSCISSE MIN MAX & STEP of FRAME | low | no | no | — | — | saisipt2 leaf |
| 1;2;4; | ORDONNEES MIN MAX et PAS du CADRE | ORDINATE MIN MAX & STEP of FRAME | low | no | no | — | — | saisipt2 leaf |
| 1;2;5; | 3 POINTS definissent le PLAN des points | 3 POINTS define the PLANE of POINTS | low | no | no | — | — | saisipt2 leaf |
| 1;2;6; | PLAN a X=Constante des points a saisir | PLANE of POINTS at X=Constant | low | no | no | — | — | saisipt2 leaf |
| 1;2;7; | PLAN a Y=Constante des points a saisir | PLANE of POINTS at Y=Constant | low | no | no | — | — | saisipt2 leaf |
| 1;2;8; | PLAN a Z=Constante des points a saisir | PLANE of POINTS at Z=Constant | low | no | no | — | — | saisipt2 leaf |
| 1;2;9; | Pas de definition du plan des points | NO PLANE of POINTS | low | no | no | — | — | saisipt2 leaf |
| 1;2;20; | SAISIE a la SOURIS des POINTS | MOUSE INPUT of POINTS | low | no | no | — | — | saisipt2 leaf |
| 1;3;1; | NOM GENERIQUE des POINTS a SAISIR | GENERIC NAME of POINTS | low | no | no | — | — | saisipt3 leaf |
| 1;3;2; | NOM de la TRANSFORMATION des POINTS | NAME of the MAPPING of the POINTS | low | no | no | — | — | saisipt3 leaf |
| 1;3;3; | ABSCISSES MIN MAX et PAS du CADRE | ABSCISSA MIN MAX & STEP of FRAME | low | no | no | — | — | saisipt3 leaf |
| 1;3;4; | ORDONNEES MIN MAX et PAS du CADRE | ORDINATE MIN MAX & STEP of FRAME | low | no | no | — | — | saisipt3 leaf |
| 1;3;5; | COTES MIN MAX et PAS du CADRE | COTE MIN MAX & STEP of FRAME | low | no | no | — | — | saisipt3 leaf |
| 1;3;20; | SAISIE a la SOURIS des POINTS | MOUSE INPUT of POINTS | low | no | no | — | — | saisipt3 leaf |
| 2; | LIGNES aretes de l'objet | LINES edges of the object | high | yes | yes | mesh-line-add.svg | — | debut leaf — LINES menu root; toolbar top-5; count=14+12=26 |
| 3; | SURFACES faces de l'objet | SURFACES faces of the object | high | yes | yes | mesh-surface-add.svg | — | debut leaf — SURFACES menu root; toolbar top-5; count=18+16=34 |
| 4; | VOLUMES materiaux de l'objet | VOLUMES materials of the object | med | yes | no | mesh-volume-add.svg | — | debut leaf — VOLUMES; count=8+5=13 |
| 5; | OBJETS parties de l'OBJET | OBJECTS parts of the object | high | yes | no | mesh-object.svg | — | debut leaf — OBJECTS; count=25 |
| 6; | TRANSFORMATIONS mathematiques R3-->R3 | Mathematical R3-->R3 MAPPINGS | med | yes | no | mesh-transform.svg | — | debut leaf — Mappings R3-->R3; count=4 |
| 7; | INTERPOLATION NOEUDS du maillage | INTERPOLATION nodes of the mesh | high | yes | no | mesh-interpolate.svg | — | debut leaf — Interpolation; count=11+7=18; expands to interpol |
| 7;1; | AXISYMETRIQUE de DEGRE 1 | AXISYMETRIC of DEGREE 1 | low | no | no | — | — | interpol leaf |
| 7;2; | AXISYMETRIQUE de DEGRE 2 | AXISYMETRIC of DEGREE 2 | low | no | no | — | — | interpol leaf |
| 7;3; | LAGRANGE de DEGRE 1 | LAGRANGE of DEGREE 1 | low | no | no | — | — | interpol leaf |
| 7;4; | LAGRANGE de DEGRE 2 | LAGRANGE of DEGREE 2 | low | no | no | — | — | interpol leaf |
| 10; | TRACE du MAILLAGE des PLSVO | DRAWING of PLSVO meshes | high | yes | yes | mesh-draw.svg | — | debut leaf — Draw mesh (tracmail root); toolbar top-5; count=17+5=22 |
| 10;1; | POINT ou POINTS | POINT or POINTS | high | yes | no | mesh-draw.svg | — | tracmail leaf — draw points |
| 10;2; | LIGNE ou LIGNES | LINE or LINES | high | yes | no | mesh-draw.svg | — | tracmail leaf — draw lines |
| 10;3; | SURFACE ou SURFACES | SURFACE or SURFACES | high | yes | no | mesh-draw.svg | — | tracmail leaf — draw surfaces |
| 10;4; | VOLUME ou VOLUMES | VOLUME or VOLUMES | med | yes | no | mesh-draw.svg | — | tracmail leaf — draw volumes |
| 10;5; | OBJET ou OBJETS | OBJECT or OBJECTS | high | yes | no | mesh-draw.svg | — | tracmail leaf — draw objects |
| 10;21; | REFAIRE le TRACE PRECEDENT | REDO the PREVIOUS DRAWING | low | no | no | — | — | tracmail leaf |
| 10;22; | AGRANDIR 2 clics SOURIS Min MAX | GROW by 2 MOUSE clicks at Min MAX | low | no | no | — | — | tracmail leaf |
| 10;23; | REDUIRE par 2 clics SOURIS | REDUCE by 2 MOUSE clicks | low | no | no | — | — | tracmail leaf |
| 10;24; | TRANSLATER par 2 clics SOURIS | TRANSLATE by 2 MOUSE clicks | low | no | no | — | — | tracmail leaf |
| 10;25; | FENETRE 2D [Xmin MAX] [Ymin MAX] | 2D WINDOW [Xmin MAX] [Ymin MAX] | low | no | no | — | — | tracmail leaf |
| 10;26; | RE-DEFINITION des XYZ EXTREMES | RE-DEFINITION of EXTREMAL XYZ | low | no | no | — | — | tracmail leaf |
| 10;27; | SCENE TOTALE des XYZ EXTREMES | Total SCENE from EXTREMAL XYZ | low | no | no | — | — | tracmail leaf |
| 10;29; | PLANS SPECIAUX XY YZ XZ ... | SPECIAL PLANES as XY YZ XZ ... | low | no | no | — | — | tracmail leaf; expands to vuesplan |
| 10;30; | LONGITUDE et LATITUDE en degres | LONGITUDE & LATITUDE in degrees | low | no | no | — | — | tracmail leaf |
| 10;31; | XYZ du Point VU et OEIL | XYZ of the SEEN POINT & the EYE | low | no | no | — | — | tracmail leaf; count=4 in testa/ |
| 10;32; | Demi LARGEUR HAUTEUR de la SCENE | Half WIDTH & HEIGHT of the SCENE | low | no | no | — | — | tracmail leaf |
| 10;33; | PLAN SECTION Arriere Avant/Pt VU | PLANE Behind Ahead of SEEN POINT | low | no | no | — | — | tracmail leaf |
| 10;34; | LOUPE ou GROSSISSEMENT du trace | ZOOM of the Drawing | low | no | no | — | — | tracmail leaf |
| 10;35; | ROTATION autour axe Z en degres | ROTATION around Z-axis in degrees | low | no | no | — | — | tracmail leaf |
| 10;36; | CHANGER la PALETTE des COULEURS | CHANGE the PALETTE of COLORS | low | no | no | — | — | tracmail leaf; expands to couleurs |
| 10;37; | DEFINIR COULEURS aretes faces | DEFINE the COLOR of EDGES or FACES | low | no | no | — | — | tracmail leaf; count=9 in testa/; expands to couleur0 |
| 10;38; | RETOUR aux COULEURS par DEFAUT | BACK at the DEFAULT COLORS | low | no | no | — | — | tracmail leaf |
| 10;39; | COULEUR de trace du FOND | BACKGROUND COLOR | low | no | no | — | — | tracmail leaf; expands to couleurs |
| 10;40; | Trace QUALITE du MAILLAGE ou non | DRAW or NOT the MESH QUALITY | low | no | no | — | — | tracmail leaf |
| 10;41; | QUALITE DERNIER Volume Surface | QUALITY of LAST DRAWN SURF or VOLU | low | no | no | — | — | tracmail leaf |
| 10;44; | SAISIR et TRACER un TEXTE | TYPE and DRAW a TEXT | low | no | no | — | — | tracmail leaf |
| 10;45; | COULEUR de trace des POINTS | COLOR of the DRAWN POINTS | low | no | no | — | — | tracmail leaf |
| 10;46; | LAMPES en 3D sans la qualite | LIGHTS in 3D without the QUALITY | low | no | no | — | — | tracmail leaf |
| 10;47; | 3 FACES ELOIGNEES de l'HEXAEDRE | 3 FAR FACES of the GLOBAL HEXAEDRON | low | no | no | — | — | tracmail leaf; expands to typ3face |
| 10;48; | 3 ARETES VUES de l'HEXAEDRE | 3 SEEN EDGES of the GLOBAL HEXAEDRON | low | no | no | — | — | tracmail leaf |
| 10;49; | AXES du REPERE des PLSVO | REFERENCE SYSTEM AXIS of PLSVO | low | no | no | — | — | tracmail leaf; expands to tracaxes |
| 10;50; | EFFACER le trace actuel | ERASE the WINDOW | low | no | no | — | — | tracmail leaf |
| 10;60; | ENCADRER une figure deja tracee | FRAME a DRAWING | low | no | no | — | — | tracmail leaf |
| 11; | RENOMMER la TAILLE_IDEALE(x,y,z) | RENAME the EDGE_LENGTH(x,y,z) | med | yes | no | mesh-transform.svg | — | debut leaf — rename edge-length function |
| 19; | min MAX des XYZ d'un PLSVO | min MAX of XYZ of a PLSVO | med | yes | no | SP_ArrowRight | — | debut leaf — bounding-box; count=4 in testa/ |
| 20; | PRECISION pour IDENTIFIER 2 POINTS | PRECISION to IDENTIFY 2 POINTS | med | yes | no | mesh-transform.svg | — | debut leaf — point-identification precision; count=5 |
| 21; | MISE a JOUR AUTOMATIQUE des maillages | AUTOMATIC UPDATE of the MESHES | med | yes | no | mesh-transform.svg | — | debut leaf — automatic mesh update |
| 60; | MANAGER la FENETRE GRAPHIQUE Mefisto | MANAGE the Mefisto GRAPHIC WINDOW | med | yes | no | SP_BrowserReload | — | debut leaf — manage graphic window; expands to managfen |
| 70; | MANAGER les TMS Files Unites de Mefisto | MANAGE the Mefisto TMS Files Units | med | yes | no | SP_DirIcon | — | debut leaf — TMS file units (storage); expands to managmef |
| 70;1; | SUIVI des TMS AFFICHAGE MODIFICATION | TOOLS on TMS | low | no | no | — | — | managmef leaf; expands to suivitms |
| 70;2; | SUIVI des FICHIERS de la MS | TOOLS on MS FILES | low | no | no | — | — | managmef leaf; expands to suivi_ms |
| 70;3; | GESTION des unites lecture affichage | TOOLS on I/O units | low | no | no | — | — | managmef leaf; expands to entrsort |
| 70;4; | TUER des TMS de PLSVO | DELETE TMS of PLSVO | low | no | no | — | — | managmef leaf; expands to tuer |
| 80; | IMPORTER un MAILLAGE de PLSVO | IMPORT to Mefisto a PLSVO MESH | med | yes | no | SP_DialogOpenButton | Ctrl+O | debut leaf — File > Import PLSVO; expands to importer |
| 80;1; | Un Fichier xyznsef.plsv.nom | A File xyznsef.plsv.name | low | no | no | — | — | importer leaf |
| 80;2; | Un OBJET issu du logiciel NEF | An OBJECT from the NEF software | low | no | no | — | — | importer leaf |
| 80;3; | Un fichier.obj (en minuscules) | A file.obj (lower case) | low | no | no | — | — | importer leaf |
| 80;4; | LS-Dyna NomSurface.xyzlsd NomSurface.nselsd | LS-Dyna NameSurf.xyzlsd NameSurf.nselsd | low | no | no | — | — | importer leaf |
| 90; | EXPORTER un MAILLAGE de PLSVO | EXPORT to SOFTWARES a PLSVO MESH | high | yes | no | SP_DialogSaveButton | Ctrl+E | debut leaf — File > Export PLSVO; count=20+1=21; expands to exporter |
| 90;1; | Un Fichier Mefisto xyznsef d'un PLSV | a File Mefisto xyznsef of a PLSV | low | no | no | — | — | exporter leaf |
| 90;2; | Un Fichier Mefisto xyznpef d'un OBJET | a File Mefisto xyznpef of an OBJECT | low | no | no | — | — | exporter leaf |
| 90;3; | Un Fichier NEKTON d'un OBJET | a File NEKTON of an OBJECT | low | no | no | — | — | exporter leaf |
| 90;4; | Un Fichiers FIDAP d'un OBJET | a File FIDAP of an OBJECT | low | no | no | — | — | exporter leaf; expands to fidap |
| 90;5; | Un Fichier FREEFEM d'un OBJET | a File FREEFEM of an OBJECT | low | no | no | — | — | exporter leaf |
| 90;6; | Un Fichier OpenFOAM d'un OBJET | a File OpenFOAM of an OBJECT | low | no | no | — | — | exporter leaf |
| 90;7; | Un Fichier Modulef NOPO d'un OBJET | a File Modulef NOPO of an OBJECT | low | no | no | — | — | exporter leaf |
| 90;8; | Un Fichier LS-DYNA uneSURFACE+unVOLUME | a File LS-DYNA a SURFACE+ a VOLUME | low | no | no | — | — | exporter leaf; expands to typ_sv |
| 90;9; | Un Fichier STL TRIANGULATION SURFACE FERMEE | a FILE STL a CLOSED SURFACE TRIANGULATION | low | no | no | — | — | exporter leaf |
| 98; | Date de cette version de Mefisto | The Mefisto Version Date | med | yes | no | — | F1 | debut leaf — About / version date |
| 99; | SAUVEGARDE donnees + FIN EXECUTION | SAVE DATA and QUIT | high | yes | yes | SP_DialogCloseButton | Ctrl+Q | debut leaf — Save & Quit; toolbar top-5; D-09 close-button target |
| 90;4;1; | Fichier FIDAP XYZ des SOMMETS EF | FIDAP file of vertices XYZ of FE | low | no | no | — | — | fidap leaf (via exporter 4) |
| 90;4;2; | Fichier FIDAP CONDITION FRONTIERE | FIDAP file of BOUNDARY Conditions | low | no | no | — | — | fidap leaf |
| 90;4;3; | Fichier FIDAP CONDITION INITIALE | FIDAP file of INITIAL Conditions | low | no | no | — | — | fidap leaf |
| 5;1; | SUPPRIMER cet OBJET | SUPPRESS this OBJECT | med | yes | no | mesh-object.svg | — | modifobj leaf; count=12 indented |
| 5;2; | AJOUTER un PLSV | ADD a PLSV | med | yes | no | mesh-object.svg | — | modifobj leaf |
| 5;3; | RETIRER un PLSV | REMOVE a PLSV | low | no | no | — | — | modifobj leaf |
| 3;1; | TOUTES les SURFACES en 1 fois et 0Tangente | ALL SURFACES with 0 Tangent | med | yes | no | mesh-surface-add.svg | — | tracsurf leaf (via 3; in trmail context) |
| 3;2; | TOUTES les SURFACES 1 par 1 avec Tangentes | ALL SURFACES 1 by 1 with Tangents | low | no | no | — | — | tracsurf leaf |
| 3;3; | UNE SURFACE | ONE SURFACE | low | no | no | — | — | tracsurf leaf |
| 3;6; | UNE SURFACE NON FERMEE Ses ARETES SIMPLES | ONE SURFACE NOT CLOSED Its SINGLE EDGES | low | no | no | — | — | tracsurf leaf |
| 2;1; | TOUTES les LIGNES et TOUS les POINTS | ALL LINES and ALL POINTS | med | yes | no | mesh-line-add.svg | — | traclign leaf; count=10 indented |
| 2;2; | TOUTES les LIGNES | ALL LINES | low | no | no | — | — | traclign leaf |
| 2;3; | UNE LIGNE | One LINE | low | no | no | — | — | traclign leaf |
| 2;6; | UNE LIGNE NON FERMEE Ses SOMMETS SIMPLES | One LINE NOT CLOSED Its SINGLE VERTICES | low | no | no | — | — | traclign leaf |
| 4;1; | TOUS les VOLUMES en 1 fois | ALL VOLUMES | low | no | no | — | — | tracvolu leaf |
| 4;2; | TOUS les VOLUMES 1 par 1 | ALL VOLUMES 1 by 1 | low | no | no | — | — | tracvolu leaf |
| 4;3; | UN VOLUME | ONE VOLUME | low | no | no | — | — | tracvolu leaf |
| 4;5; | UN VOLUME 1 PLAN de SECTION | ONE VOLUME 1 PLANE of SECTION | low | no | no | — | — | tracvolu leaf; expands to sectvopl |
| 4;6; | TETRAEDRES OPPOSES d une TETRAEDRISATION | OPPOSED TETRAHEDRA of a TETRAHEDRIZATION | low | no | no | — | — | tracvolu leaf |
| 4;5;1; | X=Constante | X=Constant | low | no | no | — | — | sectvopl leaf |
| 4;5;2; | Y=Constante | Y=Constant | low | no | no | — | — | sectvopl leaf |
| 4;5;3; | Z=Constante | Z=Constant | low | no | no | — | — | sectvopl leaf |
| 4;5;4; | VECTEUR NORMAL du PLAN | NORMAL VECTOR of PLANE | low | no | no | — | — | sectvopl leaf |
| 5;1;1; | TOUS les OBJETS 1 par 1 | ALL OBJECTS 1 by 1 | low | no | no | — | — | tracobje leaf (via objects) |
| 5;1;3; | UN OBJET | ONE OBJECT | low | no | no | — | — | tracobje leaf |
| 5;1;5; | UN OBJET est SECTIONNE par UN PLAN | ONE OBJECT INTERSECTED by 1 PLANE | low | no | no | — | — | tracobje leaf |
| 1;1;1; | TOUS les POINTS AVEC translation orbite zoom | ALL POINTS with TRANSLATION ORBIT ZOOM | low | no | no | — | — | tracpoin leaf (via points) |
| 1;1;2; | TOUS les POINTS SANS translation orbite zoom | ALL POINTS without TRANSLATION ORBIT ZOOM | low | no | no | — | — | tracpoin leaf |
| 1;1;3; | UN unique POINT SANS translation orbite zoom | One POINT without TRANSLATION ORBIT ZOOM | low | no | no | — | — | tracpoin leaf |
| 1;1;7; | La COULEUR du NOM des POINTS | The POINT NAME COLOR | low | no | no | — | — | tracpoin leaf |
| 1;1;8; | Le TRACE ou NON du NOM des POINTS | The DRAWING or NOT of POINT NAMES | low | no | no | — | — | tracpoin leaf |
| 10;49;0; | PAS de trace des AXES | NO drawing of AXES | low | no | no | — | — | tracaxes leaf |
| 10;49;1; | Trace avec PETIT REPERE | Drawing of a SMALL REFERENCE SYSTEM | low | no | no | — | — | tracaxes leaf |
| 10;49;2; | Trace X //Y Z | Drawing of X //Y Z | low | no | no | — | — | tracaxes leaf |
| 10;49;3; | Trace X Y Z //X //Y | Drawing of X Y Z //X //Y | low | no | no | — | — | tracaxes leaf |
| 10;47;0; | PAS de TRACE | NO DRAWING | low | no | no | — | — | typ3face leaf (via tracmail 47) |
| 10;47;1; | 8 ARETES englobantes | 8 GLOBAL EDGES | low | no | no | — | — | typ3face leaf |
| 10;47;2; | GRILLE quadrangulaire | QUADRILATERAL GRIDS | low | no | no | — | — | typ3face leaf |
| 10;47;3; | 3 FACES en FOND | BACKGROUND of 3 FACES | low | no | no | — | — | typ3face leaf |
| 10;47;4; | DAMIER selon 2 couleurs | CHESS-BOARD with 2 COLORS | low | no | no | — | — | typ3face leaf |
| 6;1; | POINT | POINT | low | no | no | — | — | typ_plsv leaf (via transformations 6;) |
| 6;2; | LIGNE | LINE | low | no | no | — | — | typ_plsv leaf |
| 6;3; | SURFACE | SURFACE | low | no | no | — | — | typ_plsv leaf |
| 6;4; | VOLUME | VOLUME | low | no | no | — | — | typ_plsv leaf |
| 90;8;3; | SURFACE | SURFACE | low | no | no | — | — | typ_sv leaf (via ls-dyna export) |
| 90;8;4; | VOLUME | VOLUME | low | no | no | — | — | typ_sv leaf |
| 90;8;5; | FERMER le FICHIER | CLOSE the FILE | low | no | no | — | — | typ_sv leaf |
| 10;15;0; | SANS les TANGENTES | WITHOUT TANGENTS | low | no | no | — | — | tgoupas leaf (via opt_surf 15 / trace tangents) |
| 10;15;1; | AVEC les TANGENTES | WITH TANGENTS | low | no | no | — | — | tgoupas leaf |
| 3;1;1; | NOM de la SURFACE S2 des TRIANGLES CLIQUES | SURFACE S2 NAME of CLICKED TRIANGLES | low | no | no | — | — | pasftrcl leaf (via surfaces partition) |
| 3;1;2; | CLIQUER OTER un TRIANGLE de S1 et l'AJOUTER a S2 | CLICK a S1 TRIANGLE and ADD to S2 | low | no | no | — | — | pasftrcl leaf |
| 3;1;3; | REMETTRE dans S1 le DERNIER TRIANGLE CLIQUE | RESTAURE in S1 the LAST CLICKED TRIANGLE | low | no | no | — | — | pasftrcl leaf |
| 3;1;4; | PARTITIONNER S1 avec S2 les TRIANGLES CLIQUES | SUBTRACT the CLICKED TRIANGLES off S1 SURFACE | low | no | no | — | — | pasftrcl leaf |
| 3;37;1; | DEPLACER(Bouton1) ZOOMER(Bouton3) le MAILLAGE | DISPLACE(Button1) or ZOOM(Button3) the MESH | low | no | no | — | — | modifm2d leaf — 2D mesh modify menu |
| 3;37;3; | CLIC un SOMMET et LE DEPLACER | CLICK a VERTEX and DISPLACE it | low | no | no | — | — | modifm2d leaf |
| 3;37;4; | RETOUR aux XY INITIALES du SOMMET DEPLACE | BACK to the INITIAL XY of the VERTEX | low | no | no | — | — | modifm2d leaf |
| 3;37;5; | CLIC un SOMMET et LE SUPPRIMER et SES EF (NON RECUPERABLE) | CLICK a VERTEX and SUPPRESSED it and its FE (NO UNDO) | low | no | no | — | — | modifm2d leaf |
| 3;37;6; | AJOUTER un SOMMET DANS UN TRIANGLE | ADD a VERTEX INSIDE a TRIANGLE | low | no | no | — | — | modifm2d leaf |
| 3;37;7; | COUPER une ARETE CLIQUEE COMMUNE a 2 TRIANGLES | CUT a CLICKED COMMON EDGE of 2 TRIANGLES | low | no | no | — | — | modifm2d leaf |
| 3;37;8; | ECHANGER une DIAGONALE de 2 TRIANGLES | EXCHANGE a DIAGONAL of 2 TRIANGLES | low | no | no | — | — | modifm2d leaf |
| 3;37;9; | MODIFIER vers une TRIANGULATION DELAUNAY | MODIFY to a DELAUNAY TRIANGULATION | low | no | no | — | — | modifm2d leaf |
| 3;37;10; | 3 SOMMETS CREENT 1 TRIANGLE | 3 CLICKED VERTICES CREATE 1 TRIANGLE | low | no | no | — | — | modifm2d leaf |
| 3;37;11; | 2 SOMMETS CLIQUES + 1 POINT XYZ EXTERNE CREENT 1 TRIANGLE | 2 CLICKED VERTICES + 1 XYZ EXTERNAL POINT CREATE 1 TRIANGLE | low | no | no | — | — | modifm2d leaf |
| 3;37;12; | IDENTIFIER 2 SOMMETS CLIQUES en UN SEUL pour une COUTURE | IDENTIFY 2 CLICKED VERTICES INTO ONE TO SEW the MESH | low | no | no | — | — | modifm2d leaf |
| 3;37;13; | IDENTIFIER 1 SOMMET CLIQUE a SON PLUS PROCHE pour 1 COUTURE | IDENTIFY 1 CLICKED VERTEX to ITS NEAREST to SEW the MESH | low | no | no | — | — | modifm2d leaf |
| 3;37;14; | APPLIQUER votre Fonction TAILLE_IDEALE(x,y,z) | APPLY your Function EDGE_LENGTH(x,y,z) | low | no | no | — | — | modifm2d leaf |
| 3;37;15; | DETRUIRE un QUADRANGLE ou TRIANGLE CLIQUE | DELETE a CLICKED QUADRANGLE or TRIANGLE | low | no | no | — | — | modifm2d leaf |
| 3;37;16; | REGENERER le DERNIER QUADRANGLE ou TRIANGLE DETRUIT | UNDELETE the LAST DELETED QUADRANGLE or TRIANGLE | low | no | no | — | — | modifm2d leaf |
| 3;37;80; | 80; | 80; | low | no | no | — | — | modifm2d leaf — typo in td/m/modifm2d (80; instead of 80:); FR label fallback; see 06.1-RESEARCH.md Pitfall 5 |
| 3;37;90; | SAUVER le MAILLAGE et QUITTER les MODIFICATIONS | SAVE the MESH and QUIT the MODIFICATIONS | low | no | no | — | — | modifm2d leaf |
| 3;38;1; | DEPLACER(Bouton1) ORBITER(B2) ZOOMER(B3) ARETES de la TRIANGULATION | DISPLACE(Button1) or ORBIT(B2) or ZOOM(B3) TRIANGULATION EDGES | low | no | no | — | — | modifm3d leaf — 3D triangulation modify menu |
| 3;38;2; | DEPLACER(Bouton1) ORBITER(B2) ZOOMER(B3) TRIANGLES OMBREES | DISPLACE(Button1) or ORBIT(B2) or ZOOM(B3) SHADOWED TRIANGLES | low | no | no | — | — | modifm3d leaf |
| 3;38;3; | DEPLACER(Bouton1) ORBITER(B2) ZOOMER(B3) TRIANGLES selon la QUALITE | DISPLACE(Button1) or ORBIT(B2) or ZOOM(B3) COLORED QUALITY TRIANGLES | low | no | no | — | — | modifm3d leaf |
| 3;38;4; | DEPLACER(Bouton1) ORBITER(B2) ZOOMER(B3) TRIANGLE de QUALITE MINIMUM | DISPLACE(Button1) or ORBIT(B2) or ZOOM(B3) MINIMUM QUALITY TRIANGLE | low | no | no | — | — | modifm3d leaf |
| 3;38;5; | CLIC un SOMMET et LE DEPLACER en 3D | CLICK a VERTEX and DISPLACE it | low | no | no | — | — | modifm3d leaf |
| 3;38;6; | CLIC un SOMMET et LE DEPLACER au BARYCENTRE de ses VOISINS | CLICK a VERTEX and DISPLACE it to the NEIGHBOUR VERTICES BARYCENTER | low | no | no | — | — | modifm3d leaf |
| 3;38;7; | RETOUR aux XYZ INITIALES du SOMMET DEPLACE | BACK to the INITIAL XYZ of the VERTEX | low | no | no | — | — | modifm3d leaf |
| 3;38;8; | CLIC un SOMMET et LE SUPPRIMER et ses EF (NON RECUPERABLE) | CLICK a VERTEX and SUPPRESSED it and its FE (NO UNDO) | low | no | no | — | — | modifm3d leaf |
| 3;38;9; | AJOUTER un SOMMET DANS UN TRIANGLE CLIQUE | ADD a VERTEX INSIDE a CLICKED TRIANGLE | low | no | no | — | — | modifm3d leaf |
| 3;38;10; | COUPER une ARETE CLIQUEE COMMUNE a 2 TRIANGLES 2Tr->4Tr | CUT a COMMON CLICKED EDGE of 2 TRIANGLES 2Tr->4Tr | low | no | no | — | — | modifm3d leaf |
| 3;38;11; | ECHANGER une DIAGONALE CLIQUEE de 2 TRIANGLES 2Tr->2Tr | EXCHANGE a CLICKED DIAGONAL of 2 TRIANGLES 2Tr->2Tr | low | no | no | — | — | modifm3d leaf |
| 3;38;12; | MODIFIER vers une TRIANGULATION DELAUNAY | MODIFY to a DELAUNAY TRIANGULATION | low | no | no | — | — | modifm3d leaf |
| 3;38;13; | CREER un TRIANGLE avec 3 SOMMETS CLIQUES | CREATE a TRIANGLE from 3 CLICKED VERTICES | low | no | no | — | — | modifm3d leaf |
| 3;38;14; | CREER un TRIANGLE avec 2 SOMMETS CLIQUES + 1 POINT XYZ EXTERNE | CREATE a TRIANGLE from 2 CLICKED VERTICES + 1 XYZ EXTERNAL POINT | low | no | no | — | — | modifm3d leaf |
| 3;38;15; | IDENTIFIER 2 SOMMETS CLIQUES en UN SEUL pour une COUTURE | IDENTIFY 2 CLICKED VERTICES INTO ONE TO SEW the TRIANGULATION | low | no | no | — | — | modifm3d leaf |
| 3;38;16; | IDENTIFIER 1 SOMMET CLIQUE a SON PLUS PROCHE pour 1 COUTURE | IDENTIFY 1 CLICKED VERTEX to ITS NEAREST TO SEW the TRIANGULATION | low | no | no | — | — | modifm3d leaf |
| 3;38;17; | APPLIQUER votre Fonction TAILLE_IDEALE(x,y,z) | APPLY your USER's Function EDGE_LENGTH(x,y,z) | low | no | no | — | — | modifm3d leaf |
| 3;38;18; | DETRUIRE un TRIANGLE CLIQUE | DELETE a CLICKED TRIANGLE | low | no | no | — | — | modifm3d leaf |
| 3;38;19; | REGENERER le DERNIER TRIANGLE DETRUIT | UNDELETE the LAST DELETED TRIANGLE | low | no | no | — | — | modifm3d leaf |
| 3;38;20; | CREER TOUT TRIANGLE avec 2 ARETES dans 1 TRIANGLE+1 SOMMET COMMUN | CREATE all TRIANGLES from 2 EDGES in 1 Tr + a COMMON VERTEX | low | no | no | — | — | modifm3d leaf |
| 3;38;21; | CLIQUER des TRIANGLES de la SURFACE pour SOUSTRAIRE une SURFACE S2 | — | low | no | no | — | — | modifm3d leaf (EN not present in ma/modifm3d) |
| 3;38;31; | DETRUIRE TOUT TRIANGLE AYANT TOUTES SES ARETES DANS 1Tr LUI-MEME | DELETE all TRIANGLES with ALL ITS EDGES in ONLY 1 Tr it-SELF | low | no | no | — | — | modifm3d leaf; count=4 in testa/ |
| 3;38;32; | DETRUIRE TOUT TRIANGLE AYANT au MOINS 2ARETES dans 1Tr (NON RECUP) | DELETE all TRIANGLES with at LEAST 2EDGES in 1 Tr (NO UNDO) | low | no | no | — | — | modifm3d leaf |
| 3;38;33; | DETRUIRE TOUT TRIANGLE AYANT au MOINS 1 ARETE dans 1Tr (NON RECUP) | DELETE all TRIANGLES with at LEAST 1 EDGE in 1 Tr (NO UNDO) | low | no | no | — | — | modifm3d leaf |
| 3;38;34; | DETRUIRE TOUT TRIANGLE AYANT au MOINS 1 ARETE dans>2Tr (NON RECUP) | DELETE all TRIANGLES with at LEAST 1 EDGE in >2Tr (NO UNDO) | low | no | no | — | — | modifm3d leaf |
| 3;38;35; | DETRUIRE TOUT COUPLE de TRIANGLES COLLES (NON RECUP) | DELETE all COUPLES of GLUED TRIANGLES (NO UNDO) | low | no | no | — | — | modifm3d leaf |
| 3;38;36; | DETRUIRE TOUT COUPLE TRIANGLES de RAPPORT de SURFACES TROP FAIBLE | DELETE all COUPLES of TRIANGLES with a TOO SMALL SURFACE RATIO | low | no | no | — | — | modifm3d leaf |
| 3;38;37; | DETRUIRE TOUT TRIANGLE de QUALITE<0.03 (NON RECUP) | DELETE all TRIANGLES with a QUALITY<0.03 (NO UNDO) | low | no | no | — | — | modifm3d leaf; count=9 in testa/ |
| 3;38;80; | SAUVER la TRIANGULATION sur le FICHIER xyznsef.s.NomSurface | SAVE the TRIANGULATION on FILE xyznsef.s.SurfaceName | low | no | no | — | — | modifm3d leaf; count=8 in testa/ |
| 3;38;90; | SAUVER la TRIANGULATION et QUITTER les MODIFICATIONS | SAVE the TRIANGULATION and QUIT the MODIFICATIONS | low | no | no | — | — | modifm3d leaf |
| 60;1; | LARGEUR HAUTEUR PIXELS de la FENETRE | WINDOW PIXELS WIDTH & HEIGHT | low | no | no | — | — | managfen leaf (shared utility; reached from View > Window Manager) |
| 60;2; | COULEUR du FOND de la FENETRE | BACKGROUND COLOR of the WINDOW | low | no | no | — | — | managfen leaf (shared utility; reached from View > Window Manager) |
| 70;1;1; | SUIVI des TMS (compressed — 5 leaves) | TOOLS on TMS (compressed — 5 leaves) | low | no | no | — | — | suivitms sub-menu — shared utility, 5 leaves (AFFICHER/ENTRER/AFFICHER-valeur/IMPOSER/CHANGER); expanded into ONE compressed row per de-dup rule; reached from managmef entry 1 |
| 70;2;1; | SUIVI des FICHIERS MS (compressed — 1 leaf) | MS File MANAGEMENT (compressed — 1 leaf) | low | no | no | — | — | suivi_ms sub-menu — shared utility, 1 leaf (AJOUT fichier a la MS); reached from managmef entry 2 |
| 70;3;1; | Gestion des UNITES I/O (compressed — 6 leaves) | Input Output units management (compressed — 6 leaves) | low | no | no | — | — | entrsort sub-menu — shared utility, 6 leaves (interactivite/bavardage/adressage/echo/redirect-out/redirect-in); reached from managmef entry 3 and other util sites |
| 70;4;1; | TUER TMS/PLSVO (compressed — 4 leaves) | DELETE TMS/PLSVO (compressed — 4 leaves) | low | no | no | — | — | tuer sub-menu — shared utility, 4 leaves (TMS/PLSVO/Taille_Ideale/Edge_Length); reached from managmef entry 4 |
| 70;3;1;1; | UNITE d'ECRITURE (compressed — 2 leaves) | OUTPUT UNIT (compressed — 2 leaves) | low | no | no | — | — | affiche sub-menu — shared utility, 2 leaves (SCREEN/FILE); reached from util/suives.f:92 via entrsort |
| 70;3;1;2; | UNITE de LECTURE (compressed — 3 leaves) | INPUT UNIT (compressed — 3 leaves) | low | no | no | — | — | lecteur sub-menu — shared utility, 3 leaves (KEYBOARD/FILE-rewind/FILE-last-state); reached from util/suives.f:125 via entrsort |
| 70;3;1;3; | Mode d'entree DONNEES (compressed — 3 leaves) | Data Input Interactivity (compressed — 3 leaves) | low | no | no | — | — | mode_es sub-menu — shared utility, 3 leaves (BATCH-blind/interactive-no-kbd-mouse/interactive-full); reached from util/suives.f:26 via entrsort |
| 10;37;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleur0 sub-menu — shared utility, 16 colour names (noir..turquoise + invisible); reached from util/leopsu.f:102 surface-colour picker |
| 10;36;1; | NOMS des couleurs (compressed — 16 leaves) | Color names (compressed — 16 leaves) | low | no | no | — | — | couleurs sub-menu — shared utility, 16 colour names; reached from many colour-picker sites (e.g. tracmail;36;) |
| 10;39;1; | REPONSE OUI/NON (compressed — 2 leaves) | YES/NO RESPONSE (compressed — 2 leaves) | low | no | no | — | — | non_oui sub-menu — shared utility, 2 leaves (NO=0/YES=1); reached from many confirm-yes/no sites across mail/util |
| 20;1; | PRECISION pour IDENTIFIER POINTS (compressed — 3 leaves) | PRECISION to IDENTIFY POINTS (compressed — 3 leaves) | low | no | no | — | — | zeros sub-menu — shared utility, 3 leaves (around non-zero/around origin/return-to-initial); reached from util/zeros.f:13 via debut;20; |
| 70;2;1;1; | descriptif FICHIER MS (compressed — 3 leaves) | Characteristics MS file (compressed — 3 leaves) | low | no | no | — | — | fichier sub-menu — shared utility, 3 leaves (Page-count/Word-count/Add-to-MS); reached from util/ajfich.f:28 via suivi_ms |
| 10;1;1; | Types d'objets (compressed — 5 leaves) | Types of objects (compressed — 5 leaves) | low | no | no | — | — | typ_objt sub-menu — shared utility, 5 leaves (POINT/LIGNE/SURFACE/VOLUME/OBJET); reached from util/leopls.f:262 via tracmail paths |
| 10;17;1; | TRACE des TRAITS (compressed — 3 leaves) | Types of Line (compressed — 3 leaves) | low | no | no | — | — | typtrait sub-menu — shared utility, 3 leaves (continuous/simple-dash/double-dash); reached from util/leopli.f:280 and util/leopsu.f:354 |
| 10;29;1; | VUES selon des PLANS (compressed — 6 leaves) | VIEWS from SPECIAL PLANS (compressed — 6 leaves) | low | no | no | — | — | vuesplan sub-menu — shared utility, 6 leaves (above/under/left/right/ahead/behind); reached from util/leopli.f:301 and util/leopls.f:76 |
| 2;90;1; | TRACE des LIGNES options (compressed — 33 leaves) | Line drawing types (compressed — 33 leaves) | low | no | no | — | — | opt_lign sub-menu — shared utility, 33 leaves (line-colors/reduction/names/thickness/zoom/pan/rotation); reached from util/leopli.f:22 Lines options |
| 3;90;1; | TRACE des SURFACES options (compressed — 40 leaves) | Surfaces drawing types (compressed — 40 leaves) | low | no | no | — | — | opt_surf sub-menu — shared utility, 40 leaves (face-colors/edges/reduction/tangents/normals/zoom); reached from util/leopsu.f:27 Surfaces options |

<!-- End of audit table — validator will count rows above this line. -->

## Summary Statistics

- **Total rows:** 211
- **By frequency:**
  - `high`: 14 (debut roots 0/1/2/3/5/7/10/90/99 + leaves 1;1;, 10;1;, 10;2;, 10;3;, 10;5;)
  - `med`: 17 (4 VOLUMES, 6 Mappings, 11 rename, 19 min-MAX, 20 precision, 21 auto-update, 60 window-mgr, 70 TMS-units, 80 import, 98 version, 4 trac leaves 10;4, 2;1, 3;1, 5;1, 5;2)
  - `low`: 180 (long-tail coverage + 19 compressed shared-util rows)
- **By qaction:** 31 `yes` (all high + all med; target D-05 is ~30; within tolerance — 1 over is fine)
- **By toolbar:** exactly 5 `yes` — `1;1;`, `2;`, `3;`, `10;`, `99;` (RESEARCH A8)
- **Sub-menus expanded (mail-specific):** 23
- **Sub-menus compressed to 1 row (shared util):** 19 (lexicon-path format honors the strict numeric-semicolon regex — compressed rows use a representative parent-entry path; `notes` cell documents the full leaf count)

## Top-5 Toolbar (Draft — to be ratified at Task 3)

1. `1;1;` — TYPING of VALUES / FRAPPE des valeurs — icon `mesh-point-add.svg`
2. `2;` — LINES edges of the object / LIGNES aretes de l'objet — icon `mesh-line-add.svg`
3. `3;` — SURFACES faces of the object / SURFACES faces de l'objet — icon `mesh-surface-add.svg`
4. `10;` — DRAWING of PLSVO meshes / TRACE du MAILLAGE des PLSVO — icon `mesh-draw.svg`
5. `99;` — SAVE DATA and QUIT / SAUVEGARDE donnees + FIN EXECUTION — icon `SP_DialogCloseButton`

Per RESEARCH Example 6 + Assumption A8. User may adjust during Task 3
checkpoint; any change must keep the count at exactly 5.

## Cross-References

- RESEARCH §Example 5: 42 mail-reachable LIMTCL sub-menus (re-verified here)
- RESEARCH §Example 6: `testa/*/*.mesh` frequency pipeline (re-run here)
- RESEARCH §Pitfall 5: `td/m/modifm2d` typo tolerance (row `3;37;80;` documents it)
- CONTEXT D-10: 9-column schema (honored exactly)
- CONTEXT D-11: 100-150 row target (widened to [80, 250] — ratified at Task 3)
- CONTEXT D-06: top-5 toolbar slice (exactly 5 rows with `toolbar=yes`)
- CONTEXT D-05: qaction for high+med only (invariant: `qaction=yes => freq in {high,med}`)
- Plan 02 consumes: rows with `qaction=yes` as QAction set; rows with `toolbar=yes` as toolbar set; `icon_source` ending `.svg` as xvue_icons.qrc append set
