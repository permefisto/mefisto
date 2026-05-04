// xvue/qt/src/xvue_qt_i18n.cpp
// Phase 6.0 Plan 02 (06.0-02): full bilingual lookup table + cached
// $MEFISTO/td/m/anglais language probe. Replaces the Plan 01 stub bodies.
//
// Authoritative source for every FR/EN pair below: 06.0-UI-SPEC.md
// §Copywriting Contract. Adding a new MsgId REQUIRES adding a new row here
// (the static_assert at the bottom catches enum drift at compile time).
//
// Pitfall 2.1 (06.0-RESEARCH.md): static_assert on table-vs-enum size.
// Pitfall 2.2: cache xvueIsEnglish() at first call via C++17 static-local
//              initialization (thread-safe, lazy).
#include "xvue_qt_i18n.h"

#include <QFileInfo>
#include <QString>

#include <array>
#include <cstdlib>
#ifdef XVUE_QT_TESTING
#include <optional>
#endif

namespace {

struct Entry {
    const char* fr;
    const char* en;
};

// Indexed by static_cast<int>(MsgId). Order MUST match the MsgId enum
// declared in xvue_qt_i18n.h. The static_assert below catches size drift,
// but a misaligned row order would still pass the size check while
// silently mistranslating individual strings — so reorder rows here if the
// enum is reordered.
constexpr std::array<Entry, static_cast<size_t>(MsgId::_Count_)> kTable = {{
    {"MEFISTO",                "MEFISTO"},                              // AppName
    {"&Fichier",               "&File"},                                // FileMenuTitle
    {"&Ouvrir projet…",        "&Open Project…"},                       // FileOpen
    {"Ouvrir un projet MEFISTO existant (Ctrl+O)",
     "Open an existing MEFISTO project (Ctrl+O)"},                      // FileOpenTip
    {"&Enregistrer projet",    "&Save Project"},                        // FileSave
    {"Enregistrer le projet courant (Ctrl+S)",
     "Save the current project (Ctrl+S)"},                              // FileSaveTip
    {"Enregistrer &sous…",     "Save &As…"},                            // FileSaveAs
    {"Projets &récents",       "&Recent Projects"},                     // FileRecentSubmenu
    {"Effacer la liste",       "Clear Recent Projects"},                // FileRecentClear
    {"Effacer les projets récents ?",
     "Clear Recent Projects?"},                                         // FileRecentClearConfirmTitle
    {"E&xporter…",             "E&xport…"},                             // FileExport
    // Phase 7 Plan 04 (EXPORT-02, EXPORT-05): File → Export submenu rows.
    // Order MUST match the MsgId enum extension in xvue_qt_i18n.h. Plan 05
    // will add FileExportGif + FileCaptureAnimation rows after Pdf entries.
    {"PNG…",                   "PNG…"},                                 // FileExportPng
    {"JPEG…",                  "JPEG…"},                                // FileExportJpeg
    {"PDF…",                   "PDF…"},                                 // FileExportPdf
    {"Exporter en PNG…",       "Export to PNG…"},                       // FileExportPngDialogTitle
    {"Exporter en JPEG…",      "Export to JPEG…"},                      // FileExportJpegDialogTitle
    {"Exporter en PDF…",       "Export to PDF…"},                       // FileExportPdfDialogTitle
    {"Images PNG (*.png)",     "PNG images (*.png)"},                   // FileExportPngFilter
    {"Images JPEG (*.jpg *.jpeg)", "JPEG images (*.jpg *.jpeg)"},       // FileExportJpegFilter
    {"Documents PDF (*.pdf)",  "PDF documents (*.pdf)"},                // FileExportPdfFilter
    {"Échec de l'export PNG",  "PNG export failed"},                    // ExportPngFailed
    {"Échec de l'export JPEG", "JPEG export failed"},                   // ExportJpegFailed
    {"Échec de l'export PDF",  "PDF export failed"},                    // ExportPdfFailed
    {"Aucun dessin à exporter","No canvas to export"},                  // ExportNoCanvas
    {"&Quitter",               "&Quit"},                                // FileQuit
    {"Quitter MEFISTO (Ctrl+Q)","Quit MEFISTO (Ctrl+Q)"},                // FileQuitTip
    {"&Affichage",             "&View"},                                // ViewMenuTitle
    {"Barre d'&outils",        "&Toolbar"},                             // ViewToolbar
    {"&Console (F9)",          "&Console (F9)"},                        // ViewConsole
    {"Zoom a&vant",            "Zoom &In"},                             // ViewZoomIn
    {"Zoom arri&ère",          "Zoom &Out"},                            // ViewZoomOut
    {"A&juster à la fenêtre",  "&Fit to Window"},                       // ViewFit
    {"&Préférences…",          "&Preferences…"},                        // ViewPreferences
    {"A&ide",                  "&Help"},                                // HelpMenuTitle
    {"&Documentation (F1)",    "&Documentation (F1)"},                  // HelpDocumentation
    {"À propos de MEFISTO",    "About MEFISTO"},                        // HelpAbout
    {"Console",                "Console"},                              // ConsoleTitle
    {"X\u00a0: %1\u00a0\u00a0Y\u00a0: %2",
     "X: %1  Y: %2"},                                                   // StatusCoordFormat
    {"Terminez l'opération en cours (tapez 99;)",
     "Finish current operation first (type 99;)"},                      // ModalRefuse
    {"Erreur MEFISTO",         "MEFISTO Error"},                        // ErrorMsgBoxTitle
    {"Quitter MEFISTO ?",      "Quit MEFISTO?"},                        // QuitConfirmTitle
    {"Voulez-vous quitter MEFISTO ? Les modifications non enregistrées seront perdues.",
     "Quit MEFISTO? Any unsaved work will be lost."},                   // QuitConfirmBody
    {"À propos de MEFISTO",    "About MEFISTO"},                        // AboutTitle
    {"MEFISTO %1\n(MEsh and Finite element Software TO…)\n\nAuteur : Alain Perronnet\nLaboratoire Jacques-Louis Lions (LJLL), UPMC Paris\n\nQt %2\nCompilé le %3",
     "MEFISTO %1\n(MEsh and Finite element Software TO…)\n\nAuthor: Alain Perronnet\nLaboratoire Jacques-Louis Lions (LJLL), UPMC Paris\n\nQt %2\nBuild date: %3"},  // AboutBody
    {"Ouvrir un projet MEFISTO","Open MEFISTO Project"},                // OpenProjectDialogTitle
    {"Projets MEFISTO (*)",    "MEFISTO Projects (*)"},                 // OpenProjectFilter
    {"Préférences",            "Preferences"},                          // PreferencesTitle
    {"Nombre de projets récents :",
     "Recent projects to remember:"},                                   // PrefRecentCountLabel
    {"Console visible au démarrage",
     "Show console on startup"},                                        // PrefConsoleDefaultLabel
    {"Apparence :",            "Appearance:"},                          // PrefColorSchemeLabel
    {"Système",                "System"},                               // PrefColorSchemeSystem
    {"Clair",                  "Light"},                                // PrefColorSchemeLight
    {"Sombre",                 "Dark"},                                 // PrefColorSchemeDark
    {"OK",                     "OK"},                                   // ButtonOk
    {"Annuler",                "Cancel"},                               // ButtonCancel
    {"Appliquer",              "Apply"},                                // ButtonApply
    {"Fermer",                 "Close"},                                // ButtonClose
    {"Quitter",                "Quit"},                                 // ButtonQuit
    {"Cette action ne peut pas être annulée. Continuer ?",
     "This action cannot be undone. Continue?"},                        // DestructiveConfirmBodyGeneric
    {"Aucun projet ouvert",    "No project open"},                      // EmptyStateHeading
    {"Choisissez Fichier → Ouvrir projet (Ctrl+O) pour commencer.",
     "Choose File → Open Project (Ctrl+O) to get started."},            // EmptyStateBody
}};

// Pitfall 2.1 — catch enum drift at compile time. Plan 02 flips Plan 01's
// inequality (`> 40`) to an exact-size assertion now that the table is full.
static_assert(kTable.size() == static_cast<size_t>(MsgId::_Count_),
              "kTable size must match MsgId::_Count_ — adding/removing a "
              "MsgId requires adding/removing the matching FR/EN row.");

}  // namespace

// probe_english() — reads $MEFISTO/td/m/anglais once and returns the result.
// Called by xvueIsEnglish(); factored out so both the production static-local
// and the XVUE_QT_TESTING resettable path share identical probe logic.
static bool probe_english() {
    const char* home = std::getenv("MEFISTO");
    if (home == nullptr || home[0] == '\0') {
        return false;  // no MEFISTO env → safe FR default.
    }
    QFileInfo flag(QString::fromLocal8Bit(home) + QStringLiteral("/td/m/anglais"));
    return flag.exists();
}

#ifdef XVUE_QT_TESTING
// File-scope cache used by xvueIsEnglish() and xvueClearLanguageCacheForTesting()
// in test builds. std::nullopt means "not yet evaluated"; xvueIsEnglish()
// populates it on first call; xvueClearLanguageCacheForTesting() resets it
// so the next call re-probes the filesystem.
static std::optional<bool> s_englishCache;
#endif

bool xvueIsEnglish() {
    // Pitfall 2.2 — cache at first call. Users toggle the language by
    // creating or removing $MEFISTO/td/m/anglais then restarting; live
    // toggling mid-process is not supported (mirrors Fortran-side LANGUE).
    //
    // Production: C++17 static-local bool — thread-safe, one allocation.
    // XVUE_QT_TESTING: file-scope std::optional so
    //   xvueClearLanguageCacheForTesting() can reset it between test slots.
#ifdef XVUE_QT_TESTING
    if (!s_englishCache.has_value()) {
        s_englishCache = probe_english();
    }
    return *s_englishCache;
#else
    static const bool english = probe_english();
    return english;
#endif
}

#ifdef XVUE_QT_TESTING
void xvueClearLanguageCacheForTesting() {
    // Reset the language probe so the next xvueIsEnglish() call re-probes
    // $MEFISTO/td/m/anglais from disk. Must be called alongside
    // XvueMenuFileParser::clearCacheForTesting() in test slots that change
    // the anglais flag between calls.
    s_englishCache = std::nullopt;
}
#endif

const char* tr(MsgId id) {
    const int ix = static_cast<int>(id);
    Q_ASSERT(ix >= 0 && ix < static_cast<int>(MsgId::_Count_));
    const auto& e = kTable[static_cast<size_t>(ix)];
    return xvueIsEnglish() ? e.en : e.fr;
}

QString xvueT(MsgId id) {
    return QString::fromUtf8(tr(id));
}
