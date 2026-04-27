// xvue/qt/src/xvue_qt_i18n.h
// Phase 6.0 Plan 01 (scaffold): MsgId enum (locked enum width per UI-SPEC
// §Copywriting Contract) + xvueT()/tr() declarations. Plan 02 wires the real
// FR/EN lookup table and the $MEFISTO/td/m/anglais language probe.
#pragma once
#include <QString>

enum class MsgId : int {
    AppName,
    FileMenuTitle, FileOpen, FileOpenTip, FileSave, FileSaveTip,
    FileSaveAs, FileRecentSubmenu, FileRecentClear, FileRecentClearConfirmTitle,
    FileExport, FileQuit, FileQuitTip,
    ViewMenuTitle, ViewToolbar, ViewConsole, ViewZoomIn, ViewZoomOut, ViewFit, ViewPreferences,
    HelpMenuTitle, HelpDocumentation, HelpAbout,
    ConsoleTitle, StatusCoordFormat, ModalRefuse, ErrorMsgBoxTitle,
    QuitConfirmTitle, QuitConfirmBody,
    AboutTitle, AboutBody,
    OpenProjectDialogTitle, OpenProjectFilter,
    PreferencesTitle, PrefRecentCountLabel, PrefConsoleDefaultLabel,
    PrefColorSchemeLabel, PrefColorSchemeSystem, PrefColorSchemeLight, PrefColorSchemeDark,
    ButtonOk, ButtonCancel, ButtonApply, ButtonClose, ButtonQuit,
    DestructiveConfirmBodyGeneric,
    EmptyStateHeading, EmptyStateBody,
    _Count_
};

// True when $MEFISTO/td/m/anglais selects the English locale. Plan 02 wires
// the actual probe; Plan 01 stub returns false (FR default).
bool xvueIsEnglish();

// Test-only helper — resets the cached language probe so the next call to
// xvueIsEnglish() re-probes $MEFISTO/td/m/anglais from disk. Compiled only
// when XVUE_QT_TESTING is defined. Must be called alongside
// XvueMenuFileParser::clearCacheForTesting() in tests that toggle the
// anglais flag between test slots.
#ifdef XVUE_QT_TESTING
void xvueClearLanguageCacheForTesting();
#endif

// Returns the localized C-string for id (UTF-8 encoded). Plan 02 fills the
// 47-row lookup table; Plan 01 stub returns "".
const char* tr(MsgId id);

// Convenience wrapper returning a QString. Plan 01 stub forwards tr().
QString xvueT(MsgId id);
