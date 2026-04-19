// xvue/qt/src/xvue_qt_i18n.cpp
// Phase 6.0 Plan 01 (scaffold): stub bodies. Plan 02 fills the FR/EN table
// and the $MEFISTO/td/m/anglais runtime probe.
#include "xvue_qt_i18n.h"

// Plan 02 will replace this with the exact-size assertion once the table
// is filled. Catches enum truncation at compile time today.
static_assert(static_cast<int>(MsgId::_Count_) > 40,
              "MsgId enum must declare all UI-SPEC entries");

bool xvueIsEnglish() {
    // Plan 02 wires the real $MEFISTO/td/m/anglais detection; FR default.
    return false;
}

const char* tr(MsgId id) {
    // Plan 02 fills the lookup; scaffold returns empty string for every id.
    (void)id;
    return "";
}

QString xvueT(MsgId id) {
    return QString::fromUtf8(tr(id));
}
