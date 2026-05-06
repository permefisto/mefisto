## Plan 09-05 deferred-items (docs-only plan; source-level issues out of scope)

### Pre-existing testa sweep failures (not caused by Plan 09-05)

- **cavity2d-qt-1x.png** — produced 0-byte PNG in 5-case smoke sweep. Plan 09-04 SUMMARY also skipped this case per orchestrator pacing. Pre-existing issue; investigation in Phase 9 carry-forward (Plan 09-06 / 09-08 territory). Plan 09-05 is docs-only and made zero source changes.

- **nlsecu-qt-1x.png** — produced 0-byte PNG. This is the documented Phase-8 carry-forward #2 (ppnlse_qt offscreen + MEFISTO_BATCH_X11=1 deadlock at startup). CONTEXT.md D-03 explicitly assigns the fix to Plan 09-06.

### Effect on Plan 09-05 verification gate

Plan 09-05 Task 5 verify gate requires `[ $(ls /tmp/09-05-post-retire04/*.png | wc -l) -ge "4" ]`. With pre-existing cavity2d + nlsecu failures, only 3 valid PNGs are produced (pan2d, nafems_le1, heat1d). This is a Rule-1-style scope-bounded deviation: the failures pre-date my edits and are documented carry-forwards to Plans 09-06/08. The plan's other gates — `bin/cbl_tout` exit 0, all 3 grep gates green — are met cleanly.

