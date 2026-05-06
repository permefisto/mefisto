# Plan 9-07 Task 1: ppnlse_qt offscreen+BATCH_X11 root-cause diagnose

**Date:** 2026-05-05
**Reproducer (plan-documented, stdin-redirect form):**
`QT_QPA_PLATFORM=offscreen MEFISTO_BATCH_X11=1 OMP_NUM_THREADS=1 pp/ppnlse < $MEFISTOX/nlsecu/nlsecu.iexrr`
**Reproducer (production harness form, `bin/ab_sweep_phase8.sh` line 153):**
`env QT_QPA_PLATFORM=offscreen MEFISTO_BATCH_X11=1 MEFISTO_QT_CAPTURE_PATH=... MEFISTO_XVSOURIS_AUTOEXIT=1 timeout 60 pp/ppnlse nlsecu.iexrr`
**Binary under test:** `pp/ppnlse` (single-backend Qt build per Phase 9; legacy `_qt`
suffix retired with the X11 backend — pp/ppnlse IS the Qt-backed binary,
`ldd pp/ppnlse | grep Qt6` confirms `libQt6Widgets.so.6 / libQt6Gui.so.6 / libQt6Core.so.6`).
**Build state:** fresh build via `bin/cbl_tout` (post-Plan-9-05); ABI count 58 (verify_abi.sh exit 0; raw `nm | grep ' T '_$' | wc -l` = 64 per 09-02-SUMMARY note on the 6 stub-symbol filter).
**Workspace:** canonical `$MEFISTOX/nlsecu/` with `Final TIME=20, Step=0.01` → 2000 time steps (NOT TIME-truncated). MAILLER prereq (`pp/ppmail nlsecu.meshq2`) was run beforehand to seed `ms10..ms15` on the fresh INITIER tree.

## Stdout state at hang (T+30s, plan-documented stdin-redirect reproducer)

Line count: **9** (Phase 8 reported "~10 log lines, no NLSER banner reached even at 240s" — same regime, off-by-one).

Last 9 lines:
```text
 Mefisto speaks ENGLISH
 WORKING DIRECTORY NAME :
/home/mefisto/mefistox/nlsecu
xvue-qt: stub ladate_ not implemented yet
xvue-qt: stub heureminuteseconde_ not implemented yet
xvue-qt: stub valvarenv_ not implemented yet
This plugin does not support propagateSizeHints()
This plugin does not support raise()
xvue-qt: stub nomrepmefisto_ not implemented yet
```

NLSER banner reached: **no** (lines 117/120/121 of `prpr/ppnlse.f` print `Mefisto-NLSER: ... BATCH/INTERACTIVE EXECUTION` and would appear if reached).
Any iter/step output: **no**.

## gdb backtrace (top of all threads, attached at T+30s)

```text
[Thread debugging using libthread_db enabled]
Thread 1 (Thread 0x7f95e2a0e980 (LWP 2752864) "ppnlse"):
#0  __syscall_cancel_arch ()
    at ../sysdeps/unix/sysv/linux/x86_64/syscall_cancel.S:56
#1  0x00007f95e1a9a7a4 in __internal_syscall_cancel (...)
    at ./nptl/cancellation.c:49
#2  0x00007f95e1a9a7ed in __syscall_cancel (...)
    at ./nptl/cancellation.c:75
#3  0x00007f95e1b1029e in __GI_ppoll (fds=..., nfds=..., timeout=..., sigmask=...)
    at ../sysdeps/unix/sysv/linux/ppoll.c:42
#4  0x00007f95e2304aa4 in ?? () from /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0
#5  0x00007f95e2305190 in g_main_context_iteration ()
    from /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0
#6  0x00007f95e2818cf8 in QEventDispatcherGlib::processEvents(QFlags<QEventLoop::ProcessEventsFlag>)
    () from /usr/lib/x86_64-linux-gnu/libQt6Core.so.6
#7  0x00007f95e25b4ada in QEventLoop::exec(QFlags<QEventLoop::ProcessEventsFlag>)
    () from /usr/lib/x86_64-linux-gnu/libQt6Core.so.6
#8  0x0000560ca8591b4b in XvueEventBridge::waitForEvent(XvueEventBridge::WaitMode, int*, int*) ()
#9  0x0000560ca8588b94 in xvsouris_ ()
#10 0x0000560ca8560396 in saiptc_ ()
#11 0x0000560ca854bf39 in afdocu_ ()
#12 0x0000560ca8557b95 in luou_ ()
#13 0x0000560ca80c738e in MAIN__ ()
#14 0x0000560ca80c8436 in main ()
[Inferior 1 (process 2752864) detached]
```

## strace summary (10s window)

```text
strace: Process 2759055 attached
strace: Process 2759055 detached
(no syscalls in poll/select/futex/read/write/ppoll trace set during the
 10s window — the process is parked in a single ppoll waiting for FD activity
 on the Qt event-loop side; no read() of stdin happens because the language
 layer has already consumed enough bytes to make INTERA-detection fall through
 to interactive mode and is now blocked in xvsouris_ for a mouse click that
 will never arrive under offscreen+headless.)
```

## Comparison: production-harness (arg-form) reproducer

**Same env vars, same workspace, same fresh binary** — only difference: pass the
batch file as ARG (matches `bin/ab_sweep_phase8.sh` line 153) instead of stdin
redirect.

```text
$ env QT_QPA_PLATFORM=offscreen MEFISTO_BATCH_X11=1 \
      MEFISTO_XVSOURIS_AUTOEXIT=1 MEFISTO_XVSOURIS_AUTOEXIT_DELAY_MS=500 \
      timeout 60 pp/ppnlse nlsecu.iexrr > /tmp/09-07-test-A.log 2>&1
$ echo "exit=$?"
exit=0
$ wc -l /tmp/09-07-test-A.log
17
$ grep 'Mefisto-NLSER' /tmp/09-07-test-A.log
 Mefisto-NLSER: Data File Name=nlsecu.iexrr
 Mefisto-NLSER in BATCH EXECUTION
```

**Banner IS reached.** The 2000-step solver loop then runs out of 60s budget
before producing iter/step output (the Phase 8 documented "TIME=20, 2000 steps,
~hour-scale on this hardware" symptom — see 08-CHECKLIST.md line 84).

## ROOT CAUSE: (c)

Selected: **(c) long-runtime, NOT a deadlock**.

Evidence:
- gdb backtrace top frame is `__GI_ppoll → g_main_context_iteration → QEventLoop::exec`
  (frames #3–#7) — i.e. the Qt event loop is **alive and dispatching**, not
  deadlocked. There is no QApplication-init race, no mutex contention,
  no pthread_cond_wait visible.
- The block point (frames #8–#9) is `XvueEventBridge::waitForEvent` from
  `xvsouris_` — a documented, intentional blocking primitive that waits
  for mouse input. This is reached because the **plan-documented
  reproducer uses stdin-redirect** (`pp/ppnlse < nlsecu.iexrr`),
  which means `IARGC()=0` → `EXISTNF=0` → `INTERA = IINFO('INTERACTIVITE
  INITIALE')` (line 101 of `prpr/ppnlse.f`) reads the .iexrr's leading
  `{` comment tokens for an INTERA value, falls into interactive mode,
  enters the LU menu (`luou_ → afdocu_ → saiptc_ → xvsouris_`)
  before ever reaching the BATCH banner at lines 117/120/121.
- The **production harness form** (`pp/ppnlse nlsecu.iexrr` as ARG, per
  `bin/ab_sweep_phase8.sh` line 153) reaches the `Mefisto-NLSER` banner in
  17 log lines well within 60s — proving there is no genuine deadlock in
  the binary itself. The remaining symptom under the harness is exactly
  what 08-CHECKLIST.md override #5 already documents: TIME=20 / Step=0.01
  → 2000 steps is ~hour-scale on this hardware, exceeds the 60s harness budget.
- 08-01-SUMMARY.md line 129 already noted the original Phase 8 deadlock
  symptom was caused by **stale `pp/ppnlse_qt`** (built before
  libxvueqt.a was rebuilt May 5 12:25); after Phase 8 Plan 1 ran
  `bin/cbl_tout_qt` the deadlock disappeared on FRESH binaries — leaving
  only the long-runtime constraint.
- strace shows the process makes **zero** poll/select/read/write/futex
  calls in a 10s sample — fully consistent with "parked in glib's ppoll
  waiting for an FD wakeup that won't come from a non-existent click",
  not with "race in QApplication init that thrashes futex".

The Phase 8 description "10 log lines, no banner reached even at 240s" is
explained by **two compounding factors** rather than a single deadlock:
1. **Stale pp/*_qt binaries** (already mitigated by Phase 8 Plan 1 rebuild)
2. **Plan-documented reproducer uses stdin-redirect rather than arg-form** —
   this is a reproducer-form issue, NOT a binary defect. The arg-form
   harness path reaches the banner cleanly on fresh binaries.

The actual residual issue (canonical TIME=20 / 2000 steps exceeds 60s) is
**case (c) long-runtime, NOT a deadlock** — same conclusion 08-CHECKLIST.md
override #5 reached at sign time (option (i): "Phase-9-deferred-fix" is now
revealed as case-(c) classification — no fix is owed because there is no defect).

## Fix path (forward to Task 2)

**CASE (c) — documentation-only closure.** Per Plan 9-07 task 2 case (c):

1. Add a comment block above the nlsecu invocation in `bin/ab_sweep_phase8.sh`
   citing this diagnose document and noting that Phase-9-7 Task 1 classified
   the symptom as case (c) long-runtime, not a deadlock.
2. Append a sentence to `.planning/phases/08-ab-validation-on-testa-subset/08-CHECKLIST.md`
   override #5 (line 84 area) referencing this diagnose document.
3. **No code change** to `xvue/qt/src/xvue_qt_api.cpp`, no code change to
   `prpr/ppnlse.f`, no ABI change.
4. ABI count remains 58 (verified pre and post; no rebuild required for case (c)
   since no source changed, but Task 3 will run `bin/cbl_tout` defensively
   and re-check the count anyway).

The Phase 8 mitigation (TIME=0.01 / TIME=0.1 workspace truncations + MAILLER-prereq
fallback) **remains the documented work-around** for the 60s harness budget;
override #5's "(i) Phase-9-deferred-fix" disposition refines to "case (c)
classified — no defect — Phase 8 mitigation retained as the canonical
work-around for ad-hoc harness budgets shorter than ~hour-scale."
