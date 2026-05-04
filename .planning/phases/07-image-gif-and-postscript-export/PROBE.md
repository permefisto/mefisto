# Phase 7 — QImageWriter probe (EXPORT-01)

Probe run: 2026-05-04T16:35:52Z
Host: Linux 6.19.11+deb14-amd64 x86_64
Source: `xvue/qt/probes/qimagewriter_probe.cpp` (Plan 01)

## Raw probe output

```
qt_version=6.10.2
supported_write_formats=bmp cur icns ico jfif jpeg jpg pbm pgm png ppm tif tiff wbmp webp xbm xpm 
gif_write_supported=0
gif_animation_supported=0
```

## Interpretation

* If `gif_write_supported=0` then Phase 7 Plan 05 takes the
  ffmpeg fallback path (CONTEXT.md D-11).
* If `gif_write_supported=1` then Phase 7 Plan 05 takes the
  native QImageWriter multi-frame path (CONTEXT.md D-10).
* In either case Phase 7 forbids ImageMagick `convert` shell-outs
  inside `xvue/qt/` (CONTEXT.md D-16, EXPORT-06 grep gate).

## ffmpeg availability

```
ffmpeg version 8.1-3+b1 Copyright (c) 2000-2026 the FFmpeg developers
```
