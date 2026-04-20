## Test environments

- Local macOS Tahoe 26.4, R 4.5.3
- Local source build and installed-package smoke on the same platform

## R CMD check results

0 errors | 0 warnings | 2 NOTEs

## Notes

1. CRAN incoming feasibility reported:
   - `New submission`
   - `Version contains large components (0.0.0.9000)`
   The version will be bumped to a release form before any CRAN submission.
2. HTML validation was skipped locally because the installed `tidy` binary is not a
   recent enough HTML Tidy release.

## Additional comments

- `logistf` is listed in `Suggests` only. The package includes a fallback
  initialization path when `logistf` is unavailable, so checks were run with
  `_R_CHECK_FORCE_SUGGESTS_=false`.
- The package includes adapted GPL-compatible code from the archived `binarySimCLF`
  and `geefirthr` codebases. Attribution and provenance are recorded in
  `LICENSE.note` and `inst/CREDITS.md`.
- This file is a draft for the first submission cycle, not a final CRAN cover note.
