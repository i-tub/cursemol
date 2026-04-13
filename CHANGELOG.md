# Changelog

## [4.1.1] - 2026-04-13

### Fixed

- Crash when resizing the window with an empty canvas.
- Some bonds not shown when there are atoms outside the canvas.

### Changed

- Major refactoring: split the monolithic package into modules.
- Other refactorings: removed unnecessary `y_offset`.
- Simplified redraw logic in `main_loop()`.
- Added a `Coords` class.
- Added unit tests (still a work in progress).
- Added type annotations.
- Added GitHub CI workflow.

## [4.1.0] - 2026-04-05

Initial release.
