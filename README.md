## CurseMol - Molecular sketcher for the terminally committed

Have you ever been working on a terminal with programs that use SMILES for
input or output, and switching to your GUI sketcher is too annoying? (Not
to mention reading or writing the SMILES by hand.)

Enter CurseMol, the curses-based molecular editor, with keybindings inspired
by vi! The molecules are depicted using ASCII art (actually, not exactly; a few
non-ASCII Unicode characters are used, as well.)

You can read in a SMILES string or start from a blank canvas, and when you
exit, the final SMILES is printed to stdout, which can be useful for piping
directly into other programs.

[Screenshot](screenshot.png)

```
Controls:
  h, j, k, l       - Move cursor left, down, up, right
  H, J, K, L       - Move cursor faster (10 cells horizontal, 4 cells vertical)
  Space            - Snap cursor to nearest atom
  m                - Enter move mode (hjkl moves molecule, Esc to exit)
  s                - Enter a SMILES string to replace the current molecule
  S                - Toggle SMILES display
  i                - Insert/modify atom at cursor position
  a                - Append atoms from SMILES to atom or bond under cursor
                     (appending to a bond forms a ring by connecting the bond
                     atoms to the first and last atoms from the SMILES)
  c, n, o          - Insert/modify carbon/nitrogen/oxygen atom
  x                - Delete atom or bond
  X                - Area delete (select rectangle, Enter to delete, Esc to cancel)
  +, -             - Increase/decrease formal charge on atom
  <, >             - Zoom out/in
  1, 2, 3          - Add bond or change bond (order 1/2/3) between nearest atoms
  w, d             - Add/change to wedge or dash bond (press again to reverse)
  @                - Clear canvas (reset to blank slate)
  u, r             - Undo/redo
  Ctrl-L           - Clean up (regenerate coordinates)
  ?                - Show this help
  q                - Quit and print SMILES to stdout
```

### Requirements

- Python (tested with 3.11)
- RDKit (tested with 2025.09.6)
