# Bioinformatics Utilities

A small collection of Python scripts for quick chemoinformatics tasks—no bells and whistles, just simple tools you can run from the command line.

---

## Included Scripts

### 1. `ascii_molecule_renderer.py`

Parses a SMILES string, computes 2D coordinates with RDKit, and renders a rough ASCII‐art sketch of the molecule.

Example: Dabrafenib
$ python ascii_molecule_renderer.py --molecule Dabrafenib
       ==N==     --N
    C==     ==C--
    -          -
    C==      --N
       ===C--           ---C=              F--
          -          C--     ==               --C----C==   
          -          =         =C      O        -       ==C
        --C--     --C-         -        =      -         = 
     S--     --C--    --    ===C--    ---S-----C==    ---C 
      =       =         -C==      -N--    =       =C--     
       C-----N           -                 O       -       
C--   --                 F                         F        
   --C--
   --   --C
  C
