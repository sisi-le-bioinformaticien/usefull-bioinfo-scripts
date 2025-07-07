# Bioinformatics Utilities

A small collection of Python scripts for quick chemoinformatics tasks - simple tools you can run from the command line.

---

## Included Scripts

### 1. `smiles_to_ascii_art.py`

Parses a SMILES string, computes 2D coordinates and renders a rough ASCII‐art sketch of the molecule.

Example: Dabrafenib
```
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
```


2. compute_molecular_weight.py (not yet implemented)

Description:
Calculate molecular weight from a SMILES string.

Example:

$ python compute_molecular_weight.py --smiles "CC(=O)OC1=CC=CC=C1C(=O)O"
