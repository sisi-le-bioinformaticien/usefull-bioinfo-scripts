from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
import random


EXAMPLES_SMILES = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Resveratrol": "C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O",
    "Dabrafenib": "CC(C)(C)C1=NC(=C(S1)C2=NC(=NC=C2)N)C3=C(C(=CC=C3)NS(=O)(=O)C4=C(C=CC=C4F)F)F",
    "Vitamin D": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]\\2[C@@]1(CCC/C2=C\C=C/3\C[C@H](CCC3=C)O)C",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Diazepam": "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3"
}

# This funtion start by creating empty grid, then start placing atoms inside, then it places bonds between atoms and prints the final result
def molecule_to_ascii(mol, width, height, positions):
    grid = [[' ' for _ in range(width)] for _ in range(height)]
    atom_coords = []
    for idx, (x, y) in enumerate(positions):
        gx = int(round(x))
        gy = int(round(y))
        atom_coords.append((gx, gy))
        atom = mol.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()[0]
        grid[height - gy - 1][gx] = symbol

    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        x0, y0 = atom_coords[start]
        x1, y1 = atom_coords[end]
        steps = max(abs(x1 - x0), abs(y1 - y0))
        rand_char = random.choice(['-','='])
        for i in range(steps + 1):
            xi = int(round(x0 + (x1 - x0) * i / steps))
            yi = int(round(y0 + (y1 - y0) * i / steps))
            if grid[height - yi - 1][xi] == ' ':
                if bond_type == rdchem.BondType.AROMATIC: 
                    grid[height - yi - 1][xi] = rand_char
                elif bond_type == rdchem.BondType.DOUBLE:  
                    grid[height - yi - 1][xi] = '='
                else:
                    grid[height - yi - 1][xi] = '-'
    for row in grid:
        print(''.join(row))


# This funtion initiates the grid depending on the molecule size, it also centers atom coordinates
def determine_grid_size(positions, scale_x = 4, scale_y = 1.5):
    x_pos, y_pos = [],[]
    for p in positions:
        x_pos.append(p.x) 
        y_pos.append(p.y) 
    min_x = min(x_pos)
    min_y = min(y_pos)
    centered_positions = [((x - min_x) * scale_x, (y - min_y) * scale_y) for x, y in zip(x_pos, y_pos)]
    max_x = max(x for x, _ in centered_positions)
    max_y = max(y for _, y in centered_positions)

    width = int(round(max_x)) + 1
    height = int(round(max_y)) + 1
    return width, height, centered_positions

# This funtion parses the molecule and returns all atoms positions
def molecule_setup(smile):
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        print("Invalid SMILE. Please enter a valid SMILE representation for your molecule")
        exit(0)

    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    atom_positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    return atom_positions, mol

if __name__ == "__main__":

    # ENTER YOUR CHOICE HERE (for example: "Aspirin")
    smile = EXAMPLES_SMILES["Dabrafenib"] 

    atom_positions, mol = molecule_setup(smile)
    width, height, centered_positions = determine_grid_size(atom_positions)
    molecule_to_ascii(mol, width, height, centered_positions)
