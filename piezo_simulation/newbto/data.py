import numpy as np

a = 4.00  # Lattice constant
nx, ny, nz = 2, 2, 2  # Supercell size

core_shell_pairs = {
    "Ba": {"core_type": 1, "shell_type": 2, "core_charge": 1.5, "shell_charge": 0.5},
    "Ti": {"core_type": 3, "shell_type": 4, "core_charge": 2.0, "shell_charge": 2.0},
    "O":  {"core_type": 5, "shell_type": 6, "core_charge": -1.0, "shell_charge": -1.0}
}

# Basis positions (relative, perovskite)
basis = [
    ("Ba", [0.0, 0.0, 0.0]),
    ("Ti", [0.5, 0.5, 0.5]),
    ("O",  [0.5, 0.5, 0.0]),
    ("O",  [0.5, 0.0, 0.5]),
    ("O",  [0.0, 0.5, 0.5]),
]

atoms = []
bonds = []
atom_id = 1
bond_id = 1
mol_id = 1
dr = 0.01  # Core-shell displacement offset

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            for atom_name, rel_pos in basis:
                x = (i + rel_pos[0]) * a
                y = (j + rel_pos[1]) * a
                z = (k + rel_pos[2]) * a

                # Get types and charges
                cs = core_shell_pairs[atom_name]
                core_type = cs["core_type"]
                shell_type = cs["shell_type"]
                q_core = cs["core_charge"]
                q_shell = cs["shell_charge"]

                # Add core atom
                atoms.append((atom_id, mol_id, core_type, q_core, x, y, z))
                id_core = atom_id
                atom_id += 1

                # Add shell atom (displaced slightly along x for numerical stability)
                atoms.append((atom_id, mol_id, shell_type, q_shell, x + dr, y, z))
                id_shell = atom_id
                atom_id += 1

                # Add bond between core and shell
                bonds.append((bond_id, 1, id_core, id_shell))
                bond_id += 1
                mol_id += 1

# Box dimensions
xlo, xhi = 0, nx * a
ylo, yhi = 0, ny * a
zlo, zhi = 0, nz * a

with open("batio3.data", "w") as f:
    f.write("LAMMPS data file for BaTiO3 core-shell\n\n")
    f.write(f"{len(atoms)} atoms\n")
    f.write(f"{len(bonds)} bonds\n")
    f.write("6 atom types\n")
    f.write("1 bond types\n\n")
    f.write(f"{xlo:.3f} {xhi:.3f} xlo xhi\n")
    f.write(f"{ylo:.3f} {yhi:.3f} ylo yhi\n")
    f.write(f"{zlo:.3f} {zhi:.3f} zlo zhi\n\n")

    f.write("Masses\n\n")
    f.write("1 137.33\n2 0.001\n")  # Ba core/shell
    f.write("3 47.867\n4 0.001\n")  # Ti core/shell
    f.write("5 15.999\n6 0.001\n\n")  # O core/shell

    f.write("Atoms\n\n")
    for atom in atoms:
        f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]:.3f} {atom[4]:.3f} {atom[5]:.3f} {atom[6]:.3f}\n")

    f.write("\nBonds\n\n")
    for bond in bonds:
        f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")

print("Core-shell BaTiO3 data file written as 'batio3_core_shell.data'")

with open("ovito.data", "w") as f:
    f.write("LAMMPS data file for BaTiO3 core-shell (OVITO-compatible)\n\n")
    f.write(f"{len(atoms)} atoms\n\n")
    f.write("6 atom types\n")
    f.write("Atoms\n\n")
    for atom in atoms:
        f.write(f"{atom[0]} {atom[2]} {atom[4]:.4f} {atom[5]:.4f} {atom[6]:.4f}\n")

print("VITO-compatible data file saved as 'batio3_core_shell_ovito.data'")
