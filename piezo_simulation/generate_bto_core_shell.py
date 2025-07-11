from ase.spacegroup import crystal
from ase import Atoms
from ase.io import write
import numpy as np

# Parameters
a = 4.0
rep = (1, 1, 1)
offset = 0.01

types = {'Ba': (1, 2), 'Ti': (3, 4), 'O': (5, 6)}
charges = {1: 2.0, 2: -2.0, 3: 4.0, 4: -4.0, 5: -2.6667, 6: 2.6667}
masses = {
    1: 137.33,    2: 0.00001,    # Ba
    3: 47.867,    4: 0.00001,    # Ti
    5: 15.999,    6: 0.00001     # O
}

# Create unit cell
bto = crystal(
    symbols=['Ba', 'Ti', 'O'],
    basis=[(0, 0, 0), (0.5, 0.5, 0.5),
           (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)],
    spacegroup=221,
    cellpar=[a, a, a, 90, 90, 90]
)

# Build supercell
bto *= rep

# Data storage
positions, atom_types, charges_list = [], [], []
atom_lines, bond_lines = [], []
bond_type_map = {
    1: 1,  # Ba core → bond type 1
    3: 2,  # Ti core → bond type 2
    5: 3   # O core  → bond type 3
}

dump_lines, xyz_lines = [], []
bond_id, atom_id, molecule_id = 1, 1, 1

for atom in bto:
    sym = atom.symbol
    core_type, shell_type = types[sym]
    charge_core = charges[core_type]
    charge_shell = charges[shell_type]
    core_pos = atom.position
    shell_pos = core_pos + np.array([offset, offset, offset])

    # Core
    atom_lines.append(f"{atom_id} {molecule_id} {core_type} {charge_core:.4f} {core_pos[0]:.6f} {core_pos[1]:.6f} {core_pos[2]:.6f}")
    dump_lines.append(f"{atom_id} {core_type} {core_pos[0]:.6f} {core_pos[1]:.6f} {core_pos[2]:.6f}")
    xyz_lines.append(f"{sym}_core {core_pos[0]:.6f} {core_pos[1]:.6f} {core_pos[2]:.6f}")
    atom_id += 1

    # Shell
    atom_lines.append(f"{atom_id} {molecule_id} {shell_type} {charge_shell:.4f} {shell_pos[0]:.6f} {shell_pos[1]:.6f} {shell_pos[2]:.6f}")
    dump_lines.append(f"{atom_id} {shell_type} {shell_pos[0]:.6f} {shell_pos[1]:.6f} {shell_pos[2]:.6f}")
    xyz_lines.append(f"{sym}_shell {shell_pos[0]:.6f} {shell_pos[1]:.6f} {shell_pos[2]:.6f}")

    #bond_lines.append(f"{bond_id} {core_type} {atom_id - 1} {atom_id}")
    core_type, shell_type = types[sym]
    bond_type = bond_type_map[core_type]
    bond_lines.append(f"{bond_id} {bond_type} {atom_id - 1} {atom_id}")

    bond_id += 1
    atom_id += 1
    molecule_id += 1

# Cell box
nx, ny, nz = rep
xhi, yhi, zhi = nx*a, ny*a, nz*a

# === Write .data file ===
with open("BTO_core_shell.data", "w") as f:
    f.write("LAMMPS data file for BaTiO3 core-shell model\n\n")
    f.write(f"{len(atom_lines)} atoms\n")
    f.write(f"{len(bond_lines)} bonds\n")
    f.write(f"{len(masses)} atom types\n")
    f.write("3 bond types\n\n")
    f.write(f"0.0 {xhi:.4f} xlo xhi\n0.0 {yhi:.4f} ylo yhi\n0.0 {zhi:.4f} zlo zhi\n\n")

    f.write("Masses\n\n")
    for i in sorted(masses):
        label = [k for k, v in types.items() if i in v][0]
        role = "core" if i % 2 == 1 else "shell"
        f.write(f"{i} {masses[i]:.5f}  # {label}_{role}\n")
    f.write("\nAtoms # full\n\n")
    for line in atom_lines:
        f.write(line + "\n")
    f.write("\nBonds\n\n")
    for line in bond_lines:
        f.write(line + "\n")

print("✅ LAMMPS .data file saved: BTO_core_shell.data")

# === Write .xyz file ===
with open("BTO_preview.xyz", "w") as f:
    f.write(f"{len(xyz_lines)}\nBaTiO3 core-shell preview\n")
    for line in xyz_lines:
        f.write(line + "\n")
print("✅ .xyz file for OVITO saved: BTO_preview.xyz")

# === Write .dump file ===
with open("BTO_dump.lammpstrj", "w") as f:
    f.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
    f.write(f"{len(dump_lines)}\n")
    f.write("ITEM: BOX BOUNDS pp pp pp\n")
    f.write(f"0.0 {xhi:.4f}\n0.0 {yhi:.4f}\n0.0 {zhi:.4f}\n")
    f.write("ITEM: ATOMS id type x y z\n")
    for line in dump_lines:
        f.write(line + "\n")
print("✅ LAMMPS .dump file saved: BTO_dump.lammpstrj")

# === Write LAMMPS groups ===
with open("group_definitions.in", "w") as f:
    f.write("# Group definitions for core-shell BaTiO3\n")
    for label, (core, shell) in types.items():
        f.write(f"group {label} type {core} {shell}\n")
print("✅ LAMMPS group definitions written: group_definitions.in")
