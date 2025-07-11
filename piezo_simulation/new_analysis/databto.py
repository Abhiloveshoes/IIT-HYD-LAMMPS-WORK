# Generates a core-shell LAMMPS data file for BaTiO3 (tetragonal unit cell)
# Author:

# Lattice parameters
a = 4.00  # Å
c = 4.04  # Å

# Atom types
atom_types = {
    "Ba_core": 1,
    "Ba_shell": 2,
    "Ti_core": 3,
    "Ti_shell": 4,
    "O_core": 5,
    "O_shell": 6
}

# Charges
charges = {
    "Ba_core": 1.9,
    "Ba_shell": -1.9,
    "Ti_core": 2.1,
    "Ti_shell": -2.1,
    "O_core": 0.7,
    "O_shell": -0.7
}

# Masses
masses = {
    1: 137.327,   # Ba_core
    2: 0.001,     # Ba_shell
    3: 47.867,    # Ti_core
    4: 0.001,     # Ti_shell
    5: 15.999,    # O_core
    6: 0.001      # O_shell
}

# Core-shell bond offset
offset = 0.01

# Atom list: (element_label, x, y, z)
atoms = [
    # Ba at (0,0,0)
    ("Ba_core", 0.0, 0.0, 0.0),
    ("Ba_shell", offset, offset, offset),

    # Ti at (0.5, 0.5, 0.5)
    ("Ti_core", a/2, a/2, c/2),
    ("Ti_shell", a/2 + offset, a/2 + offset, c/2 + offset),

    # O1 at (0.5, 0.5, 0)
    ("O_core", a/2, a/2, 0.0),
    ("O_shell", a/2 + offset, a/2 + offset, offset),

    # O2 at (0.5, 0, 0.5)
    ("O_core", a/2, 0.0, c/2),
    ("O_shell", a/2 + offset, offset, c/2 + offset),

    # O3 at (0, 0.5, 0.5)
    ("O_core", 0.0, a/2, c/2),
    ("O_shell", offset, a/2 + offset, c/2 + offset),
]

# Bond list: (bond-ID, bond-type, atom1-ID, atom2-ID)
bonds = []
for i in range(0, len(atoms), 2):
    bonds.append((i//2 + 1, 1, i + 1, i + 2))  # bond-ID, type=1, atom1, atom2

# Write to file
with open("batio3.data", "w") as f:
    f.write("LAMMPS data file for BaTiO3 with core-shell model\n\n")
    f.write(f"{len(atoms)} atoms\n")
    f.write(f"{len(bonds)} bonds\n\n")
    f.write("6 atom types\n")
    f.write("1 bond types\n\n")
    f.write(f"0.0 {a:.3f} xlo xhi\n")
    f.write(f"0.0 {a:.3f} ylo yhi\n")
    f.write(f"0.0 {c:.3f} zlo zhi\n\n")

    # Masses
    f.write("Masses\n\n")
    for type_id in sorted(masses):
        f.write(f"{type_id} {masses[type_id]:.4f}\n")
    f.write("\n")

    # Atoms
    f.write("Atoms # ID mol-ID type charge x y z\n\n")
    for idx, (label, x, y, z) in enumerate(atoms, start=1):
        mol_id = (idx + 1) // 2  # one molecule per core-shell pair
        type_id = atom_types[label]
        charge = charges[label]
        f.write(f"{idx} {mol_id} {type_id} {charge:.3f} {x:.4f} {y:.4f} {z:.4f}\n")
    f.write("\n")

    # Bonds
    f.write("Bonds\n\n")
    for bond_id, bond_type, a1, a2 in bonds:
        f.write(f"{bond_id} {bond_type} {a1} {a2}\n")
