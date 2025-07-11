import numpy as np

# Parameters
a = 4.00  # Lattice constant (Å)
nx, ny, nz = 2, 2, 2  # Supercell size

# Atom type mapping
# 1: Ba_core, 2: Ba_shell
# 3: Ti_core, 4: Ti_shell
# 5: O_core, 6: O_shell

def write_data(filename):
    atom_id = 1
    bond_id = 1
    atoms = []
    bonds = []

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x0, y0, z0 = i * a, j * a, k * a

                # Ba
                atoms.append((atom_id, 1, 1.5, x0, y0, z0))
                atom_id += 1
                atoms.append((atom_id, 2, 0.5, x0 + 0.01, y0, z0))
                bonds.append((bond_id, 1, atom_id - 1, atom_id))
                atom_id += 1
                bond_id += 1

                # Ti
                atoms.append((atom_id, 3, 2.0, x0 + a/2, y0 + a/2, z0 + a/2))
                atom_id += 1
                atoms.append((atom_id, 4, 2.0, x0 + a/2 + 0.01, y0 + a/2, z0 + a/2))
                bonds.append((bond_id, 1, atom_id - 1, atom_id))
                atom_id += 1
                bond_id += 1

                # O (3 directions)
                for dx, dy, dz in [(a/2, 0, 0), (0, a/2, 0), (0, 0, a/2)]:
                    atoms.append((atom_id, 5, -1.0, x0 + dx, y0 + dy, z0 + dz))
                    atom_id += 1
                    atoms.append((atom_id, 6, -1.0, x0 + dx + 0.01, y0 + dy, z0 + dz))
                    bonds.append((bond_id, 1, atom_id - 1, atom_id))
                    atom_id += 1
                    bond_id += 1

    with open(filename, 'w') as f:
        f.write("LAMMPS data file for BaTiO3 core-shell 4x4x4\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(bonds)} bonds\n\n")
        f.write("6 atom types\n1 bond types\n\n")
        f.write(f"0.0 {nx*a:.3f} xlo xhi\n")
        f.write(f"0.0 {ny*a:.3f} ylo yhi\n")
        f.write(f"0.0 {nz*a:.3f} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write("1 137.327\n2 0.0001\n3 47.867\n4 0.0001\n5 15.999\n6 0.0001\n\n")

        f.write("Atoms\n\n")
        for idx, (aid, atype, q, x, y, z) in enumerate(atoms):
            f.write(f"{aid} {atype} {atype} {q:.4f} {x:.4f} {y:.4f} {z:.4f}\n")

        f.write("\nBonds\n\n")
        for bid, btype, a1, a2 in bonds:
            f.write(f"{bid} {btype} {a1} {a2}\n")

    print(f" Data file written: {filename}")

# Run it
write_data("batio3_core_shell_4x4x4.data")
