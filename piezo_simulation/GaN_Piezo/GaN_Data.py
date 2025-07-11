import numpy as np

# ----- GaN Wurtzite Parameters -----
a = 3.189  # Angstrom
c = 5.185  # Angstrom
u = 0.377  # internal parameter for N

# ----- Supercell Size -----
nx, ny, nz = 5, 5, 5  # you can change this

# ----- Atom Types -----
# 1 = Ga core, 2 = Ga shell, 3 = N core, 4 = N shell

# ----- Basis Positions -----
basis = [
    ("Ga", 0.0000, 0.0000, 0.0000),
    ("Ga", 1/3, 2/3, 0.5),
    ("N", 0.0000, 0.0000, u),
    ("N", 1/3, 2/3, 0.5 + u)
]

# ----- Lattice Vectors (orthorhombic for LAMMPS box) -----
lx = a * nx
ly = a * np.sqrt(3) * ny
lz = c * nz

# ----- Output Containers -----
atoms = []
bonds = []
atom_id = 1
bond_id = 1

def frac_to_cart(xf, yf, zf, i, j, k):
    x = (xf + i) * a
    y = (yf + j) * a * np.sqrt(3)
    z = (zf + k) * c
    return x, y, z

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            for element, xf, yf, zf in basis:
                x, y, z = frac_to_cart(xf, yf, zf, i, j, k)

                if element == "Ga":
                    core_type = 1
                    shell_type = 2
                    core_charge = 3.0
                else:
                    core_type = 3
                    shell_type = 4
                    core_charge = -3.0

                shell_offset = 0.02  # slight offset for shell
                atoms.append((atom_id, core_type, core_charge, x, y, z))  # core
                core_id = atom_id
                atom_id += 1

                atoms.append((atom_id, shell_type, 0.0, x + shell_offset, y, z))  # shell
                shell_id = atom_id
                atom_id += 1

                bonds.append((bond_id, 1, core_id, shell_id))
                bond_id += 1

# ----- Write .data File -----
with open("gan_core_shell.data", "w") as f:
    f.write("LAMMPS data file for GaN core-shell model\n\n")
    f.write(f"{len(atoms)} atoms\n")
    f.write(f"{len(bonds)} bonds\n")
    f.write("4 atom types\n1 bond types\n\n")
    f.write(f"0.0 {lx:.4f} xlo xhi\n")
    f.write(f"0.0 {ly:.4f} ylo yhi\n")
    f.write(f"0.0 {lz:.4f} zlo zhi\n\n")

    f.write("Masses\n\n")
    f.write("1 69.723\n")
    f.write("2 0.0\n")
    f.write("3 14.007\n")
    f.write("4 0.0\n\n")

    f.write("Atoms\n\n")
    for atom in atoms:
        f.write(f"{atom[0]} {atom[1]} {atom[2]:.4f} {atom[3]:.4f} {atom[4]:.4f} {atom[5]:.4f}\n")

    f.write("\nBonds\n\n")
    for bond in bonds:
        f.write(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}\n")

print(f"[✔] Generated GaN core-shell data file for {nx}x{ny}x{nz} unit cells.")
