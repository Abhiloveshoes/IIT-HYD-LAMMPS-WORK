input_file = "BTO_10x5x5.txt"
output_file = "BTO_10x5x5_fixed.txt"

inside_atoms = False

with open(input_file, "r") as fin, open(output_file, "w") as fout:
    for line in fin:
        stripped = line.strip()

        if stripped.lower().startswith("atoms"):
            inside_atoms = True
            fout.write(line)
            continue

        # Skip blank lines
        if stripped == "":
            fout.write(line)
            continue

        # Modify atoms section
        if inside_atoms and stripped[0].isdigit():
            # Split at comment if any
            parts = line.split("#")[0].split()
            comment = "#" + line.split("#")[1] if "#" in line else ""

            if len(parts) == 6:
                atom_id, atom_type, charge, x, y, z = parts
                new_line = f"{atom_id} 0 {atom_type} {charge} {x} {y} {z} {comment}".strip() + "\n"
                fout.write(new_line)
            else:
                fout.write(line)
        else:
            fout.write(line)
