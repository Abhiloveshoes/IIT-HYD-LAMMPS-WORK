log_filename = "log.lammps"
output_filename = "ke_data.txt"

with open(log_filename, "r") as f:
    lines = f.readlines()

data_lines = []
header_found = False

for line in lines:
    if not header_found:
        if line.strip().startswith("Step"):
            header_found = True
            continue
    else:
        if line.strip() == "":
            header_found = False
            continue
        tokens = line.strip().split()
        if len(tokens) < 4:
            continue
        step = tokens[0]
        kineng = tokens[3]  # 4th column (KinEng)
        data_lines.append(f"{step} {kineng}\n")

with open(output_filename, "w") as out:
    out.writelines(data_lines)

print(f"✅ Extracted {len(data_lines)} lines to '{output_filename}'")
