import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import os
import re

# -------------------- PARAMETERS --------------------
timestep_ps = 0.001
data_dirs = ['efield_-2.0Vnm', 'efield_-1.0Vnm', 'efield_0Vnm', 'efield_1.0Vnm', 'efield_2.0Vnm']
field_values, stress_values, frequencies, q_factors = [], [], [], []

# -------------------- MAIN LOOP --------------------
for folder in data_dirs:
    print(f"\n📂 Processing: {folder}")
    match = re.search(r'efield_([-0-9.]+)Vnm', folder)
    if not match:
        continue
    efield = float(match.group(1))

    # --- Stress from log.lammps ---
    stress_bar = None
    with open(os.path.join(folder, 'log.lammps')) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if 'Step' in line and 'v_stress_z' in line:
            header = line.strip().split()
            stress_col = header.index('v_stress_z')
            break
    else:
        continue

    # Read last value from stress column
    for line in reversed(lines):
        if re.match(r'^\s*\d+', line):
            try:
                stress_bar = float(line.strip().split()[stress_col])
                break
            except:
                continue
    if stress_bar is None:
        continue

    stress_gpa = stress_bar * 0.00001  # bar → GPa (Corrected conversion!)
    field_values.append(efield)
    stress_values.append(stress_gpa)

    # --- Kinetic Energy and FFT ---
    ke_path = os.path.join(folder, 'ke_data.txt')
    data = np.loadtxt(ke_path)
    steps, ke = data[:, 0], data[:, 1]
    time_ps = steps * timestep_ps
    ke_detrended = ke - np.mean(ke)

    freq = np.fft.fftfreq(len(time_ps), d=timestep_ps)
    fft_vals = np.fft.fft(ke_detrended)
    freq_pos = freq[freq > 0]
    amp = np.abs(fft_vals[freq > 0])
    amp_norm = amp / np.max(amp) * 10  # Match amplitude scaling from paper
    peak_freq = freq_pos[np.argmax(amp)] / 2  # THz
    freq_ghz = peak_freq   # GHz

    frequencies.append(freq_ghz)

    # --- Individual FFT Plot ---
    plt.figure()
    plt.plot(freq_pos , amp_norm, label=f'{efield} V/nm')
    plt.axvline(peak_freq , color='r', linestyle='--', label=f'Peak ≈ {freq_ghz:.1f} GHz')
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Amplitude (a.u.)")
    plt.xticks([0, 100, 200, 300])
    plt.yticks([0, 2, 4, 6, 8, 10])
    plt.xlim(0, 300)
    plt.ylim(0, 10.5)
    plt.grid(True)
    plt.legend()
    plt.title(f"FFT Spectrum for {efield} V/nm")
    plt.tight_layout()
    plt.savefig(f"fft_field_{efield:+.1f}Vnm.png", dpi=300)
    plt.close()

    # --- Q-Factor via peak decay ---
    peaks, _ = find_peaks(ke, distance=20)
    t_peaks, ke_peaks = time_ps[peaks], ke[peaks]
    def decay(t, A, Q, f): return A * np.exp(-np.pi * f * t / Q)
    try:
        popt, _ = curve_fit(decay, t_peaks, ke_peaks, p0=[max(ke_peaks), 1000, peak_freq])
        q_factors.append(popt[1])
    except:
        q_factors.append(0)

# -------------------- CONVERT ARRAYS --------------------
fields = np.array(field_values)
stresses = np.array(stress_values)
freqs = np.array(frequencies)
qs = np.array(q_factors)

# -------------------- PIEZOELECTRIC COEFFICIENT e33 --------------------
coeff = np.polyfit(fields, stresses, 1)
slope = coeff[0]
e33_sim = -slope  # σ = -e33 * E
print(f"\nComputed piezoelectric coefficient from MD data: e33 ≈ {e33_sim:.2f} C/m²")

# -------------------- PLOTS --------------------

# 1. Frequency vs Electric Field
plt.figure()
plt.plot(fields, freqs, 'o-', color='green', label='MD simulations')
plt.xlabel(r"$E_3$ (V/nm)")
plt.ylabel("f (GHz)")
plt.title("Resonant Frequency vs Electric Field")
plt.grid(True)
plt.tight_layout()
plt.savefig("freq_vs_field.png", dpi=300)

# 2. Q-Factor vs Electric Field
plt.figure()
plt.plot(fields, qs, 'o-', color='darkred', label='MD simulations')
q_fit = np.polyfit(fields, qs, 1)
plt.plot(fields, np.poly1d(q_fit)(fields), 'r-', label='Linear fit')
plt.xlabel(r"$E_3$ (V/nm)")
plt.ylabel("Quality Factor")
plt.title("Quality Factor vs Electric Field")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("q_vs_field.png", dpi=300)

# 3. Stress vs Electric Field
plt.figure()
plt.plot(fields, stresses, 'o-', label=f"MD (e33 ≈ {e33_sim:.2f} C/m²)")
plt.plot(fields, -1.65 * fields, 'r--', label="Paper Fit (e33 = -1.65 C/m²)")
plt.xlabel(r"$E_3$ (V/nm)")
plt.ylabel(r"Stress $\sigma_{33}$ (GPa)")
plt.title("Stress vs Electric Field")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("stress_vs_field.png", dpi=300)

# 4. f² vs Stress
f_squared = freqs ** 2
plt.figure()
plt.plot(stresses, f_squared, 'o-', label='MD Data')
fit = np.polyfit(stresses, f_squared, 1)
plt.plot(stresses, np.poly1d(fit)(stresses), 'r--', label='Linear Fit')
plt.xlabel(r"Stress $\sigma_{33}$ (GPa)")
plt.ylabel(r"$f^2$ (GHz²)")
plt.title(r"$f^2$ vs Stress (Beam Theory Validation)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("f_squared_vs_stress.png", dpi=300)

# 5. Combined FFT Plot
colors = {"-2.0": 'red', "-1.0": 'orange', "0.0": 'blue', "1.0": 'purple', "2.0": 'green'}
plt.figure(figsize=(6, 4))
for folder in data_dirs:
    match = re.search(r'efield_([-0-9.]+)Vnm', folder)
    if not match:
        continue
    field = match.group(1)
    ke_path = os.path.join(folder, 'ke_data.txt')
    try:
        data = np.loadtxt(ke_path)
        steps, ke = data[:, 0], data[:, 1]
        ke_detrended = ke - np.mean(ke)
        freq = np.fft.fftfreq(len(steps), d=timestep_ps)
        fft_vals = np.fft.fft(ke_detrended)
        freq_pos = freq[freq > 0]
        amp = np.abs(fft_vals[freq > 0])
        amp_norm = amp / np.max(amp) * 10
        plt.plot(freq_pos * 1000, amp_norm, label=f"{field} V/nm", color=colors.get(field, 'black'))
    except:
        continue

plt.xlim(0, 300)
plt.ylim(0, 10.5)
plt.xticks([0, 100, 200, 300])
plt.yticks([0, 2, 4, 6, 8, 10])
plt.xlabel("Frequency (GHz)")
plt.ylabel("Amplitude (a.u.)")
plt.title("FFT Spectrum vs Electric Field")
plt.legend(title="E-field")
plt.grid(True)
plt.tight_layout()
plt.savefig("fft_combined_all_fields.png", dpi=300)

print("\n✅ All graphs saved successfully.")
