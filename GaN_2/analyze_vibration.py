# ----------- Final Corrected Python Script (analyze_piezo.py) -----------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import os
import re

# --- PARAMETERS ---
timestep_ps = 0.001  # LAMMPS timestep is 0.001 ps (1 fs)
data_dirs = ['efield_-2.0Vnm', 'efield_-1.0Vnm', 'efield_0Vnm', 'efield_1.0Vnm', 'efield_2.0Vnm']

field_values = []
stress_values = []
frequencies = []
q_factors = []

# --- LOOP THROUGH SIMULATION FOLDERS ---
for folder in data_dirs:
    print(f"Processing: {folder}")

    # Extract electric field from folder name
    match = re.search(r'efield_([-0-9.]+)Vnm', folder)
    if not match:
        print(f"Skipping folder {folder}, cannot parse electric field.")
        continue
    efield = float(match.group(1))

    # --- Parse stress from log.lammps ---
    log_path = os.path.join(folder, 'log.lammps')
    stress_bar = None

    with open(log_path, 'r') as log_file:
        for line in log_file:
            if 'v_stress_z' in line:
                break  # Found thermo column header

        for line in log_file:
            if not line.strip() or not re.match(r'^\s*\d+', line):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                stress_bar = float(parts[2])  # v_stress_z is the 3rd column
            except ValueError:
                continue

    if stress_bar is None:
        print(f"Warning: No valid stress found in {folder}, skipping this field.")
        continue

    stress_gpa = stress_bar * 0.00001 # convert from bar to GPa

    # Append electric field and stress only if stress is valid
    field_values.append(efield)
    stress_values.append(stress_gpa)

    # --- Load kinetic energy data ---
    ke_data_path = os.path.join(folder, 'ke_data.txt')
    data = np.loadtxt(ke_data_path)
    steps = data[:, 0]
    ke = data[:, 1]
    time_ps = steps * timestep_ps

    # --- FFT to Find Fundamental Frequency ---
    ke_detrended = ke - np.mean(ke)
    freq = np.fft.fftfreq(len(time_ps), d=timestep_ps)  # in THz
    fft_vals = np.fft.fft(ke_detrended)

    freq_pos = freq[freq > 0]
    amp = np.abs(fft_vals[freq > 0])

    peak_freq_thz = freq_pos[np.argmax(amp)] / 2  # divide by 2 due to KE oscillations at 2f
    freq_ghz = peak_freq_thz * 1000
    frequencies.append(freq_ghz)

    # --- Extract Envelope (Peaks) for Q-Factor Fit ---
    peak_indices, _ = find_peaks(ke, distance=20)
    t_peaks = time_ps[peak_indices]
    ke_peaks = ke[peak_indices]

    # Fit exponential decay to peaks
    def decay(t, A, Q, f):
        return A * np.exp(-np.pi * f * t / Q)

    try:
        p0 = [np.max(ke_peaks), 1000.0, peak_freq_thz]
        popt, _ = curve_fit(decay, t_peaks, ke_peaks, p0=p0)
        q_factors.append(popt[1])
    except Exception as e:
        print(f"Warning: Q-factor fitting failed for {folder}: {e}")
        q_factors.append(0)

# --- Check Data Lengths ---
print(f"Fields: {field_values}")
print(f"Stress: {stress_values}")

# --- PLOTTING ---

# Resonant Frequency vs Electric Field
plt.figure()
plt.plot(field_values, frequencies, 'o-', label='Resonant Frequency')
plt.xlabel('Electric Field (V/nm)')
plt.ylabel('Resonant Frequency (GHz)')
plt.title('Frequency vs Electric Field')
plt.grid(True)
plt.tight_layout()
plt.savefig("Resonant Frequency vs Electric Field.png", dpi=300)
plt.show()
# Q-Factor vs Electric Field
plt.figure()
plt.plot(field_values, q_factors, 'o-', label='Q Factor')
plt.xlabel('Electric Field (V/nm)')
plt.ylabel('Q-Factor')
plt.title('Q-Factor vs Electric Field')
plt.grid(True)
plt.tight_layout()

plt.savefig("# Q-Factor vs Electric Field.png", dpi=300)
plt.show()
# Stress vs Electric Field → Calculate e33
# --- Stress vs Electric Field with paper's e33 = 1.65 C/m² ---
field_arr = np.array(field_values)
stress_arr = np.array(stress_values)
coeff = np.polyfit(field_arr, stress_arr, 1)
e33_sim = -coeff[0]

plt.figure()
plt.plot(field_arr, stress_arr, 'o-', color='green', label=f'MD Stress (e33 ≈ {e33_sim:.2f} C/m²)')

# Fit your MD data to see what e33 you're getting (optional)


# Paper's reference line using e33 = 1.65 C/m² → slope = -1.65 GPa / (V/nm)
stress_ref = -1.65 * field_arr  # Linear relation from paper
plt.plot(field_arr, stress_ref, 'r--', label='Paper Fit (e33 = -1.65 C/m²)')

plt.xlabel('Electric Field (V/nm)')
plt.ylabel('Stress σ₃₃ (GPa)')
plt.title('Stress vs Electric Field')
plt.legend()
plt.grid(True)
plt.tight_layout()


plt.savefig("Stress vs Electric Field → Calculate e33.png", dpi=300)
print(f"\n🔍 Your simulation e33 ≈ {e33_sim:.2f} C/m²")
print(f"📊 Paper reference line drawn with e33 = 1.65 C/m²")
plt.show()
# --- FFT spectrum plot ---
plt.figure()
plt.plot(freq_pos, amp, label=f'Field = {efield} V/nm')
plt.axvline(peak_freq_thz, color='r', linestyle='--', label=f'Peak = {peak_freq_thz*1000:.1f} GHz')
plt.xlabel('Frequency (THz)')
plt.ylabel('Amplitude')
plt.title(f'FFT - Field {efield} V/nm')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f'fft_field_{efield:+.1f}Vnm.png', dpi=300)
plt.close()

# Configuration
timestep_ps = 0.001
data_dirs = ['efield_-2.0Vnm', 'efield_-1.0Vnm', 'efield_0Vnm', 'efield_1.0Vnm', 'efield_2.0Vnm']
colors = {
    "-2.0": 'red',
    "-1.0": 'orange',
    "0.0": 'blue',
    "1.0": 'purple',
    "2.0": 'green'
}

plt.figure(figsize=(5, 4))

for folder in data_dirs:
    match = re.search(r'efield_([-0-9.]+)Vnm', folder)
    if not match:
        continue
    field = match.group(1)

    ke_path = os.path.join(folder, 'ke_data.txt')
    try:
        data = np.loadtxt(ke_path)
    except:
        continue

    steps = data[:, 0]
    ke = data[:, 1]
    time_ps = steps * timestep_ps

    # FFT
    ke_detrended = ke - np.mean(ke)
    freq = np.fft.fftfreq(len(time_ps), d=timestep_ps)
    fft_vals = np.fft.fft(ke_detrended)
    freq_pos = freq[freq > 0]
    amp = np.abs(fft_vals[freq > 0])
    amp_norm = amp / np.max(amp)

    freq_ghz = freq_pos * 1000  # THz → GHz

    plt.plot(freq_ghz, amp_norm, label=f'{field} V/nm', color=colors.get(field, 'black'))

# Formatting like the paper
plt.xlim(0, 300)
plt.ylim(0, 1.05)
plt.xlabel("Frequency (GHz)")
plt.ylabel("Amplitude (a.u.)")
plt.title("FFT Spectrum vs Electric Field")
plt.legend(title="E-field")
plt.grid(True)
plt.tight_layout()
plt.savefig("fft_combined_all_fields.png", dpi=300)
plt.show()