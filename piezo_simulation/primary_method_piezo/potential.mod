# analysis_piezo.py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def analyze_piezoelectric_response():
    """Analyze LAMMPS piezoelectric simulation results"""
    
    # Read polarization data
    polar_data = np.loadtxt('polarization_time.dat', skiprows=2)
    time_steps = polar_data[:, 0]
    px = polar_data[:, 1]
    py = polar_data[:, 2] 
    pz = polar_data[:, 3]
    
    # Read stress data
    stress_data = np.loadtxt('stress_evolution.dat', skiprows=2)
    stress_zz = stress_data[:, 3]  # Normal stress in z-direction
    
    # Calculate piezoelectric coefficient d33
    # d33 = ΔP3/Δσ3 (polarization change per stress change)
    if len(stress_zz) > 100:
        initial_stress = np.mean(stress_zz[:10])
        final_stress = np.mean(stress_zz[-10:])
        initial_polar = np.mean(pz[:10])
        final_polar = np.mean(pz[-10:])
        
        d33 = (final_polar - initial_polar) / (final_stress - initial_stress)
        print(f"Calculated d33 coefficient: {d33:.2e} C/N")
    
    # Plot results
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Polarization vs time
    ax1.plot(time_steps, pz, 'b-', linewidth=2)
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel('Polarization Pz (e·Å)')
    ax1.set_title('Polarization Evolution')
    ax1.grid(True)
    
    # Stress vs time
    ax2.plot(time_steps[:len(stress_zz)], stress_zz, 'r-', linewidth=2)
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('Stress σ33 (GPa)')
    ax2.set_title('Applied Stress')
    ax2.grid(True)
    
    # Polarization components
    ax3.plot(time_steps, px, 'r-', label='Px', linewidth=2)
    ax3.plot(time_steps, py, 'g-', label='Py', linewidth=2)
    ax3.plot(time_steps, pz, 'b-', label='Pz', linewidth=2)
    ax3.set_xlabel('Time Step')
    ax3.set_ylabel('Polarization (e·Å)')
    ax3.set_title('All Polarization Components')
    ax3.legend()
    ax3.grid(True)
    
    # Stress-Polarization relationship
    if len(stress_zz) == len(pz):
        ax4.scatter(stress_zz, pz, c=time_steps[:len(stress_zz)], cmap='viridis')
        ax4.set_xlabel('Stress σ33 (GPa)')
        ax4.set_ylabel('Polarization Pz (e·Å)')
        ax4.set_title('Piezoelectric Response')
        ax4.grid(True)
        cbar = plt.colorbar(ax4.collections[0], ax=ax4)
        cbar.set_label('Time Step')
    
    plt.tight_layout()
    plt.savefig('piezoelectric_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    analyze_piezoelectric_response()