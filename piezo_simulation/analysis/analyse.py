#!/usr/bin/env python3
"""
Complete Analysis and Visualization Tool for BaTiO₃ Piezoelectric Simulation
This script analyzes LAMMPS output files and creates publication-quality plots
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from scipy import stats
import os

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class PiezoelectricAnalyzer:
    def __init__(self, data_dir='.'):
        self.data_dir = data_dir
        self.polarization_data = None
        self.strain_data = None
        
    def load_data(self):
        """Load all simulation data files"""
        try:
            # Load polarization vs time data
            polar_file = os.path.join(self.data_dir, 'polarization.dat')
            if os.path.exists(polar_file):
                self.polarization_data = self.read_lammps_data(polar_file)
                print(f"Loaded polarization data: {len(self.polarization_data)} time steps")
            
            # Load strain vs polarization data
            strain_file = os.path.join(self.data_dir, 'strain_polar.dat')
            if os.path.exists(strain_file):
                self.strain_data = self.read_strain_data(strain_file)
                print(f"Loaded strain data: {len(self.strain_data)} strain points")
                
        except Exception as e:
            print(f"Error loading data: {e}")
    
    def read_lammps_data(self, filename):
        """Read LAMMPS thermo-style output data"""
        data = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Skip header lines and read data
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            try:
                values = [float(x) for x in line.split()]
                if len(values) >= 4:  # timestep, Px, Py, Pz
                    data.append(values)
            except ValueError:
                continue
                
        return np.array(data)
    
    def read_strain_data(self, filename):
        """Read strain vs polarization data"""
        data = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            try:
                values = [float(x) for x in line.split()]
                if len(values) >= 2:  # strain, polarization
                    data.append(values)
            except ValueError:
                continue
                
        return np.array(data)
    
    def calculate_piezoelectric_coefficient(self):
        """Calculate piezoelectric coefficient d33"""
        if self.strain_data is None or len(self.strain_data) < 2:
            return None
            
        strain = self.strain_data[:, 0]
        polarization = self.strain_data[:, 1]
        
        # Linear fit to get d33 coefficient
        slope, intercept, r_value, p_value, std_err = stats.linregress(strain, polarization)
        
        # d33 is the slope of P vs strain relationship
        d33 = slope  # in units consistent with your simulation
        
        return {
            'd33': d33,
            'r_squared': r_value**2,
            'p_value': p_value,
            'std_error': std_err
        }
    
    def create_comprehensive_plots(self):
        """Create all analysis plots"""
        fig = plt.figure(figsize=(16, 12))
        
        # Create subplot layout
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # Plot 1: Polarization evolution over time
        ax1 = fig.add_subplot(gs[0, :2])
        self.plot_polarization_time(ax1)
        
        # Plot 2: Polarization components
        ax2 = fig.add_subplot(gs[0, 2])
        self.plot_polarization_components(ax2)
        
        # Plot 3: Strain-Polarization relationship
        ax3 = fig.add_subplot(gs[1, :2])
        self.plot_strain_polarization(ax3)
        
        # Plot 4: Piezoelectric coefficient visualization
        ax4 = fig.add_subplot(gs[1, 2])
        self.plot_piezo_coefficient(ax4)
        
        # Plot 5: 3D polarization trajectory
        ax5 = fig.add_subplot(gs[2, :], projection='3d')
        self.plot_3d_polarization(ax5)
        
        plt.suptitle('BaTiO₃ Piezoelectric Analysis', fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.savefig('piezoelectric_analysis_complete.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Create individual detailed plots
        self.create_detailed_plots()
    
    def plot_polarization_time(self, ax):
        """Plot polarization evolution over time"""
        if self.polarization_data is None:
            ax.text(0.5, 0.5, 'No polarization data available', 
                   ha='center', va='center', transform=ax.transAxes)
            return
            
        time_steps = self.polarization_data[:, 0]
        px = self.polarization_data[:, 1]
        py = self.polarization_data[:, 2]
        pz = self.polarization_data[:, 3]
        
        ax.plot(time_steps, px, 'r-', linewidth=2, label='Px', alpha=0.8)
        ax.plot(time_steps, py, 'g-', linewidth=2, label='Py', alpha=0.8)
        ax.plot(time_steps, pz, 'b-', linewidth=2, label='Pz', alpha=0.8)
        
        ax.set_xlabel('Time Step')
        ax.set_ylabel('Polarization (e·Å)')
        ax.set_title('Polarization Evolution During Stress Application')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def plot_polarization_components(self, ax):
        """Plot final polarization components as bar chart"""
        if self.polarization_data is None:
            return
            
        # Take final values
        final_px = self.polarization_data[-1, 1]
        final_py = self.polarization_data[-1, 2]
        final_pz = self.polarization_data[-1, 3]
        
        components = ['Px', 'Py', 'Pz']
        values = [final_px, final_py, final_pz]
        colors = ['red', 'green', 'blue']
        
        bars = ax.bar(components, values, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Final Polarization (e·Å)')
        ax.set_title('Final Polarization Components')
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                   f'{value:.3f}', ha='center', va='bottom', fontweight='bold')
    
    def plot_strain_polarization(self, ax):
        """Plot strain vs polarization relationship"""
        if self.strain_data is None:
            ax.text(0.5, 0.5, 'No strain data available', 
                   ha='center', va='center', transform=ax.transAxes)
            return
            
        strain = self.strain_data[:, 0]
        polarization = self.strain_data[:, 1]
        
        # Plot data points
        ax.scatter(strain * 100, polarization, c='blue', s=60, alpha=0.7, edgecolors='navy')
        
        # Fit and plot trend line
        if len(strain) > 1:
            z = np.polyfit(strain, polarization, 1)
            p = np.poly1d(z)
            ax.plot(strain * 100, p(strain), 'r--', linewidth=2, alpha=0.8, label=f'Linear fit: slope = {z[0]:.2f}')
            ax.legend()
        
        ax.set_xlabel('Strain (%)')
        ax.set_ylabel('Polarization (e·Å)')
        ax.set_title('Piezoelectric Response: Strain vs Polarization')
        ax.grid(True, alpha=0.3)
    
    def plot_piezo_coefficient(self, ax):
        """Visualize piezoelectric coefficient"""
        coeff_data = self.calculate_piezoelectric_coefficient()
        
        if coeff_data is None:
            ax.text(0.5, 0.5, 'Insufficient data\nfor coefficient calculation', 
                   ha='center', va='center', transform=ax.transAxes)
            return
        
        # Create a gauge-like visualization
        d33 = coeff_data['d33']
        r_squared = coeff_data['r_squared']
        
        # Typical BaTiO3 d33 ranges from 85-190 pC/N
        # Scale our simulation units appropriately
        ax.text(0.5, 0.7, 'Piezoelectric Coefficient', ha='center', va='center', 
               transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.5, 0.5, f'd₃₃ = {d33:.3f}\n(simulation units)', ha='center', va='center',
               transform=ax.transAxes, fontsize=11)
        ax.text(0.5, 0.3, f'R² = {r_squared:.3f}', ha='center', va='center',
               transform=ax.transAxes, fontsize=10)
        
        # Create colored background based on R² quality
        if r_squared > 0.9:
            color = 'lightgreen'
            quality = 'Excellent fit'
        elif r_squared > 0.7:
            color = 'lightyellow' 
            quality = 'Good fit'
        else:
            color = 'lightcoral'
            quality = 'Poor fit'
            
        bbox = FancyBboxPatch((0.1, 0.1), 0.8, 0.8, boxstyle="round,pad=0.1", 
                             facecolor=color, alpha=0.3, transform=ax.transAxes)
        ax.add_patch(bbox)
        ax.text(0.5, 0.1, quality, ha='center', va='center',
               transform=ax.transAxes, fontsize=9, style='italic')
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    
    def plot_3d_polarization(self, ax):
        """3D visualization of polarization vector evolution"""
        if self.polarization_data is None:
            return
            
        px = self.polarization_data[:, 1]
        py = self.polarization_data[:, 2]
        pz = self.polarization_data[:, 3]
        
        # Color code by time
        colors = plt.cm.viridis(np.linspace(0, 1, len(px)))
        
        ax.scatter(px, py, pz, c=colors, s=20, alpha=0.6)
        
        # Draw arrow from origin to final polarization
        ax.quiver(0, 0, 0, px[-1], py[-1], pz[-1], 
                 color='red', arrow_length_ratio=0.1, linewidth=3, alpha=0.8)
        
        ax.set_xlabel('Px (e·Å)')
        ax.set_ylabel('Py (e·Å)')
        ax.set_zlabel('Pz (e·Å)')
        ax.set_title('3D Polarization Vector Evolution')
    
    def create_detailed_plots(self):
        """Create individual detailed plots"""
        
        # High-resolution strain-polarization plot
        if self.strain_data is not None:
            fig, ax = plt.subplots(figsize=(10, 6))
            strain = self.strain_data[:, 0]
            polarization = self.strain_data[:, 1]
            
            ax.scatter(strain * 100, polarization, c='darkblue', s=100, alpha=0.7, 
                      edgecolors='navy', linewidth=2, label='Simulation Data')
            
            # Fit curve
            if len(strain) > 1:
                z = np.polyfit(strain, polarization, 1)
                p = np.poly1d(z)
                strain_fit = np.linspace(strain.min(), strain.max(), 100)
                ax.plot(strain_fit * 100, p(strain_fit), 'r-', linewidth=3, alpha=0.8,
                       label=f'Linear Fit: d₃₃ = {z[0]:.3f}')
            
            ax.set_xlabel('Applied Strain (%)', fontsize=14)
            ax.set_ylabel('Induced Polarization (e·Å)', fontsize=14)
            ax.set_title('BaTiO₃ Piezoelectric Response Curve', fontsize=16, fontweight='bold')
            ax.legend(fontsize=12)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig('strain_polarization_detailed.png', dpi=300, bbox_inches='tight')
            plt.show()
    
    def generate_report(self):
        """Generate analysis report"""
        report = []
        report.append("="*50)
        report.append("BaTiO₃ PIEZOELECTRIC SIMULATION ANALYSIS REPORT")
        report.append("="*50)
        
        if self.polarization_data is not None:
            report.append(f"\nPolarization Data Summary:")
            report.append(f"- Total time steps: {len(self.polarization_data)}")
            report.append(f"- Initial polarization (Px, Py, Pz): ({self.polarization_data[0,1]:.3f}, {self.polarization_data[0,2]:.3f}, {self.polarization_data[0,3]:.3f})")
            report.append(f"- Final polarization (Px, Py, Pz): ({self.polarization_data[-1,1]:.3f}, {self.polarization_data[-1,2]:.3f}, {self.polarization_data[-1,3]:.3f})")
            
            # Calculate changes
            dpx = self.polarization_data[-1,1] - self.polarization_data[0,1]
            dpy = self.polarization_data[-1,2] - self.polarization_data[0,2]  
            dpz = self.polarization_data[-1,3] - self.polarization_data[0,3]
            report.append(f"- Polarization changes (ΔPx, ΔPy, ΔPz): ({dpx:.3f}, {dpy:.3f}, {dpz:.3f})")
        
        if self.strain_data is not None:
            report.append(f"\nStrain Analysis:")
            report.append(f"- Strain points analyzed: {len(self.strain_data)}")
            report.append(f"- Maximum strain: {self.strain_data[-1,0]*100:.2f}%")
            
            coeff_data = self.calculate_piezoelectric_coefficient()
            if coeff_data:
                report.append(f"- Piezoelectric coefficient d₃₃: {coeff_data['d33']:.6f} (simulation units)")
                report.append(f"- Linear fit quality (R²): {coeff_data['r_squared']:.4f}")
                if coeff_data['r_squared'] > 0.9:
                    report.append("- Fit quality: EXCELLENT")
                elif coeff_data['r_squared'] > 0.7:
                    report.append("- Fit quality: GOOD")
                else:
                    report.append("- Fit quality: POOR - Consider longer simulation or different parameters")
        
        report.append(f"\nFiles Generated:")
        report.append(f"- piezoelectric_analysis_complete.png: Comprehensive analysis")
        report.append(f"- strain_polarization_detailed.png: Detailed strain response")
        report.append(f"- analysis_report.txt: This report")
        
        report.append(f"\nRecommendations:")
        report.append(f"- For VMD visualization: Load trajectory.dump and color atoms by type")
        report.append(f"- For quantitative analysis: Use the strain-polarization relationship")
        report.append(f"- For comparison: Literature d₃₃ for BaTiO₃ is typically 85-190 pC/N")
        
        report_text = "\n".join(report)
        
        # Save report
        with open('analysis_report.txt', 'w') as f:
            f.write(report_text)
        
        print(report_text)
        return report_text

def main():
    """Main analysis function"""
    print("Starting BaTiO₃ Piezoelectric Analysis...")
    
    # Initialize analyzer
    analyzer = PiezoelectricAnalyzer()
    
    # Load data
    analyzer.load_data()
    
    # Create visualizations
    analyzer.create_comprehensive_plots()
    
    # Generate report
    analyzer.generate_report()
    
    print("\nAnalysis complete! Check the generated PNG files and analysis_report.txt")
    print("\nFor VMD visualization:")
    print("1. Open VMD")
    print("2. Load trajectory.dump")
    print("3. Graphics -> Representations")
    print("4. Color by 'Type' to distinguish Ba (red), Ti (green), O (blue)")
    print("5. Add polarization vectors if desired")

if __name__ == "__main__":
    main()