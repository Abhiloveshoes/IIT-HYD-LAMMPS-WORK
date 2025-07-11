# funcs.py

import math

def myPy(t):
    """
    LAMMPS will call this function every timestep to compute the electric field.
    t = simulation step (e.g., 100000 to 200000)
    Returns: electric field amplitude in V/Å
    """
    E0 = 0.01  # Peak electric field (V/Å)
    freq = 0.001  # Oscillation frequency (arbitrary units)

    # Sinusoidal electric field: E(t) = E0 * sin(2πft)
    return E0 * math.sin(2 * math.pi * freq * t)
