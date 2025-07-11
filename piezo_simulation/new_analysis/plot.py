import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('polarization.dat', comments='#')
timesteps = data[:, 0]
px = data[:, 1]
py = data[:, 2]
pz = data[:, 3]

plt.plot(timesteps, pz, label='Pz')
plt.plot(timesteps, px, label='Px')
plt.plot(timesteps, py, label='Py')
plt.xlabel('Timestep')
plt.ylabel('Polarization (e·Å)')
plt.legend()
plt.title('Polarization vs Time')
plt.grid(True)
plt.show()
