import numpy as np
import ase.io
import matplotlib.pyplot as plt

print("Start")

def calculate_density(position, grid ,grid_spacing=1.5):
    val = 0
    g = grid
    for coordinates in position:
        #a =  np.linalg.norm(abs(g-coordinates))
        
        if np.linalg.norm(abs(g - coordinates))/grid_spacing <= 1:
            num = np.dot(g-coordinates,g-coordinates)
            e = np.exp(-1/(1-(num/grid_spacing**2)))
            val += e
        else:
            val += 0

    a= 1/(grid_spacing**(3))
    return val*(grid_spacing**(-3))


traj = ase.io.read("C:/Users/sande/Desktop/lammps/flow/rand.dump" , format='lammps-dump-text',index=':')
x_den = []
y_den = []
z_den = []
grid_spacing = 1.5
x_range = np.arange(0, 25.0, grid_spacing)
y_range = np.arange(0, 30.0, grid_spacing)
z_range = np.arange(0, 5.0, grid_spacing)
grid = np.array(np.meshgrid(x_range,y_range,z_range)).T.reshape(-1,3)
grid_val = np.array(np.meshgrid(x_range,y_range,z_range)).T.reshape(-1,3)
grid_val[:] = 0
i = 0

for frame in traj:
    print(i)
    i = i+1
    positions = frame.get_scaled_positions()

    for g, gval in zip(grid,grid_val):
        gval += calculate_density(positions,g)

    x = sum(grid_val[:, 0])/len(grid_val[:, 0])
    y = sum(grid_val[:, 1])/len(grid_val[:, 1])
    z = sum(grid_val[:, 2])/len(grid_val[:, 2])
    x_den.append(x)
    y_den.append(y)
    z_den.append(z)


avg_den = [a + b + c for a, b, c in zip(x_den, y_den, z_den)]
plt.figure(1)
plt.plot(x_den)
plt.xlabel("Time step")
plt.ylabel("X axis density")
plt.show()

plt.figure(2)
plt.plot(y_den)
plt.xlabel("Time step")
plt.ylabel("Y axis density")
plt.show()

plt.figure(3)
plt.plot(z_den)
plt.xlabel("Time step")
plt.ylabel("Z axis density")
plt.show()

plt.figure(4)
plt.plot(avg_den)
plt.ylabel("total density")
plt.xlabel("Time steps")
plt.show()
