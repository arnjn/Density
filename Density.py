import numpy as np
import ase.io
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

print("Start")
#Center difference 1 sided Gradient

def calculate_stress_tensor(positions, forces, box_volume):
    stress_tensor = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            stress_tensor[i, j] = np.sum(positions[:, i] * forces[:, j]) + np.sum(positions[:, j] * forces[:, i])
            if i == j:
                stress_tensor[i, j] += np.sum([np.dot(f, r) for f, r in zip(forces, positions)])
            stress_tensor[i, j] /= box_volume
            stress_tensor[j, i] = stress_tensor[i, j] # stress tensor is symmetric
    return stress_tensor

def calculate_density(position, grid ,grid_spacing=1.5):
    val = 0
    g = grid
    for coordinates in position:
        #a =  np.linalg.norm(abs(g-coordinates))
        if np.linalg.norm(abs(g - coordinates))/grid_spacing < 1:
            num = np.dot(g-coordinates,g-coordinates)
            a = (np.sqrt(num)/grid_spacing)
            e = np.exp(-1/(1-a**2))
            val += e
        else:
            val += 0

    a= 1/(grid_spacing**(3))
    return val*(grid_spacing**(-3))

def calculate_reynolds_stress_tensor(positions, forces, velocities, box_volume):
    reynolds_stress_tensor = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            reynolds_stress_tensor[i, j] = np.sum(positions[:, i] * forces[:, j]) + np.sum(positions[:, j] * forces[:, i])
            reynolds_stress_tensor[i, j] /= box_volume
            reynolds_stress_tensor[j, i] = reynolds_stress_tensor[i, j]  # reynolds stress tensor is symmetric
    
    for i in range(3):
        for j in range(i, 3):
            reynolds_stress_tensor[i, j] -= np.sum(velocities[:, i] * velocities[:, j])
            reynolds_stress_tensor[j, i] = reynolds_stress_tensor[i, j]  # reynolds stress tensor is symmetric
    
    return reynolds_stress_tensor

traj = ase.io.read("C:/Users/sande/Desktop/lammps/flow/sphere.dump" , format='lammps-dump-text',index=':')
grid_spacing = 1.5
#30,20,20
den = np.zeros((30,20,20))
stress = []
reynolds_stress = []
d = []
i = 0
# frame = traj[0]
# positions = frame.get_scaled_positions()
# for z in range(0,2):
#     for x in range(0,10):
#         for y in range(0,12):
#             print(y)
#             grid = [(x)*grid_spacing, (y)*grid_spacing,(z)*grid_spacing]
#             den[x,y,z] += calculate_density(positions, grid ,grid_spacing)
xden = []
yden = []
zden = []
for j in range(0,1):
    print(i)
    frame = traj[0]
    i = i+1
    positions = frame.get_positions()
    forces = (frame.get_forces())
    velocities = frame.get_velocities()
    stress_tensor = calculate_stress_tensor(positions, forces, 12000)
    stress.append(stress_tensor)
    reynolds_stress_tensor = calculate_reynolds_stress_tensor(positions, forces, velocities, 12000)
    reynolds_stress.append(reynolds_stress_tensor)

    # for z in range(0,20):
    #     for x in range(0,30):
    #         for y in range(0,20):
    #             grid = [x*grid_spacing, y*grid_spacing,z*grid_spacing]
    #             den[x,y,z] = calculate_density(positions, grid ,grid_spacing)

    #             #stress = calculate_stress_tensor(positions,forces,12000)
    
    #     # Compute gradient along each axis
    # # Apply Gaussian smoothing
    # den_smoothed = gaussian_filter(den, sigma=4)
    
    # # Compute gradient along each axis
    # grad_x = np.gradient(den_smoothed, axis=0, edge_order=1)
    # grad_y = np.gradient(den_smoothed, axis=1, edge_order=1)
    # grad_z = np.gradient(den_smoothed, axis=2, edge_order=1)
    
    # # grad_y = np.gradient(den, axis=1, edge_order=2)
    # # grad_z = np.gradient(den, axis=2, edge_order=2)
    
    # # Compute sum along each axis
    # sum_x = np.sum(grad_x, axis=(1, 2))
    # sum_y = np.sum(grad_y, axis=(0, 2))
    # sum_z = np.sum(grad_z, axis=(0, 1))
    # xden.append(np.sum(sum_x))
    # yden.append(np.sum(sum_y))
    # zden.append(np.sum(sum_z))
    # integrated_density = np.trapz(np.trapz(np.trapz(den, axis=0), axis=0), axis=0)
    # d.append(integrated_density)

timesteps = range(0,50)

# Plot the contour of stress tensor components along different axes
print((stress))

# for i in range(3):
#     plt.figure(figsize=(8, 6))
#     plt.plot(timesteps, [stress_tensor[i, i] for stress_tensor in stress])
#     plt.xlabel('Timestep')
#     plt.ylabel('Stress Component')
#     plt.title(f'Stress Component along {["X", "Y", "Z"][i]} Axis')
    
# plt.show()

# plt.figure(figsize=(8, 6))
# plt.plot(timesteps, xden, label='X Axis')
# plt.plot(timesteps, yden, label='Y Axis')
# plt.plot(timesteps, zden, label='Z Axis',linestyle='--')
# plt.title('Sum of Gradient along Each Axis over Time')
# plt.xlabel('Timestep')
# plt.ylabel('Sum of Gradient')
# plt.legend()
# plt.show()

# Plot contour of x-z plane at y=0 from the d matrix
#-----------------Extra Plot/Contour
# sum_z_values = np.sum(den, axis=2)


# # Plot the contour of the x-y projection
# plt.figure(figsize=(8, 6))
# plt.contourf(np.arange(sum_z_values.shape[0]), np.arange(sum_z_values.shape[1]), sum_z_values.T, cmap='viridis')
# plt.colorbar()
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Contour Plot of Density (x-y Projection)')
# # Create a grid of x and y values
# # Sum the stress tensors along the 0-axis
s = np.sum(np.array(stress), axis=0)

# # Reshape the 3x3 array into a 1D array
zp = s.reshape(1, -1)

# Create meshgrids for x and y
x = np.arange(s.shape[1])
y = np.arange(s.shape[0])
X, Y = np.meshgrid(x, y)
zp2 = zp.reshape(X.shape)
# Plot the contour
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, zp2, cmap='viridis')
plt.colorbar(label='Stress')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Stress Tensor')
#plt.show()

print("--------------------------------")
subtracted_list = [x - y for x, y in zip(stress, reynolds_stress)]
print(subtracted_list)

s = np.sum(np.array(subtracted_list), axis=0)

# # Reshape the 3x3 array into a 1D array
zp = s.reshape(1, -1)

# Create meshgrids for x and y
x = np.arange(s.shape[1])
y = np.arange(s.shape[0])
X, Y = np.meshgrid(x, y)
zp2 = zp.reshape(X.shape)
# Plot the contour
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, zp2, cmap='viridis')
plt.colorbar(label='subtracted_list')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Stress Tensor')
plt.show()
