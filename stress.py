import numpy as np
import matplotlib.pyplot as plt
import ase.io

# Placeholder function for calculating stress tensor and Reynolds stress tensor
def calculate_stress_tensor(positions, forces, box_volume):
    num_particles = len(positions)
    stress_tensor = np.zeros((nx, ny, nz, 3, 3))  # 3D stress tensor based on grid points

    for i in range(3):
        for j in range(i, 3):
            stress_tensor[:, :, :, i, j] = np.sum(positions[:, i] * forces[:, j]) + np.sum(positions[:, j] * forces[:, i])
            if i == j:
                stress_tensor[:, :, :, i, j] += np.sum([np.dot(f, r) for f, r in zip(forces, positions)])
            stress_tensor[:, :, :, i, j] /= box_volume
            stress_tensor[:, :, :, j, i] = stress_tensor[:, :, :, i, j]  # stress tensor is symmetric

    return stress_tensor

def calculate_reynolds_stress_tensor(positions, forces, velocities, box_volume):
    num_particles = len(positions)
    reynolds_stress_tensor = np.zeros((nx, ny, nz, 3, 3))  # 3D Reynolds stress tensor based on grid points

    for i in range(3):
        for j in range(i, 3):
            reynolds_stress_tensor[:, :, :, i, j] = np.sum(positions[:, i] * forces[:, j]) + np.sum(positions[:, j] * forces[:, i])
            reynolds_stress_tensor[:, :, :, i, j] /= box_volume
            reynolds_stress_tensor[:, :, :, j, i] = reynolds_stress_tensor[:, :, :, i, j]  # reynolds stress tensor is symmetric
    
    for i in range(3):
        for j in range(i, 3):
            reynolds_stress_tensor[:, :, :, i, j] -= np.sum(velocities[:, i] * velocities[:, j])
            reynolds_stress_tensor[:, :, :, j, i] = reynolds_stress_tensor[:, :, :, i, j]  # reynolds stress tensor is symmetric

    return reynolds_stress_tensor

# Placeholder function for interactive analysis of eigenvalues at selected grid points
def interactive_analysis(grid_points, stress_tensor, reynolds_stress_tensor):
    for grid_point in grid_points:
        stress_tensor_point = stress_tensor[grid_point[0], grid_point[1], grid_point[2]]
        reynolds_stress_tensor_point = reynolds_stress_tensor[grid_point[0], grid_point[1], grid_point[2]]
        diff_tensor = stress_tensor_point - reynolds_stress_tensor_point
        eigenvalues = np.linalg.eigvals(diff_tensor)
        print(f"Eigenvalues at grid point {grid_point}: {eigenvalues}")

# Example usage:
nx, ny, nz = 10, 10, 10  # Number of grid points along x, y, z axes

# Assume positions, forces, velocities, and box_volume are obtained from simulation data
traj = ase.io.read("C:/Users/sande/Desktop/lammps/flow/10size.dump", format='lammps-dump-text', index=':')

# Replace these placeholders with actual data from your simulation
x_e = []
y_e = []
z_e = []
d = 5
positions = traj[d].get_positions()  # Example positions of 200 particles in 3D space
forces = traj[d].get_forces()        # Example forces acting on the particles
velocities = traj[d].get_velocities()  # Example velocities of the particles
box_volume = 1093500  # Example box volume

# Calculate stress tensor and Reynolds stress tensor
stress_tensor = calculate_stress_tensor(positions, forces, box_volume)
reynolds_stress_tensor = calculate_reynolds_stress_tensor(positions, forces, velocities, box_volume)

# Calculate stress tensor along each axis
stress_tensor_x = np.sum(stress_tensor, axis=(1, 2, 3))
stress_tensor_y = np.sum(stress_tensor, axis=(0, 2, 3))
stress_tensor_z = np.sum(stress_tensor, axis=(0, 1, 3))

print("-------------x----------------")
print(stress_tensor_x[0][0])
print("-------------Y----------------")
print(stress_tensor_y[0][1])
print("-------------Z----------------")
print(stress_tensor_z[0][2])

# Choose one or more grid points (replace with actual grid indices)
grid_points = [(0, 1, 0), (3, 3, 3), (5, 5, 5)]

# Perform interactive analysis
interactive_analysis(grid_points, stress_tensor, reynolds_stress_tensor)

