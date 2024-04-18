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
box_volume = 1093500  # Example box volume

# Iterate through each frame in traj
for frame in traj:
    positions = frame.get_positions()  # Example positions of 200 particles in 3D space
    forces = frame.get_forces()        # Example forces acting on the particles
    velocities = frame.get_velocities()  # Example velocities of the particles
    
    # Calculate stress tensor and Reynolds stress tensor for the current frame
    stress_tensor = calculate_stress_tensor(positions, forces, box_volume)
    reynolds_stress_tensor = calculate_reynolds_stress_tensor(positions, forces, velocities, box_volume)

    # Calculate stress tensor along each axis for the current frame
    stress_tensor_x = np.sum(stress_tensor, axis=(1, 2, 3))
    stress_tensor_y = np.sum(stress_tensor, axis=(0, 2, 3))
    stress_tensor_z = np.sum(stress_tensor, axis=(0, 1, 3))

    # Append x_e, y_e, and z_e values for the current frame
    x_e.append(stress_tensor_x[0][0])  # Adjust indices based on your stress tensor shape
    y_e.append(stress_tensor_y[0][1])
    z_e.append(stress_tensor_z[0][2])

# Plotting x_e, y_e, and z_e values over simulation time
time_steps = np.arange(len(traj))
plt.plot(time_steps, x_e, label='Stress Tensor along X-axis')
plt.plot(time_steps, y_e, label='Stress Tensor along Y-axis')
plt.plot(time_steps, z_e, label='Stress Tensor along Z-axis')
plt.xlabel('Simulation Time')
plt.ylabel('Stress Tensor Values')
plt.title('Stress Tensor Variation Over Time')
plt.legend()
plt.show()


# Placeholder function for interactive analysis of eigenvalues at selected grid points
def interactive_analysis(grid_points, stress_tensor, reynolds_stress_tensor):
    time_steps = np.arange(len(traj))  # Define time steps for the entire simulation
    
    for grid_point in grid_points:
        eigenvalues_over_time = []  # Store eigenvalues for each time step
        for frame in traj:
            positions = frame.get_positions()  # Example positions of 200 particles in 3D space
            forces = frame.get_forces()        # Example forces acting on the particles
            velocities = frame.get_velocities()  # Example velocities of the particles
            
            # Calculate stress tensor and Reynolds stress tensor for the current frame
            stress_tensor = calculate_stress_tensor(positions, forces, box_volume)
            reynolds_stress_tensor = calculate_reynolds_stress_tensor(positions, forces, velocities, box_volume)

            # Calculate eigenvalues for the current grid point and frame
            stress_tensor_point = stress_tensor[grid_point[0], grid_point[1], grid_point[2]]
            reynolds_stress_tensor_point = reynolds_stress_tensor[grid_point[0], grid_point[1], grid_point[2]]
            diff_tensor = stress_tensor_point - reynolds_stress_tensor_point
            eigenvalues = np.linalg.eigvals(diff_tensor)
            
            eigenvalues_over_time.append(eigenvalues)

        # Plot eigenvalues over simulation time for the current grid point
        eigenvalues_over_time = np.array(eigenvalues_over_time)
        plt.figure(1)
        plt.plot(time_steps, eigenvalues_over_time[:, 0], label='Eigenvalue 1',c="blue")
        plt.xlabel('Simulation Time')
        plt.ylabel('Eigenvalues')
        plt.title(f'Eigenvalues Variation Over Time for Grid Point {grid_point}')
        plt.legend()
        plt.figure(2)
        plt.plot(time_steps, eigenvalues_over_time[:, 1], label='Eigenvalue 2',c="green")
        plt.xlabel('Simulation Time')
        plt.ylabel('Eigenvalues')
        plt.title(f'Eigenvalues Variation Over Time for Grid Point {grid_point}')
        plt.legend()
        plt.figure(3)
        plt.plot(time_steps, eigenvalues_over_time[:, 2], label='Eigenvalue 3',c="red")
        plt.xlabel('Simulation Time')
        plt.ylabel('Eigenvalues')
        plt.title(f'Eigenvalues Variation Over Time for Grid Point {grid_point}')
        plt.legend()
        plt.show()

# Choose one or more grid points (replace with actual grid indices)
grid_points = [(1, 1, 0)]

# Perform interactive analysis
interactive_analysis(grid_points, stress_tensor, reynolds_stress_tensor)
