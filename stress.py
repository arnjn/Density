import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
import ase.io
from scipy.ndimage import gaussian_filter


# Read LAMMPS trajectory using ASE
traj = ase.io.read("C:/Users/sande/Desktop/lammps/flow/10size.dump", format='lammps-dump-text', index=':')

# Define grid parameters and initialize grid arrays
nx, ny, nz = 100, 100, 100  # Number of grid points along x, y, z axes
stress_tensor = np.zeros((nx, ny, nz, 3, 3))  # Assuming a 3x3 stress tensor
reynolds_stress_tensor = np.zeros((nx, ny, nz, 3, 3))  # Assuming a 3x3 Reynolds stress tensor

# Assuming 'grid_point' is the grid point you are interested in (e.g., grid_point = (i, j, k))

# Initialize lists to store eigenvalues and time steps
eigenvalues_list = []
time_steps = []

# Print available keys in the trajectory to check for stress tensor data
print("Available keys in trajectory arrays:")
for key in traj[0].arrays.keys():
    print(key)

# Loop through each frame in the trajectory
for idx, frame in enumerate(traj):
    # Check if 'virial' array is available in the frame
    if 'virial' in frame.arrays:
        # Get stress tensor components from the frame
        stress_components = frame.arrays['virial']

        # Calculate stress tensor at each grid point based on atomic positions and stress components
        for atom_pos, stress_comp in zip(frame.get_positions(), stress_components):
            # Calculate grid indices based on atomic positions
            i, j, k = int(atom_pos[0] / (lx / nx)), int(atom_pos[1] / (ly / ny)), int(atom_pos[2] / (lz / nz))

            # Assign stress tensor components to the corresponding grid point
            stress_tensor[i, j, k] += stress_comp.reshape((3, 3))

# Optionally, perform smoothing on the stress tensor and Reynolds stress tensor data
smoothed_stress_tensor = gaussian_filter(stress_tensor, sigma=1)

# Loop through each frame in the trajectory again to calculate Reynolds stress tensor
for idx, frame in enumerate(traj):
    if 'virial' in frame.arrays:
        # Calculate Reynolds stress tensor based on stress tensor components
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    reynolds_stress_tensor[i, j, k] += stress_tensor[i, j, k]  # Example: Assigning stress tensor for Reynolds stress (adjust as per your calculation)

# Check if stress tensor and Reynolds stress tensor data are available
if np.any(stress_tensor) and np.any(reynolds_stress_tensor):
    # Loop through each frame in the trajectory again to calculate eigenvalues at the specified grid point
    for idx, frame in enumerate(traj):
        # Get stress tensor and Reynolds stress tensor values at the specified grid point
        stress_tensor_point = smoothed_stress_tensor[grid_point[0], grid_point[1], grid_point[2]]
        reynolds_stress_tensor_point = smoothed_reynolds_stress_tensor[grid_point[0], grid_point[1], grid_point[2]]

        # Calculate the difference between stress tensor and Reynolds stress tensor
        diff_tensor = stress_tensor_point - reynolds_stress_tensor_point

        # Calculate eigenvalues of the difference tensor
        eigenvalues = eigvals(diff_tensor)

        # Append eigenvalues to the list
        eigenvalues_list.append(eigenvalues)
        
        # Append time step to the list (assuming each frame corresponds to a time step)
        time_steps.append(idx)

    # Convert lists to NumPy arrays for plotting
    eigenvalues_array = np.array(eigenvalues_list)
    time_steps_array = np.array(time_steps)

    # Plot eigenvalues through time
    plt.figure()
    for i in range(len(eigenvalues_array[0])):
        plt.plot(time_steps_array, eigenvalues_array[:, i], label=f'Eigenvalue {i+1}')

    # Set labels and title
    plt.xlabel('Time Step')
    plt.ylabel('Eigenvalues')
    plt.title('Eigenvalues of Stress Tensor - Reynolds Stress Tensor')
    plt.legend()
    plt.show()
else:
    print("Stress tensor or Reynolds stress tensor data not available in the trajectory.")
