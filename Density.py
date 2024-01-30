import numpy as np
import ase.io
import matplotlib.pyplot as plt

print("Start")

def calculate_density(position, grid ,grid_spacing=2.5):
    val = 0
    g = grid
    for coordinates in position:
        #a =  np.linalg.norm(abs(g-coordinates))
        if np.linalg.norm(abs(g - coordinates))/grid_spacing <= 1:
            num = np.dot(g-coordinates,g-coordinates)
            a = (num/grid_spacing)
            e = np.exp(-1/(1-a**2))
            val += e
        else:
            val += 0

    a= 1/(grid_spacing**(3))
    return val*(grid_spacing**(-3))


traj = ase.io.read("C:/Users/sande/Desktop/lammps/flow/rand.dump" , format='lammps-dump-text',index=':')
grid_spacing = 2.5
den = np.zeros((10,12,2))
d = []
i = 0
frame = traj[0]
positions = frame.get_scaled_positions()
for z in range(0,2):
    for x in range(0,10):
        for y in range(0,12):
            print(y)
            grid = [(x)*grid_spacing, (y)*grid_spacing,(z)*grid_spacing]
            den[x,y,z] += calculate_density(positions, grid ,grid_spacing)
# for frame in traj:
#     print(i)
#     #frame = traj[frame]
#     i = i+1
#     positions = frame.get_scaled_positions()
#     for z in range(0,2):
#         for x in range(0,10):
#             for y in range(0,12):
#                 grid = [x*grid_spacing, y*grid_spacing,z*grid_spacing]
#                 den[x,y,z] = calculate_density(positions, grid ,grid_spacing)
    
#     d.append(np.sum(den))

x = np.arange(2.5, 25.1, 2.5)

# Arrange y from 0 to 30 by 2.5
y = np.arange(2.5, 30.1, 2.5)

# Create a meshgrid for x and y
X, Y = np.meshgrid(x, y)

# Plot a 2D contour plot with z-values always 1
plt.contour(X, Y, den[:, :, 1].T, levels=[1], colors='k')

# Set labels and title
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Contour Plot with z=1')

# # Show the plot
plt.show()




    