import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt("particle path.dat", unpack=True)
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, color='blue', linewidth=1.5, alpha=0.8)
ax.scatter(x[0], y[0], z[0], color='green', s=100, marker='o', edgecolors='black')
ax.scatter(x[-1], y[-1], z[-1], color='red', s=100, marker='o', edgecolors='black')
ax.set_xlabel('x [M]', fontsize=12, labelpad=10)
ax.set_ylabel('y [M]', fontsize=12, labelpad=10)
ax.set_zlabel('z [M]', fontsize=12, labelpad=10)
max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
mid_x = (x.max()+x.min()) * 0.5
mid_y = (y.max()+y.min()) * 0.5
mid_z = (z.max()+z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)
ax.view_init(elev=20, azim=-120)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('particle orbit.png', dpi=200, bbox_inches='tight')
plt.show()