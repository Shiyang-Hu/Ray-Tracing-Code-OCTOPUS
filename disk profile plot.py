import numpy as np
import matplotlib.pyplot as plt
x, y = np.loadtxt("disk profile.dat", unpack=True)
x_min, x_max = -30, 30
y_min, y_max = -30, 30  
plt.figure(figsize=(10, 10))
ax = plt.gca()
plt.scatter(x, y, s=8, alpha=0.7, color='black', edgecolors='none')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.gca().set_aspect('equal') 
plt.grid(True, color='gray', linestyle='--', alpha=0.3)
plt.xlabel('x [M]', fontsize=25)
plt.ylabel('y [M]', fontsize=25)
ax.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()
plt.savefig('disk contour.png', dpi=600, bbox_inches='tight')
plt.show()