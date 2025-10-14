import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
x, y, z = np.loadtxt("image86.dat", unpack=True)
#x, y, z = np.loadtxt("image230.dat", unpack=True)
#x, y, z = np.loadtxt("hot-spot image.dat", unpack=True)
x_min, x_max = -15, 15
y_min, y_max = -15, 15
resolution = 1000
x_edges = np.linspace(x_min, x_max, resolution)
y_edges = np.linspace(y_min, y_max, resolution)
X, Y = np.meshgrid(x_edges, y_edges)
Z_grid = np.full_like(X, np.nan, dtype=float)
for i, (x_val, y_val, z_val) in enumerate(zip(x, y, z)):
    x_index = np.argmin(np.abs(x_edges - x_val))
    y_index = np.argmin(np.abs(y_edges - y_val))
    Z_grid[y_index, x_index] = z_val
eht_colors = [
    (0.00, '#000000'),
    (0.15, '#330000'),  
    (0.30, '#660000'),  
    (0.45, '#990000'),  
    (0.60, '#CC0000'),  
    (0.75, '#FF3300'),  
    (0.85, '#FF6600'),  
    (0.92, '#FF9900'),  
    (0.97, '#FFCC00'),  
    (1.00, '#FFFFFF')   
]
eht_cmap = LinearSegmentedColormap.from_list('eht_black_hole', eht_colors)
plt.figure(figsize=(10, 10), facecolor='black')
Z_masked = np.ma.masked_invalid(Z_grid)
eht_cmap.set_bad(color='black')
plt.imshow(Z_masked, 
           cmap=eht_cmap, 
           origin='lower', 
           extent=[x_min, x_max, y_min, y_max], 
           interpolation='none', 
           vmin=0, 
           vmax=0.5, 
           aspect='auto')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.gca().set_aspect('equal')
plt.axis('off')
#plt.text(-14, 12, '(a)', fontsize=60, fontname='Times New Roman', color='white')
plt.gca().set_facecolor('black')
plt.gcf().set_facecolor('black')
plt.tight_layout()
plt.savefig('image.png', dpi=300, bbox_inches='tight', facecolor='black')
plt.show()