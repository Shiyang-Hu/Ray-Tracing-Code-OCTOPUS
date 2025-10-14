import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt("hamiltonian error.dat", unpack=True)
z_max_data = np.max(z)
print(f"max error: {z_max_data:.6f}")
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
plt.figure(figsize=(10, 8))
#cmap = plt.cm.rainbow  
#cmap = plt.cm.jet        
#cmap = plt.cm.viridis    
#cmap = plt.cm.plasma     
#cmap = plt.cm.inferno    
cmap = plt.cm.turbo      
Z_masked = np.ma.masked_invalid(Z_grid)
cmap.set_bad(color='black')
im = plt.imshow(Z_masked, cmap=cmap, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                interpolation='none', vmin=-16, vmax=-6, aspect='auto')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.gca().set_aspect('equal')
plt.axis('off')
#plt.text(-14, 12, '(c)', fontsize=50, fontname='Times New Roman', color='white')
cbar = plt.colorbar(im, shrink=1.0, pad=0.01)
cbar.set_label(r'$\log_{10}\Delta$', fontsize=25, fontname='Times New Roman')
cbar.ax.tick_params(labelsize=14)
for label in cbar.ax.get_yticklabels():
    label.set_fontname('Times New Roman')
plt.tight_layout()
plt.savefig('error.jpeg', dpi=100, bbox_inches='tight')
plt.show()