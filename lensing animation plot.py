import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import time
import os

start_time = time.time()
data = np.loadtxt("hot-spot image.dat")
num_frames = 200    
resolution = 500  
points_per_frame = resolution * resolution
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
eht_cmap = LinearSegmentedColormap.from_list('eht_style', eht_colors)

fig = plt.figure(figsize=(10, 10), facecolor='black')
ax = fig.add_subplot(111)
ax.set_facecolor('black')

x_min, x_max = np.min(data[:, 0]), np.max(data[:, 0])
y_min, y_max = np.min(data[:, 1]), np.max(data[:, 1])
z_min, z_max = np.min(data[:, 2]), np.max(data[:, 2])

frame_data = data[:points_per_frame]
x, y, z = frame_data[:, 0], frame_data[:, 1], frame_data[:, 2]

image_data = np.full((resolution, resolution), np.nan)
x_edges = np.linspace(x_min, x_max, resolution)
y_edges = np.linspace(y_min, y_max, resolution)

for i in range(len(x)):
    x_idx = np.argmin(np.abs(x_edges - x[i]))
    y_idx = np.argmin(np.abs(y_edges - y[i]))
    image_data[y_idx, x_idx] = z[i]

im = ax.imshow(image_data, cmap=eht_cmap, origin='lower', 
               extent=[x_min, x_max, y_min, y_max], 
               vmin=z_min, vmax=z_max, aspect='auto')
im.set_clim(z_min, z_max)
ax.axis('off')
plt.tight_layout()
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white', 
                   fontsize=14, fontweight='bold')
output_dir = os.path.dirname(os.path.abspath(__file__))
os.makedirs(output_dir, exist_ok=True)
def update(frame):
    print(f"Processing {frame+1}/{num_frames} frame...")
    start_idx = frame * points_per_frame
    end_idx = start_idx + points_per_frame
    frame_data = data[start_idx:end_idx]
    x, y, z = frame_data[:, 0], frame_data[:, 1], frame_data[:, 2]
    
    image_data.fill(np.nan)
    for i in range(len(x)):
        x_idx = np.argmin(np.abs(x_edges - x[i]))
        y_idx = np.argmin(np.abs(y_edges - y[i]))
        image_data[y_idx, x_idx] = z[i]
    
    im.set_array(image_data)
    time_text.set_text(f'Frame: {frame+1}')
    
    
    frame_filename = os.path.join(output_dir, f"frame_{frame+1:03d}.png")
    plt.savefig(frame_filename, dpi=300, facecolor='black', bbox_inches='tight', pad_inches=0.1)

    return im, time_text

animation = FuncAnimation(fig, update, frames=num_frames, 
                         interval=200, blit=True)

animation.save('hot-spot.gif', writer='pillow', 
              fps=5, dpi=600, savefig_kwargs={'facecolor':'black'})

plt.show()