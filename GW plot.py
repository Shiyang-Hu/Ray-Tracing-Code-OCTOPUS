import numpy as np
import matplotlib.pyplot as plt
t, h1, h2 = np.loadtxt("GW emission.dat", unpack=True)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
ax1.plot(t, h1, color='blue', linewidth=1.5)
ax1.set_ylabel('plus', fontsize=12)
ax1.grid(True, alpha=0.3)
ax2.plot(t, h2, color='red', linewidth=1.5)
ax2.set_xlabel('t', fontsize=12)
ax2.set_ylabel('cross', fontsize=12)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('GW.png', dpi=200, bbox_inches='tight')
plt.show()
