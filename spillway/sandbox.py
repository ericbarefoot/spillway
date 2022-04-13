import numpy as np
from matplotlib import pyplot as plt

import spillway as sp
from importlib import reload

plt.ion() 
# reload(sp)

self = sp.Spillway()
self.set_defaults()
self.set_x(np.arange(-4.5, 6.5, 0.01))
self.set_eta(0.1, 100, 0.0005)

plt.plot(self.x, self.eta)
plt.ylim([0, 0.16])


pip install -e .