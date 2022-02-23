import os
import sys
import time

import matplotlib.pyplot as plt

from plotter.plotter import Plotter

dir_path = os.path.dirname(os.path.realpath("isochronplotter"))
csv_path = os.path.join(dir_path, "csvs", "CIA dense 197'.csv")

t0 = time.time()
iso = Plotter(csv_path, omit=None, verbose=False, iverbose=False)
iso.plot()
tf = time.time()
print("\nrun time:", tf-t0)

plt.show()
