import os
from plotter.plotter import Plotter

file_name = "CIA dense 197'.csv"
dir_path = os.path.dirname(os.path.realpath("isochronplotter"))
csv_path = os.path.join(dir_path, "csvs", file_name)

iso = Plotter(csv_path, omit=[], verbose=False, iverbose=False)
iso.plot()
