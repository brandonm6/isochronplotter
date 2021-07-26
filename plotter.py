import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from variables import *
from helpers import *
from locate_plateaus import find_plateaus

# for debugging
import sys
import time

csv_dir = os.path.join(os.getcwd(), "csvs")

class Isochron:
    class SectionSteps:
        def __init__(self, df, start=None, stop=None, has_slope=False, every=None):
            self.df = df
            if every:
                self.start = 0
                self.stop = None
            else:
                self.start = self.try_steps(start, 1)
                self.stop = self.try_steps(stop, -1) + 1
            self.has_slope = has_slope

        def try_steps(self, pos, sign):
            result = None
            while result is None:
                try:
                    result = self.df.index[self.df["Step"] == pos].tolist()[0]
                except IndexError:
                    pos += sign
            return result

    def __init__(self, csv, x_ticks=np.arange(0, .35, step=.05), y_ticks=np.arange(0, .005, step=0.001)):
        # x_tics and y_ticks used for cheating purposes
        self.orig_df = pd.read_csv(csv)
        self.j_val = float(self.orig_df.columns[-1])
        self.name = csv[csv.rfind("/") + 1:csv.find(".csv")]
        self.df = self.organize()
        self.removed_steps = []

        self.x_ticks = x_ticks
        self.y_ticks = y_ticks

    def organize(self):
        self.orig_df["Step"] = np.arange(len(self.orig_df)) + 1
        self.orig_df["39Ar/40Ar"] = 1 / self.orig_df["40Ar/39Ar"]
        self.orig_df["36Ar/40Ar"] = self.orig_df["36Ar/39Ar"] / self.orig_df["40Ar/39Ar"]
        self.orig_df["36Ar/40Ar"] = self.orig_df["36Ar/39Ar"] / self.orig_df["40Ar/39Ar"]
        new_df = self.orig_df[["Step", "36Ar/40Ar", "39Ar/40Ar", "36Ar/39Ar", "%39Ark", "1SD"]]

        return new_df

    def plot(self):

        plt.plot()
        ax = plt.gca()

        # plateaus = find_plateaus(self.df, self.j_val)
        """
        add filtering mechanism to choose which plateau is the best for each trapped argon level
        """

        for subset in subsets:

            x, y = lambda df: df["39Ar/40Ar"], lambda df: df["36Ar/40Ar"]

            # outliers = locate_outliers_resid(x(df), y(df))
            # self.removed_steps += df["Step"].iloc[outliers].tolist()
            # df = df.drop(outliers).reset_index(drop=True)

            xs, ys = x(df), y(df)
            ax.scatter(xs, ys)
            for i in df["Step"]:
                pos = df.index[df["Step"] == i].tolist()[0]
                plt.annotate(i, (xs[pos], ys[pos]), horizontalalignment="center")

            [m, b] = np.polyfit(xs, ys, 1)
            age = (1 / total_decay_40k) * np.log(((1 / (-b / m)) * self.j_val) + 1) / 1000000  # age in Ma
            plt.plot([0, (-b / m)], [b, 0], label=(
                    "$^{40}$Ar/$^{36}$Ar = " + str(round(1 / b)) + "\n Age = " + str(round(age, 1)) + " Ma"))
            plt.legend(loc="upper right")

        plt.margins(x=0, y=0)
        plt.title(self.name)
        plt.annotate("Air", (0, 1 / atm_argon), horizontalalignment='right')
        plt.xlabel("$^{39}$Ar/$^{40}$Ar")
        plt.ylabel("$^{36}$Ar/$^{40}$Ar")
        plt.xticks(self.x_ticks)
        plt.yticks(self.y_ticks)
        plt.show()



skip = []  #"P K-feldspar.csv","B K-feldspar.csv", "Buckskin Peak Biotite.csv"
for file in os.listdir(csv_dir):
    if file not in skip:  # for debugging purposes

        # title = file[:file.find(".csv")]
        # orig_stdout = sys.stdout
        # f = open((os.getcwd() + "/outputs/" + title + ".txt"), 'w')
        # sys.stdout = f

        print(file[:file.find(".csv")])
        csv_path = os.path.join(csv_dir, file)
        Isochron(csv_path).plot()
        print("\n\n")

        # sys.stdout = orig_stdout
        # f.close()


