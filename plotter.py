import os
import matplotlib.pyplot as plt
import pandas as pd
from helpers import *

csvs_dir = os.path.join(os.getcwd(), "csvs")


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
        self.x_ticks = x_ticks
        self.y_ticks = y_ticks

    def cheat_steps(self):
        df = self.df
        if self.name == "P K-feldspar":
            df_1 = Isochron.SectionSteps(df, 1, 10, False)
            df_2 = Isochron.SectionSteps(df, 11, 18, True)
            separated = [df_1, df_2]
        elif self.name == "B K-feldspar":
            df_1 = Isochron.SectionSteps(df, 1, 10, True)  # 2, 10
            df_2 = Isochron.SectionSteps(df, 11, 18, True)
            separated = [df_1, df_2]
        else:
            if self.name == "Buckskin Peak Biotite":
                self.x_ticks = np.arange(0.053, .062, step=.001)
                self.y_ticks = np.arange(0, .00012, step=.00002)
            df_1 = Isochron.SectionSteps(df, has_slope=False, every=True)
            separated = [df_1]
        return separated

    def organize(self):
        self.orig_df["Step"] = np.arange(len(self.orig_df)) + 1
        self.orig_df["39Ar/40Ar"] = 1 / self.orig_df["40Ar/39Ar"]
        self.orig_df["36Ar/40Ar"] = self.orig_df["36Ar/39Ar"] / self.orig_df["40Ar/39Ar"]
        self.orig_df["36Ar/40Ar"] = self.orig_df["36Ar/39Ar"] / self.orig_df["40Ar/39Ar"]
        new_df = self.orig_df[["Step", "36Ar/40Ar", "39Ar/40Ar", "%39Ark", "Age", "1SD"]]

        return new_df

    def plot(self):
        atm_argon = 295.5  # 298.6
        total_decay_40k = 5.54e-10

        plt.plot()
        ax = plt.gca()

        sections = self.cheat_steps()
        for section in sections:
            df = self.df.iloc[section.start:section.stop]
            x, y = lambda df: df["39Ar/40Ar"], lambda df: df["36Ar/40Ar"]

            x_outliers = locate_outliers(x(df))
            y_outliers = locate_outliers(y(df))
            # if section.has_slope:  # so initial nonhomogeneous scattered data not discarded
            df = df.drop((x_outliers + y_outliers)).reset_index(drop=True)

            infl = locate_influential(x(df), y(df))
            # df = df.drop(check_removal(infl))

            xs, ys = x(df), y(df)
            ax.scatter(xs, ys)
            for i in df["Step"]:
                pos = df.index[df["Step"] == i].tolist()[0]
                plt.annotate(i, (xs[pos], ys[pos]), horizontalalignment="center")
            if section.has_slope:
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


skip = ["P K-feldspar.csv", "Buckskin Peak Biotite.csv"]
for file in os.listdir(csvs_dir):
    if file not in skip:  # for debugging purposes
        print(file[:file.find(".csv")])
        csv_path = os.path.join(csvs_dir, file)
        Isochron(csv_path).plot()
        print("\n\n")

