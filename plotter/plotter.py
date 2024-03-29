import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.offsetbox import AnchoredText
import numpy as np
import itertools

from constants import *

from builddataframe.builddataframe import BuildDataframe
from locateplateaus.locateplateaus import LocatePlateaus

from plotter.mahon import Mahon
from plotter.helpers import *


class Plotter:
    def __init__(self, csv_path, omit=None, save=True, verbose=True, iverbose=False):
        """

        :param str csv_path: path to full raw run data file (from Mass Spec) in csv format
        :param list omit: list of letters or numbers to remove all of (ex. all As, all 2s, etc.)
        :param save: saves plots as pngs to isochronplotter/outputs
        :param verbose: prints three consecutive steps and plateaus found by locate_plateaus
        :param iverbose: prints dataframe for every incremented trapped argon along with groups of overlapping steps
        """
        self.output_path = os.path.join(os.path.dirname(os.path.realpath("isochronplotter")), "outputs")
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
        if omit is None:
            omit = []
        self.save = save
        self.verbose = verbose
        self.iverbose = iverbose
        self.name = csv_path[csv_path.rfind("/") + 1:csv_path.find(".csv")]
        unpacked = BuildDataframe(csv_path, omit)

        # Some variables to call
        self.force_removed = unpacked.force_removed  # removed run ids that had NaN values
        self.j_val = unpacked.j_val
        self.df = unpacked.main_df
        self.splits = unpacked.splits  # sorted list of run id#'s

    def make_plots(self, name):
        ax = plt.gca()

        # plot settings to match Turrin's
        set_size(9, 5.25, ax)

        ax.set_xticks(np.arange(0, .09 + .015, step=.015))
        ax.set_yticks(np.arange(-.0030, .0055, step=0.0005))
        ax.set_xlim(left=-.005, right=.09)
        ax.set_ylim(top=.0055)

        ax.annotate("Air", (ax.get_xlim()[0], 1 / ATM_ARGON), horizontalalignment='right')
        ax.set_xlabel("$^{39}$Ar/$^{40}$Ar")
        ax.set_ylabel("$^{36}$Ar/$^{40}$Ar")

        df = self.df.iloc[np.where((self.df["Run ID"].str[:-1] == name))[0]].reset_index(drop=True)

        aspect_ratio = get_aspect(ax)
        df["covar xy"] = df["corr 36/39"] * df["39Ar/40Ar er"] * df["36Ar/40Ar er"]
        df["ang"] = 0.5 * (np.arctan(
            (1 / aspect_ratio) * ((2 * df["covar xy"]) / ((df["39Ar/40Ar er"] ** 2) - (df["36Ar/40Ar er"] ** 2))))) * (
                            180 / np.pi)

        ok_df = df.iloc[np.where(df["Status"] == "OK")[0]].reset_index(drop=True)

        if ok_df.empty:  # may be because all the run ids were ommitted
            print(f"{name} is empty.")
            plt.clf()
            return

        plateaus = LocatePlateaus(ok_df, self.j_val, verbose=self.verbose, iverbose=self.iverbose).plateaus
        if self.verbose:
            print("Plateaus:", plateaus)

        # condense plateaus to create subsets (removes repeats and complete overlaps)
        if len(plateaus):  # if there are plateaus
            plat_status = "Plateaus found"
            condensed = list(set(itertools.chain.from_iterable(plateaus.values())))
            subsets = remove_complete_overlaps(sorted(condensed))
        else:
            plat_status = "No plateaus"
            subsets = [(ok_df["Step"].iloc[0], ok_df["Step"].iloc[-1])]

        for subset in subsets:
            # +1 for stop index b/c stop for .iloc is length-1
            sub_df = ok_df.iloc[step2index(ok_df, subset[0]):step2index(ok_df, subset[1]) + 1].reset_index(drop=True)

            # draw ellipse
            for i in np.arange(len(df)):
                covar_mat = np.array([[df["39Ar/40Ar er"].iloc[i] ** 2, df["covar xy"].iloc[i]],
                                      [df["covar xy"].iloc[i], df["36Ar/40Ar er"].iloc[i] ** 2]])
                eigenvalues = np.linalg.eigvals(covar_mat)
                if df["39Ar/40Ar er"].iloc[i] > df["36Ar/40Ar er"].iloc[i]:
                    x = np.sqrt(np.max(eigenvalues))
                    y = np.sqrt(np.min(eigenvalues))
                else:
                    x = np.sqrt(np.min(eigenvalues))
                    y = np.sqrt(np.max(eigenvalues))

                if df["Step"].iloc[i] in sub_df["Step"].to_numpy():
                    color = "C0"
                else:
                    color = "C4"

                ax.add_artist(
                    Ellipse(xy=(df["39Ar/40Ar"].iloc[i], df["36Ar/40Ar"].iloc[i]), width=(x * 2), height=(y * 2), 
                            angle=df["ang"].iloc[i], color=color))
                ax.annotate(df["Step"][i], (df["39Ar/40Ar"][i], df["36Ar/40Ar"][i]), horizontalalignment="center")

            # create tmp df that is readable by Mahon
            tmp = sub_df[["39Ar/40Ar", "39Ar/40Ar er", "36Ar/40Ar", "36Ar/40Ar er", "corr 36/39"]]
            tmp.columns = ["x", "x unc", "y", "y unc", "corr"]

            york = Mahon(tmp)
            xinter, xinterunc = york.xinter, york.xinterunc
            yinter, yinterunc = york.yinter, york.yinterunc
            mswd = york.mswd

            # convert yinterunc of 36Ar/40Ar to 40Ar/36Ar
            yinterunc = (yinterunc / yinter) * (1 / yinter)

            age = (1 / TOTAL_DECAY_40K) * np.log(((1 / xinter) * self.j_val) + 1) / 1000000  # age in Ma
            ax.plot([ax.get_xlim()[0], xinter], [yinter, 0], label=(
                    "Age = " + str(round(age, 3)) + " ± " + str(round(xinterunc, 3)) + " Ma" +
                    "\n$^{40}$Ar/$^{36}$Ar = " + str(round((1/yinter), 1)) + " ± " + f'{float(f"{yinterunc:.2g}"):g}' +
                    "\nMSWD = " + str(round(mswd, 1)) +
                    ", Steps: " + str(subset[0]) + "-" + str(subset[1])))
            ax.add_artist(ax.legend(loc='lower left'))

        ax.add_artist(AnchoredText(plat_status, loc="lower right", frameon=False, pad=0))
        ax.set_title(name)

        if self.save:
            plt.savefig(self.output_path + "/" + name + '.png')

    def plot(self):
        # create plots for splits individually
        for split in self.splits:
            plt.figure(self.splits.index(split))
            if self.verbose:
                print("\n\nRun ID :", split)
            self.make_plots(split)
