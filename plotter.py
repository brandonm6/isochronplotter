import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.offsetbox import AnchoredText
import numpy as np
import itertools

from variables import *
from builddataframe import BuildDataframe
from locateplateaus import LocatePlateaus
from mahon import Mahon


def set_size(w, h, ax):
    """ w, h: width, height in inches """
    left = ax.figure.subplotpars.left
    right = ax.figure.subplotpars.right
    top = ax.figure.subplotpars.top
    btm = ax.figure.subplotpars.bottom
    figw = float(w) / (right - left)
    figh = float(h) / (top - btm)
    ax.figure.set_size_inches(figw, figh)


def get_aspect(ax):
    ylower, yupper = ax.get_ylim()
    xlower, xupper = ax.get_xlim()
    data_ratio = (yupper - ylower) / (xupper - xlower)

    # total figure size
    fig_w, fig_h = ax.get_figure().get_size_inches()
    # axis size on figure
    _, _, w, h = ax.get_position().bounds
    disp_ratio = (fig_h * h) / (fig_w * w)

    return data_ratio / disp_ratio


def step2index(df, step):
    return df.index[df['Step'] == step].tolist()[0]


def remove_complete_overlaps(lst):
    """
    returns list with complete overlaps removed
    lst: sorted(list)
    """
    if len(lst) <= 1:
        return lst
    # first interval falls inside second
    elif (lst[0][0] >= lst[1][0]) and (lst[0][1] <= lst[1][1]):
        return remove_complete_overlaps(lst[1:])
    # second interval falls inside first
    elif (lst[0][0] <= lst[1][0]) and (lst[0][1] >= lst[1][1]):
        return remove_complete_overlaps([lst[0]] + lst[2:])
    else:
        return [lst[0]] + remove_complete_overlaps(lst[1:])


class Plotter:
    def __init__(self, csv_path, omit=None, verbose=True, iverbose=False):
        """
        csv_path: path to full raw run data file
        omit: list of letters to remove (removes all  ex. all As, all Bs, etc.)
        """
        if omit is None:
            omit = []
        self.verbose = verbose
        self.iverbose = iverbose
        self.name = csv_path[csv_path.rfind("/") + 1:csv_path.find(".csv")]
        unpacked = BuildDataframe(csv_path, omit)
        self.force_removed = unpacked.force_removed  # removed run ids that had NaN values
        self.j_val = unpacked.j_val
        self.step_key = unpacked.step_key
        self.df = unpacked.main_df
        self.splits = unpacked.splits

    def make_plots(self, ax, name):
        # plot settings to match Turrin's
        ax.set_xticks(np.arange(0, .09 + .015, step=.015))
        ax.set_yticks(np.arange(-.0030, .0055, step=0.0005))
        ax.set_xlim(left=-.005, right=.09)
        ax.set_ylim(top=.0055)

        ax.annotate("Air", (ax.get_xlim()[0], 1 / atm_argon), horizontalalignment='right')
        ax.set_xlabel("$^{39}$Ar/$^{40}$Ar")
        ax.set_ylabel("$^{36}$Ar/$^{40}$Ar")

        aspect_ratio = get_aspect(ax)

        df = self.df.iloc[np.where((self.df["Run ID"].str[:-1] == name))[0]].reset_index(drop=True)

        df["covar xy"] = df["corr 36/39"] * df["39Ar/40Ar er"] * df["36Ar/40Ar er"]
        df["ang"] = 0.5 * (np.arctan(
            (1 / aspect_ratio) * ((2 * df["covar xy"]) / ((df["39Ar/40Ar er"] ** 2) - (df["36Ar/40Ar er"] ** 2))))) * (
                            180 / np.pi)

        ok_df = df.iloc[np.where(df["Status"] == "OK")[0]].reset_index(drop=True)

        plateaus = LocatePlateaus(ok_df, self.j_val, verbose=self.verbose, iverbose=self.iverbose).plateaus
        if self.verbose:
            print("Plateaus:", plateaus)

        # condense plateaus to create subsets (removes repeats and complete overlaps)
        subsets = []
        if len(plateaus):  # if there are plateaus
            plat_status = "Plateaus found"
            condensed = list(set(itertools.chain.from_iterable(plateaus.values())))
            condensed = remove_complete_overlaps(sorted(condensed))
            for section in condensed:
                subsets.append((section[0], section[1]))
        else:
            plat_status = "No plateaus"
            subsets.append((ok_df["Step"].iloc[0], ok_df["Step"].iloc[-1]))

        for subset in subsets:
            # +1 for stop index b/c stop for .iloc is length-1
            sub_df = ok_df.iloc[step2index(ok_df, subset[0]):step2index(ok_df, subset[1]) + 1].reset_index(drop=True)

            def draw_ellipse(i):
                covar_mat = np.array([[df["39Ar/40Ar er"].iloc[i] ** 2, df["covar xy"].iloc[i]],
                                      [df["covar xy"].iloc[i], df["36Ar/40Ar er"].iloc[i] ** 2]])
                eigenvalues = np.linalg.eigvals(covar_mat)
                if df["39Ar/40Ar er"].iloc[i] > df["36Ar/40Ar er"].iloc[i]:
                    x = np.sqrt(max(eigenvalues))
                    y = np.sqrt(min(eigenvalues))
                else:
                    x = np.sqrt(min(eigenvalues))
                    y = np.sqrt(max(eigenvalues))

                if df["Step"].iloc[i] in sub_df["Step"].to_numpy():
                    color = "C0"
                else:
                    color = "C4"

                ax.add_artist(
                    Ellipse((df["39Ar/40Ar"].iloc[i], df["36Ar/40Ar"].iloc[i]), (x * 2), (y * 2), df["ang"].iloc[i],
                            color=color))
                ax.annotate(df["Step"][i], (df["39Ar/40Ar"][i], df["36Ar/40Ar"][i]), horizontalalignment="center")

            vdraw_ellipse = np.vectorize(draw_ellipse)
            vdraw_ellipse(np.arange(len(df)))

            # create tmp df that is readable by Mahon
            tmp = sub_df[["39Ar/40Ar", "39Ar/40Ar er", "36Ar/40Ar", "36Ar/40Ar er", "corr 36/39"]]
            tmp.columns = ["x", "x unc", "y", "y unc", "corr"]

            york = Mahon(tmp)
            xinter, xinterunc = york.xinter, york.xinterunc
            yinter, yinterunc = york.yinter, york.yinterunc
            mswd = york.mswd

            age = (1 / total_decay_40k) * np.log(((1 / xinter) * self.j_val) + 1) / 1000000  # age in Ma
            ax.plot([ax.get_xlim()[0], xinter], [yinter, 0], label=(
                    "Age = " + str(round(age, 3)) + " ± " + str(round(xinterunc, 3)) + " Ma" +
                    "\n$^{40}$Ar/$^{36}$Ar = " + str(round((1 / yinter), 1)) + " ± " + str(round(yinterunc, 1)) +
                    "\nMSWD = " + str(round(mswd, 1)) +
                    ", Steps: " + str(subset[0]) + "-" + str(subset[1])))
            ax.add_artist(ax.legend(loc='lower left'))

        ax.add_artist(AnchoredText(plat_status, loc="lower right", frameon=False, pad=0))

        ax.set_title(name)

    def plot(self):
        fig, axs = plt.subplots(len(self.splits), 1)

        # create plots for splits individually
        for split in self.splits:
            if self.verbose:
                print("\nRun ID :", split)
            ax = axs[self.splits.index(split)]
            self.make_plots(ax, split)

        set_size(22.5 / 2.54, ((13 / 2.54) * len(self.splits)), ax)
        plt.show()
