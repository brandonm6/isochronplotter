import numpy as np
import pandas as pd

from variables import *

pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('expand_frame_repr', False)


def index2step(indexes, df):
    # indexes: array of start and stop indexes ex. [[2 5] [8 12]]
    # df: dataframe to compare to
    return [(df["Step"].iloc[i[0]], df["Step"].iloc[i[1]]) for i in indexes]


def check_overlap(a0, af, b0, bf):
    # returns True if there is overlap
    return (a0 <= bf) & (af >= b0)


def check_three(df, iverbose=False):
    df["btm"] = df["Age"] - df["2SD"]
    df["top"] = df["Age"] + df["2SD"]

    df["cum btm"] = df["btm"]
    df["cum top"] = df["top"]

    for i in np.arange(len(df)):
        if check_overlap(df["cum btm"].iloc[i], df["cum top"].iloc[i], np.roll(df["cum btm"], 1)[i],
                         np.roll(df["cum top"], 1)[i]):
            df["cum btm"].iloc[i] = np.maximum(df["cum btm"].iloc[i], np.roll(df["cum btm"], 1)[i])
            df["cum top"].iloc[i] = np.minimum(df["cum top"].iloc[i], np.roll(df["cum top"], 1)[i])

    df["cum btm"], df["cum top"] = np.roll(df["cum btm"], 1), np.roll(df["cum top"], 1)

    # overlap = 1 where there is overlap with following x
    df["overlap"] = np.roll(np.where(check_overlap(df["btm"], df["top"], df["cum btm"], df["cum top"]), 1, 0), -1)
    df.iloc[-1, df.columns.get_loc("overlap")] = 0  # to prevent last step from overlapping with first step

    if iverbose:
        print(df)

    indexes = np.nonzero(df["overlap"].to_numpy())[0]

    try:
        # groups indexes that overlap  ex. [[2, 6],] - index 2 through 6 (inclusive) are overlapping
        groups = np.array([(s[0], s[-1] + 1) for s in np.split(indexes, np.where(np.diff(indexes) != 1)[0] + 1)])
        if iverbose:
            print(groups)
        return np.array([group for group in groups if (group[1] - group[0]) + 1 >= 3])
    except IndexError:  # indexes is empty (no overlaps at all)
        return indexes  # empty numpy array


class LocatePlateaus:
    def __init__(self, df, j_val, verbose=True, iverbose=False, level=0.01, start=(1 / atm_argon), stop=None):
        self.verbose = verbose
        self.iverbose = iverbose  # incremental verbose - prints dataframe for every incremented trapped argon
        self.df = df
        self.level = level
        self.start = start
        if stop is None:
            stop = self.df["36Ar/40Ar"].min()
        self.stop = stop
        self.j_val = j_val

        # some variables you can call
        self.plateaus = None
        self.find_plateaus()
        """
        returns dictionary of trapped argon values as keys and values that are
        lists of tuples where each tuple is in the format of (start_step, end_step)
        """

    def find_plateaus(self):
        num_incs = int(np.log(self.stop / self.start) / np.log((1 - self.level)))
        incs = self.start * (1 / (1 - self.level)) ** np.arange(0, -num_incs, -1)

        to_check = {}  # have at least three steps, format is trapped:
        for i in incs:
            tmp = self.df[["Step"]]
            ar40star_div_ar39 = ((1 / self.df["39Ar/40Ar"]) - ((1 / i) * self.df["36Ar/39Ar"]))
            tmp["Age"] = ((1 / total_decay_40k) * np.log((ar40star_div_ar39 * self.j_val) + 1)) / 1000000

            # using error for atmospheric argon since we dont have measurements to calculate the actual
            # uncertainty for each age on the incremented level of excess argon
            tmp["2SD"] = self.df["Age er"] * 2
            if self.iverbose:
                print("Trapped:", (1 / i))

            three_consec = check_three(tmp, self.iverbose)
            if three_consec.size:
                to_check[1 / i] = three_consec

        if self.verbose:
            print("At least three consecutive steps:")
            # convert index to steps for easier reading
            print("".join([f"{trap}: {index2step(to_check.get(trap), self.df)}\n" for trap in to_check]))

        plateaus = {}  # have three consecutive steps and 50% Ar39
        for trap in to_check:
            least_fifty = np.array([i for i in to_check.get(trap) if
                                    (self.df["cum %39Ark"].iloc[i[1]] - self.df["cum %39Ark"].iloc[i[0]]) >= 50])
            if least_fifty.size:  # if there are any plateaus
                # convert index units to step units
                plateaus[trap] = index2step(least_fifty, self.df)

        self.plateaus = plateaus
