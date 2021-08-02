import numpy as np
import pandas as pd

from variables import *


pd.options.mode.chained_assignment = None  # default='warn'


def check_overlap(a0, af, b0, bf):
    # returns True if there is overlap
    return (a0 <= bf) & (af >= b0)


def build_intervals(df):
    df["btm"] = df["x"] - df["err"]
    df["top"] = df["x"] + df["err"]

    df["cur btm"] = df["btm"]
    df["cur top"] = df["top"]

    def create_intervals(i):
        if check_overlap(df["cur btm"].iloc[i], df["cur top"].iloc[i], np.roll(df["cur btm"], 1)[i],
                         np.roll(df["cur top"], 1)[i]):
            # check of overlap with previous
            df["cur btm"].iloc[i] = np.maximum(df["cur btm"].iloc[i], np.roll(df["cur btm"], 1)[i])
            df["cur top"].iloc[i] = np.minimum(df["cur top"].iloc[i], np.roll(df["cur top"], 1)[i])

    vcreate_intervals = np.vectorize(create_intervals)
    vcreate_intervals(np.arange(len(df)))

    df["cur btm"], df["cur top"] = np.roll(df["cur btm"], 1), np.roll(df["cur top"], 1)

    return df


def group_overlaps(df):
    # returns groups of indexes that overlap
    # overlap = 1 where there is overlap with following x
    df["overlap"] = np.roll(np.where(check_overlap(df["btm"], df["top"], df["cur btm"], df["cur top"]), 1, 0), -1)

    df.iloc[-1, df.columns.get_loc("overlap")] = 0  # because last step has no following step

    index = np.nonzero(df["overlap"].to_numpy())[0]

    if not len(index):  # if no overlaps at all (needed b/c grouped.apply fails - try to optimize w/o grouped.apply)
        return index

    grp = np.split(index, np.where(np.diff(index) != 1)[0] + 1)
    grouped = pd.DataFrame(grp)
    mod = grouped.apply(lambda x: sorted(x, key=pd.notnull), 1)
    mod = pd.DataFrame.from_dict(dict(zip(mod.index, mod.values))).T
    grouped["next"] = mod.iloc[:, -1:] + 1

    grouped = grouped[[grouped.columns[0], "next"]].to_numpy()

    return grouped


def check_three(df):
    df = build_intervals(df)
    groups = group_overlaps(df)
    try:
        ind = pd.DataFrame(groups, columns=["first", "last"])
    except ValueError:
        return []

    least_three = ind.loc[np.where((ind["last"] - ind["first"] + 1) >= 3)[0]].to_numpy()

    return least_three


class LocatePlateaus:
    def __init__(self, df, j_val, verbose=True, level=None, start=None, stop=None):
        self.verbose = verbose
        self.df = df
        if level is None:
            level = 0.01
        if start is None:
            start = 1 / atm_argon
        if stop is None:
            stop = self.df["36Ar/40Ar"].min()
        self.level = level
        self.start = start
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

        to_check = {}  # have at least three steps
        for i in incs:
            tmp = self.df[["Step"]]
            ar40star_div_ar39 = ((1 / self.df["39Ar/40Ar"]) - ((1 / i) * self.df["36Ar/39Ar"]))
            tmp["Age"] = ((1 / total_decay_40k) * np.log((ar40star_div_ar39 * self.j_val) + 1)) / 1000000

            # using 1% since we dont have measurements to calculate uncertainty of age
            tmp["1SD"] = tmp["Age"].multiply(.01)

            tmp.columns = ["Step", "x", "err"]
            least_three = check_three(tmp)
            if len(least_three):
                to_check[1 / i] = least_three

        if self.verbose:
            print("At least three steps:")
            for trapped in to_check.keys():
                print("%s: %s" % (trapped, to_check.get(trapped)))

        plateaus = {}  # have three steps and 50% Ar39
        for pot_trap in to_check.keys():
            indexes = to_check.get(pot_trap)
            tmp = pd.DataFrame(indexes)

            tmp["last %39"] = self.df["%39Ark"].iloc[tmp.loc[:, tmp.columns[1]]].to_numpy()
            tmp["first %39"] = self.df["%39Ark"].iloc[tmp.loc[:, tmp.columns[0]]].to_numpy()
            tmp["%39ark diff"] = tmp["last %39"] - tmp["first %39"]

            ind = np.where((tmp["%39ark diff"] > 50))[0]

            if len(ind):  # if there are any plateaus
                keep = tmp.iloc[ind]

                # convert index units to step
                keep["step0"] = self.df["Step"].iloc[keep.loc[:, keep.columns[0]]].to_numpy()
                keep["stepf"] = self.df["Step"].iloc[keep.loc[:, keep.columns[1]]].to_numpy()
                steps = list(keep[["step0", "stepf"]].itertuples(index=False, name=None))

                plateaus[pot_trap] = steps

        self.plateaus = plateaus
