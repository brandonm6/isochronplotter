import numpy as np
import pandas as pd

from variables import *

pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('expand_frame_repr', False)


def index2step(indexes, df):
    """
    indexes: array of start and stop indexes ex. [[2 5] [8 12]]
    df: dataframe to compare to
    """
    tmp = pd.DataFrame(indexes)
    tmp["step0"] = df["Step"].iloc[tmp.loc[:, tmp.columns[0]]].to_numpy()
    tmp["stepf"] = df["Step"].iloc[tmp.loc[:, tmp.columns[1]]].to_numpy()

    return list(tmp[["step0", "stepf"]].itertuples(index=False, name=None))


def check_overlap(a0, af, b0, bf):
    # returns True if there is overlap
    return (a0 <= bf) & (af >= b0)


def build_intervals(df):
    df["btm"] = df["Age"] - df["2SD"]
    df["top"] = df["Age"] + df["2SD"]

    df["cum btm"] = df["btm"]
    df["cum top"] = df["top"]

    def create_intervals(i):
        if check_overlap(df["cum btm"].iloc[i], df["cum top"].iloc[i], np.roll(df["cum btm"], 1)[i],
                         np.roll(df["cum top"], 1)[i]):
            # check of overlap with previous
            df["cum btm"].iloc[i] = np.maximum(df["cum btm"].iloc[i], np.roll(df["cum btm"], 1)[i])
            df["cum top"].iloc[i] = np.minimum(df["cum top"].iloc[i], np.roll(df["cum top"], 1)[i])

    vcreate_intervals = np.vectorize(create_intervals)
    vcreate_intervals(np.arange(len(df)))

    df["cum btm"], df["cum top"] = np.roll(df["cum btm"], 1), np.roll(df["cum top"], 1)

    return df


def group_overlaps(df):
    # returns groups of indexes that overlap
    # overlap = 1 where there is overlap with following x
    df["overlap"] = np.roll(np.where(check_overlap(df["btm"], df["top"], df["cum btm"], df["cum top"]), 1, 0), -1)

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


def check_three(df, iverbose):
    df = build_intervals(df)
    groups = group_overlaps(df)
    if iverbose:
        print(df)
        print(groups)
    try:
        ind = pd.DataFrame(groups, columns=["first", "last"])
    except ValueError:
        return []

    least_three = ind.loc[np.where((ind["last"] - ind["first"] + 1) >= 3)[0]].to_numpy()

    return least_three


class LocatePlateaus:
    def __init__(self, df, j_val, verbose=True, iverbose=False, level=None, start=None, stop=None):
        self.verbose = verbose
        self.iverbose = iverbose  # incremental verbose - prints dataframe for every incremented trapped argon
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

            least_three = check_three(tmp, self.iverbose)
            if len(least_three):
                to_check[1 / i] = least_three

        if self.verbose:
            print("At least three steps:")
            for trapped in to_check.keys():
                # convert index to steps for clarity purposes
                steps = index2step(to_check.get(trapped), self.df)
                print("%s: %s" % (trapped, steps))

        plateaus = {}  # have three steps and 50% Ar39
        for pot_trap in to_check.keys():
            indexes = to_check.get(pot_trap)
            tmp = pd.DataFrame(indexes)

            tmp["last %39"] = self.df["cum %39Ark"].iloc[tmp.loc[:, tmp.columns[1]]].to_numpy()
            tmp["first %39"] = self.df["cum %39Ark"].iloc[tmp.loc[:, tmp.columns[0]]].to_numpy()
            tmp["cum %39ark diff"] = tmp["last %39"] - tmp["first %39"]

            ind = np.where((tmp["cum %39ark diff"] > 50))[0]

            if len(ind):  # if there are any plateaus
                # convert index units to step units
                keep = tmp.iloc[ind].to_numpy()
                steps = index2step(keep, self.df)
                plateaus[pot_trap] = steps

        self.plateaus = plateaus
