import numpy as np
import pandas as pd

from variables import *


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
    '''
    returns groups of STEPS (not indexes) that overlap
    '''
    # overlap = 1 where there is overlap with the following x
    df["overlap"] = np.roll(np.where(check_overlap(df["btm"], df["top"], df["cur btm"], df["cur top"]), 1, 0), -1)

    df.iloc[-1, df.columns.get_loc("overlap")] = 0  # because last step has no following step

    index = np.nonzero(df["overlap"].to_numpy())[0]

    # make function return groups of steps that overlap rather than their indexes
    index += 1

    if not len(index):  # if no overlaps at all (needed b/c grouped.apply fails - try to stop using grouped.apply)
        return index

    grp = np.split(index, np.where(np.diff(index) != 1)[0] + 1)
    grouped = pd.DataFrame(grp)
    grouped = grouped.apply(lambda x: sorted(x, key=pd.notnull), 1)
    grouped = pd.DataFrame.from_dict(dict(zip(grouped.index, grouped.values))).T
    grouped["next"] = grouped.iloc[:, -1:] + 1

    grouped = grouped.to_numpy()

    return grouped

def check_three(x, err):
    df = build_intervals(pd.concat([x, err], axis=1, keys=["x", "err"]))
    groups = group_overlaps(df)

    ind = pd.DataFrame(groups)
    least_three = ind.loc[np.where(ind.count(axis=1) >= 3)[0]].to_numpy()

    return least_three


def find_plateaus(orig_df, j_val):
    level = 0.01
    start = 1 / atm_argon
    stop = orig_df["36Ar/40Ar"].min()

    num_incs = int(np.log(stop / start) / np.log((1 - level)))
    incs = start * (1 / (1 - level)) ** np.arange(0, -num_incs, -1)

    to_check = {}  # have at least three steps
    for i in incs:
        df = orig_df.copy()
        ar40star_div_ar39 = ((1 / df["39Ar/40Ar"]) - ((1/i) * df["36Ar/39Ar"]))
        df["Age"] = ((1 / total_decay_40k) * np.log((ar40star_div_ar39 * j_val) + 1)) / 1000000

        # because dont have measurements to calculate uncertainty of age
        df["1SD"] = df["1SD"].multiply(1.75)

        least_three = check_three(df["Age"], df["1SD"])
        if len(least_three):
            to_check[1 / i] = least_three

    plateaus = {}  # have three steps and 50% Ar39
    for pot_trap in to_check.keys():
        steps = to_check.get(pot_trap)
        tmp = pd.DataFrame(steps)

        # subtract 1 in .iloc because using steps not indexes (indexes = steps - 1)
        last_ind = orig_df["%39Ark"].iloc[tmp.loc[tmp.index, tmp.columns[-1]] - 1].to_numpy()
        first_ind = orig_df["%39Ark"].iloc[tmp.loc[tmp.index, tmp.columns[0]]-1].to_numpy()
        tmp["diff"] = last_ind - first_ind

        keep = np.where((tmp["diff"]) > 50)[0]
        if len(keep):
            plateaus[pot_trap] = steps[keep]

    return plateaus