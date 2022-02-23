import numpy as np
import pandas as pd


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
        # group indexes that overlap  ex. [[2, 6],] - index 2 through 6 (inclusive) are overlapping
        groups = np.array([(s[0], s[-1] + 1) for s in np.split(indexes, np.where(np.diff(indexes) != 1)[0] + 1)])
        if iverbose:
            print(groups)
        return np.array([group for group in groups if (group[1] - group[0]) + 1 >= 3])
    except IndexError:  # indexes is empty (no overlaps at all)
        return indexes  # empty numpy array
