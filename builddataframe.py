import pandas as pd
import numpy as np


class BuildDataframe:
    # unpacks and builds dataframe

    def __init__(self, csv_path, omit=None):
        self.orig_df = pd.read_csv(csv_path)

        # run ids that were removed from the data set (contained NaN values)
        self.force_removed = None
        self.j_val = None
        # main df with corrections
        self.main_df = self.build_main()
        # sorted list of Run ID#s
        self.splits = self.find_splits()

        if omit is not None:
            self.omit = set(np.char.upper(omit))
            self.update_omit()

    def build_main(self):
        start = np.where(self.orig_df[self.orig_df.columns[0]] == "Run_ID")[0][0]
        # drop empty rows before the first run id
        self.orig_df = self.orig_df.rename(columns=self.orig_df.iloc[start]).drop(
            self.orig_df.index[:start + 1]).reset_index(drop=True)
        # drop rows after the last run id
        self.orig_df = self.orig_df.iloc[:np.where(self.orig_df["Run_ID"].fillna("") == "")[0][0]].reset_index(
            drop=True)

        self.j_val = float(self.orig_df["J"].iloc[0])
        self.orig_df["Step"] = np.arange(len(self.orig_df)) + 1
        self.orig_df["Status"] = "OK"

        # grab select columns
        grabbed = self.orig_df[
            ["Step", "Run_ID", "Status", "Isoch_39_Over_40", "Pct_i39_Over_40_Er", "Isoch_36_Over_40",
             "Pct_i36_Over_40_Er", "Correl_36_Over_39", "Ar36_Over_Ar39", "Ar39_Moles", "Age_Er"]]

        # drop NaN values
        self.force_removed = list(grabbed["Run_ID"][set(np.where(grabbed.isna())[0])])
        grabbed = grabbed.dropna().reset_index(drop=True)

        # convert string values to floats
        grabbed = pd.concat([grabbed.iloc[:, :3], grabbed.iloc[:, 3:].astype(float)], axis=1)

        grabbed["Pct_i39_Over_40_Er"] = (grabbed["Pct_i39_Over_40_Er"] / 100) * grabbed["Isoch_39_Over_40"]
        grabbed["Pct_i36_Over_40_Er"] = (grabbed["Pct_i36_Over_40_Er"] / 100) * grabbed["Isoch_36_Over_40"]

        # rename columns for convenience
        grabbed.columns = ["Step", "Run ID", "Status", "39Ar/40Ar", "39Ar/40Ar er", "36Ar/40Ar", "36Ar/40Ar er",
                           "corr 36/39", "36Ar/39Ar", "Ar39 moles", "Age er"]

        self.step_key = grabbed[["Step", "Run ID"]]

        return grabbed

    def find_splits(self):
        # find splits in the sample
        id_nums = sorted({run_id[:run_id.find(next(filter(str.isalpha, run_id)))] for run_id in self.main_df["Run ID"]})

        pct_39ark = np.array([])
        for num in id_nums:
            # build splits dataframes to get cum %39Ar
            steps = np.where(self.main_df["Run ID"].str.contains(num))[0]
            tmp = self.main_df.iloc[int(steps[0]):int(steps[-1]) + 1]
            pct_39ark = np.append(pct_39ark, ((tmp["Ar39 moles"] / tmp["Ar39 moles"].sum()) * 100).cumsum().to_numpy())

        # updates main dataframe with cum %39Ark column
        self.main_df["cum %39Ark"] = pct_39ark

        return list(id_nums)

    def update_omit(self):
        # updates status of omitted letters in main_df
        for i in self.main_df.index:
            # check if run id letter is in omit
            if self.main_df["Run ID"].iloc[i][-1] in self.omit:
                self.main_df["Status"].iloc[i] = "Omitted"
