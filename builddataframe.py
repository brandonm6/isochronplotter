import pandas as pd
import numpy as np


class BuildDataframe:
    # unpacks and builds dataframe

    def __init__(self, csv_path, omit):
        self.omit = omit

        # run ids removed from the data set (NaN values)
        self.force_removed = None
        self.j_val = None
        # dataframe of steps that correspond with their run ids
        self.step_key = None

        self.orig_df = pd.read_csv(csv_path)

        # main df with corrected values
        self.main_df = self.build_main()
        # list of Run ID#s
        self.splits = self.find_splits()

        self.update_omit()

    def build_main(self):
        # remove empty rows and correct column names
        start = np.where(self.orig_df[self.orig_df.columns[0]] == "Run_ID")[0][0]
        self.orig_df = self.orig_df.rename(columns=self.orig_df.iloc[start]).drop(
            self.orig_df.index[:start + 1]).reset_index(drop=True)
        stop = np.where(self.orig_df["Run_ID"].fillna("") == "")[0][0]
        self.orig_df = self.orig_df.iloc[:stop].reset_index(drop=True)

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
        id_nums = np.sort(
            list(set([run_id[:run_id.find(next(filter(str.isalpha, run_id)))] for run_id in self.main_df["Run ID"]])))

        pct_39ark = np.array([])
        for num in id_nums:
            # builds splits dataframes to get cum %39Ar
            steps = np.where(self.main_df["Run ID"].str.contains(num))[0]
            tmp = self.main_df.iloc[int(steps[0]):int(steps[-1]) + 1]
            pct_39ark = np.append(pct_39ark, ((tmp["Ar39 moles"] / tmp["Ar39 moles"].sum()) * 100).cumsum().to_numpy())

        # updates main dataframe with cum %39Ark column
        self.main_df["cum %39Ark"] = pct_39ark

        return list(id_nums)

    def update_omit(self):
        # updates status of omitted letters in main_df
        for num in self.main_df["Run ID"]:
            # check if run id letter is in omit
            if num[-1] in self.omit:
                self.main_df["Status"].iloc[self.main_df["Run ID"][self.main_df["Run ID"] == num].index] = "Omitted"
