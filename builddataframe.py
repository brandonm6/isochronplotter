import pandas as pd
import numpy as np


class BuildDataframe:
    # unpacks and builds dataframe

    def __init__(self, csv_path, remove, verbose=True):
        self.verbose = verbose
        self.removed = remove
        self.j_val = None
        self.step_key = None
        self.orig_df = pd.read_csv(csv_path)
        # main df with corrected values
        self.main_df = self.build_main()
        # dictionary of Run ID#s and dataframes for each #
        self.splits = self.find_splits()

    def build_main(self):
        # remove empty rows and correct column names
        start = np.where(self.orig_df[self.orig_df.columns[0]] == "Run_ID")[0][0]
        self.orig_df = self.orig_df.rename(columns=self.orig_df.iloc[start]).drop(
            self.orig_df.index[:start + 1]).reset_index(drop=True)
        stop = np.where(self.orig_df["Run_ID"].fillna("") == "")[0][0]
        self.orig_df = self.orig_df.iloc[:stop].reset_index(drop=True)

        self.j_val = float(self.orig_df["J"].iloc[0])
        self.orig_df["Step"] = np.arange(len(self.orig_df)) + 1
        self.orig_df["Ltr"] = self.orig_df["Run_ID"].str[-1]

        # grab select columns
        grabbed = self.orig_df[
            ["Step", "Run_ID", "Ltr", "Isoch_39_Over_40", "Pct_i39_Over_40_Er", "Isoch_36_Over_40",
             "Pct_i36_Over_40_Er",
             "Correl_36_Over_39", "Ar36_Over_Ar39", "Ar39_Moles"]]

        # drop self.removed values and NaN values
        grabbed = grabbed[~grabbed["Ltr"].isin(self.removed)].reset_index(drop=True)
        self.removed += list(grabbed["Run_ID"][set(np.where(grabbed.isna())[0])])
        grabbed = grabbed.dropna().reset_index(drop=True)

        # convert string values to floats
        grabbed = pd.concat([grabbed.iloc[:, :3], grabbed.iloc[:, 3:].astype(float)], axis=1)

        grabbed["Pct_i39_Over_40_Er"] = (grabbed["Pct_i39_Over_40_Er"] / 100) * grabbed["Isoch_39_Over_40"]
        grabbed["Pct_i36_Over_40_Er"] = (grabbed["Pct_i36_Over_40_Er"] / 100) * grabbed["Isoch_36_Over_40"]

        # rename columns for convenience
        grabbed.columns = ["Step", "Run ID", "Ltr", "39Ar/40Ar", "39Ar/40Ar er", "36Ar/40Ar", "36Ar/40Ar er",
                           "corr 36/39", "36Ar/39Ar", "Ar39 moles"]

        self.step_key = grabbed[["Step", "Run ID"]]

        return grabbed

    def update_main(self, splits):
        """
        updates main dataframe by adding 39Ar cumulative percentage column (%39Ark)
        if there are splits in the sample
        """
        tmp = pd.DataFrame()
        for split in np.sort(list(splits.keys())):
            tmp = pd.concat([tmp, splits.get(split)["%39Ark"]])
        self.main_df["%39Ark"] = tmp.reset_index(drop=True)

    def find_splits(self):
        # check if there are splits in the sample
        id_nums = set([run_id[:run_id.find(next(filter(str.isalpha, run_id)))] for run_id in self.main_df["Run ID"]])
        if len(id_nums) == len(self.main_df["Run ID"]):
            if self.verbose:
                print("The sample contains no splits.")
            return

        splits = {}
        for num in id_nums:
            # builds splits dataframes
            steps = np.where(self.main_df["Run ID"].str.contains(num))[0]
            df = self.main_df.iloc[int(steps[0]):int(steps[-1]) + 1]

            df["%39Ark"] = ((df["Ar39 moles"] / df["Ar39 moles"].sum()) * 100).cumsum()
            splits[num] = df.drop("Ar39 moles", axis=1).reset_index(drop=True)

        self.update_main(splits)

        return splits
