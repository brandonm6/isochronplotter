import numpy as np
import pandas as pd

from constants import *
from locateplateaus.helpers import *


pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('expand_frame_repr', False)


class LocatePlateaus:
    def __init__(self, df, j_val, verbose=True, iverbose=False, start=(1 / ATM_ARGON), stop=None):
        self.verbose = verbose
        self.iverbose = iverbose  # incremental verbose - prints dataframe for every incremented trapped argon
        self.df = df
        self.start = start
        if stop is None:
            stop = self.df["36Ar/40Ar"].min()
        self.stop = stop
        self.j_val = j_val

        # variable to call get found plateaus
        self.plateaus = self.find_plateaus()

    def find_plateaus(self):
        """Finds plateaus where a plateau has at least three consecutive steps and at least 50% cumulative 39Ar

        :return: dictionary where keys are trapped argon values and values are lists of tuples with each tuple in the
            format of (start_step, end_step)
        """
        # gets final exponent of the increments assuming base of 1-level (default=.99)
        final_exp = int((np.log(self.stop / self.start) / np.log(MULTIPLIER)) + 0.5)
        # builds array of values to increment over where each is MULTIPLIER times the previous
        incs = self.start * (MULTIPLIER ** np.arange(0, final_exp + 1, 1))

        to_check = {}  # have at least three steps, format is trapped:
        for i in incs:
            tmp = self.df[["Step"]]
            ar40star_div_ar39 = ((1 / self.df["39Ar/40Ar"]) - ((1 / i) * self.df["36Ar/39Ar"]))
            tmp["Age"] = ((1 / TOTAL_DECAY_40K) * np.log((ar40star_div_ar39 * self.j_val) + 1)) / 1000000

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

        plateaus = {}  # have three consecutive steps and 50% 39Ar
        for trap in to_check:
            least_fifty = np.array([i for i in to_check.get(trap) if
                                    (self.df["cum %39Ark"].iloc[i[1]] - self.df["cum %39Ark"].iloc[i[0]]) >= 50])
            if least_fifty.size:  # if there are any plateaus
                # convert index units to step units
                plateaus[trap] = index2step(least_fifty, self.df)

        return plateaus
