import numpy as np
from statsmodels.api import OLS
from statsmodels.stats.outliers_influence import OLSInfluence
from scipy import stats
# import sklearn.cluster as sc


def locate_outliers_resid(x, y):
    """
    Finds outliers based on external studentized residuals > 3
    :param x: independent variable column of dataframe
    :param y: dependent variable column of dataframe
    :return: list of indexes of steps to be removed
    """
    y = y.to_frame()
    x = x.to_frame()
    x.insert(loc=0, column="const", value=1.0)

    model = OLS(y, x)
    stdnt_resid = OLSInfluence(model.fit()).resid_studentized_external
    outliers = stdnt_resid[abs(stdnt_resid) > 3].index.values.tolist()

    return outliers


# def locate_outliers(x):
#     '''
#     :param x: dataframe column to check
#     :return: list of indexes of outliers
#     '''
#     q1 = x.quantile(.25)
#     q3 = x.quantile(.75)
#     iqr = q3 - q1
#     outliers = np.where((x < (q1 - 1.5 * iqr)) | (x > (q3 + 1.5 * iqr)))[0]
#
#     return outliers.tolist()
#
# def locate_influential(x, y):
#     """
#     Finds influential points using Difference in Fits (DFFITS) then checks to remove
#     them by seeing if data becomes linear (where previously non-linear) without them
#     :param x: independent variable column of dataframe
#     :param y: dependent variable column of dataframe
#     :return: list of indexes of influential points to be removed
#     """
#     y = y.to_frame()
#     x = x.to_frame()
#     x.insert(loc=0, column="const", value=1.0)
#
#     model = OLS(y, x)
#     dffits, _ = OLSInfluence(model.fit()).dffits
#     k, n = 1, len(dffits)
#     threshold_1 = 2 * np.sqrt((k + 1) / n)
#     threshold_2 = 2 * np.sqrt((k + 2) / (n - k - 2))
#     threshold_3 = 2 * np.sqrt(k / n)
#     avg_threshold = np.mean([threshold_1, threshold_2, threshold_3])
#     infl = dffits[abs(dffits) > avg_threshold].index.values.tolist()
#
#     x_new, y_new = x.copy(), y.copy()
#     remove = []
#     for val in infl:
#         x_without, y_without = x_new.drop([val]), y_new.drop([val])
#         # Calculate Pearson's correlation coefficient and p-value
#         r_with, p_with = stats.pearsonr(x_new, y_new)  # including val
#         r_without, p_without = stats.pearsonr(x_without, y_without)  # excluding val
#
#         # if (not linear with val) and (linear without val)
#         if (p_with > .05) and (p_without < .05):
#             # print("test 1 remove:", val)
#             x_new, y_new = x_without, y_without
#             remove.append(val)
#
#     return remove


def identify_nonhomogeneity():
    """
    nonhomogeneous trapped argon can be disregarded
    identify through scatter in initial increments
    maybe use cluster identification algorithm
    min requirements for plateau are:
     - at least 3 steps and 50% cumulative %39Ar
    """
