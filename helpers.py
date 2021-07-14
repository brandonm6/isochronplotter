import numpy as np
from statsmodels.api import OLS
from statsmodels.stats.outliers_influence import OLSInfluence
from scipy import stats
# import sklearn.cluster as sc


def check_removal(lst, x, y):
    '''
    Returns list of values to be removed from dataframe [finish writing]
    :param lst: list of values to check for removal
    :param x: independent variable column of dataframe
    :param y: dependent variable column of dataframe
    :return:
    '''
    x_new, y_new = x.copy(), y.copy()
    remove = []
    for val in lst:
        '''
        Check with significance test of correlation coefficients if data 
        becomes linear (as opposed to previously nonlinear) by removing val 
        '''
        x_without, y_without = x_new.drop([val]), y_new.drop([val])
        # Calculate Pearson's correlation coefficient and p-value
        r_with, p_with = stats.pearsonr(x_new, y_new)  # including val
        r_without, p_without = stats.pearsonr(x_without, y_without)  # excluding val

        # if (not linear with val) and (linear without val)
        if (p_with > .05) and (p_without < .05):
            # print("test 1 remove:", val)
            x_new, y_new = x_without, y_without
            remove.append(val)

        '''
        Check for outliers in the residuals
        '''

    return remove


def locate_outliers(x):
    '''
    :param x: dataframe column to check
    :return: list of indexes of outliers
    '''
    q1 = x.quantile(.25)
    q3 = x.quantile(.75)
    iqr = q3 - q1
    outliers = np.where((x < (q1 - 1.5 * iqr)) | (x > (q3 + 1.5 * iqr)))[0]

    return outliers.tolist()


def locate_influential(x, y):
    """
    Finds influential points using Difference in Fits (DFFITS) - diffits quantifies the
    number of standard deviations that the fitted value changes when the ith point is omitted
    :param x: independent variable column of dataframe
    :param y: dependent variable column of dataframe
    :return: list of indexes of influential points
    """
    y = y.to_frame()
    x = x.to_frame()
    x.insert(loc=0, column="const", value=1.0)

    model = OLS(y, x)
    dffits = OLSInfluence(model.fit()).summary_frame().filter(["dffits"])

    k, n = 1, len(dffits)
    threshold_1 = 2 * np.sqrt((k + 1) / n)
    threshold_2 = 2 * np.sqrt((k + 2) / (n - k - 2))
    threshold_3 = 2 * np.sqrt(k / n)
    avg_threshold = np.mean([threshold_1, threshold_2, threshold_3])

    infl = dffits[abs(dffits["dffits"]) > avg_threshold].index.values.tolist()

    return infl


def identify_nonhomogeneity():
    """
    nonhomogeneous trapped argon can be disregarded
    identify through scatter in initial increments
    maybe use cluster identification algorithm
    min requirements for plateau are:
     - at least 3 steps and 50% cumulative %39Ar
    """
    # x = self.df[["39Ar/40Ar", "36Ar/40Ar"]].to_numpy()
    # db = sc.DBSCAN(eps=.025).fit(x)  # (eps=3, min_samples=2
