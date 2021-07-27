from plotter import *

def cheat_steps(name, df):
    if name == "P K-feldspar":
        df_1 = Isochron.SectionSteps(df, 1, 10, False)
        df_2 = Isochron.SectionSteps(df, 11, 18, True)
        separated = [df_1, df_2]
    elif name == "B K-feldspar":
        df_1 = Isochron.SectionSteps(df, 1, 10, True)  # 2, 10
        df_2 = Isochron.SectionSteps(df, 11, 18, True)
        separated = [df_1, df_2]
    else:
        if name == "Buckskin Peak Biotite":
            x_ticks = np.arange(0.053, .062, step=.001)
            y_ticks = np.arange(0, .00012, step=.00002)
        df_1 = Isochron.SectionSteps(df, has_slope=False, every=True)
        separated = [df_1]
    return separated