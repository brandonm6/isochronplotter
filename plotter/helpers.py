
import matplotlib.pyplot as plt


def set_size(w, h, ax):
    # w, h: width, height in inches
    left = ax.figure.subplotpars.left
    right = ax.figure.subplotpars.right
    top = ax.figure.subplotpars.top
    btm = ax.figure.subplotpars.bottom
    figw = float(w) / (right - left)
    figh = float(h) / (top - btm)
    ax.figure.set_size_inches(figw, figh)


def get_aspect(ax):
    ylower, yupper = ax.get_ylim()
    xlower, xupper = ax.get_xlim()
    data_ratio = (yupper - ylower) / (xupper - xlower)

    # total figure size
    fig_w, fig_h = ax.get_figure().get_size_inches()
    # axis size on figure
    _, _, w, h = ax.get_position().bounds
    disp_ratio = (fig_h * h) / (fig_w * w)

    return data_ratio / disp_ratio


def step2index(df, step):
    return df.index[df['Step'] == step].tolist()[0]


def remove_complete_overlaps(lst):
    """Removes complete overlaps in a list of intervals

    - ex. [[3, 7], [3, 6], [5, 9]] --> [[3, 7], [5, 9]]
    :param lst: sorted(list) where list holds intervals
    :return: list with complete overlaps removed
    """
    if len(lst) <= 1:
        return lst
    # first interval falls inside second
    elif (lst[0][0] >= lst[1][0]) and (lst[0][1] <= lst[1][1]):
        return remove_complete_overlaps(lst[1:])
    # second interval falls inside first
    elif (lst[0][0] <= lst[1][0]) and (lst[0][1] >= lst[1][1]):
        return remove_complete_overlaps([lst[0]] + lst[2:])
    else:
        return [lst[0]] + remove_complete_overlaps(lst[1:])
