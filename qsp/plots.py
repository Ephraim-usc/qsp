import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import pandas as pd

norm = clr.FuncNorm((lambda x: x/(x+10), lambda y: 10*y/(1-y)), vmin=0, vmax=400)
cmap = clr.LinearSegmentedColormap.from_list('custom',
                                             [(0,      'black'),
                                              (0.25,   'tab:blue'),
                                              (0.5,    'tab:green'),
                                              (0.75,   'gold'),
                                              (1,      'tab:orange')], N=256)

def heatmap(df, xcol, ycol, xlabel, ylabel, zcols, zlabels, norm = norm, cmap = cmap, output = None):
  fig, axs = plt.subplots(1, len(zcols))
  for ax, zcol, zlabel in zip(axs, zcols, zlabels):
    data = df.pivot(index = ycol, columns = xcol, values = zcol)
    im = ax.imshow(data, origin = "lower", norm = norm, cmap = cmap)
    for i in range(data.shape[1]):
      for j in range(data.shape[0]):
        text = ax.text(i, j, np.round(data.iloc[j, i], 2), ha="center", va="center", color="w", fontsize = 6)
    
    ax.set_xticks(np.arange(data.shape[1]), labels = data.columns, fontsize = 6); ax.set_xlabel(xlabel)
    ax.set_yticks(np.arange(data.shape[0]), labels = data.index, fontsize = 6); ax.set_ylabel(ylabel)
  
  fig.tight_layout()
  plt.gcf().set_size_inches(7 * len(zcols), 7)
  if output is None:
    fig.savefig(f"tmp.png", dpi = 300)
  else:
    fig.savefig(output, dpi = 300)
