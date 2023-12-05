import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import pandas as pd

def heatmap(df, x, y, z, xlabel, ylabel, zlabel, midpoint = 10, maxpoint = 50, out = "tmp"):
  norm = clr.FuncNorm((lambda x: x/(x+midpoint), lambda y: midpoint*y/(1-y)), vmin=0, vmax=maxpoint)
  cmap = clr.LinearSegmentedColormap.from_list('custom',
                                             [(0,      'black'),
                                              (0.25,   'tab:blue'),
                                              (0.5,    'tab:green'),
                                              (0.75,   'gold'),
                                              (1,      'tab:orange')], N=256)
  
  
  Z = df.pivot(index = y, columns = x, values = z)
  if Z.columns.dtype != "O":
    X = np.round(Z.columns, 2)
  else:
    X = Z.columns
  if Z.index.dtype != "O":
    Y = np.round(Z.index, 2)
  else:
    Y = Z.index
  
  fig, ax = plt.subplots()
  im = ax.imshow(Z, origin = "lower", norm = norm, cmap = cmap)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  cbar = plt.colorbar(im, cax = cax)
  cbar.ax.set_ylabel(zlabel, rotation=-90, va="bottom")
  cbar.set_ticks([0, midpoint, maxpoint])
  cbar.set_ticklabels([0, midpoint, maxpoint])
  
  for i in range(len(X)):
    for j in range(len(Y)):
      text = ax.text(i, j, np.round(Z.iloc[j, i], 2), ha="center", va="center", color="w", fontsize = 6)
  
  ax.set_xticks(np.arange(len(X)), labels = X, fontsize = 6); ax.set_xlabel(xlabel)
  ax.set_yticks(np.arange(len(Y)), labels = Y, fontsize = 6); ax.set_ylabel(ylabel)
  
  fig.tight_layout()
  plt.gcf().set_size_inches(8, 7)
  fig.savefig(f"{out}.png", dpi = 300)
  plt.close(fig)
