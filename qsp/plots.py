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
    ax.set_title(zlabel)
  
  #fig.tight_layout()
  plt.gcf().set_size_inches(7 * len(zcols), 7)
  if output is None:
    fig.savefig(f"tmp.png", dpi = 300)
  else:
    fig.savefig(output, dpi = 300)

def compare_systems(systems, labels, analytes, compartments = None, colors = None, linthresh = 1e-3, output = None):
  indices = [systems[0].analytes.index(analyte) for analyte in analytes]
  
  if compartments is None:
    compartments = systems[0].compartments
  compartments = [systems[0].compartments.index(compartment) for compartment in compartments]
  
  if labels is None:
    labels = [f"{i}" for i, system in enumerate(systems)]
  
  if colors is None:
    colors = list(mcolors.TABLEAU_COLORS.values())
  
  Xmax = max([t for t, x in systems[0].history])
  Ymax = max([x[indices, compartment].sum() for system in systems for t, x in system.history for compartment in compartments])
  Ymax = 10**np.ceil(np.log10(Ymax))
  
  fig, axs = plt.subplots(nrows = 1, ncols = len(compartments), figsize = (4*len(compartments), 3), squeeze = False)
  axs = axs.ravel().tolist()
  for ax, compartment in zip(axs, compartments):
    for system, label, color in zip(systems, labels, colors):
      X = [t for t, x in system.history]
      Y = [x[indices, compartment].sum() for t, x in system.history]
      AVG = np.trapz(Y, X) / (X[-1] - X[0])
      ax.plot(X, Y, label = f"{label}, avg={AVG:.3}nM", color = color)
      if Xmax > 100:
        ax.set_xticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
      else:
        ax.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    
    ax.set_xlim(0, Xmax)
    ax.set_yscale('symlog', linthresh = linthresh)
    ax.set_yticks([y for y in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0, 1, 10, 100, 1000, 10000, 1e5, 1e6] if y >= linthresh])
    ax.set_ylim(0, Ymax)
    ax.grid(axis = "y", color = "grey", linewidth = 1)
    ax.set_title(systems[0].compartments[compartment])
    ax.legend(loc = "upper right", prop={'size': 6})
    
    if output is None:
      fig.show()
    else:
      fig.savefig(output, dpi = 300)
      plt.close(fig)
