import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import fsolve


def parse_time(time_str):
  if "-" in time_str:
    dose_str, time_str = time_str.split("-")
    dose = (int(dose_str[4:]) - 1) * 168.0
  else:
    dose = 0.0
  
  if time_str[-1] == "h":
    time = int(time_str[:-1])
  elif time_str[-3:] == "min":
    time = int(time_str[:-3]) / 60
  
  return dose + time

def hyperbolic(x, baseline, emax, aff):
  return baseline + emax * x / (x + aff)

def hyperbolic_inverse(x, baseline, emax, aff):
  return aff / (emax/(x - baseline) - 1)

# ax: the ax object to plot on, if None then do not plot but just return fitted parameters
# drugs: the vector of drug forms
# concs: the vector of concentrations if known, use NaN for observed samples
def fit_standard_curves(ax, drugs, concs, values, cutoff):
  idx = np.logical_and(values < 3.5, values.notnull())
  concs_train, values_train = concs[idx], values[idx]
  popt, pcov = curve_fit(hyperbolic, concs_train, values_train)
  
  if ax:
    x_fit = np.power(10.0, np.arange(-2, 5, 0.1))
    y_fit = [hyperbolic(_, *popt) for _ in x_fit]
    ax.plot(x_fit, y_fit)





