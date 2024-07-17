import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from scipy.optimize import curve_fit
from scipy.optimize import fsolve



# T211-Dose2-168h(20X)
def parse_sample(sample):
  name, dilution_str = sample.split("(")
  dilution = float(dilution_str[:-2])
  
  id, dose_str, time_str = name.split("-")
  dose = (int(dose_str[4:]) - 1) * 168.0
  
  if time_str[-1] == "h":
    time = int(time_str[:-1])
  elif time_str[-3:] == "min":
    time = int(time_str[:-3]) / 60
  
  return id, dose + time, dilution

def parse_samples(samples):
  ids = []
  times = []
  dilutions = []
  for sample in samples:
    if sample is np.nan:
      ids.append(np.nan)
      times.append(np.nan)
      dilutions.append(np.nan)
    else:
      id, time, dilution = parse_sample(sample)
      ids.append(id)
      times.append(time)
      dilutions.append(dilution)
  return ids, times, dilutions

def hyperbolic(x, baseline, emax, aff):
  return baseline + emax * x / (x + aff)

def hyperbolic_inverse(x, baseline, emax, aff):
  return aff / (emax/(x - baseline) - 1)

def bivariate_hyperbolic(x, baseline, emax, aff1, aff2):
  x1, x2 = x
  return baseline + emax * (x1/aff1 + x2/aff2) / (1 + x1/aff1 + x2/aff2)

def fit(concs, values, cutoff = 3.5):
  idx = np.logical_and(values < 3.5, ~np.isnan(values))
  concs, values = concs[idx], values[idx]
  popt, pcov = curve_fit(hyperbolic, concs, values, p0 = [0.01, 1.0, 10.0], bounds = ([0.0, 0.0, 0.0], [10.0, 10.0, np.inf]))
  baseline, emax, aff = popt
  return baseline, emax, aff

# concs should be of shape (2, M)
def bivariate_fit(concs, values, cutoff = 3.5):
  idx = np.logical_and(values < 3.5, ~np.isnan(values))
  concs, values = concs[:, idx], values[idx]
  popt, pcov = curve_fit(bivariate_hyperbolic, concs, values, p0 = [0.01, 1.0, 10.0, 10.0], bounds = ([0.0, 0.0, 0.0, 0.0], [10.0, 10.0, np.inf, np.inf]))
  baseline, emax, aff1, aff2 = popt
  return baseline, emax, aff1, aff2


def solve(baseline, emax, aff1, aff2, total):
  func = lambda x: bivariate_hyperbolic(x, total - x, baseline, emax, aff1, aff2)
  x1 = fsolve(func, x0 = 100)
  return x1


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





