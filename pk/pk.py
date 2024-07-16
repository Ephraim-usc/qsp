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


