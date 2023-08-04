import math
import numpy as np
import pandas as pd

from scipy.linalg import expm
from tqdm import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt

import unum
import unum.units as units

units.l = unum.new_unit('l', 1e-3 * units.m ** 3)
units.ml = unum.new_unit('ml', 1e-3 * units.l)
units.pl = unum.new_unit('pl', 1e-12 * units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)
units.avagadro = unum.new_unit('avagadro', 6.0221415e23 / units.mol)

np.set_printoptions(suppress=True)

def dict2array(x, names):
  buffer = np.zeros(len(names), dtype = object)
  for key, value in x.items():
    buffer[names.index(key)] += value
  return buffer

def array2dict(x, names):
  buffer = {name:x_ for name, x_ in zip(names, x)}
  return buffer

def array_number(x, unit):
  buffer = np.array([x_.number(unit) for x_ in x])
  return buffer


class System:
  def __init__(self, analytes, compartments, variables = None):
    self.analytes = analytes
    self.n_analytes = len(analytes)
    self.compartments = compartments
    self.n_compartments = len(compartments)
    
    variables = [] if variables is None else variables
    self.variables = variables
    self.n_variables = len(variables)
    
    self.volumes = np.zeros([self.n_analytes, self.n_compartments])
    self.flows = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments])
    self.Qs = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments])
    self.reactions = []
    
    self.t = 0
    self.x = np.zeros([self.n_analytes, self.n_compartments])
    self.z = np.zeros(self.n_variables)
    self.history = None
  
  def clear_x(self):
    self.x = np.zeros([self.n_analytes, self.n_compartments])
  
  def clear_t(self):
    self.t = 0
    self.history = []
  
  def set_x(self, analyte, compartment, concentration):
    concentration = concentration.number(units.nM)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.x[analyte, compartment] = concentration
  
  def set_z(self, variable):
    variable = self.variables.index(variable)
    self.z[variable] = variable
  
  def run(self, t_end, t_step = 1/60 * units.h, t_record = 1 * units.h):
    t_end = t_end.number(units.h)
    t_step = t_step.number(units.h)
    t_record = t_record.number(units.h)
    flowing_analytes = [analyte for analyte in range(self.n_analytes) if self.Qs[analyte].any()]
    
    pbar = tqdm(total = t_end, unit = "h", bar_format = "{desc}: {percentage:3.0f}%|{bar}| {n:.2f}/{total_fmt} [{elapsed}<{remaining},  {rate_fmt}{postfix}]")
    pbar.update(self.t)
    while True:
      t_ = self.t
      self.t = min(self.t + t_step, t_end)
      for analyte in flowing_analytes:
        self.x[analyte] = np.dot(self.x[analyte], expm((self.t - t_) * self.Qs[analyte]))
      for compartment, reaction in self.reactions:
        x = array2dict(self.x[:, compartment] * units.nM, self.analytes)
        z = array2dict(self.z[variable], self.variables)
        delta = dict2array(reaction(x, self.t * units.h), self.analytes) * (t_step * units.h)
        self.x[:, compartment] += array_number(delta, units.nM)
      if math.floor(self.t / t_record) > math.floor(t_ / t_record):
        self.history.append((self.t, self.x.copy()))
      pbar.update(self.t - t_)
      if math.isclose(self.t, t_end, rel_tol = 0, abs_tol = 1e-9):
        break
    pbar.close()
  
  def plot(self, compartments):
    compartments = [self.compartments.index(compartment) for compartment in compartments]
    fig, axs = plt.subplots(nrows = 1, ncols = len(compartments), squeeze = False)
    axs = axs.ravel().tolist()
    for ax, compartment in zip(axs, compartments):
      for analyte in range(len(self.analytes)):
        X = [t for t, x in self.history]
        Y = [x[analyte, compartment] for t, x in self.history]
        AVG = np.trapz(Y, X) / (X[-1] - X[0])
        ax.plot(X, Y, label = f"{self.analytes[analyte]}, avg={AVG:.2f}nM")
      ax.set_yscale('symlog', linthresh = 1)
      ax.set_yticks([0, 1, 10, 100, 1000, 10000])
      ax.set_title(self.compartments[compartment])
      ax.legend(prop={'size': 6})
    fig.show()
  
  # if volumes is a vector, we assume all analytes share the same volume in each department
  def set_volumes(self, volumes):
    if volumes.ndim == 1:
      volumes = np.tile(volumes, reps = [self.n_analytes, 1])
    volumes = np.array([[(i.number(units.ml) if type(i) is unum.Unum else i) for i in row] for row in volumes])
    self.volumes = volumes
  
  def set_volume(self, analyte, compartment, volume):
    volume = volume.number(units.ml)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.volumes[analyte, compartment] = volume
  
  def get_volume(self, analyte, compartment):
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    return self.volumes[analyte, compartment] * units.ml
  
  # set compartment_dest as None if it is a clearance
  def add_flow(self, analyte, compartment_source, compartment_dest, rate):
    rate = rate.number(units.ml/units.h)
    analyte = self.analytes.index(analyte)
    compartment_source = self.compartments.index(compartment_source)
    self.flows[analyte, compartment_source, compartment_source] -= rate
    self.Qs[analyte, compartment_source, compartment_source] -= rate / self.volumes[analyte, compartment_source]
    if compartment_dest is not None:
      compartment_dest = self.compartments.index(compartment_dest)
      self.flows[analyte, compartment_source, compartment_dest] += rate
      self.Qs[analyte, compartment_source, compartment_dest] += rate / self.volumes[analyte, compartment_dest]
  
  # set compartment as None if it happens everywhere
  # powers: dict of powers for each analyte
  # deltas: dict of concentration changes of the analytes per reaction
  def add_reaction(self, compartment, reaction):
    compartment = self.compartments.index(compartment)
    self.reactions.append([compartment, reaction])
  
  def str_reaction(self, reaction):
    compartment, coefficient, powers, deltas = reaction
    str = f"[{compartment}] "
    str += "[rate = " + " * ".join([f"{coefficient:.6f}"] + [f"{self.analytes[analyte]}^{power}" for analyte,power in enumerate(powers) if power > 0]) + "] "
    
    inputs = [f"{-delta}*{self.analytes[analyte]}" for analyte,delta in enumerate(deltas) if delta < 0]
    outputs = [f"{delta}*{self.analytes[analyte]}" for analyte,delta in enumerate(deltas) if delta > 0]
    str += "[" + " + ".join(inputs) + " → " + " + ".join(outputs) + "]"
    return str
  
  def print(self):
    volumes = pd.DataFrame(self.volumes, index = self.analytes, columns = self.compartments)
    print("<volumes>", flush = True)
    print(volumes, flush = True)
    print(" ", flush = True)
    for i, analyte in enumerate(self.analytes):
      flows = pd.DataFrame(self.flows[i,:,:], index = self.compartments, columns = self.compartments)
      if not flows.values.any():
        continue
      print(f"<flow of {analyte}>", flush = True)
      print(flows, flush = True)
      print(" ", flush = True)
    print("<reactions>", flush = True)
    for reaction in self.reactions:
      print(self.str_reaction(reaction), flush = True)
