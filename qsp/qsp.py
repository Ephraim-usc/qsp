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
units.microliter = unum.new_unit('microliter', 1e-6 * units.l)
units.pl = unum.new_unit('pl', 1e-12 * units.l)

units.micromolar = unum.new_unit('micromolar', 1e-6 * units.mol / units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)

units.avagadro = unum.new_unit('avagadro', 6.0221415e23 / units.mol)

np.set_printoptions(suppress=True)

def dict2array(x, names):
  buffer = np.zeros(len(names), dtype = object)
  for key, value in x.items():
    buffer[names.index(key)] += value
  return buffer

def array2dict(x, names, trim = False):
  if trim:
    buffer = {name:x_ for name, x_ in zip(names, x) if x_ != 0}
  else:
    buffer = {name:x_ for name, x_ in zip(names, x)}
  return buffer

def array_number(x, unit):
  buffer = np.array([x_.number(unit) for x_ in x])
  return buffer


class Reaction:
  def __init__(self, compartment, reactants, products, forward, backward = None, side_compartment = None, side_products = None):
    if side_compartment is not None:
      assert side_products is not None, "side compartment is given, but products not provided!"
    if side_products is not None:
      assert side_compartment is not None, "side products are given, but compartment not specified!"

    self.compartment = compartment
    self.reactants = reactants
    self.products = products
    self.forward = forward
    self.side_compartment = side_compartment
    self.side_products = side_products
  
  def apply(self, system, t):
    rate = np.power(system.x[:, compartment], self.reactants).prod() * self.forward - np.power(system.x[:, compartment], self.products).prod() * self.backward
    delta = (self.products - self.reactants) * rate * t
    delta_side = self.side_products * rate * t
    
    system.x[:, compartment] += delta
    system.x[:, side_compartment] += delta_side


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
    self.processes = []
    
    self.t = 0
    self.x = np.zeros([self.n_analytes, self.n_compartments])
    self.z = np.zeros(self.n_variables, dtype = float)
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
  
  def set_z(self, variable, value):
    variable = self.variables.index(variable)
    self.z[variable] = value
  
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
        z = array2dict(self.z, self.variables)
        delta = dict2array(reaction(x, z), self.analytes) * (t_step * units.h)
        self.x[:, compartment] += array_number(delta, units.nM)
      for process in self.processes:
        z = array2dict(self.z, self.variables)
        delta = dict2array(process(z), self.variables) * (t_step * units.h)
        self.z += delta
      
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
        ax.axhline(y = 1, color = 'lightgrey', lw = 1)
        ax.axhline(y = 10, color = 'lightgrey', lw = 1)
        ax.axhline(y = 100, color = 'lightgrey', lw = 1)
        ax.axhline(y = 1000, color = 'lightgrey', lw = 1)
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

  def add_reaction(self, compartment, reactants, products, forward, backward = None, side_compartment = None, side_products = None):
    compartment = self.compartments.index(compartment)
    reactants = dict2array(reactants, self.analytes)
    products = dict2array(products, self.analytes)
    forward = forward.number(units.nM / units.h / units.nM**(reactants.sum().astype(int)))
    if backward is not None:
      backward = backward.number(units.nM / units.h / units.nM**(products.sum().astype(int)))
    if side_compartment is not None:
      assert side_products is not None, "side compartment is given, but products not provided!"
      side_compartment = self.compartments.index(side_compartment)
    if side_products is not None:
      assert side_compartment is not None, "side products are given, but compartment not specified!"
      side_products = dict2array(side_products, self.analytes)
    
    reaction = Reaction(compartment, reactants, products, forward, backward, side_compartment, side_products)
    self.reactions.append(reaction)
  
  def add_reaction(self, compartment, reaction):
    compartment = self.compartments.index(compartment)
    self.reactions.append([compartment, reaction])
  
  def add_process(self, process):
    self.processes.append(process)
  
  def print(self):
    volumes = pd.DataFrame(self.volumes, index = self.analytes, columns = self.compartments)
    print("<volumes>", flush = True)
    print(volumes, flush = True)
    print(" ", flush = True)
    for i, analyte in enumerate(self.analytes):
      flows = self.flows[i,:,:]
      if not flows.any():
        continue
      for j, compartment in enumerate(self.compartments):
        flow = array2dict(np.round(flows[j,:], 6), self.compartments, trim = True)
        print(f"{compartment} {flow}", flush = True)
      print(" ", flush = True)
    print(f"<{len(self.reactions)} reactions>", flush = True)
    print(f"<{len(self.processes)} processes>", flush = True)
