import math
import numpy as np
import pandas as pd
import functools

from scipy.linalg import expm
from tqdm import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt

import unum
import unum.units as units

units.microgram = unum.new_unit('microgram', 1e-6 * units.g)

units.l = unum.new_unit('l', 1e-3 * units.m ** 3)
units.ml = unum.new_unit('ml', 1e-3 * units.l)
units.microliter = unum.new_unit('microliter', 1e-6 * units.l)
units.pl = unum.new_unit('pl', 1e-12 * units.l)

units.micromolar = unum.new_unit('micromolar', 1e-6 * units.mol / units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)

units.avagadro = unum.new_unit('avagadro', 6.0221415e23 / units.mol)

np.set_printoptions(suppress=True)

def dict2array(x, names, dtype = None):
  if dtype is None:
    buffer = np.zeros(len(names), dtype = object)
  else:
    buffer = np.zeros(len(names), dtype = dtype)
  for key, value in x.items():
    buffer[names.index(key)] += value
  return buffer

def array2dict(x, names, trim = False):
  if trim:
    buffer = {name:x_ for name, x_ in zip(names, x) if x_ != 0}
  else:
    buffer = {name:x_ for name, x_ in zip(names, x)}
  return buffer


def reaction_general(system, compartment, reactants, products, forward, backward, side_compartment, side_products, t):
  rate = np.power(system.x[:, compartment], reactants).prod() * forward
  if backward is not None:
    rate -= np.power(system.x[:, compartment], products).prod() * backward
  
  delta = (products - reactants) * rate * t
  system.x[:, compartment] += delta
  
  if side_compartment is not None:
    delta_side = side_products * rate * t
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
    
    self.V = np.zeros([self.n_analytes, self.n_compartments], dtype = float) # in units.ml
    self.Q = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments], dtype = float) # in 1/units.h
    self.reactions = []
    self.processes = []
    
    self.t = 0
    self.x = np.zeros([self.n_analytes, self.n_compartments], dtype = float) # in units.nM
    self.z = np.zeros(self.n_variables, dtype = float) # any object
    self.history = []
  
  def get_volume(self, analyte, compartment):
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    return self.V[analyte, compartment] * units.ml
  
  def set_volume(self, analyte, compartment, value):
    value = value.number(units.ml)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.V[analyte, compartment] = value
  
  # set compartment_dest as None if it is a clearance
  def add_flow(self, analyte, compartment_source, compartment_dest, rate):
    rate = rate.number(units.ml/units.h)
    analyte = self.analytes.index(analyte)
    compartment_source = self.compartments.index(compartment_source)
    self.Q[analyte, compartment_source, compartment_source] -= rate / self.V[analyte, compartment_source]
    if compartment_dest is not None:
      compartment_dest = self.compartments.index(compartment_dest)
      self.Q[analyte, compartment_source, compartment_dest] += rate / self.V[analyte, compartment_dest]
  
  def add_reaction(self, compartment, reactants, products, forward, backward = None, side_compartment = None, side_products = None):
    compartment = self.compartments.index(compartment)
    reactants = dict2array(reactants, self.analytes, dtype = int)
    products = dict2array(products, self.analytes, dtype = int)
    forward = forward.number(units.nM / units.h / units.nM**(reactants.sum()))
    if backward is not None:
      backward = backward.number(units.nM / units.h / units.nM**(products.sum()))
    if side_compartment is not None:
      assert side_products is not None, "side compartment is given, but products not provided!"
      side_compartment = self.compartments.index(side_compartment)
    if side_products is not None:
      assert side_compartment is not None, "side products are given, but compartment not specified!"
      side_products = dict2array(side_products, self.analytes, dtype = int)
    
    reaction = functools.partial(reaction_general, self, compartment, reactants, products, forward, backward, side_compartment, side_products)
    self.reactions.append(reaction)
  
  def add_process(self, process_func):
    process = functools.partial(process_func, self)
    self.processes.append(process)
  
  
  def get_x(self, analyte, compartment):
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    return self.x[analyte, compartment] * units.nM
  
  def set_x(self, analyte, compartment, value):
    value = value.number(units.nM)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.x[analyte, compartment] = value
  
  def add_x(self, analyte, compartment, value):
    value = value.number(units.nM)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.x[analyte, compartment] += value
  
  def clear_x(self):
    self.x = np.zeros([self.n_analytes, self.n_compartments], dtype = float)
  
  def get_z(self, variable):
    variable = self.variables.index(variable)
    return self.z[variable]
  
  def set_z(self, variable, value):
    variable = self.variables.index(variable)
    self.z[variable] = value
  
  def print(self):
    V = pd.DataFrame(self.V, index = self.analytes, columns = self.compartments)
    print("<volumes>", flush = True)
    print(V, flush = True)
    print(" ", flush = True)
    
    for i, analyte in enumerate(self.analytes):
      Q = self.Q[i,:,:]
      if not Q.any():
        continue
      Q = pd.DataFrame(Q, index = self.compartments, columns = self.compartments)
      print(f"<Q matrix for {analyte}>", flush = True)
      print(Q, flush = True)
      print(" ", flush = True)
      #for j, compartment in enumerate(self.compartments):
      #  Q = array2dict(np.round(flows[j,:], 6), self.compartments, trim = True)
      #  print(f"{compartment} {flow}", flush = True)
      #print(" ", flush = True)
    
    print(f"<{len(self.reactions)} reactions>", flush = True)
    print(f"<{len(self.processes)} processes>", flush = True)
    print(" ", flush = True)
    
    x = pd.DataFrame(self.x, index = self.analytes, columns = self.compartments)
    print("<x>", flush = True)
    print(x, flush = True)
    print(" ", flush = True)
  
  
  def clear_t(self):
    self.t = 0
    self.history = []
  
  def run(self, t_end, t_step = 1/60 * units.h, t_record = 1 * units.h):
    t_end = t_end.number(units.h)
    t_step = t_step.number(units.h)
    t_record = t_record.number(units.h)
    flowing_analytes = [analyte for analyte in range(self.n_analytes) if self.Q[analyte].any()]
    
    pbar = tqdm(total = t_end, unit = "h", bar_format = "{desc}: {percentage:3.0f}%|{bar}| {n:.2f}/{total_fmt} [{elapsed}<{remaining},  {rate_fmt}{postfix}]")
    pbar.update(self.t)
    while True:
      t_ = self.t
      self.t = min(self.t + t_step, t_end)
      for analyte in flowing_analytes:
        self.x[analyte] = np.dot(self.x[analyte], expm((self.t - t_) * self.Q[analyte]))
      for reaction in self.reactions:
        reaction(self.t - t_)
      for process in self.processes:
        process((self.t - t_) * units.h)
      
      if math.floor(self.t / t_record) > math.floor(t_ / t_record):
        self.history.append((self.t, self.x.copy()))
      pbar.update(self.t - t_)
      if math.isclose(self.t, t_end, rel_tol = 0, abs_tol = 1e-9):
        break
    pbar.close()
  
  def plot(self, compartments, output = None):
    compartments = [self.compartments.index(compartment) for compartment in compartments]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    Xmax = max([t for t, x in self.history])
    Ymax = max([x.max() for t, x in self.history]); Ymax = 10**np.ceil(np.log10(Ymax))
    
    fig, axs = plt.subplots(nrows = 1, ncols = len(compartments), figsize = (4*len(compartments), 3), squeeze = False)
    axs = axs.ravel().tolist()
    for ax, compartment in zip(axs, compartments):
      for analyte in range(len(self.analytes)):
        X = [t for t, x in self.history]
        Y = [x[analyte, compartment] for t, x in self.history]
        AVG = np.trapz(Y, X) / (X[-1] - X[0])
        if AVG > 0:
          ax.plot(X, Y, label = f"{self.analytes[analyte]}, avg={AVG:.3}nM", color = colors[analyte])
      ax.set_xticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
      ax.set_xlim(0, Xmax)
      ax.set_yscale('symlog', linthresh = 1e-2)
      ax.set_yticks([1e-2, 1e-1, 0, 1, 10, 100, 1000, 10000, 1e5, 1e6])
      ax.set_ylim(0, Ymax)
      ax.grid(axis = "y", color = "grey", linewidth = 1)
      ax.set_title(self.compartments[compartment])
      ax.legend(loc = "upper right", prop={'size': 6})
    
    if output is None:
      fig.show()
    else:
      fig.savefig(output, dpi = 300)
