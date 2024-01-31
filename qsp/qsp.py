import math
import numpy as np
import pandas as pd
import functools

from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

from tqdm import tqdm
from time import time as tt

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import unum
import unum.units as units

units.l = unum.new_unit('l', 1e-3 * units.m ** 3)
units.ml = unum.new_unit('ml', 1e-3 * units.l)
units.ul = unum.new_unit('ul', 1e-6 * units.l)
units.pl = unum.new_unit('pl', 1e-12 * units.l)

units.M = unum.new_unit('M', 1 * units.mol / units.l)
units.uM = unum.new_unit('uM', 1e-6 * units.mol / units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)

units.kDa = unum.new_unit('kDa', units.kg / units.mol)
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


class RS: # reaction system
  def __init__(self, n_analytes):
    self.active = False
    self.n = n_analytes
    self.Q = np.zeros([n_analytes, n_analytes]) # linear term coefficients
    self.QQ = np.zeros([n_analytes, n_analytes, n_analytes]) # quadratic term coefficients
  
  def refresh(self):
    self.linear_i, self.linear_o = np.where(self.Q != 0)
    self.linear_k = self.Q[self.Q != 0]
    self.quadratic_i, self.quadratic_j, self.quadratic_o = np.where(self.QQ != 0)
    self.quadratic_k = self.QQ[self.QQ != 0]
  
  def add_simple_forward(self, reactants, products, forward): # reactants can be of length 1 or 2, products can be of any length
    self.active = True
    if len(reactants) == 1:
      self.Q[reactants, reactants] -= forward
      self.Q[reactants, products] += forward
    else:
      self.QQ[reactants[0], reactants[1], reactants] -= forward
      self.QQ[reactants[0], reactants[1], products] += forward
  
  def add_simple(self, reactants, products, forward, backward):
    self.active = True
    self.add_simple_forward(reactants, products, forward)
    if backward > 0:
      self.add_simple_forward(products, reactants, backward)
  
  def rate(self, _, x):
    buffer = np.zeros(self.n)
    np.add.at(buffer, self.linear_o, self.linear_k * x[self.linear_i])
    np.add.at(buffer, self.quadratic_o, self.quadratic_k * x[self.quadratic_i] * x[self.quadratic_j])
    return buffer
  
  def jac(self, _, x):
    buffer = np.zeros((self.n, self.n))
    np.add.at(buffer, (self.linear_o, self.linear_i), self.linear_k)
    np.add.at(buffer, (self.quadratic_o, self.quadratic_i), self.quadratic_k * x[self.quadratic_j])
    np.add.at(buffer, (self.quadratic_o, self.quadratic_j), self.quadratic_k * x[self.quadratic_i])
    return buffer
  
  def __call__(self, x, t):
    return solve_ivp(self.rate, jac = self.jac, t_span = (0, t), y0 = x, t_eval=[t], method = "BDF", rtol = 1e-10, atol = 1e-8).y[:,0]


class System:
  def __init__(self, analytes, compartments, variables = None, cells = None):
    variables = [] if variables is None else variables
    cells = [] if cells is None else cells
    
    self.analytes = analytes
    self.n_analytes = len(analytes)
    
    self.compartments = compartments
    self.n_compartments = len(compartments)
    
    self.variables = variables
    self.n_variables = len(variables)
    
    self.cells = cells
    self.n_cells = len(cells)
    
    self.V = np.zeros([self.n_analytes, self.n_compartments], dtype = float) # in units.ml
    self.Q = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments], dtype = float) # in 1/units.h
    self.reactions = []
    self.RS = [RS(self.n_analytes) for compartment in self.compartments]
    self.processes = []
    
    self.t = 0
    self.x = np.zeros([self.n_analytes, self.n_compartments], dtype = float) # in units.nM
    self.z = np.zeros(self.n_variables, dtype = float) # any object
    self.c = np.zeros([self.n_cells, self.n_compartments], dtype = float) # in units.nM
    
    self.history = []
    self.history_cells = []
  
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

  # this function is for backward compatibility only
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
  
  def add_simple(self, compartment, reactants, products, forward, backward = None):
    compartment = self.compartments.index(compartment)
    reactants = [self.analytes.index(reactant) for reactant in reactants]
    products = [self.analytes.index(product) for product in products]
    forward = forward.number(units.nM / units.h / units.nM**(len(reactants)))
    if backward is None:
      backward = 0.0
    else:
      backward = backward.number(units.nM / units.h / units.nM**(len(products)))
    self.RS[compartment].add_simple(reactants, products, forward, backward)
  
  def add_process(self, process):
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
  
  def add_z(self, variable, value):
    variable = self.variables.index(variable)
    self.z[variable] += value
  
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
    print(f"<{len([compartment for compartment in range(self.n_compartments) if self.RS[compartment].active])} reactive compartments>", flush = True)
    print(f"<{len(self.processes)} processes>", flush = True)
    print(" ", flush = True)
    
    x = pd.DataFrame(self.x, index = self.analytes, columns = self.compartments)
    print("<x>", flush = True)
    print(x, flush = True)
    print(" ", flush = True)
  
  
  def clear_t(self):
    self.t = 0
    self.history = []
  
  def run_reactions(self, t):
    t = t.number(units.h)
    reacting_compartments = [compartment for compartment in range(self.n_compartments) if self.RS[compartment].active]
    for compartment in reacting_compartments:
      self.RS[compartment].refresh()
    for compartment in reacting_compartments:
      self.x[:, compartment] = self.RS[compartment](self.x[:, compartment], t)
    
    self.t = self.t + t
    self.history_cells.append((self.t, self.c.copy()))
    self.history.append((self.t, self.x.copy()))
  
  def run(self, t, t_step = 1/60 * units.h, t_record = 1 * units.h):
    t_begin = self.t
    t = t.number(units.h)
    t_step = t_step.number(units.h)
    t_record = t_record.number(units.h)
    flowing_analytes = [analyte for analyte in range(self.n_analytes) if self.Q[analyte].any()]
    reacting_compartments = [compartment for compartment in range(self.n_compartments) if self.RS[compartment].active]
    for compartment in reacting_compartments:
      self.RS[compartment].refresh()
    
    pbar = tqdm(total = t, unit = "h", bar_format = "{desc}: {percentage:3.0f}%|{bar}| {n:.2f}/{total_fmt} [{elapsed}<{remaining},  {rate_fmt}{postfix}]")
    pbar.update(0.0); A, B, C, D = 0.0, 0.0, 0.0, 0.0
    while True:
      t_prev = self.t
      self.t = min(self.t + t_step, t_begin + t)
      t_delta = self.t - t_prev
      for analyte in flowing_analytes:
        A -= tt()
        self.x[analyte] = np.dot(self.x[analyte], expm(t_delta * self.Q[analyte]))
        A += tt()
      for reaction in self.reactions:
        B -= tt()
        reaction(t_delta)
        B += tt()
      for compartment in reacting_compartments:
        C -= tt()
        self.x[:, compartment] = self.RS[compartment](self.x[:, compartment], t_delta)
        C += tt()
      for process in self.processes:
        D -= tt()
        process(self, t_delta * units.h)
        D += tt()
      
      if math.floor(self.t / t_record) > math.floor(t_prev / t_record):
        self.history_cells.append((self.t, self.c.copy()))
        self.history.append((self.t, self.x.copy()))
      pbar.update(t_delta)
      if math.isclose(self.t, t_begin + t, rel_tol = 0, abs_tol = 1e-9):
        break
    pbar.close()
    print(f"time in computing flows: {A:.8f}s\ntime in computing reactions: {B:.8f}s\ntime in computing reactions: {C:.8f}s\ntime in computing processes: {D:.8f}s\n", flush = True)
  
  def plot(self, compartments = None, groups = None, labels = None, colors = None, linestyles = None, linthresh = 1e-3, output = None):
    if compartments is None:
      compartments = self.compartments
    compartments = [self.compartments.index(compartment) for compartment in compartments]
    
    if groups is None:
      groups = [[i] for i in range(self.n_analytes)]
    else:
      groups = [[self.analytes.index(analyte) for analyte in group] for group in groups]
    
    if labels is None:
      labels = [" + ".join([self.analytes[analyte] for analyte in group]) for group in groups]
    
    if colors is None:
      colors = list(mcolors.TABLEAU_COLORS.values())
    if linestyles is None:
      linestyles = ["solid"] * 10
    
    Xmax = max([t for t, x in self.history])
    Ymax = max([x[group, compartment].sum() for t, x in self.history for group in groups for compartment in compartments])
    Ymax = 10**np.ceil(np.log10(Ymax))
    
    fig, axs = plt.subplots(nrows = 1, ncols = len(compartments), figsize = (4*len(compartments), 3), squeeze = False)
    axs = axs.ravel().tolist()
    for ax, compartment in zip(axs, compartments):
      for group, label, color, linestyle in zip(groups, labels, colors, linestyles):
        X = [t for t, x in self.history]
        Y = [x[group, compartment].sum() for t, x in self.history]
        AVG = np.trapz(Y, X) / (X[-1] - X[0])
        if AVG > 0:
          ax.plot(X, Y, label = f"{label}, avg={AVG:.3}nM", color = color, linestyle = linestyle)
      if Xmax > 100:
        ax.set_xticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
      else:
        ax.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
      ax.set_xlim(0, Xmax)
      ax.set_yscale('symlog', linthresh = linthresh)
      ax.set_yticks([y for y in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0, 1, 10, 100, 1000, 10000, 1e5, 1e6] if y >= linthresh])
      ax.set_ylim(0, Ymax)
      ax.grid(axis = "y", color = "grey", linewidth = 1)
      ax.set_title(self.compartments[compartment])
      ax.legend(loc = "upper right", prop={'size': 6})
    
    if output is None:
      fig.show()
    else:
      fig.savefig(output, dpi = 300)
      plt.close(fig)
  
  def summary(self, analytes):
    analytes = [self.analytes.index(analyte) for analyte in analytes]
    avgs = []; maxs = []; hfws = []
    for compartment in range(self.n_compartments):
      X = np.array([t for t, x in self.history])
      Y = np.array([x[analytes, compartment].sum(axis = 0) for t, x in self.history])
      steps = X[1:] - X[:-1]
      widths = (np.append(0, steps) + np.append(steps, 0))/2
      idx = np.argsort(Y)[::-1]
      cumsums = np.cumsum((Y * widths)[idx])
      tmp = np.where(cumsums >= cumsums[-1]/2)[0].min() # minimum number of intervals needed to have 50% of the AUC
      halfwidth = widths[idx[:tmp]].sum()
      
      avgs.append((widths*Y).sum() / (X[-1] - X[0]))
      maxs.append(Y.max())
      hfws.append(halfwidth)
    
    buffer = pd.DataFrame({"average":avgs, "maximum":maxs, "halfwidth":hfws}, index = self.compartments)
    return buffer

