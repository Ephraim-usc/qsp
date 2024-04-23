import math
import numpy as np
import pandas as pd
import functools

from scipy.linalg import expm
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



class RS: # linear and quadratic reaction system
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
  def __init__(self, compartments, analytes, cells = None):
    self.compartments = compartments
    self.n_compartments = len(compartments)
    
    self.analytes = analytes
    self.n_analytes = len(analytes)
    
    cells = [] if cells is None else cells
    self.cells = cells
    self.n_cells = len(cells)
    
    # list of lists of analyte indices for each cell index
    self.ligands = [[] for _ in self.cells]
    for i, cell in enumerate(self.cells):
      for j, analyte in enumerate(self.analytes):
        if analyte.startswith(f"[{cell}]"):
          self.ligands[i].append(j)
    
    self.V = np.zeros(self.n_compartments, dtype = float) # volume of each compartment, in units.ml
    self.Q = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments], dtype = float) # flow matrix of analytes, in 1/units.h
    self.M = np.zeros([self.n_cells, self.n_compartments, self.n_compartments], dtype = float) # migration matrix of cells, in 1/units.h
    self.RS = [RS(self.n_analytes) for compartment in self.compartments]
    self.processes = []
    
    self.t = 0
    self.x = np.zeros([self.n_analytes, self.n_compartments], dtype = float) # concentration of analytes, in units.nM
    self.c = np.zeros([self.n_cells, self.n_compartments], dtype = float) # concentration of cells, in 1/units.ml
    
    self.history = []
    self.history_cells = []
  
  def get_volume(self, compartment):
    compartment = self.compartments.index(compartment)
    return self.V[compartment] * units.ml
  
  def set_volume(self, compartment, value):
    value = value.number(units.ml)
    compartment = self.compartments.index(compartment)
    self.V[compartment] = value
  
  # set compartment_dest as None if it is a clearance
  def add_flow(self, analyte, compartment_source, compartment_dest, rate):
    rate = rate.number(units.ml/units.h)
    analyte = self.analytes.index(analyte)
    compartment_source = self.compartments.index(compartment_source)
    self.Q[analyte, compartment_source, compartment_source] -= rate / self.V[compartment_source]
    if compartment_dest is not None:
      compartment_dest = self.compartments.index(compartment_dest)
      self.Q[analyte, compartment_source, compartment_dest] += rate / self.V[compartment_dest]
  
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
  
  def get_x(self, compartment, analyte):
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    return self.x[analyte, compartment] * units.nM
  
  def set_x(self, compartment, analyte, value):
    value = value.number(units.nM)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.x[analyte, compartment] = value
  
  def add_x(self, compartment, analyte, value):
    value = value.number(units.nM)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.x[analyte, compartment] += value
  
  def get_c(self, compartment, cell):
    cell = self.cells.index(cell)
    compartment = self.compartments.index(compartment)
    return self.c[cell, compartment] * 1/units.ml

  # adding a type of cell with certain ligands
  def add_c(self, compartment, cell, value, ligands, copys):
    for ligand, copy in zip(ligands, copy):
      analyte = f"{cell}-{ligand}"
      self.add_x(analyte, compartment, value * copy / units.avagadro)
    
    value = value.number(1/units.ml)
    cell = self.cells.index(cell)
    compartment = self.compartments.index(compartment)
    self.c[cell, compartment] += value
  
  def decay_c(self, compartment, cell, value):
    ligands = self.ligands[cell]
    self.x[ligands, compartment] *= 1 - value
    
    value = value.number(1/units.ml)
    cell = self.cells.index(cell)
    compartment = self.compartments.index(compartment)
    self.c[cell, compartment] *= 1 - value
  
  ### system running functions
  def run_flows(self, t):
    t = t.number(units.h)
    flowing_analytes = [analyte for analyte in range(self.n_analytes) if self.Q[analyte].any()]
    for analyte in flowing_analytes:
      self.x[analyte] = np.dot(self.x[analyte], expm(t * self.Q[analyte]))
    migrating_cells = [cell for cell in range(self.n_cells) if self.M[cell].any()]
    for cell in migrating_cells:
      ligands = self.ligands[cell]
      self.x[ligands] = np.dot(self.x[ligands],  expm(t * self.M[cell]))
      self.c[cell] = np.dot(self.c[cell], expm(t * self.M[cell]))
    
    self.t = self.t + t
    self.history_cells.append((self.t, self.x.copy(), self.c.copy()))
  
  def run_reactions(self, t):
    t = t.number(units.h)
    reacting_compartments = [compartment for compartment in range(self.n_compartments) if self.RS[compartment].active]
    for compartment in reacting_compartments:
      self.RS[compartment].refresh()
    for compartment in reacting_compartments:
      self.x[:, compartment] = self.RS[compartment](self.x[:, compartment], t)
    
    self.t = self.t + t
    self.history_cells.append((self.t, self.x.copy(), self.c.copy()))
  
  def run_processes(self, t):
    t = t.number(units.h)
    for process in self.processes:
      process(self, t * units.h)
    
    self.t = self.t + t
    self.history_cells.append((self.t, self.x.copy(), self.c.copy()))
  
  def run(self, t, t_step = 1/60 * units.h, t_record = 1 * units.h):
    t = t.number(units.h)
    t_start = self.t
    t_end = t_start + t
    t_step = t_step.number(units.h)
    t_record = t_record.number(units.h)
    flowing_analytes = [analyte for analyte in range(self.n_analytes) if self.Q[analyte].any()]
    migrating_cells = [cell for cell in range(self.n_cells) if self.M[cell].any()]
    reacting_compartments = [compartment for compartment in range(self.n_compartments) if self.RS[compartment].active]
    for compartment in reacting_compartments:
      self.RS[compartment].refresh()
    
    pbar = tqdm(total = t, unit = "h", bar_format = "{desc}: {percentage:3.0f}%|{bar}| {n:.2f}/{total_fmt} [{elapsed}<{remaining},  {rate_fmt}{postfix}]")
    pbar.update(0.0); A, B, C, D = 0.0, 0.0, 0.0, 0.0
    while True:
      t_prev = self.t
      self.t = min(self.t + t_step, t_end)
      t_delta = self.t - t_prev
      for analyte in flowing_analytes:
        A -= tt()
        self.x[analyte] = np.dot(self.x[analyte], expm(t_delta * self.Q[analyte]))
        A += tt()
      for cell in migrating_cells:
        A -= tt()
        ligands = self.ligands[cell]
        self.x[ligands] = np.dot(self.x[ligands], expm(t_delta * self.M[cell]))
        self.c[cell] = np.dot(self.c[cell], expm(t_delta * self.M[cell]))
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
        self.history_cells.append((self.t, self.x.copy(), self.c.copy()))
      pbar.update(t_delta)
      if math.isclose(self.t, t_end, rel_tol = 0, abs_tol = 1e-9):
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
  
  def summary_cells(self, cells):
    cells = [self.cells.index(cell) for cell in cells]
    avgs = []; maxs = []; hfws = []
    for compartment in range(self.n_compartments):
      X = np.array([t for t, x in self.history_cells])
      Y = np.array([x[cells, compartment].sum(axis = 0) for t, x in self.history_cells])
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
