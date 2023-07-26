import math
import numpy as np
import pandas as pd

np.set_printoptions(suppress=True)

import matplotlib as mpl
import matplotlib.pyplot as plt

import unum
import unum.units as units

units.l = unum.new_unit('l', 1e-3 * units.m ** 3)
units.ml = unum.new_unit('ml', 1e-3 * units.l)
units.pl = unum.new_unit('pl', 1e-12 * units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)
units.avagadro = unum.new_unit('avagadro', 6.0221415e23 / units.mol)



class System:
  def __init__(self, analytes, compartments, volumes = None):
    self.analytes = analytes
    self.n_analytes = len(analytes)
    self.compartments = compartments
    self.n_compartments = len(compartments)
    
    self.x = np.zeros([self.n_analytes, self.n_compartments])
    self.volumes = np.zeros([self.n_analytes, self.n_compartments])
    self.flows = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments])
    self.flowing_analytes = []
    self.reactions = []
  
  def clear_x(self):
    self.x = np.zeros([self.n_analytes, self.n_compartments])
  
  def set_x(self, analyte, compartment, concentration):
    concentration = concentration.number(units.nM)
    analyte = self.analytes.index(analyte)
    compartment = self.compartments.index(compartment)
    self.x[analyte, compartment] = concentration
  
  def step(self, t):
    t = t.number(units.h)
    for analyte in self.flowing_analytes:
      x = self.x[analyte, :]
      flows = self.flows[analyte, :, :]
      self.x[analyte, :] += t * np.dot(x, flows)
  
  def run(self, t, t_step = 0.001 * units.h, t_record = 1 * units.h):
    t = t.number(units.h)
    t_step = t_step.number(units.h)
    t_record = t_record.number(units.h)
    
    records = [self.x.copy()]
    for t_ in np.arange(t_step, t, t_step):
      self.step(t_step * units.h)
      if t_ / t_record >= len(records):
        records.append(self.x.copy())
    
    records = np.dstack(records) # 3d array of size n_analytes x n_compartments x n_ts
    return records
  
  def plot(self, analyte, compartments, records):
    analyte = self.analytes.index(analyte)
    compartments = [self.compartments.index(compartment) for compartment in compartments]
    fig, axs = plt.subplots(nrows = 1, ncols = len(compartments), squeeze = False)
    axs = axs.ravel().tolist()
    for ax, compartment in zip(axs, compartments):
      ax.plot(records[analyte,compartment,:])
      ax.set_yscale('symlog', linthresh = 1)
      ax.set_yticks([0, 1, 10, 100, 1000, 10000])
      ax.set_title(self.compartments[compartment])
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
    if analyte not in self.flowing_analytes:
      self.flowing_analytes.append(analyte)
      self.flowing_analytes.sort()
    self.flows[analyte, compartment_source, compartment_source] -= rate
    if compartment_dest is not None:
      compartment_dest = self.compartments.index(compartment_dest)
      self.flows[analyte, compartment_source, compartment_dest] += rate
  
  # set compartment as None if it happens everywhere
  # powers: dict of powers for each analyte
  # deltas: dict of concentration changes of the analytes per reaction
  def add_reaction(self, compartment, coefficient, powers, deltas):
    powers_ = np.zeros([self.n_analytes])
    for analyte, power in powers.items():
      powers_[self.analytes.index(analyte)] += power
    n_powers = powers_.sum()
    
    deltas_ = np.zeros([self.n_analytes])
    for analyte, delta in deltas.items():
      deltas_[self.analytes.index(analyte)] += delta
    
    if type(coefficient) is unum.Unum:
      coefficient = coefficient.number(1/units.h/(units.nM ** (n_powers-1)))
    
    if compartment is not None:
      compartment = self.compartments.index(compartment)
    self.reactions.append([compartment, coefficient, powers_, deltas_])
  
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
