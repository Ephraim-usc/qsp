import math
import numpy as np
import pandas as pd

import unum
import unum.units as units

units.l = unum.new_unit('l', 1e-3 * units.m ** 3)
units.pl = unum.new_unit('pl', 1e-12 * units.l)
units.nM = unum.new_unit('nM', 1e-9 * units.mol / units.l)
units.avagadro = unum.new_unit('avagadro', 6.0221415e23 / units.mol)



class System:
  def __init__(self, analytes, compartments, volumes = None):
    self.analytes = analytes
    self.n_analytes = len(analytes)
    self.compartments = compartments
    self.n_compartments = len(compartments)
    
    self.volumes = np.zeros([self.n_analytes, self.n_compartments])
    self.flows = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments])
    self.reactions = []
  
  def set_volumes(self, volumes):
    volumes = np.array([[(i.number(units.l) if type(i) is unum.Unum else i) for i in row] for row in volumes])
    self.volumes = volumes
  
  def set_volume(self, analyte, compartment, volume):
    if type(volume) is unum.Unum:
      volume = volume.number(units.l)
    self.volumes[self.analytes.index(analyte), self.compartments.index(compartment)] = volume
  
  # set compartment_dest as None if it is a clearance
  def add_flow(self, analyte, compartment_source, compartment_dest, coefficient):
    if type(coefficient) is unum.Unum:
      coefficient = coefficient.number(units.l/units.h)
    
    self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_source), ] -= coefficient
    if compartment_dest is None:
      return
    self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_dest)] += coefficient
  
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
    str += "[" + " + ".join(inputs) + " → " + " + ".join(inputs) + "]"
    return str
  
  def print(self):
    volumes = pd.DataFrame(self.volumes, index = self.analytes, columns = self.compartments)
    print("<volumes>", flush = True)
    print(volumes, flush = True)
    print(" ", flush = True)
    for i, analyte in enumerate(self.analytes):
      flows = pd.DataFrame(self.flows[i,:,:], index = self.compartments, columns = self.compartments)
      print(f"<flow of {analyte}>", flush = True)
      print(flows, flush = True)
      print(" ", flush = True)
    print("<reactions>", flush = True)
    for reaction in self.reactions:
      print(self.str_reaction(reaction), flush = True)


compartments = ["central", "peripheral", "extracellular", "intracellular"]
analytes = ["adc", "drug"]
volumes = np.array([[0.084 * units.l, 0.051 * units.l, np.nan, np.nan], [0.136 * units.l, 0.523 * units.l, np.nan, np.nan]])

system = System(analytes, compartments)
system.set_volumes(volumes)

system.add_flow("adc", "central", None, 0.033 * (units.l/units.d))
system.add_flow("drug", "central", None, 18.4 * (units.l/units.d))
system.add_flow("adc", "central", "peripheral", 0.0585 * (units.l/units.d))

system.add_reaction(None, 0.323 * (1/units.d), {"adc":1}, {"adc":-1, "drug":1})



