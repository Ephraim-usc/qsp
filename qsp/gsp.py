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
  def __init__(self, compartments, analytes):
    self.compartments = compartments
    self.n_compartments = len(compartments)
    self.analytes = analytes
    self.n_analytes = len(analytes)
    
    self.volumes = np.full([self.n_compartments, self.n_analytes], np.nan)
    self.flows = np.zeros([self.n_compartments, self.n_compartments], np.nan)
  
  def set_volume(self, compartment, volumes):
    self.compartments.append(compartment)
    self.volumes.append(volumes)
  
  def add_flow(self, analyte, compartment_a, compartment_b, coefficient):
    self.flows.append()
  
  # powers: dict of powers for each analyte
  # deltas: dict of concentration changes of the analytes per reaction
  def add_reaction(self, compartment, coefficient, powers, deltas):
    self.reactions.append()
  
  def print(self):
    volumes = pd.DataFrame(self.volumes, index = self.compartments, columns = self.analytes)
    print(volumes, flush = True)


def antigen_on(params):
  return 5 * params["adc"] * params["drug"]


compartments = ["central", "peripheral", "extracellular", "intracellular"]
analytes = ["adc", "drug"]

system = System(compartments, analytes)

