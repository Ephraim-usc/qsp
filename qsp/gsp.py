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
    
    self.volumes = np.full([self.n_analytes, self.n_compartments], np.nan)
    self.flows = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments])
    self.reactions = []
  
  def set_volumes(self, volumes):
    volumes = np.array([[(i.number(units.l) if type(i) is unum.Unum else i) for i in row] for row in volumes])
    self.volumes = volumes
  
  def set_volume(self, compartment, analyte, volume):
    volume = volume.number(units.l) if type(volume) is unum.Unum else volume
    self.volumes[self.analytes.index(analyte), self.compartments.index(compartment)] = volume
  
  # set compartment_dest as None for clearance
  def add_flow(self, compartment_source, compartment_dest, analyte, coefficient):
    self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_source), ] -= coefficient
    if compartment_dest is None:
      return
    self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_dest)] += coefficient
  
  # powers: dict of powers for each analyte
  # deltas: dict of concentration changes of the analytes per reaction
  def add_reaction(self, compartment, coefficient, powers, deltas):
    self.reactions.append()
  
  def print(self):
    volumes = pd.DataFrame(self.volumes, index = self.compartments, columns = self.analytes)
    print("<volumes>", flush = True)
    print(volumes, flush = True)
    print(" ", flush = True)
    for i, analyte in enumerate(self.analytes):
      flows = pd.DataFrame(self.flows[:,:,i], index = self.compartments, columns = self.compartments)
      print(f"<flow of {analyte}>", flush = True)
      print(flows, flush = True)
      print(" ", flush = True)


compartments = ["central", "peripheral", "extracellular", "intracellular"]
analytes = ["adc", "drug"]
volumes = np.array([[0.084 * units.l, 0.136 * units.l], [0.051 * units.l, 0.523 * units.l], [np.nan, np.nan], [np.nan, np.nan]])

volumes = np.array([[0.084 * units.l, 0.136 * units.l], [0.051 * units.l, 0.523 * units.l], [np.nan, np.nan], [np.nan, np.nan]])

system = System(compartments, analytes)
system.set_volumes(volumes)
system.add_flow("central", "peripheral", "drug", 3.2)
