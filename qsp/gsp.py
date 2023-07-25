import math
import numpy as np
import pandas as pd

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
    
    self.volumes = np.zeros([self.n_analytes, self.n_compartments])
    self.flows = np.zeros([self.n_analytes, self.n_compartments, self.n_compartments])
    self.reactions = []
  
  # if volumes is a vector, we assume all analytes share the same volume in each department
  def set_volumes(self, volumes):
    if volumes.ndim == 1:
      volumes = np.tile(volumes, reps = [self.n_analytes, 1])
    volumes = np.array([[(i.number(units.ml) if type(i) is unum.Unum else i) for i in row] for row in volumes])
    self.volumes = volumes
  
  def set_volume(self, analyte, compartment, volume):
    volume = volume.number(units.ml)
    self.volumes[self.analytes.index(analyte), self.compartments.index(compartment)] = volume
  
  def get_volume(self, analyte, compartment):
    return self.volumes[self.analytes.index(analyte), self.compartments.index(compartment)] * units.ml
  
  # set compartment_dest as None if it is a clearance
  def add_flow(self, analyte, compartment_source, compartment_dest, rate):
    rate = rate.number(units.ml/units.h)
    self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_source)] -= rate
    if compartment_dest is not None:
      self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_dest)] += rate
  
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


# two-compartment model
compartments = ["central", "peripheral", "extracellular", "intracellular"]
analytes = ["adc", "drug", "antigen", "substrate"]
system = System(analytes, compartments)

system.set_volume("adc", "central", 0.084 * units.l)
system.set_volume("adc", "peripheral", 0.051 * units.l)

system.set_volume("drug", "central", 0.136 * units.l)
system.set_volume("drug", "peripheral", 0.523 * units.l)

system.add_flow("adc", "central", None, 0.033 * (units.l/units.d))
system.add_flow("drug", "central", None, 18.4 * (units.l/units.d))
system.add_flow("adc", "central", "peripheral", 0.0585 * (units.l/units.d))

system.add_reaction(None, 0.323 * (1/units.d), {"adc":1}, {"adc":-1, "drug":1})


# Dhaval et al. model
organs = ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"]
plasma_flows = {"heart": 36.5 * units.ml/units.h, 
                "lung": 373 * units.ml/units.h,
                "muscle": 86.1 * units.ml/units.h,
                "skin": 27.8 * units.ml/units.h,
                "adipose": 13.4 * units.ml/units.h,
                "bone": 15.2 * units.ml/units.h,
                "brain": 11.8 * units.ml/units.h,
                "kidney": 68.5 * units.ml/units.h,
                "liver": 10.3 * units.ml/units.h,
                "SI": 58.1 * units.ml/units.h,
                "LI": 17.3 * units.ml/units.h,
                "pancreas": 6.24 * units.ml/units.h,
                "thymus": 1.19 * units.ml/units.h,
                "spleen": 8.18 * units.ml/units.h,
                "other": 10.9 * units.ml/units.h}
volumes = np.array([0.00585, 0.00479, 0.0217, 0.000760, 0.119,
                    0.0295, 0.0241, 0.0384, 0.00102, 0.111,
                    0.249, 0.204, 1.47, 0.0566, 9.34,
                    0.188, 0.154, 1.66, 0.0251, 3.00,
                    0.0218, 0.0178, 0.337, 0.00991, 1.60,
                    0.0621, 0.0508, 0.525, 0.0141, 2.17,
                    0.0107, 0.00873, 0.0873, 0.00243, 0.376,
                    0.0289, 0.0236, 0.0788, 0.00263, 0.391,
                    0.164, 0.134, 0.385, 0.00963, 1.23,
                    0.0116, 0.00950, 0.127, 0.00364, 0.577,
                    0.0050, 0.00409, 0.0545, 0.00157, 0.248,
                    0.00534, 0.00437, 0.0169, 0.000485, 0.0699,
                    0.0005, 0.000405, 0.00153, 0.00005, 0.00653, 
                    0.0154, 0.0126, 0.0254, 0.000635, 0.0730,
                    0.0195, 0.0160, 0.0797, 0.00233, 0.348,
                    0.944, 0.773, 0.113
                   ]) * units.ml
vascular_reflection_coefficients = {"heart": 0.95,
                                   "lung": 0.95,
                                   "muscle": 0.95,
                                   "skin": 0.95,
                                   "adipose": 0.95,
                                   "bone": 0.85,
                                   "brain": 0.99,
                                   "kidney": 0.9,
                                   "liver": 0.85,
                                   "SI": 0.9,
                                   "LI": 0.95,
                                   "pancreas": 0.9,
                                   "thymus": 0.9,
                                   "spleen": 0.85,
                                   "other": 0.95}
lymphatic_reflection_coefficient = 0.2
endosomal_pinocytosis_rate = 3.66e-2 / units.h
endosomal_degradation_rate = 42.9 / units.h

compartments = [f"{organ}_{tissue}" for organ in organs for tissue in ["plasma", "BC", "interstitial", "endosomal", "cellular"]] + ["plasma", "BC", "lymph"]
analytes = ["T-vc-MMAE", "MMAE", "HER2", "tubulin"]
system = System(analytes, compartments)
system.set_volumes(volumes)


# plasma circle
system.add_flow("T-vc-MMAE", "plasma", "lung_plasma", 373 * units.ml / units.h)

for organ in [organ for organ in organs if organ not in ["lung"]]:
  system.add_flow("T-vc-MMAE", "lung_plasma", f"{organ}_plasma", plasma_flows[organ])

for organ in ["SI", "LI", "spleen", "pancreas"]:
  system.add_flow("T-vc-MMAE", f"{organ}_plasma", "liver_plasma", plasma_flows[organ] * (499/500))
  system.add_flow("T-vc-MMAE", "liver_plasma", "plasma", plasma_flows[organ] * (499/500))

for organ in [organ for organ in organs if organ not in ["lung"] + ["SI", "LI", "spleen", "pancreas"]]:
  system.add_flow("T-vc-MMAE", f"{organ}_plasma", "plasma", plasma_flows[organ] * (499/500))


# lymph circle
for organ in organs:
  system.add_flow("T-vc-MMAE", f"{organ}_plasma", f"{organ}_interstitial", plasma_flows[organ] * (1/500) * (1 - vascular_reflection_coefficients[organ]))
  system.add_flow("T-vc-MMAE", f"{organ}_interstitial", "lymph", plasma_flows[organ] * (1/500) * (1 - lymphatic_reflection_coefficient))

system.add_flow("T-vc-MMAE", "lymph", "plasma", 373*(1/9.1) * units.ml / units.h)


# endosomal degradation
for organ in organs:
  system.add_flow("T-vc-MMAE", f"{organ}_plasma", f"{organ}_endosomal", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal"))
  system.add_flow("T-vc-MMAE", f"{organ}_endosomal", f"{organ}_plasma", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal") * 0.715)
  
  system.add_flow("T-vc-MMAE", f"{organ}_interstitial", f"{organ}_endosomal", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal"))
  system.add_flow("T-vc-MMAE", f"{organ}_endosomal", f"{organ}_interstitial", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal") * (1-0.715))
  
  system.add_flow("T-vc-MMAE", f"{organ}_endosomal", None, endosomal_degradation_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal"))




