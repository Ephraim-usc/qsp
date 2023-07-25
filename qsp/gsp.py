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
  
  def set_volumes(self, volumes):
    volumes = np.array([[(i.number(units.ml) if type(i) is unum.Unum else i) for i in row] for row in volumes])
    self.volumes = volumes
  
  def set_volume(self, analyte, compartment, volume):
    if type(volume) is unum.Unum:
      volume = volume.number(units.ml)
    self.volumes[self.analytes.index(analyte), self.compartments.index(compartment)] = volume
  
  # set compartment_dest as None if it is a clearance
  def add_flow(self, analyte, compartment_source, compartment_dest, coefficient):
    if type(coefficient) is unum.Unum:
      coefficient = coefficient.number(units.ml/units.h)
    
    self.flows[self.analytes.index(analyte), self.compartments.index(compartment_source), self.compartments.index(compartment_source)] -= coefficient
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
compartments = [f"{organ}_{tissue}" for organ in ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"] for tissue in ["plasma", "BC", "interstitial", "endosomal", "cellular"]]
compartments += ["plasma", "BC", "lymph"]
analytes = ["T-vc-MMAE", "MMAE", "FcRn", "HER2", "tubulin"]
system = System(analytes, compartments)

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
volumes = np.tile(volumes, reps = [len(analytes),1])
system.set_volumes(volumes)

# plasma to lung plasma flow
system.add_flow("T-vc-MMAE", "plasma", "lung_plasma", 373 * units.ml / units.h)

# lung plasma to tissue plasma flow
system.add_flow("T-vc-MMAE", "lung_plasma", "heart_plasma", 36.5 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "muscle_plasma", 86.1 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "skin_plasma", 27.9 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "adipose_plasma", 13.4 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "bone_plasma", 15.2 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "brain_plasma", 11.8 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "kidney_plasma", 68.5 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "liver_plasma", 10.3 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "thymus_plasma", 1.19 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "other_plasma", 10.9 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "SI_plasma", 58.1 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "LI_plasma", 17.3 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "pancreas_plasma", 6.24 * units.ml / units.h)
system.add_flow("T-vc-MMAE", "lung_plasma", "spleen_plasma", 8.18 * units.ml / units.h)

# SLSP plasma to liver plasma flow
system.add_flow("T-vc-MMAE", "SI_plasma", "liver_plasma", 58.1*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "LI_plasma", "liver_plasma", 17.3*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "spleen_plasma", "liver_plasma", 8.18*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "pancreas_plasma", "liver_plasma", 6.24*(499/500) * units.ml / units.h)

# other tissue plasma to plasma flow
system.add_flow("T-vc-MMAE", "heart_plasma", "plasma", 36.5*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "muscle_plasma", "plasma", 86.1*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "skin_plasma", "plasma", 27.8*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "adipose_plasma", "plasma", 13.4*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "bone_plasma", "plasma", 15.2*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "brain_plasma", "plasma", 11.8*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "kidney_plasma", "plasma", 68.5*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "liver_plasma", "plasma", (10.3 + 58.1 + 17.3 + 8.18 + 6.24)*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "thymus_plasma", "plasma", 1.19*(499/500) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "other_plasma", "plasma", 10.9*(499/500) * units.ml / units.h)

# tissue interstitial to lymph flow
system.add_flow("T-vc-MMAE", "lung_interstitial", "lymph", 373*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "heart_interstitial", "lymph", 36.5*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "muscle_interstitial", "lymph", 86.1*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "skin_interstitial", "lymph", 27.8*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "adipose_interstitial", "lymph", 13.4*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "bone_interstitial", "lymph", 15.2*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "brain_interstitial", "lymph", 11.8*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "kidney_interstitial", "lymph", 68.5*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "liver_interstitial", "lymph", 10.3*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "thymus_interstitial", "lymph", 1.19*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "SI_interstitial", "lymph", 58.1*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "LI_interstitial", "lymph", 17.3*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "pancreas_interstitial", "lymph", 6.24*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "spleen_interstitial", "lymph", 8.18*(1/500)*(1-0.2) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "other_interstitial", "lymph", 10.9*(1/500)*(1-0.2) * units.ml / units.h)

# lymph to plasma flow (should be what the authors mean, but note that this is much faster than the input into lymph)
system.add_flow("T-vc-MMAE", "lymph", "plasma", 373*(1/9.1) * units.ml / units.h)

# tissue plasma to interstitial flow
system.add_flow("T-vc-MMAE", "lung_plasma", "lung_interstitial", 373*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "heartg_plasma", "heart_interstitial", 36.5*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "muscle_plasma", "muscle_interstitial", 86.1*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "skin_plasma", "skin_interstitial", 27.8*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "adipose_plasma", "adipose_interstitial", 13.4*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "bone_plasma", "bone_interstitial", 15.2*(1/500)*(1-0.85) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "brain_plasma", "brain_interstitial", 11.8*(1/500)*(1-0.99) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "kidney_plasma", "kidney_interstitial", 68.5*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "liver_plasma", "liver_interstitial", 10.3*(1/500)*(1-0.85) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "thymus_plasma", "thymus_interstitial", 1.19*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "SI_plasma", "SI_interstitial", 58.1*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "LI_plasma", "LI_interstitial", 17.3*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "pancreas_plasma", "pancreas_interstitial", 6.24*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "spleen_plasma", "spleen_interstitial", 8.18*(1/500)*(1-0.85) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "other_plasma", "other_interstitial", 10.9*(1/500)*(1-0.95) * units.ml / units.h)

# endosomal take-up
system.add_flow("T-vc-MMAE", "lung_plasma", "lung_interstitial", 373*(1/500)*(1-0.95) * units.ml / units.l)
system.add_flow("T-vc-MMAE", "heartg_plasma", "heart_interstitial", 36.5*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "muscle_plasma", "muscle_interstitial", 86.1*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "skin_plasma", "skin_interstitial", 27.8*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "adipose_plasma", "adipose_interstitial", 13.4*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "bone_plasma", "bone_interstitial", 15.2*(1/500)*(1-0.85) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "brain_plasma", "brain_interstitial", 11.8*(1/500)*(1-0.99) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "kidney_plasma", "kidney_interstitial", 68.5*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "liver_plasma", "liver_interstitial", 10.3*(1/500)*(1-0.85) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "thymus_plasma", "thymus_interstitial", 1.19*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "SI_plasma", "SI_interstitial", 58.1*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "LI_plasma", "LI_interstitial", 17.3*(1/500)*(1-0.95) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "pancreas_plasma", "pancreas_interstitial", 6.24*(1/500)*(1-0.9) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "spleen_plasma", "spleen_interstitial", 8.18*(1/500)*(1-0.85) * units.ml / units.h)
system.add_flow("T-vc-MMAE", "other_plasma", "other_interstitial", 10.9*(1/500)*(1-0.95) * units.ml / units.h)



