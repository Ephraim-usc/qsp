from .qsp import *

### this model is mostly from (A Lindauer et al. 2016) and (Mark Stroh et al. 2019)

analytes = ["bimasked", "monomasked", "unmasked", "target", "target-bimasked", "target-monomasked", "target-unmasked", "FcRn", "FcRn-bimasked", "FcRn-monomasked", "FcRn-unmasked"]
compartments = ["central", "peripheral", "tumor_plasma", "tumor_endosomal", "tumor_interstitial"]


def nonlinear_clearance(system, t):
  molecular_weight = 150000 * units.g/units.mol
  MAX = 0.518 * units.microgram/units.d / molecular_weight
  EC50 = 0.366 * units.microgram/units.ml / molecular_weight
  
  for analyte in ["bimasked", "monomasked", "unmasked"]:
    x = system.get_x(analyte, "central")
    rate = - MAX * x / (x + EC50) / system.get_volume(analyte, "central")
    system.add_x(analyte, "central", rate * t)

mouse = {}
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"distribution": 4.82 * units.ml/units.d})
mouse.update({"clearance": 0.334 * units.ml/units.d})
mouse.update({"nonlinear_clearance": nonlinear_clearance})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
mouse.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
mouse.update({"FcRn": 49.8 * units.micromolar, "FcRn-on": 0.0806 * 1/units.nM/units.d, "FcRn-off": 6.55 / units.h})


def PD1_dynamics(system, t):
  death = 0.02 / units.d
  MAX = 94.7
  EC50 = 1.46 * units.nM
  
  for compartment in ["central", "tumor_interstitial"]:
    x_target = system.get_x("target", compartment)
    x_complex = system.get_x("target-bimasked", compartment) + system.get_x("target-monomasked", compartment) + system.get_x("target-unmasked", compartment)
    rate = death * PD1[compartment] * (1 + MAX * x_complex / (EC50 + x_complex)) - death * x_target
    system.add_x("target", compartment, rate * t)


PD1 = {}
PD1.update({"central": (1e4 * 1000/units.microliter) / units.avagadro}); PD1.update({"tumor_interstitial": PD1["central"] * 4.3})
PD1.update({"dynamics": PD1_dynamics})
PD1.update({"on": 0.34 * 1/units.nM/units.d, "off": 0.106 / units.h, "internalization": 0.0194/ units.h})


Tx = {}
Tx.update({"foldchange": 1/57, "cleavage_central": 0.0527 / units.d, "cleavage_tumor": 0.1783 / units.d})


MC38 = {}
MC38.update({"volume": 170 * units.microliter})
MC38.update({"volume_plasma_proportion": 0.07, "volume_endosomal_proportion": 0.005, "volume_interstitial_proportion": 0.55})
MC38.update({"plasma_flow_density": 12.7 / units.h})
MC38.update({"lymphatic_flow_ratio": 0.002})


def model(host, target, mask, tumor):
  system = System(analytes, compartments)
  
  for analyte in analytes:
    system.set_volume(analyte, "central", host["volume_central"])
    system.set_volume(analyte, "peripheral", host["volume_peripheral"])
    system.set_volume(analyte, "tumor_plasma", tumor["volume"] * tumor["volume_plasma_proportion"])
    system.set_volume(analyte, "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"])
    system.set_volume(analyte, "tumor_interstitial", tumor["volume"] * tumor["volume_interstitial_proportion"])
  
  # distribution and clearance
  for analyte in ["bimasked", "monomasked", "unmasked"]:
    system.add_flow(analyte, "central", None, host["clearance"])
    system.add_flow(analyte, "central", "peripheral", host["distribution"])
    system.add_flow(analyte, "peripheral", "central", host["distribution"])
  
  system.add_process(host["nonlinear_clearance"])
  
  # tumor plasma and lymphatic flow
  for analyte in ["bimasked", "monomasked", "unmasked"]:
    system.add_flow(analyte, "central", "tumor_plasma", tumor["volume"] * tumor["plasma_flow_density"])
    system.add_flow(analyte, "tumor_plasma", "central", tumor["volume"] * tumor["plasma_flow_density"] * (1 - tumor["lymphatic_flow_ratio"]))
    
    system.add_flow(analyte, "tumor_plasma", "tumor_interstitial", tumor["volume"] * tumor["plasma_flow_density"] * tumor["lymphatic_flow_ratio"] * (1 - host["vascular_reflection"]))
    system.add_flow(analyte, "tumor_interstitial", "central", tumor["volume"] * tumor["plasma_flow_density"] * tumor["lymphatic_flow_ratio"] * (1 - host["lymphatic_reflection"]))
  
  # endosomal take-up, degradation, and recycle
  for analyte in ["bimasked", "monomasked", "unmasked"]:
    system.add_flow(analyte, "tumor_plasma", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
    system.add_flow(analyte, "tumor_interstitial", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
    
    system.add_flow(analyte, "tumor_endosomal", None, tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_degradation"])
    
    system.add_reaction("tumor_endosomal", {f"{analyte}":1, "FcRn":1}, {f"FcRn-{analyte}":1}, host["FcRn-on"], backward = host["FcRn-off"])
    system.add_reaction("tumor_endosomal", {f"FcRn-{analyte}":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * host["vascular_recycle"], side_compartment = "tumor_plasma", side_products = {f"{analyte}":1})
    system.add_reaction("tumor_endosomal", {f"FcRn-{analyte}":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * (1 - host["vascular_recycle"]), side_compartment = "tumor_interstitial", side_products = {f"{analyte}":1})
  
  # cleavage of mask
  system.add_reaction("central", {"bimasked":1}, {"monomasked":1}, 2 * mask["cleavage_central"])
  system.add_reaction("central", {"monomasked":1}, {"unmasked":1}, mask["cleavage_central"])
  
  system.add_reaction("tumor_interstitial", {"bimasked":1}, {"monomasked":1}, 2 * mask["cleavage_tumor"])
  system.add_reaction("tumor_interstitial", {"monomasked":1}, {"unmasked":1}, mask["cleavage_tumor"])
  
  # PD1 dynamics
  system.add_process(target["dynamics"])
  
  # PD1 association
  for compartment in ["central", "tumor_interstitial"]:
    system.add_reaction(compartment, {"bimasked":1, "target":1}, {"target-bimasked":1}, target["on"] * mask["foldchange"], target["off"])
    system.add_reaction(compartment, {"monomasked":1, "target":1}, {"target-monomasked":1}, target["on"] * (0.5 + 0.5 * mask["foldchange"]), target["off"])
    system.add_reaction(compartment, {"unmasked":1, "target":1}, {"target-unmasked":1}, target["on"], target["off"])
  
  # complex internalization
  for analyte in ["target-bimasked", "target-monomasked", "target-unmasked"]:
    system.add_flow(analyte, "central", None, host["volume_central"] * target["internalization"])
    system.add_flow(analyte, "tumor_interstitial", None, tumor["volume"] * tumor["volume_interstitial_proportion"] * target["internalization"])
  
  # initial concentrations
  system.set_x("target", "central", target["central"])
  system.set_x("target", "tumor_interstitial", target["tumor_interstitial"])
  system.set_x("FcRn", "tumor_endosomal", host["FcRn"])
  
  return system
