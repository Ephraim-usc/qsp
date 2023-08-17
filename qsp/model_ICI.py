from .qsp import *

### this model is mostly from (A Lindauer et al. 2016) and (Mark Stroh et al. 2019)

analytes = ["masked", "unmasked", "target", "target-masked", "target-unmasked", "FcRn", "FcRn-masked", "FcRn-unmasked"]
compartments = ["central", "peripheral", "tumor_plasma", "tumor_endosomal", "tumor_interstitial"]

mouse = []
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"distribution": 4.82 * units.ml/units.d})
mouse.update({"clearance": 0.334 * units.ml/units.d})
mouse.update({"max_nonlinear_clearance": 0.518 * units.microgram/units.d, "max_nonlinear_clearance": 0.366 * units.microgram/units.ml})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
mouse.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
mouse.update({"FcRn": 49.8 * units.micromolar, "FcRn-on": 0.0806 1/units.nM/units.d, "FcRn-off": 6.55 / units.h})

PD1 = {}
PD1.update({"on": 0.34 * 1/units.nM/units.d, "off": 0.106 / units.h, "internalization": 0.0194/ units.h})

Tx = {}
Tx.update({"foldchange": 0.01, "cleavage_central": 0.05 / units.d, "cleavage_tumor": 0.2 / units.d})

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
  for analyte in ["masked", "unmasked"]:
    system.add_flow(analyte, "central", None, host["clearance"])
    system.add_flow(analyte, "central", "peripheral", host["distribution"])
    system.add_flow(analyte, "peripheral", "central", host["distribution"])

  # tumor plasma and lymphatic flow
  for analyte in ["masked", "unmasked"]:
    system.add_flow(analyte, "central", "tumor_plasma", tumor["volume"] * tumor["plasma_flow_density"])
    system.add_flow(analyte, "tumor_plasma", "central", tumor["volume"] * tumor["plasma_flow_density"] * (1 - tumor["lymphatic_flow_ratio"]))
    
    system.add_flow(analyte, "tumor_plasma", "tumor_interstitial", tumor["volume"] * tumor["plasma_flow_density"] * tumor["lymphatic_flow_ratio"] * (1 - tumor["vascular_reflection"]))
    system.add_flow(analyte, "tumor_interstitial", "central", tumor["volume"] * tumor["plasma_flow_density"] * tumor["lymphatic_flow_ratio"] * (1 - tumor["lymphatic_reflection"]))
  
  # endosomal take-up, degradation, and recycle
  for analyte in ["masked", "unmasked"]:
    system.add_flow(analyte, "tumor_plasma", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
    system.add_flow(analyte, "tumor_interstitial", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
    
    system.add_flow(analyte, "tumor_endosomal", None, tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_degradation"])
    
    system.add_reaction("tumor_endosomal", {f"{analyte}":1, "FcRn":1}, {f"FcRn-{analyte}":1}, host["FcRn-on"], backward = target["FcRn-off"])
    system.add_reaction("tumor_endosomal", {f"FcRn-{analyte}":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * host["vascular_recycle"], side_compartment = "tumor_plasma", side_products = {f"{analyte}":1})
    system.add_reaction("tumor_endosomal", {f"FcRn-{analyte}":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * (1 - host["vascular_recycle"]), side_compartment = "tumor_interstitial", side_products = {f"{analyte}":1})
  
  # PD1 dynamics
  def rate_func(system):
    x_target = system.get_x("target")
    x_complex = system.get_x("target-masked", "tumor_interstitial") + system.get_x("target-unmasked", "tumor_interstitial")
    rate = target["growth"] * (1 + target["TPemax"] * x_complex / (target["TPec50"] + x_complex)) - target["exhaustion"] * x_target
    return rate
  
  def change_func(system, t):
    system.x[system.analytes.index("PD1"), system.compartments.index("tumor_interstitial")] += (delta * t * units.h).number(units.nM)
  
  system.add_process(Process(rate_func, change_func))
  
  # PD1 association
  system.add_reaction("central", {"masked":1, "target":1}, {"target-masked":1}, target["on"] * mask["foldchange"], target["off"])
  system.add_reaction("central", {"unmasked":1, "target":1}, {"target-unmasked":1}, target["on"], target["off"])
  system.add_reaction("tumor_interstitial", {"masked":1, "target":1}, {"target-masked":1}, target["on"] * mask["foldchange"], target["off"])
  system.add_reaction("tumor_interstitial", {"unmasked":1, "target":1}, {"target-unmasked":1}, target["on"], target["off"])

  # complex internalization
  for analyte in ["target-masked", "target-unmasked"]:
    system.add_flow(analyte, "central", None, host["volume_central"] * host["internalization"])
    system.add_flow(analyte, "tumor_interstitial", None, tumor["volume"] * tumor["volume_interstitial_proportion"] * target["internalization"])
  
  # nonlinear clearance reaction
  
  
  return system
