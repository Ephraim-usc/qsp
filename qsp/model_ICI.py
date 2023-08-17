from .qsp import *

### this model is mostly from (A Lindauer et al. 2016) and (Mark Stroh et al. 2019)

analytes = ["masked", "unmasked", "PD1", "PD1-masked", "PD1-unmasked", "FcRn", "FcRn-masked", "FcRn-unmasked"]
compartments = ["central", "peripheral", "tumor_plasma", "tumor_endosomal", "tumor_interstitial"]

mouse = []
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"distribution": 4.82 * units.ml/units.d})
mouse.update({"clearance": 0.334 * units.ml/units.d})
mouse.update({"max_nonlinear_clearance": 0.518 * units.microgram/units.d, "max_nonlinear_clearance": 0.366 * units.microgram/units.ml})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
mouse.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
mouse.update({"internalization": 0.0194/ units.h})


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
  
  # endosomal take-up
  for analyte in ["masked", "unmasked"]:
    system.add_flow(analyte, "tumor_endosomal", None, tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_degradation"])
    
    system.add_flow(analyte, "tumor_plasma", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
    system.add_flow(analyte, "tumor_interstitial", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
    system.add_flow(analyte, "tumor_endosomal", "tumor_plasma", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"] * host["vascular_recycle"])
    system.add_flow(analyte, "tumor_endosomal", "tumor_interstitial", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"] * (1 - host["vascular_recycle"]))

  # complex internalization
  for analyte in ["PD1-masked", "PD1-unmasked"]:
    system.add_flow(analyte, "central", None, host["volume_central"] * host["internalization"])
    system.add_flow(analyte, "tumor_interstitial", None, tumor["volume"] * tumor["volume_interstitial_proportion"] * host["internalization"])
  
  # nonlinear clearance reaction
  # PD1 association reaction
  
  
  for organ in [organ for organ in organs if organ not in ["lung"]]:
    system.add_flow("adc", "lung_plasma", f"{organ}_plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("drug", "lung_plasma", f"{organ}_plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("adc", "lung_BC", f"{organ}_BC", host[f"BC_flow_{organ}"])
    system.add_flow("drug", "lung_BC", f"{organ}_BC", host[f"BC_flow_{organ}"])
  
  for organ in ["SI", "LI", "spleen", "pancreas"]:
    system.add_flow("adc", f"{organ}_plasma", "liver_plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
    system.add_flow("drug", f"{organ}_plasma", "liver_plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("adc", f"{organ}_BC", "liver_BC", host[f"BC_flow_{organ}"])
    system.add_flow("drug", f"{organ}_BC", "liver_BC", host[f"BC_flow_{organ}"])
    
    system.add_flow("adc", "liver_plasma", "plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
    system.add_flow("drug", "liver_plasma", "plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("adc", "liver_BC", "BC", host[f"BC_flow_{organ}"])
    system.add_flow("drug", "liver_BC", "BC", host[f"BC_flow_{organ}"])
  
  for organ in [organ for organ in organs if organ not in ["lung"] + ["SI", "LI", "spleen", "pancreas"]]:
    system.add_flow("adc", f"{organ}_plasma", "plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
    system.add_flow("drug", f"{organ}_plasma", "plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("adc", f"{organ}_BC", "BC", host[f"BC_flow_{organ}"])
    system.add_flow("drug", f"{organ}_BC", "BC", host[f"BC_flow_{organ}"])
  
  # lymphatic circle of adc
  for organ in organs:
    system.add_flow("adc", f"{organ}_plasma", f"{organ}_interstitial", host[f"lymphatic_flow_{organ}"] * (1 - host[f"vascular_reflection_{organ}"]))
    system.add_flow("adc", f"{organ}_interstitial", "lymph", host[f"lymphatic_flow_{organ}"] * (1 - host["lymphatic_reflection"]))
    system.add_flow("adc", "lymph", "plasma", host[f"lymphatic_flow_{organ}"]) # the article used a different measure, which I think is weird
  
  # endosomal take-up and plasma recycle
  for organ in organs:
    v = system.get_volume("adc", f"{organ}_endosomal")
    system.add_flow("adc", f"{organ}_plasma", f"{organ}_endosomal", host["endosomal_pinocytosis"] * v)
    system.add_flow("adc", f"{organ}_interstitial", f"{organ}_endosomal", host["endosomal_pinocytosis"] * v)
    system.add_flow("adc", f"{organ}_endosomal", f"{organ}_plasma", host["endosomal_pinocytosis"] * v * host["plasma_recycle"])
    system.add_flow("adc", f"{organ}_endosomal", f"{organ}_interstitial", host["endosomal_pinocytosis"] * v * (1 - host["plasma_recycle"]))
  
  # target binding and internalization
  for organ in organs:
    v = system.get_volume("adc", f"{organ}_interstitial")
    system.add_flow("adc", f"{organ}_interstitial", f"{organ}_membrane", target["on"] * (target[f"num_{organ}"] * host[f"cell_density_{organ}"] / units.avagadro) * v)
    system.add_flow("adc", f"{organ}_membrane", f"{organ}_interstitial", target["off"] * v)
    system.add_flow("adc", f"{organ}_membrane", f"{organ}_cellular", target["int"] * v)
  
  # fast drug equilibrium between tissues
  for organ in organs:
    system.add_flow("drug", f"{organ}_plasma", f"{organ}_endosomal", host[f"plasma_flow_{organ}"] * drug["unbound_plasma"] * 1000)
    system.add_flow("drug", f"{organ}_interstitial", f"{organ}_endosomal", host[f"plasma_flow_{organ}"] * 1000)
    system.add_flow("drug", f"{organ}_endosomal", f"{organ}_plasma", host[f"plasma_flow_{organ}"] * 1000)
    system.add_flow("drug", f"{organ}_endosomal", f"{organ}_interstitial", host[f"plasma_flow_{organ}"] * 1000)
  
  # BC and cellular permeation
  v = host["volume_BC"]
  system.add_flow("drug", "plasma", "BC", drug["permeability_BC"] * v * drug["unbound_plasma"])
  system.add_flow("drug", "BC", "plasma", drug["permeability_BC"] * v * drug["unbound_BC"])
  
  for organ in organs:
    v = host[f"volume_{organ}_BC"]
    system.add_flow("drug", f"{organ}_plasma", f"{organ}_BC", drug["permeability_BC"] * v * drug["unbound_plasma"])
    system.add_flow("drug", f"{organ}_BC", f"{organ}_plasma", drug["permeability_BC"] * v * drug["unbound_BC"])
  
  for organ in organs:
    v = host[f"volume_{organ}_cellular"]
    system.add_flow("drug", f"{organ}_interstitial", f"{organ}_cellular", drug[f"permeability_{organ}"] * v)
    system.add_flow("drug", f"{organ}_cellular", f"{organ}_interstitial", drug[f"permeability_{organ}"] * v * drug[f"unbound_{organ}"])
  
  # liver clearance
  v = host["volume_liver_interstitial"]
  system.add_flow("drug", "liver_interstitial", None, drug["liver_clearance"] * v)
  
  # dissociatoin and degradation
  def degradation_endosomal(x, z):
    DAR = z["DAR"]
    rate = linker["degradation_endosomal"] * x["adc"]
    return {"adc": -rate, "drug": DAR * rate}
  
  def degradation_cellular(x, z):
    DAR = z["DAR"]
    rate = linker["degradation_cellular"] * x["adc"]
    return {"adc": -rate, "drug": DAR * rate}
  
  def dissociation(x, z):
    DAR = z["DAR"]
    if callable(linker["dissociation"]):
      rate = linker["dissociation"](DAR) * x["adc"]
    else:
      rate = linker["dissociation"] * x["adc"]
    return {"adc": -rate, "drug": DAR * rate}
  
  def DAR_decay(z):
    DAR = z["DAR"]
    if callable(linker["dissociation"]):
      rate = linker["dissociation"](DAR)
    else:
      rate = linker["dissociation"]
    return {"DAR":-DAR * rate}
  
  for organ in organs:
    system.add_reaction(f"{organ}_endosomal", degradation_endosomal)
    system.add_reaction(f"{organ}_cellular", degradation_cellular)
  
  system.add_reaction("plasma", dissociation)
  system.add_process(DAR_decay)
  
  return system