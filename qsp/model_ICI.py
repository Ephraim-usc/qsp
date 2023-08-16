from .qsp import *

### this model is mostly from (A Lindauer et al. 2016) and (Mark Stroh et al. 2019)

analytes = ["masked", "unmasked", "PD1", "PD1-masked", "PD1-unmasked", "FcRn", "FcRn-masked", "FcRn-unmasked"]
compartments = ["central", "peripheral", "tumor_plasma", "tumor_endosomal", "tumor_interstitial"]

mouse = []
mouse.update({"volume_central":1.26*units.ml, "volume_peripheral":0.819*units.ml})


MC38 = {}
MC38.update({"volume":170 * units.microliter})
MC38.update({"ratio_plasma":0.07, "ratio_endosomal":0.005, "ratio_interstitial":0.55})


def model(host, target, mask, tumor):
  system = System(analytes, compartments)
  
  for analyte in analytes:
    system.set_volume(analyte, "central", host["volume_central"])
    system.set_volume(analyte, "peripheral", host["volume_peripheral"])
    system.set_volume(analyte, "tumor_plasma", tumor["ratio_plasma"])
    system.set_volume(analyte, "tumor_endosomal", tumor["ratio_endosomal"])
    system.set_volume(analyte, "tumor_interstitial", tumor["ratio_interstitial"])
  
  # plasma and BC circles of adc and drug
  system.add_flow("adc", "plasma", "lung_plasma", host["plasma_flow_lung"])
  system.add_flow("drug", "plasma", "lung_plasma", host["plasma_flow_lung"])
  system.add_flow("adc", "BC", "lung_BC", host["BC_flow_lung"])
  system.add_flow("drug", "BC", "lung_BC", host["BC_flow_lung"])
  
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
