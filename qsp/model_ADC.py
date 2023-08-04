organs = ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"]
tissues = ["plasma", "BC", "interstitial", "endosomal", "membrane", "cellular"]
centrals = ["plasma", "BC", "lymph"]
compartments = centrals + [f"{organ}_{tissue}" for organ in organs for tissue in tissues]

def model(host, target, linker, drug, DAR):
  analytes = ["adc", "drug"]
  system = System(analytes, compartments)

  for analyte in analytes:
    for compartment in compartments:
      system.set_volume(analyte, compartment, host[f"volume_{compartment}"])
  
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
  
  # blood cell up-take
  system.add_flow("drug", "plasma", "BC", drug["permeability_BC"] * unbound_proportion_plasma)
  system.add_flow("drug", "BC", "plasma", permeability_surface_area_coefficient_BC * PS * unbound_proportion_BC)
  
  for organ in organs:
    system.add_flow("MMAE", f"{organ}_plasma", f"{organ}_BC", permeability_surface_area_coefficient_BC * PS * unbound_proportion_plasma)
    system.add_flow("MMAE", f"{organ}_BC", f"{organ}_plasma", permeability_surface_area_coefficient_BC * PS * unbound_proportion_BC)
  
  
  # organ plasma to interstitial flow (here I do not follow the authors in multiplying plasma flow by 1000)
  for organ in organs:
    system.add_flow("MMAE", f"{organ}_plasma", f"{organ}_endosomal", plasma_flows[organ] * unbound_proportion_plasma * 1000)
    system.add_flow("MMAE", f"{organ}_endosomal", f"{organ}_plasma", plasma_flows[organ] * 1000)
    
    system.add_flow("MMAE", f"{organ}_endosomal", f"{organ}_interstitial", plasma_flows[organ] * 1000)
    system.add_flow("MMAE", f"{organ}_interstitial", f"{organ}_endosomal", plasma_flows[organ] * 1000)
  
  
  # cellular take-up
  for organ in organs:
    system.add_flow("MMAE", f"{organ}_interstitial", f"{organ}_cellular", permeability_surface_area_coefficients[organ] * PS)
    system.add_flow("MMAE", f"{organ}_cellular", f"{organ}_interstitial", permeability_surface_area_coefficients[organ] * PS * unbound_proportions[organ])
  
  
  # liver clearance
  system.add_flow("MMAE", "liver_interstitial", None, liver_clearance_rate * CL)

  
  # dissociatoin and degradation
  def degradation_endosomal(x, z):
    DAR = z["DAR"]
    rate = drug["degradation_endosomal"] * x["adc"]
    return {"adc": -rate, "drug": DAR * rate}
  
  def degradation_cellular(x, z):
    DAR = z["DAR"]
    rate = drug["degradation_cellular"] * x["adc"]
    return {"adc": -rate, "drug": DAR * rate}
  
  def dissociation(x, z):
    DAR = z["DAR"]
    if callable(drug["dissociation"]):
      rate = drug["dissociation"](DAR) * x["adc"]
    else:
      rate = drug["dissociation"] * x["adc"]
    return {"adc": -rate, "drug": DAR * rate}
  
  for organ in organs:
    system.add_reaction(f"{organ}_endosomal", degradation_endosomal)
    system.add_reaction(f"{organ}_cellular", degradation_cellular)
  
  system.add_reaction("plasma", dissociation)
  
