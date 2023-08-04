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
    system.add_flow("adc", "liver_plasma", "plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
    system.add_flow("drug", "liver_plasma", "plasma", host[f"plasma_flow_{organ}"])
  
  for organ in [organ for organ in organs if organ not in ["lung"] + ["SI", "LI", "spleen", "pancreas"]]:
    system.add_flow("adc", f"{organ}_plasma", "plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])

  #
  
