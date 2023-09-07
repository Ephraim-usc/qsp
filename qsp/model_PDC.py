from .qsp import *
from scipy.optimize import fsolve

analytes = ["bimasked", "monomasked", "unmasked", "payload", "target", "target-bimasked", "target-monomasked", "target-unmasked", "FcRn", "FcRn-bimasked", "FcRn-monomasked", "FcRn-unmasked"]

organs = ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"]
tissues = ["plasma", "BC", "interstitial", "endosomal", "cellular"]
centrals = ["plasma", "BC", "lymph"]
compartments = centrals + [f"{organ}_{tissue}" for organ in organs for tissue in tissues]

################### hosts ###################
volumes_mouse = np.array([0.944, 0.773, 0.113,
                          0.00585, 0.00479, 0.0217, 0.000760, 0.119,
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
                         ]) * units.ml

volumes_human = np.array([3126, 2558, 274,
                          13.1, 10.8, 48.8, 1.71, 267,
                          55.0, 45.0, 300, 5.00, 595,
                          662, 541, 3910, 150, 24815,
                          127, 104, 1125, 17.0, 2035,
                          148, 121, 2289, 67.3, 10840,
                          224, 183, 1891, 50.8, 7817,
                          31.9, 26.1, 261, 7.25, 1124,
                          18.2, 14.9, 49.8, 1.66, 247,
                          183, 149, 429, 10.7, 1371,
                          6.15, 5.03, 67.1, 1.93, 305,
                          8.74, 7.15, 95.3, 2.74, 434,
                          5.70, 4.66, 18.0, 0.518, 74.7,
                          0.353, 0.288, 1.09, 0.0321, 4.65,
                          26.8, 21.9, 44.3, 1.11, 127,
                          204, 167, 831, 24.3, 3626
                         ]) * units.ml

plasma_flows_mouse = np.array([36.5, 373, 86.1, 27.8, 13.4, 15.2, 11.8, 68.5, 10.3, 58.1, 17.3, 6.24, 1.19, 8.18, 10.9]) * units.ml/units.h
plasma_flows_human = np.array([7752, 181913, 33469, 11626, 11233, 2591, 21453, 36402, 13210, 12368, 12867, 3056, 353, 6343, 5521]) * units.ml/units.h

lymphatic_flows_mouse = plasma_flows_mouse * 1/500
lymphatic_flows_human = plasma_flows_human * 1/500

BC_flows_mouse = np.array([29.9, 305, 70.5, 22.8, 11.0, 12.4, 9.64, 56.1, 8.40, 47.5, 14.1, 5.10, 0.97, 6.70, 8.91]) * units.ml/units.h
BC_flows_human = np.array([6342, 148838, 27383, 9512, 9191, 2120, 17553, 29784, 10808, 10120, 10527, 2500, 289, 5189, 4517]) * units.ml/units.h

cell_densities = np.array([1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9]) * 1/units.ml

vascular_reflections = np.array([0.95, 0.95, 0.95, 0.95, 0.95, 0.85, 0.99, 0.9, 0.85, 0.9, 0.95, 0.9, 0.9, 0.85, 0.95])
lymphatic_reflection = 0.2

endosomal_pinocytosis = 3.66e-2 / units.h
vascular_recycle = 0.715
FcRn = 49.8 * units.micromolar
FcRn_on = 0.0806 * 1/units.nM/units.d
FcRn_off = 6.55 / units.h

mouse = {}
mouse.update({f"volume_{compartment}":x for compartment, x in zip(compartments, volumes_mouse)})
mouse.update({f"plasma_flow_{organ}":x for organ, x in zip(organs, plasma_flows_mouse)})
mouse.update({f"BC_flow_{organ}":x for organ, x in zip(organs, BC_flows_mouse)})
mouse.update({f"lymphatic_flow_{organ}":x for organ, x in zip(organs, lymphatic_flows_mouse)})
mouse.update({f"cell_density_{organ}":x for organ, x in zip(organs, cell_densities)})
mouse.update({f"vascular_reflection_{organ}":x for organ, x in zip(organs, vascular_reflections)})
mouse.update({"lymphatic_reflection":lymphatic_reflection})
mouse.update({"endosomal_pinocytosis":endosomal_pinocytosis, "vascular_recycle":vascular_recycle})
mouse.update({"FcRn": FcRn, "FcRn_on": FcRn_on, "FcRn_off": FcRn_off})

human = {}
human.update({f"volume_{compartment}":x for compartment, x in zip(compartments, volumes_human)})
human.update({f"plasma_flow_{organ}":x for organ, x in zip(organs, plasma_flows_human)})
human.update({f"BC_flow_{organ}":x for organ, x in zip(organs, BC_flows_human)})
human.update({f"lymphatic_flow_{organ}":x for organ, x in zip(organs, lymphatic_flows_human)})
human.update({f"cell_density_{organ}":x for organ, x in zip(organs, cell_densities)})
human.update({f"vascular_reflection_{organ}":x for organ, x in zip(organs, vascular_reflections)})
human.update({"lymphatic_reflection":lymphatic_reflection})
human.update({"endosomal_pinocytosis":endosomal_pinocytosis, "vascular_recycle":vascular_recycle})
human.update({"FcRn": FcRn, "FcRn_on": FcRn_on, "FcRn_off": FcRn_off})

################### targets ###################
# organs = ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"]
nums_CD71 = np.array([1e3, 1e5, 1e3, 1e4, 0, 1e5, 1e4, 1e4, 1e3, 1e4, 1e4, 1e3, 0, 0, 0])
nums_CD166 = np.array([0, 1e5, 1e3, 1e5, 0, 0, 1e5, 1e5, 1e5, 1e5, 1e3, 1e5, 1e5, 0, 0])
nums_CD166 = np.array([1e3, 1e4, 1e5, 1e3, 0, 0, 1e3, 1e4, 1e4, 1e4, 1e4, 1e3, 0, 0, 0])

nums_zero = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
nums_HER2 = np.array([1e4, 1e4, 1e4, 1e4, 0, 1e3, 0, 1e3, 1e3, 1e3, 1e3, 0, 1e3, 0, 0])
nums_CAIX = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1e5, 1e5, 0, 0, 0, 0, 0])
nums_EGFR = np.array([1e3, 1e5, 1e3, 1e5, 1e3, 0, 1e3, 1e4, 1e4, 1e3, 1e4, 1e4, 1e4, 0, 0])
nums_mesothelin = np.array([0, 1e5, 0, 1e3, 0, 0, 0, 0, 0, 0, 1e4, 0, 0, 0, 0])


# from Singh et al. 2020
on_HER2 = 0.03 / units.nM / units.h
off_HER2 = 0.014 / units.h
int_HER2 = 0.11 / units.h

zero = {}
zero.update({f"num_{organ}":num for organ, num in zip(organs, nums_zero)})
zero.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

HER2 = {f"num_{organ}":num for organ, num in zip(organs, nums_HER2)}
HER2.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

CAIX = {f"num_{organ}":num for organ, num in zip(organs, nums_CAIX)}
CAIX.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

EGFR = {f"num_{organ}":num for organ, num in zip(organs, nums_EGFR)}
EGFR.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

mesothelin = {f"num_{organ}":num for organ, num in zip(organs, nums_mesothelin)}
mesothelin.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

CD71 = {f"num_{organ}":num for organ, num in zip(organs, nums_CD71)}
CD71.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

CD166 = {f"num_{organ}":num for organ, num in zip(organs, nums_CD166)}
CD166.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})

AXL = {f"num_{organ}":num for organ, num in zip(organs, nums_AXL)}
AXL.update({"on":on_HER2, "off":off_HER2, "int":int_HER2})


################### linkers ###################
# dissociation_vc = 0.323 / units.d # from Adam P. Singh et al. 2020

def dissociation_vc(DAR):
  curve = lambda t: 1.5*math.exp(-0.15*t) + 3*math.exp(-0.012*t)
  t = fsolve(lambda t_: curve(t_) - DAR, 0)[0]
  delta = 1e-3
  derivative = (curve(t) - curve(t-delta))/delta
  return -derivative/DAR * 1/units.h

degradation_endosomal_vc = 42.9 / units.h
degradation_cellular_vc = 0.353 / units.h # from Adam P. Singh et al. 2017

vc = {}
vc.update({"dissociation":dissociation_vc, "degradation_endosomal":degradation_endosomal_vc, "degradation_cellular":degradation_cellular_vc})


################### payloads ###################
# organs = ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"]
PSs_MMAE = np.array([1.47, 2.47, 3.16, 0.681, 0.588, 0.568, 0.00825, 14.2, 49.2, 0.457, 0.457, 0.0657, 0.457, 0.457, 0.457]) * units.ml/units.h
PS_BC_MMAE = 0.105 * units.ml/units.h

permeabilities_MMAE = PSs_MMAE / np.array([mouse[f"volume_{organ}_cellular"] for organ in organs])
permeability_BC_MMAE = PS_BC_MMAE / mouse["volume_BC"]

unbound_plasma_MMAE = 0.8
unbound_BC_MMAE = 0.8 / 5.46
unbounds_MMAE = unbound_plasma_MMAE / np.array([22.8, 64.9, 1.51, 2.87, 3.01, 1.89, 0.530, 42.4, 3.80, 27.1*(0.728/0.577), 27.1*(0.314/0.248), 2.93, 27.1*(0.009/0.00653), 47.2, 27.1*(0.465/0.348)])

liver_clearance_MMAE = 137 * units.ml/units.h / mouse["volume_liver_interstitial"]

MMAE = {}
MMAE.update({f"permeability_{organ}":x for organ, x in zip(organs, permeabilities_MMAE)})
MMAE.update({"permeability_BC":permeability_BC_MMAE})
MMAE.update({f"unbound_{organ}":x for organ, x in zip(organs, unbounds_MMAE)})
MMAE.update({"unbound_plasma":unbound_plasma_MMAE, "unbound_BC":unbound_BC_MMAE})
MMAE.update({"liver_clearance":liver_clearance_MMAE})


################### masks ###################
Tx_M1 = {"foldchange": 1/220, "cleavage_central": 0.0527 / units.d, "cleavage_tumor": 0.1783 / units.d}
Tx_M2 = {"foldchange": 1/57, "cleavage_central": 0.0527 / units.d, "cleavage_tumor": 0.1783 / units.d}



################### the model ###################

def model(host, target, linker, payload, mask):
  variables = ["DAR"]
  system = System(analytes, compartments, variables)
  
  for analyte in analytes:
    for compartment in compartments:
      system.set_volume(analyte, compartment, host[f"volume_{compartment}"])
  
  for analyte in ["bimasked", "monomasked", "unmasked"]:
    # plasma and BC circles
    system.add_flow(analyte, "plasma", "lung_plasma", host["plasma_flow_lung"])
    system.add_flow(analyte, "BC", "lung_BC", host["BC_flow_lung"])
    for organ in [organ for organ in organs if organ not in ["lung"]]:
      system.add_flow(analyte, "lung_plasma", f"{organ}_plasma", host[f"plasma_flow_{organ}"])
      system.add_flow(analyte, "lung_BC", f"{organ}_BC", host[f"BC_flow_{organ}"])
    for organ in ["SI", "LI", "spleen", "pancreas"]:
      system.add_flow(analyte, f"{organ}_plasma", "liver_plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
      system.add_flow(analyte, "liver_plasma", "plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
      system.add_flow(analyte, f"{organ}_BC", "liver_BC", host[f"BC_flow_{organ}"])
      system.add_flow(analyte, "liver_BC", "BC", host[f"BC_flow_{organ}"])
    for organ in [organ for organ in organs if organ not in ["lung"] + ["SI", "LI", "spleen", "pancreas"]]:
      system.add_flow(analyte, f"{organ}_plasma", "plasma", host[f"plasma_flow_{organ}"] - host[f"lymphatic_flow_{organ}"])
      system.add_flow(analyte, f"{organ}_BC", "BC", host[f"BC_flow_{organ}"])
    
    # lymphatic circle
    for organ in organs:
      system.add_flow(analyte, f"{organ}_plasma", f"{organ}_interstitial", host[f"lymphatic_flow_{organ}"] * (1 - host[f"vascular_reflection_{organ}"]))
      system.add_flow(analyte, f"{organ}_interstitial", "lymph", host[f"lymphatic_flow_{organ}"] * (1 - host["lymphatic_reflection"]))
      system.add_flow(analyte, "lymph", "plasma", host[f"lymphatic_flow_{organ}"]) # the article used a different measure, which I think is weird
    
    # endosomal take-up and plasma recycle
    for organ in organs:
      v = system.get_volume(analyte, f"{organ}_endosomal")
      system.add_flow(analyte, f"{organ}_plasma", f"{organ}_endosomal", host["endosomal_pinocytosis"] * v)
      system.add_flow(analyte, f"{organ}_interstitial", f"{organ}_endosomal", host["endosomal_pinocytosis"] * v)
      
      system.add_reaction(f"{organ}_endosomal", {analyte:1, "FcRn":1}, {f"FcRn-{analyte}":1}, host["FcRn_on"], backward = host["FcRn_off"])
      system.add_reaction(f"{organ}_endosomal", {f"FcRn-{analyte}":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * host["vascular_recycle"], side_compartment = f"{organ}_plasma", side_products = {analyte:1})
      system.add_reaction(f"{organ}_endosomal", {f"FcRn-{analyte}":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * (1 - host["vascular_recycle"]), side_compartment = f"{organ}_interstitial", side_products = {analyte:1})

  
  # target binding and internalization
  for organ in organs:
    system.add_reaction(f"{organ}_interstitial", {"bimasked":1, "target":1}, {f"target-{analyte}":1}, target["on"] * mask["foldchange"], target["off"])
    system.add_reaction(f"{organ}_interstitial", {"bimasked":1, "target":1}, {f"target-{analyte}":1}, target["on"] * (0.5 + 0.5 * mask["foldchange"]), target["off"])
    system.add_reaction(f"{organ}_interstitial", {"bimasked":1, "target":1}, {f"target-{analyte}":1}, target["on"], target["off"])
    for analyte in ["bimasked", "monomasked", "unmasked"]:
      system.add_reaction(f"{organ}_interstitial", {f"target-{analyte}":1}, {"target":1}, target["int"], side_compartment = f"{organ}_cellular", side_products = {analyte:1})
  
  # cleavage of mask
  system.add_reaction("plasma", {"bimasked":1}, {"monomasked":1}, 2 * mask["cleavage_central"])
  system.add_reaction("plasma", {"monomasked":1}, {"unmasked":1}, mask["cleavage_central"])
  
  # dissociatoin and degradation
  def degradation(system, t):
    DAR = system.get_z("DAR")
    for analyte in ["bimasked", "monomasked", "unmasked"]:
      for organ in organs:
        for tissue in ["endosomal", "cellular"]:
          rate = linker[f"degradation_{tissue}"] * system.get_x(analyte, f"{organ}_{tissue}")
          system.add_x(analyte, f"{organ}_{tissue}", - rate * t)
          system.add_x("payload", f"{organ}_{tissue}", DAR * rate * t)
  
  def dissociation(system, t):
    DAR = system.get_z("DAR")
    if callable(linker["dissociation"]):
      rate = linker["dissociation"](DAR)
    else:
      rate = linker["dissociation"]
    for analyte in ["bimasked", "monomasked", "unmasked"]:
      for compartment in compartments:
        system.add_x(analyte, compartment, - rate * system.get_x(analyte, compartment) * t)
        system.add_x("payload", compartment, DAR * rate * system.get_x(analyte, compartment) * t)
    system.add_z("DAR", - DAR * rate * t)
  
  system.add_process(degradation)
  system.add_process(dissociation)

  
  # payload plasma and BC circles, and fast equilibrium between tissues
  system.add_flow("payload", "plasma", "lung_plasma", host["plasma_flow_lung"])
  system.add_flow("payload", "BC", "lung_BC", host["BC_flow_lung"])
  for organ in [organ for organ in organs if organ not in ["lung"]]:
    system.add_flow("payload", "lung_plasma", f"{organ}_plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("payload", "lung_BC", f"{organ}_BC", host[f"BC_flow_{organ}"])
  for organ in ["SI", "LI", "spleen", "pancreas"]:
    system.add_flow("payload", f"{organ}_plasma", "liver_plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("payload", "liver_plasma", "plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("payload", f"{organ}_BC", "liver_BC", host[f"BC_flow_{organ}"])
    system.add_flow("payload", "liver_BC", "BC", host[f"BC_flow_{organ}"])
  for organ in [organ for organ in organs if organ not in ["lung"] + ["SI", "LI", "spleen", "pancreas"]]:
    system.add_flow("payload", f"{organ}_plasma", "plasma", host[f"plasma_flow_{organ}"])
    system.add_flow("payload", f"{organ}_BC", "BC", host[f"BC_flow_{organ}"])
  for organ in organs:
    system.add_flow("payload", f"{organ}_plasma", f"{organ}_endosomal", host[f"plasma_flow_{organ}"] * payload["unbound_plasma"] * 1000)
    system.add_flow("payload", f"{organ}_interstitial", f"{organ}_endosomal", host[f"plasma_flow_{organ}"] * 1000)
    system.add_flow("payload", f"{organ}_endosomal", f"{organ}_plasma", host[f"plasma_flow_{organ}"] * 1000)
    system.add_flow("payload", f"{organ}_endosomal", f"{organ}_interstitial", host[f"plasma_flow_{organ}"] * 1000)
  
  # BC and cellular permeation
  v = host["volume_BC"]
  system.add_flow("payload", "plasma", "BC", payload["permeability_BC"] * v * payload["unbound_plasma"])
  system.add_flow("payload", "BC", "plasma", payload["permeability_BC"] * v * payload["unbound_BC"])
  for organ in organs:
    v = host[f"volume_{organ}_BC"]
    system.add_flow("payload", f"{organ}_plasma", f"{organ}_BC", payload["permeability_BC"] * v * payload["unbound_plasma"])
    system.add_flow("payload", f"{organ}_BC", f"{organ}_plasma", payload["permeability_BC"] * v * payload["unbound_BC"])
  for organ in organs:
    v = host[f"volume_{organ}_cellular"]
    system.add_flow("payload", f"{organ}_interstitial", f"{organ}_cellular", payload[f"permeability_{organ}"] * v)
    system.add_flow("payload", f"{organ}_cellular", f"{organ}_interstitial", payload[f"permeability_{organ}"] * v * payload[f"unbound_{organ}"])
  
  # liver clearance
  v = host["volume_liver_interstitial"]
  system.add_flow("payload", "liver_interstitial", None, payload["liver_clearance"] * v)
  
  
  # initial concentrations
  for organ in organs:
    system.set_x("FcRn", f"{organ}_endosomal", host["FcRn"])
  for organ in organs:
    system.set_x("target", f"{organ}_interstitial", target[f"num_{organ}"] * host[f"cell_density_{organ}"] / units.avagadro)
  
  
  system.set_z("DAR", 4.5)
  return system
