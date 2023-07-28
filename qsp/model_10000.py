from .qsp import *

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
volumes = np.array([0.00585, 0.00479, 0.0217, 0.000760, 0.0217, 0.119,
                    0.0295, 0.0241, 0.0384, 0.00102, 0.0384, 0.111,
                    0.249, 0.204, 1.47, 0.0566, 1.47, 9.34,
                    0.188, 0.154, 1.66, 0.0251, 1.66, 3.00,
                    0.0218, 0.0178, 0.337, 0.00991, 0.337, 1.60,
                    0.0621, 0.0508, 0.525, 0.0141, 0.525, 2.17,
                    0.0107, 0.00873, 0.0873, 0.00243, 0.0873, 0.376,
                    0.0289, 0.0236, 0.0788, 0.00263, 0.0788, 0.391,
                    0.164, 0.134, 0.385, 0.00963, 0.385, 1.23,
                    0.0116, 0.00950, 0.127, 0.00364, 0.127, 0.577,
                    0.0050, 0.00409, 0.0545, 0.00157, 0.0545, 0.248,
                    0.00534, 0.00437, 0.0169, 0.000485, 0.0169, 0.0699,
                    0.0005, 0.000405, 0.00153, 0.00005, 0.00153, 0.00653, 
                    0.0154, 0.0126, 0.0254, 0.000635, 0.0254, 0.0730,
                    0.0195, 0.0160, 0.0797, 0.00233, 0.0797, 0.348,
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
cellular_degradation_rate = 0.353 / units.h
dissociation_rate = 0.323 / units.d

permeability_surface_area_coefficient_BC = 0.105 * units.ml/units.h
permeability_surface_area_coefficients = {"heart": 1.47 * units.ml/units.h,
                                          "lung": 2.47 * units.ml/units.h,
                                          "muscle": 3.16 * units.ml/units.h,
                                          "skin": 0.681 * units.ml/units.h,
                                          "adipose": 0.588 * units.ml/units.h,
                                          "bone": 0.568 * units.ml/units.h,
                                          "brain": 0.00825 * units.ml/units.h,
                                          "kidney": 14.2 * units.ml/units.h,
                                          "liver": 49.2 * units.ml/units.h,
                                          "SI": 0.457 * units.ml/units.h,
                                          "LI": 0.457 * units.ml/units.h,
                                          "pancreas": 0.0657 * units.ml/units.h,
                                          "thymus": 0.457 * units.ml/units.h,
                                          "spleen": 0.457 * units.ml/units.h,
                                          "other": 0.457 * units.ml/units.h}

unbound_proportion_plasma = 0.8
unbound_proportion_BC = 0.8 / 5.46
unbound_proportions = {"heart": unbound_proportion_plasma / 22.8,
                       "lung": unbound_proportion_plasma / 64.9,
                       "muscle": unbound_proportion_plasma / 1.51,
                       "skin": unbound_proportion_plasma / 2.87,
                       "adipose": unbound_proportion_plasma / 3.01,
                       "bone": unbound_proportion_plasma / 1.89,
                       "brain": unbound_proportion_plasma / 0.530,
                       "kidney": unbound_proportion_plasma / 42.4,
                       "liver": unbound_proportion_plasma / 3.80,
                       "SI": unbound_proportion_plasma / 27.1 * (0.577 / 0.728),
                       "LI": unbound_proportion_plasma / 27.1 * (0.248 / 0.314),
                       "pancreas": unbound_proportion_plasma / 2.93,
                       "thymus": unbound_proportion_plasma / 27.1 * (0.00653 / 0.009),
                       "spleen": unbound_proportion_plasma / 47.2,
                       "other": unbound_proportion_plasma / 27.1 * (0.348/0.465)}
liver_clearance_rate = 137 * units.ml/units.h


compartments = [f"{organ}_{tissue}" for organ in organs for tissue in ["plasma", "BC", "interstitial", "endosomal", "membrane", "cellular"]] + ["plasma", "BC", "lymph"]
analytes = ["T-vc-MMAE", "MMAE"]
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


# endosomal take-up
for organ in organs:
  system.add_flow("T-vc-MMAE", f"{organ}_plasma", f"{organ}_endosomal", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal"))
  system.add_flow("T-vc-MMAE", f"{organ}_endosomal", f"{organ}_plasma", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal") * 0.715)
  
  system.add_flow("T-vc-MMAE", f"{organ}_interstitial", f"{organ}_endosomal", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal"))
  system.add_flow("T-vc-MMAE", f"{organ}_endosomal", f"{organ}_interstitial", endosomal_pinocytosis_rate * system.get_volume("T-vc-MMAE", f"{organ}_endosomal") * (1-0.715))


# target-specific membrane crossing
K_on_HER2 = 0.03 / units.nM / units.h
K_off_HER2 = 0.014 / units.h
K_int = 0.11 / units.h
N_HER2 = 1e4
cell_density = 1e9 / units.ml

for organ in organs:
  system.add_flow("T-vc-MMAE", f"{organ}_interstitial", f"{organ}_membrane", K_on_HER2 * N_HER2 * cell_density / units.avagadro * system.get_volume("T-vc-MMAE", f"{organ}_membrane"))
  system.add_flow("T-vc-MMAE", f"{organ}_membrane", f"{organ}_interstitial", K_off_HER2 * system.get_volume("T-vc-MMAE", f"{organ}_membrane"))
  system.add_flow("T-vc-MMAE", f"{organ}_membrane", f"{organ}_cellular", K_int * system.get_volume("T-vc-MMAE", f"{organ}_membrane"))


# dissociatoin and degradation
def degradation_endosomal(x, t):
  rate = endosomal_degradation_rate * x["T-vc-MMAE"]
  DAR = 1.5 * math.exp(-0.15/units.h * t) + 3 * math.exp(-0.012/units.h * t)
  return {"T-vc-MMAE": -rate, "MMAE": DAR * rate}

def degradation_cellular(x, t):
  rate = cellular_degradation_rate * x["T-vc-MMAE"]
  DAR = 1.5 * math.exp(-0.15/units.h * t) + 3 * math.exp(-0.012/units.h * t)
  return {"T-vc-MMAE": -rate, "MMAE": DAR * rate}

def dissociation(x, t):
  rate = dissociation_rate * x["T-vc-MMAE"]
  DAR = 1.5 * math.exp(-0.15/units.h * t) + 3 * math.exp(-0.012/units.h * t)
  return {"T-vc-MMAE": -rate, "MMAE": DAR * rate}

for organ in organs:
  system.add_reaction(f"{organ}_endosomal", degradation_endosomal)
  system.add_reaction(f"{organ}_cellular", degradation_cellular)

system.add_reaction("plasma", dissociation)

# plasma circle
system.add_flow("MMAE", "plasma", "lung_plasma", 373 * units.ml / units.h)

for organ in [organ for organ in organs if organ not in ["lung"]]:
  system.add_flow("MMAE", "lung_plasma", f"{organ}_plasma", plasma_flows[organ])

for organ in ["SI", "LI", "spleen", "pancreas"]:
  system.add_flow("MMAE", f"{organ}_plasma", "liver_plasma", plasma_flows[organ])
  system.add_flow("MMAE", "liver_plasma", "plasma", plasma_flows[organ])

for organ in [organ for organ in organs if organ not in ["lung"] + ["SI", "LI", "spleen", "pancreas"]]:
  system.add_flow("MMAE", f"{organ}_plasma", "plasma", plasma_flows[organ])


# blood cell up-take
system.add_flow("MMAE", "plasma", "BC", permeability_surface_area_coefficient_BC * unbound_proportion_plasma)
system.add_flow("MMAE", "BC", "plasma", permeability_surface_area_coefficient_BC * unbound_proportion_BC)

for organ in organs:
  system.add_flow("MMAE", f"{organ}_plasma", f"{organ}_BC", permeability_surface_area_coefficient_BC * unbound_proportion_plasma)
  system.add_flow("MMAE", f"{organ}_BC", f"{organ}_plasma", permeability_surface_area_coefficient_BC * unbound_proportion_BC)


# organ plasma to interstitial flow (here I do not follow the authors in multiplying plasma flow by 1000)
for organ in organs:
  system.add_flow("MMAE", f"{organ}_plasma", f"{organ}_endosomal", plasma_flows[organ] * unbound_proportion_plasma * 1000)
  system.add_flow("MMAE", f"{organ}_endosomal", f"{organ}_plasma", plasma_flows[organ] * 1000)
  
  system.add_flow("MMAE", f"{organ}_endosomal", f"{organ}_interstitial", plasma_flows[organ] * 1000)
  system.add_flow("MMAE", f"{organ}_interstitial", f"{organ}_endosomal", plasma_flows[organ] * 1000)


# cellular take-up
for organ in organs:
  system.add_flow("MMAE", f"{organ}_interstitial", f"{organ}_cellular", permeability_surface_area_coefficients[organ])
  system.add_flow("MMAE", f"{organ}_cellular", f"{organ}_interstitial", permeability_surface_area_coefficients[organ] * unbound_proportions[organ])


# liver clearance
system.add_flow("MMAE", "liver_interstitial", None, liver_clearance_rate)
