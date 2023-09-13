from .qsp import *

### this model is mostly from ...

cells = ["C", "Teff", "Treg"]
drugs = ["bimasked", "leftmasked", "rightmasked", "unmasked"]
C_conjugates = [f"{Ag}-{drug}" for drug in drugs for Ag in ["A", "B", "AA", "AB", "BB"]]
T_conjugates = [f"{drug}-{CD3}" for drug in drugs for CD3 in ["CD3eff", "CD3reg"]]
trimers = [f"{Ag}-{drug}-{CD3}" for drug in drugs for Ag in ["A", "B", "AA", "AB", "BB"] for CD3 in ["CD3eff", "CD3reg"]]
others =  [f"FcRn-{drug}" for drug in drugs]
analytes = cells + drugs + C_conjugates + T_conjugates + trimers + others

compartments = ["central", "peripheral", "organ_plasma", "organ_endosomal", "organ_interstitial", "tumor_plasma", "tumor_interstitial"]


############ constants ############

molecular_weight = 150000 * units.g/units.mol
tumor_cell_density = 3e8 / units.ml

def Hill(EMAX, EC50, coefficient, x):
  if coefficient == 1:
    return EMAX/(1 + (x/EC50))
  else:
    return EMAX/(1 + (x/EC50) ** coefficient)


############ host ############

def nonlinear_clearance_mouse(system, t):
  x = system.get_x("antibody", "central")
  rate = - Hill(mouse["nonlinear_clearance_EMAX"], mouse["nonlinear_clearance_EC50"], 1, x)
  system.add_x("antibody", "central", rate * t)

def nonlinear_clearance_human(system, t):
  x = system.get_x("antibody", "central")
  rate = - Hill(human["nonlinear_clearance_EMAX"], human["nonlinear_clearance_EC50"], 1, x)
  system.add_x("antibody", "central", rate * t)

mouse = {}
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"distribution": 4.82 * units.ml/units.d})
mouse.update({"clearance": 0.334 * units.ml/units.d})
mouse.update({"nonlinear_clearance_EMAX": (0.518 * units.microgram/units.d / molecular_weight).cast_unit(units.nM*units.ml/ units.h)})
mouse.update({"nonlinear_clearance_EC50": (0.366 * units.microgram/units.ml / molecular_weight).cast_unit(units.nM)})
mouse.update({"nonlinear_clearance": nonlinear_clearance_mouse})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
mouse.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
mouse.update({"FcRn": 49.8 * units.micromolar, "FcRn-on": 0.0806 * 1/units.nM/units.h, "FcRn-off": 6.55 / units.h})

human = {}
human.update({"volume_central": 2877 * units.ml, "volume_peripheral": 2857 * units.ml})
human.update({"distribution": 384 * units.ml/units.d})
human.update({"clearance": 167 * units.ml/units.d})
human.update({"nonlinear_clearance_EMAX": (114 * units.microgram/units.d / molecular_weight).cast_unit(units.nM*units.ml/ units.h)})
human.update({"nonlinear_clearance_EC50": (0.078 * units.microgram/units.ml / molecular_weight).cast_unit(units.nM)})
human.update({"nonlinear_clearance": nonlinear_clearance_human})
human.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
human.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
human.update({"FcRn": 49.8 * units.micromolar, "FcRn-on": 0.792 * 1/units.nM/units.h, "FcRn-off": 23.9 / units.h})


############ antigens ############

HER2 = {}
HER2.update({"num": 8e5})
HER2.update({"affinity": 1e-8 * units.molar})
HER2.update({"off": 0.106 / units.h, "internalization": 0.0 / units.h})
HER2.update({"on": HER2["off"] / HER2["affinity"]})


############ tumors ############

MC38 = {}
MC38.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
MC38.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
MC38.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
MC38.update({"diffusion": 10 * units.micrometer**2 / units.s})


############ organs ############

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_endosomal": 1.93 * units.ml, "volume_interstitial_proportion": 67.1 * units.ml})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})




############ model ############

def model(host, target, tumor, organs):
  system = System(analytes, compartments)
  
  for analyte in analytes:
    system.set_volume(analyte, "central", host["volume_central"])
    system.set_volume(analyte, "peripheral", host["volume_peripheral"])
    system.set_volume(analyte, "tumor_plasma", tumor["volume"] * tumor["volume_plasma_proportion"])
    system.set_volume(analyte, "tumor_interstitial", tumor["volume"] * tumor["volume_interstitial_proportion"])
    for organ in organs:
      system.set_volume(analyte, f"{organ["name"]}_plasma", organ["volume"] * organ["volume_plasma_proportion"])
      system.set_volume(analyte, f"{organ["name"]}_endosomal", organ["volume"] * organ["volume_endosomal_proportion"])
      system.set_volume(analyte, f"{organ["name"]}_interstitial", organ["volume"] * organ["volume_interstitial_proportion"])
  
  # distribution and clearance
  system.add_flow("antibody", "central", "peripheral", host["distribution"])
  system.add_flow("antibody", "peripheral", "central", host["distribution"])
  system.add_flow("antibody", "central", None, host["clearance"])
  system.add_process(host["nonlinear_clearance"])
  

  
  # tumor plasma and lymphatic flow
  system.add_flow("antibody", "central", "tumor_plasma", tumor["volume"] * tumor["plasma_flow_density"])
  system.add_flow("antibody", "tumor_plasma", "central", tumor["volume"] * tumor["plasma_flow_density"] * (1 - tumor["lymphatic_flow_ratio"]))
  
  system.add_flow("antibody", "tumor_plasma", "tumor_interstitial_0", tumor["volume"] * tumor["plasma_flow_density"] * tumor["lymphatic_flow_ratio"] * (1 - host["vascular_reflection"]))
  system.add_flow("antibody", "tumor_interstitial_0", "central", tumor["volume"] * tumor["plasma_flow_density"] * tumor["lymphatic_flow_ratio"] * (1 - host["lymphatic_reflection"]))
  
  system.add_flow("antibody", "tumor_plasma", "tumor_interstitial_0", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
  system.add_flow("antibody", "tumor_interstitial_0", "tumor_plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
  
  # endosomal take-up and degradation
  system.add_flow("antibody", "tumor_plasma", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
  system.add_flow("antibody", "tumor_interstitial_0", "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_pinocytosis"])
  system.add_flow("antibody", "tumor_endosomal", None, tumor["volume"] * tumor["volume_endosomal_proportion"] * host["endosomal_degradation"])
  
  # recycle by FcRn
  system.set_x("FcRn", "tumor_endosomal", host["FcRn"])
  system.add_reaction("tumor_endosomal", {"antibody":1, "FcRn":1}, {"FcRn-antibody":1}, host["FcRn-on"], backward = host["FcRn-off"])
  system.add_reaction("tumor_endosomal", {"FcRn-antibody":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * host["vascular_recycle"], side_compartment = "tumor_plasma", side_products = {"antibody":1})
  system.add_reaction("tumor_endosomal", {"FcRn-antibody":1}, {"FcRn":1}, host["endosomal_pinocytosis"] * (1 - host["vascular_recycle"]), side_compartment = "tumor_interstitial_0", side_products = {"antibody":1})

  # diffusion in tumor interstitial
  for i in range(n_layers - 1):
    system.add_flow("antibody", f"tumor_interstitial_{i}", f"tumor_interstitial_{i+1}", tumor["diffusion"] * tumor["area"] / tumor["depth_layer"])
    system.add_flow("antibody", f"tumor_interstitial_{i+1}", f"tumor_interstitial_{i}", tumor["diffusion"] * tumor["area"] / tumor["depth_layer"])
  
  # interaction with target
  for i in range(n_layers):
    system.set_x("target", f"tumor_interstitial_{i}", target["num"] * tumor_cell_density / units.avagadro)
    system.add_reaction(f"tumor_interstitial_{i}", {"antibody":1, "target":1}, {"target-antibody":1}, target["off"]/target["affinity"], target["off"])
    system.add_reaction(f"tumor_interstitial_{i}", {"target-antibody":1}, {"target":1}, target["internalization"])
  
  return system
