from .qsp import *

### this model is mostly from ...

cells = ["C", "Teff", "Treg"]
drugs = ["bimasked", "leftmasked", "rightmasked", "unmasked"]
C_conjugates = [f"{Ag}-{drug}" for drug in drugs for Ag in ["A", "B", "AA", "AB", "BB"]]
T_conjugates = [f"{drug}-{CD3}" for drug in drugs for CD3 in ["CD3eff", "CD3reg"]]
trimers = [f"{Ag}-{drug}-{CD3}" for drug in drugs for Ag in ["A", "B", "AA", "AB", "BB"] for CD3 in ["CD3eff", "CD3reg"]]
others =  [f"FcRn-{drug}" for drug in drugs]
analytes = cells + drugs + C_conjugates + T_conjugates + trimers + others

compartments = ["central", "peripheral", "organ", "tumor_plasma", "tumor_endosomal", "tumor_interstitial"]


def nonlinear_clearance_mouse(system, t):
  x = system.get_x("antibody", "central")
  rate = - mouse["nonlinear_clearance_EMAX"] * x / (x + mouse["nonlinear_clearance_EC50"]) / system.get_volume("antibody", "central")
  system.add_x("antibody", "central", rate * t)

def nonlinear_clearance_human(system, t):
  x = system.get_x("antibody", "central")
  rate = - human["nonlinear_clearance_EMAX"] * x / (x + human["nonlinear_clearance_EC50"]) / system.get_volume("antibody", "central")
  system.add_x("antibody", "central", rate * t)

molecular_weight = 150000 * units.g/units.mol

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

tumor_cell_density = 3e8 / units.ml

HER2 = {}
HER2.update({"num": 8e5})
HER2.update({"affinity": 1e-8 * units.molar})
HER2.update({"off": 0.106 / units.h, "internalization": 0.0 / units.h})
HER2.update({"on": HER2["off"] / HER2["affinity"]})

MC38 = {}
MC38.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_endosomal_proportion": 0.005})
MC38.update({"area": 1 * units.cm**2, "depth_layer": 0.01 * units.cm})
MC38.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
MC38.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
MC38.update({"diffusion": 10 * units.micrometer**2 / units.s})





def model(host, target, tumor):
  system = System(analytes, compartments)
  
  for analyte in analytes:
    system.set_volume(analyte, "central", host["volume_central"])
    system.set_volume(analyte, "peripheral", host["volume_peripheral"])
    system.set_volume(analyte, "tumor_plasma", tumor["volume"] * tumor["volume_plasma_proportion"])
    system.set_volume(analyte, "tumor_endosomal", tumor["volume"] * tumor["volume_endosomal_proportion"])
    for i in range(n_layers):
      system.set_volume(analyte, f"tumor_interstitial_{i}", tumor["area"] * tumor["depth_layer"])
  
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
