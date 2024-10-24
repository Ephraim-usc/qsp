from .qsp import *
from qsp.processes import *
from qsp.human import *
from qsp.tumors import *


def model(num_antigen, aff, off, int_rate, half_life,
          tumor_capillary_radius = 10 * units.um, tumor_capillary_permeability = 3e-7 * units.cm/units.s, tumor_diffusion = 10 * units.um**2 / units.s,
          tumor_capillary = 100*units.cm, tumor_layer_depth = 10*units.um, tumor_num_layers = 100):
  analytes = ["antigen", "drug", "antigen-drug"]
  centrals = [plasma, lymph]
  organs = [bone, lung, liver, SI, other]
  tumors = [FTC238.copy() for _ in range(10)]
  for i in range(10):
    tumors[i]["name"] = f"tumor_{i}"
  
  compartments = [organ["name"] for organ in centrals + organs + tumors]
  system = System(compartments, analytes)
  
  # define volumes
  for central in centrals:
    system.set_volume(central["name"], central["volume"])
  for organ in organs:
    system.set_volume(organ["name"], organ["volume_interstitial"])
  for tumor in tumors:
    system.set_volume(tumor["name"], tumor_surface_area * tumor_layer_depth)
  
  # organ distribution
  for organ in organs:
    system.add_flow("drug", "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
    system.add_flow("drug", organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
    system.add_flow("drug", "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
  
  # tumor distribution
  system.add_flow("drug", "plasma", "tumor_0", tumor_capillary * (2 * math.pi * tumor_capillary_radius) * tumor_capillary_permeability)
  system.add_flow("drug", "tumor_0", "plasma", tumor_capillary * (2 * math.pi * tumor_capillary_radius) * tumor_capillary_permeability)
  
  for i in range(tumor_num_layers - 1):
    radius = tumor_capillary_radius + i * tumor_layer_depth
    area = tumor_capillary * (2 * math.pi * radius)
    system.add_flow("drug", f"tumor_{i}", f"tumor_{i+1}", area/tumor_layer_depth * tumor_diffusion)
    system.add_flow("drug", f"tumor_{i+1}", f"tumor_{i}", area/tumor_layer_depth * tumor_diffusion)
  
  # bulk clearance
  for compartment in compartments:
    system.add_flow(drug, compartment, None, system.get_volume(compartment) * TCE["clearance"])
  
  # binding kinetics
  for compartment in compartments:
    system.add_simple(compartment, ["antigen", "drug"], ["antigen-drug"], aff*off, off)
    system.add_transform(compartment, "antigen-drug", ["antigen"], rate = int_rate)
  
  # initial concentrations
  for central in centrals:
    system.add_c(central["name"], "T", central["num_T"] / central["volume"], ["P"], [15000])
  for organ in organs:
    system.add_c(organ["name"], "T", organ["num_T"] / organ["volume_interstitial"], ["P"], [15000])
  for tumor in tumors:
    system.set_x("target", f"tumor_interstitial_{i}", target["num"] * tumor_cell_density / units.avagadro)
    system.add_c(tumor["name"], "T", tumor["density_T"] / tumor["volume_interstitial_proportion"], ["P"], [50000])





def model(TCE, plasma, lymph, organs, tumors):
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + organs + tumors]
  system = System(compartments, analytes, cells)
  system.centrals = [plasma, lymph]
  system.organs = organs
  system.tumors = tumors
  
  # define volumes
  for central in centrals:
    system.set_volume(central["name"], central["volume"])
  for organ in organs:
    system.set_volume(organ["name"], organ["volume_interstitial"])
  for tumor in tumors:
      system.set_volume(tumor["name"], tumor["volume"] * tumor["volume_interstitial_proportion"])
  
  # distribution
  for drug in drugs:
    for organ in organs:
      system.add_flow(drug, "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
      system.add_flow(drug, organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
      system.add_flow(drug, "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
    for tumor in tumors:
      system.add_flow(drug, "plasma", tumor["name"], tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
      system.add_flow(drug, tumor["name"], "plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
  
  # bulk clearance
  for compartment in compartments:
    for drug in drugs:
      system.add_flow(drug, compartment, None, system.get_volume(compartment) * TCE["clearance"])
  
  # small molecule clearance
  for small in TCE["smalls"]:
    system.add_flow(small, "plasma", None, system.get_volume("plasma") * math.log(2)/(45 * units.MIN))
  
  # mask cleavage
  for linker, drug_source, drug_dests in TCE["cleavages"]:
    add_cleavage(system, linker, drug_source, drug_dests, bindings)
  
  for binding_source, analyte_dests, rate in TCE["internalizations"]:
    add_internalization(system, binding_source, analyte_dests, rate, drugs)
  
  # binding kinetics
  for drug in drugs:
    off_P = TCE["off_P"]; on_P = {"n":TCE["off_P"] / TCE["affn_P"], "m":TCE["off_P"] / TCE["affm_P"]}[drug]
    for organ in centrals + organs + tumors:
      system.add_simple(organ["name"], ["[T]P", f"{drug}"], [f"[T]P-{drug}"], on_P, off_P)
  
  # initial concentrations
  for central in centrals:
    system.add_c(central["name"], "T", central["num_T"] / central["volume"], ["P"], [15000])
  for organ in organs:
    system.add_c(organ["name"], "T", organ["num_T"] / organ["volume_interstitial"], ["P"], [15000])
  for tumor in tumors:
    system.add_c(tumor["name"], "T", tumor["density_T"] / tumor["volume_interstitial_proportion"], ["P"], [50000])
  
  return system


############# demo ###############
'''
from qsp import *
from qsp.processes import *
from qsp.human import *
from qsp.tumors import *
from qsp.model_PD1 import *

system = model(X, plasma, lymph, [bone, lung, liver], [FTC238])
system.add_x("plasma", "n", 100 * units.nM)
system.run(168 * units.h, t_step = 1/6 * units.h)

coverages = np.array([x[-1, -1] / (x[-1, -1] + x[0, -1]) for t, x, c in system.history])
coverage = coverages.mean()

'''
