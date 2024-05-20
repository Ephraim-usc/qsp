from .qsp import *
from .processes import *
import itertools

### this model is mostly from ...

cells = ["T"]
areas = [200]

antigens = ["[T]P"]
bindings = ["[T]P"]
drugs = ["m", "n"]
dimers = [f"{binding}-{drug}" for binding in bindings for drug in drugs]
analytes = antigens + drugs + dimers


### drugs ###

linker_175 = [("plasma", 0.07 / units.d), 
              ("lymph", 0.07 / units.d), 
              ("bone", 0.07 / units.d), 
              ("liver", 0.07 / units.d), 
              ("lung", 0.07 / units.d), 
              ("SI", 0.07 / units.d), 
              ("gallbladder", 0.07 / units.d)]

X = {}
X.update({"off_P": 10**-4 / units.s, "affn_P": 0.1 * units.nM, "affm_P": 10 * units.nM})
X.update({"clearance": math.log(2)/(80 * units.h), "smalls": []})
X["cleavages"] = [(linker_175, "m", ["n"])]
X["internalizations"] = [("[T]P", ["[T]P"], 0.1 / units.h)]

############ model ############

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
system.run(units.h)
system.get_y("tumor", "[T]P")
'''
