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
X.update({"off_P": 10**-4 / units.s, "affn_P": 260 * units.nM, "affm_P": 26000 * units.nM})
X.update({"clearance": math.log(2)/(80 * units.h), "smalls": []})
X["cleavages"] = [(linker, ["m"], ["n"])]


############ model ############

def model(TCE, plasma, lymph, organs):
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + organs]
  system = System(compartments, analytes, cells)
  system.centrals = [plasma, lymph]
  system.organs = organs

  # define volumes
  for central in centrals:
    system.set_volume(central["name"], central["volume"])
  for organ in organs:
    system.set_volume(organ["name"], organ["volume_interstitial"])
  
  # distribution
  for drug in drugs:
    for organ in organs:
      system.add_flow(drug, "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
      system.add_flow(drug, organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
      system.add_flow(drug, "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
  
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


############# demo ###############
'''
from qsp import *
from qsp.processes import *
from qsp.human import *
from qsp.model_PD1 import *


system = model(BD, plasma, lymph, [bone, lung, liver])
'''
