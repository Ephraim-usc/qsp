from .qsp import *
import itertools

### this model is mostly from ...

cells = ["T"]
areas = [200]

antigens = ["[T]P"]
bindings = ["[T]P"]
drugs = ["m", "n"]
dimers = [f"{binding}-{drug}" for binding in bindings for drug in drugs]
analytes = antigens + drugs + dimers

X = {}
X.update({"off_P": 10**-4 / units.s, "affn_P": 260 * units.nM, "affm_P": 26000 * units.nM})
X.update({"clearance": math.log(2)/(80 * units.h), "smalls": []})

X["internalization"] = internalization(rates = [("[T]C", ["[T]C"], 0.1 / units.h),
                                                 ("[B]A", ["[B]A"], 0.1 / units.h)])

BD["cleavage"] = transform()
for a in ("m", "n"):
  BD["cleavage"].add(linker = linker, reactant = f"m{a}", products = [f"n{a}"])
for c in ("m", "n"):
  BD["cleavage"].add(linker = linker, reactant = f"{c}m", products = [f"{c}n"])


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
  
  # 
  for compartment in compartments:
    for drug in drugs:
      system.add_flow(drug, compartment, None, system.get_volume(compartment) * TCE["clearance"])
  
  # small forms plasma clearance
  for small in TCE["smalls"]:
    system.add_flow(small, "plasma", None, system.get_volume("plasma") * math.log(2)/(45 * units.MIN))
  


############# demo ###############
'''
from qsp import *
from qsp.processes import *
from qsp.human import *
from qsp.model_PD1 import *


system = model(BD, plasma, lymph, [bone, lung, liver])
'''
