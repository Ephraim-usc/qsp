from .qsp import *
import re
import itertools

### this model is mostly from ...

cells = ["T", "B"]
areas = [200, 254]

antigens = ["[T]CD3", "[B]CD19", "HA"] # H: hydroxyapatite, C: CD3, A: CD19
bindings = ["[T]CD3", "[B]CD19", "HA"]
drugs = [f"{cd3}{cd19}" for cd3 in ("m", "n") for cd19 in ("m", "n")]
dimers = [f"{binding}-{drug}" for binding in bindings for drug in drugs]
analytes = antigens + drugs + dimers


############ drugs ############

linker = [("plasma", 0.07 / units.d), 
          ("lymph", 0.07 / units.d), 
          ("bone", 0.07 / units.d), 
          ("liver", 0.07 / units.d), 
          ("lung", 0.07 / units.d), 
          ("SI", 0.07 / units.d), 
          ("gallbladder", 0.07 / units.d)]

BD = {}
BD.update({"off_CD3": 10**-4 / units.s, "aff_CD3": 260 * units.nM, "mask_CD3": 20})
BD.update({"off_CD19": 10**-4 / units.s, "aff_CD19": 1.49 * units.nM, "mask_CD19": 20})
BD.update({"off_HA": 10**-4 / units.s, "aff_HA": math.inf * units.nM})
BD.update({"clearance": math.log(2)/(80 * units.h)})
BD["smalls"] = []

BD["cleavages"] = []
for cd19 in ["m", "n"]:
  BD["cleavages"].append((linker, f"m{cd19}", [f"n{cd19}"]))
for cd3 in ["m", "n"]:
  BD["cleavages"].append((linker, f"{cd19}m", [f"{cd19}n"]))

BD["internalization"] = [("[T]CD3", ["[T]CD3"], 0.1 / units.h), ("[B]CD19", ["[B]CD19"], 0.1 / units.h)]


############ model ############

def model(TCE, plasma, lymph, organs):
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + organs]
  system = System(compartments, analytes, cells)
  system.centrals = [plasma, lymph]
  system.organs = organs
  
  for central in centrals:
    system.set_volume(central["name"], central["volume"])
  for organ in organs:
    system.set_volume(organ["name"], organ["volume_interstitial"])
  
  # whole-body clearance
  for compartment in compartments:
    for drug in drugs:
      system.add_flow(drug, compartment, None, system.get_volume(compartment) * TCE["clearance"])
  
  # small forms plasma clearance
  for small in TCE["smalls"]:
    system.add_flow(small, "plasma", None, system.get_volume("plasma") * math.log(2)/(45 * units.MIN))
  
  for drug in drugs:
    for organ in organs:
      system.add_flow(drug, "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
      system.add_flow(drug, organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
      system.add_flow(drug, "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))

  
  # killing
  for drug in drugs:
    off_C = TCE["off_C"]; on_C = {"n":TCE["off_C"] / TCE["affn_C"], "m":TCE["off_C"] / TCE["affm_C"]}[drug[0]]
    off_A = TCE["off_A"]; on_A = {"n":TCE["off_A"] / TCE["affn_A"], "m":TCE["off_A"] / TCE["affm_A"]}[drug[1]]
    off_H = TCE["off_H"]; on_H = TCE["off_H"] / TCE["aff_H"]
    
    for organ in centrals + organs:
      system.add_simple(organ["name"], ["[T]C", f"{drug}"], [f"[T]C-{drug}"], on_C, off_C)
      system.add_simple(organ["name"], ["[B]A", f"{drug}"], [f"[B]A-{drug}"], on_A, off_A)
      system.add_simple(organ["name"], ["H", f"{drug}"], [f"H-{drug}"], on_H, off_H)
      
      system.add_simple(organ["name"], ["[T]C", f"H-{drug}"], [f"[T]C-{drug}", "H"], on_C)
      system.add_simple(organ["name"], ["[B]A", f"H-{drug}"], [f"[B]A-{drug}", "H"], on_A)
  
  ligands_effector = np.array(system.analytes)[np.array(system.ligands[0])]
  ligands_target = np.array(system.analytes)[np.array(system.ligands[1])]
  on2ds = pd.DataFrame(0.0, index = ligands_effector, columns = ligands_target) # in unit of um**2/s
  for drug in drugs:
    on2d_C = {"n":TCE["on2dn_C"], "m":TCE["on2dm_C"]}[drug[0]].number(units.um**2 / units.s)
    on2d_A = {"n":TCE["on2dn_A"], "m":TCE["on2dm_C"]}[drug[1]].number(units.um**2 / units.s)
    on2ds.loc[f"[T]C-{drug}", f"[B]A"] = on2d_A
    on2ds.loc[f"[T]C", f"[B]A-{drug}"] = on2d_A
  
  system.add_process(kill(compartments, on2ds, synapse_efficiency = TCE["synapse_efficiency"]))
  
  
  # mask cleavage
  if TCE["cleavage"] is not None:
    system.add_process(TCE["cleavage"])
  
  # internalization
  if TCE["internalization"] is not None:
    system.add_process(TCE["internalization"])
  
  # initial concentrations
  for central in centrals:
    system.add_c(central["name"], "T", central["num_T"] / central["volume"], ["C"], [124000])
    system.add_c(central["name"], "B", central["num_B"] / central["volume"], ["A"], [20000])
  
  for organ in organs:
    system.add_c(organ["name"], "T", organ["num_T"] / organ["volume_interstitial"], ["C"], [124000])
    system.add_c(organ["name"], "B", organ["num_B"] / organ["volume_interstitial"], ["A"], [20000])

  system.add_x("bone", "H", 100 * units.nM)
  
  return system



############# plot #############

def plot(system, name):
  groups = [["[T]C"],
            ["[B]A"],
            ["H"],
            drugs,
            [f"{binding}-{drug}" for binding in ["[T]C"] for drug in drugs],
            [f"{binding}-{drug}" for binding in ["[B]A"] for drug in drugs],
            [f"{binding}-{drug}" for binding in ["H"] for drug in drugs]]
  labels = ["CD3", "CD19", "HA", "drug", "CD3-drug", "CD19-drug", "HA-drug"]
  colors = [ "tab:orange", "tab:green", "tab:blue", "black", "wheat", "lightgreen", "skyblue"]
  linestyles = ["solid", "solid", "solid", "solid", "solid", "solid", "solid"]
  system.plot(compartments = system.compartments, 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              linthresh = 1e-5,
              output = f"{name}_summary.png")


############# demo ###############
'''
from qsp import *
from qsp.human import *
from qsp.model_TCE_B import *

bone.update({"plasma_flow": 10000 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})

system = model(BD, plasma, lymph, [bone, lung, liver])
for _ in range(2):
  system.add_x("plasma", "nn", 0.009 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "unmasked")
system.plot_cell(output = "unmasked_9ug.png") # 0.009nM

system = model(BD, plasma, lymph, [bone, lung, liver])
for _ in range(2):
  system.add_x("plasma", "nn", 5*70/1000 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "unmasked_5upk")
system.plot_cell(output = "unmasked_5upk.png") # 0.35nM

system = model(BD, plasma, lymph, [bone, lung, liver])
for _ in range(3):
  system.add_x("plasma", "nn", 135*70/1000 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "unmasked_135upk")
system.plot_cell(output = "unmasked_135upk.png") # 9.45nM




TCE = BD.copy()
TCE.update({"off_H": 10**-4 / units.s, "aff_H": 1 * units.nM})
system = model(TCE, plasma, lymph, [bone, lung, liver])
for _ in range(7):
  system.add_x("plasma", "nn", 0.009 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "unmasked_HA")
system.plot_cell(output = "unmasked_HA.png")


system = model(BD, plasma, lymph, [bone, lung, liver])
for _ in range(7):
  system.add_x("plasma", "mm", 1 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "masked")
system.plot_cell(output = "masked_cell.png")


TCE = BD.copy()
TCE.update({"off_H": 10**-4 / units.s, "aff_H": 1 * units.nM})
system = model(TCE, plasma, lymph, [bone, lung, liver])
for _ in range(7):
  system.add_x("plasma", "mm", 1 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "masked_HA")
system.plot_cell(output = "masked_HA.png")

'''
