from .qsp import *
import re

### this model is mostly from ...

effectors = ["C"]; targets = ["A", "B", "AB"]

antigens_effector = ["C"]; antigens_target = ["A", "B"]
drugs = [f"{c}{a}" for c in ["m", "n"] for a in ["m", "n"]]
dimers_effector = [f"{effector}-{drug}" for effector in effectors for drug in drugs]
dimers_target = [f"{drug}-{target}" for drug in drugs for target in targets]
trimers = [f"{effector}-{drug}-{target}" for effector in effectors for drug in drugs for target in targets]
analytes = antigens_effector + antigens_target + drugs + dimers_effector + dimers_target + trimers



############ constants ############

molecular_weight = 150000 * units.g/units.mol

class equilibrium:
  def __init__(self, compartments, analytes):
    self.system = None
    self.compartments = compartments
    self.analytes = analytes
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments]
      self.analytes_ = [system.analytes.index(analyte) for analyte in self.analytes]
    
    for analyte_ in self.analytes_:
      x = system.x[analyte_, self.compartments_]
      volumes = system.V[analyte_, self.compartments_]
      system.x[analyte_, self.compartments_] = np.average(x, weights = volumes)



############ drug ############

class transform:
  def __init__(self, compartments, rates):
    self.system = None
    self.compartments = compartments
    
    Q = np.zeros([len(drugs), len(drugs)])
    for reactant, product, rate in rates:
      reactants = [i for i, drug in enumerate(drugs) if re.match(reactant + "$", drug)]
      products = [i for i, drug in enumerate(drugs) if re.match(product + "$", drug)]
      Q[reactants, reactants] -= rate.number(1/units.h)
      Q[reactants, products] += rate.number(1/units.h)
    self.Q = Q
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      if callable(self.compartments):
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments(system) if compartment in system.compartments]
      else:
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments if compartment in system.compartments]
      
      self.analyteses_ = []
      self.analyteses_.append([system.analytes.index(f"{drug}") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"{drug}-A") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"{drug}-B") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"{drug}-AB") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-A") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-B") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-AB") for drug in drugs])
    
    for compartment in self.compartments_:
      for analytes_ in self.analyteses_:
        system.x[analytes_, compartment] = system.x[analytes_, compartment] @ expm(self.Q * t.number(units.h))


class internalization:
  def __init__(self, compartments, rates_effector, rates_target):
    self.system = None
    self.compartments = compartments
    
    q_effector = np.zeros(len(dimers_effector))
    Q_effector = np.zeros([len(dimers_effector), len(antigens_effector)])
    for effector, antigens, rate in rates_effector:
      idx_effector = [dimers_effector.index(f"{effector}-{drug}") for drug in drugs]
      idx_antigens = [antigens_effector.index(antigen) for antigen in antigens]
      q_effector[idx_effector] -= rate.number(1/units.h)
      for i in idx_effector:
        np.add.at(Q_effector, (i, idx_antigens), 1)
    
    q_target = np.zeros(len(dimers_target))
    Q_target = np.zeros([len(dimers_target), len(antigens_target)])
    for target, antigens, rate in rates_target:
      idx_target = [dimers_target.index(f"{drug}-{target}") for drug in drugs]
      idx_antigens = [antigens_target.index(antigen) for antigen in antigens]
      q_target[idx_target] -= rate.number(1/units.h)
      for i in idx_target:
        np.add.at(Q_target, (i, idx_antigens), 1) # there may be repeated antigens
    
    self.q_effector = q_effector
    self.Q_effector = Q_effector
    self.q_target = q_target
    self.Q_target = Q_target
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      if callable(self.compartments):
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments(system) if compartment in system.compartments]
      else:
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments if compartment in system.compartments]
      
      self.idx_dimers_effector = [system.analytes.index(analyte) for analyte in dimers_effector]
      self.idx_dimers_target = [system.analytes.index(analyte) for analyte in dimers_target]
      self.idx_antigens_effector = [system.analytes.index(analyte) for analyte in antigens_effector]
      self.idx_antigens_target = [system.analytes.index(analyte) for analyte in antigens_target]
    
    for compartment in self.compartments_:
      delta_effector = system.x[self.idx_dimers_effector, compartment] * (1 - np.exp(self.q_effector * t.number(units.h)))
      system.x[self.idx_dimers_effector, compartment] -= delta_effector
      system.x[self.idx_antigens_effector, compartment] += delta_effector @ self.Q_effector
      
      delta_target = system.x[self.idx_dimers_target, compartment] * (1 - np.exp(self.q_target * t.number(units.h)))
      system.x[self.idx_dimers_target, compartment] -= delta_target
      system.x[self.idx_antigens_target, compartment] += delta_target @ self.Q_target


VIBX = {}
VIBX.update({"off_C": 10**-4 / units.s, "affn_C": 10 * units.nM, "affm_C": 200 * units.nM})
VIBX.update({"off_A": 10**-4 / units.s, "affn_A": 10 * units.nM, "affm_A": 200 * units.nM})
VIBX.update({"off_B": 10**-4 / units.s, "aff_B": 10 * units.nM})
VIBX.update({"avidity_effector": 19, "avidity_target": 19})
VIBX.update({"clearance": math.log(2)/(70 * units.h)}); VIBX["smalls"] = []
VIBX["cleavage_plasma"] = transform(compartments = lambda system: [central["name"] for central in system.centrals] + [organ["name"] for organ in system.organs], 
                                    rates = [("m.", "n.", 0.05 / units.d), (".m", ".n", 0.05 / units.d)])
VIBX["cleavage_tumor"] = transform(compartments = lambda system: [tumor["name"] for tumor in system.tumors], 
                                   rates = [("m.", "n.", 0.15 / units.d), (".m", ".n", 0.15 / units.d)])
VIBX["internalization"] = internalization(compartments = lambda system: system.compartments,
                                          rates_effector = [("C", ["C"], 0.1 / units.h)],
                                          rates_target = [("A", ["A"], 0.02 / units.h), ("B", ["B"], 0.02 / units.h), ("AB", ["A", "B"], 0.02 / units.h)])


VIBY = {}
VIBY.update({"off_C": 10**-4 / units.s, "affn_C": 10 * units.nM, "affm_C": 200 * units.nM})
VIBY.update({"off_A": 10**-4 / units.s, "affn_A": 10 * units.nM, "affm_A": 200 * units.nM})
VIBY.update({"off_B": 10**-4 / units.s, "aff_B": 10 * units.nM})
VIBY.update({"avidity_effector": 19, "avidity_target": 19})
VIBY.update({"clearance": math.log(2)/(70 * units.h)}); VIBY["smalls"] = ["mn", "nn"]
VIBY["cleavage_plasma"] = transform(compartments = lambda system: [central["name"] for central in system.centrals] + [organ["name"] for organ in system.organs], 
                                    rates = [("m.", "n.", 0.05 / units.d), (".m", ".n", 0.05 / units.d)])
VIBY["cleavage_tumor"] = transform(compartments = lambda system: [tumor["name"] for tumor in system.tumors], 
                                   rates = [("m.", "n.", 0.15 / units.d), (".m", ".n", 0.15 / units.d)])
VIBY["internalization"] = internalization(compartments = lambda system: system.compartments,
                                          rates_effector = [("C", ["C"], 0.1 / units.h)],
                                          rates_target = [("A", ["A"], 0.02 / units.h), ("B", ["B"], 0.02 / units.h), ("AB", ["A", "B"], 0.02 / units.h)])



############ tumors ############

UT44 = {"name": "tumor"}
UT44.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
UT44.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
UT44.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
UT44.update({"diffusion": 10 * units.micrometer**2 / units.s})
UT44.update({"density_cell": 3e8 * 0.44 / units.ml, "density_T": 3e8 * 0.15 / units.ml, "density_NK": 3e8 * 0.02 / units.ml})
UT44.update({"num_A": 7e5, "num_B": 1.45e6})

FTC238 = {"name": "tumor"}
FTC238.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
FTC238.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
FTC238.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
FTC238.update({"diffusion": 10 * units.micrometer**2 / units.s})
FTC238.update({"density_cell": 3e8 * 0.44 / units.ml, "density_T": 3e8 * 0.15 / units.ml, "density_NK": 3e8 * 0.02 / units.ml})
FTC238.update({"num_A": 1e5, "num_B": 1e5})


############ organs ############

plasma = {"name": "plasma"}
plasma.update({"volume": 3126 * units.ml})
plasma.update({"num_T": 7.9E+09, "num_NK": 1.6E+09})
plasma.update({"conc_A": 5 * units.ug/units.ml / (170 * units.kDa), "conc_B": 0 * units.nM})

lymph = {"name": "lymph"}
lymph.update({"volume": 274 * units.ml})
lymph.update({"num_T": 3.6E+11, "num_NK": 6.7E+08})
lymph.update({"conc_A": 0 * units.nM, "conc_B": 0 * units.nM})

bone = {"name": "bone"}
bone.update({"volume_plasma": 224 * units.ml, "volume_interstitial": 1891 * units.ml})
bone.update({"plasma_flow": 2591 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
bone.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
bone.update({"num_cell": 4.77E+09 * 0.5, "num_T": 2.1E+10, "num_NK": 3.3E+09})
bone.update({"num_A": 0, "num_B": 0})

lung = {"name": "lung"}
lung.update({"volume_plasma": 55 * units.ml, "volume_interstitial": 300 * units.ml})
lung.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
lung.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"num_cell": 2.36E+11 * 0.5, "num_T": 1.3E+10, "num_NK": 7.2E+08})
lung.update({"num_A": 133439, "num_B": 0}) # num_B = 1019 from Liyuan

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"num_cell": 7.2e11 * 0.5, "num_T": 1.8E+10, "num_NK": 8.1E+08})
SI.update({"num_A": 57075, "num_B": 39649})

other = {"name": "other"}
other.update({"volume_plasma": 1000 * units.ml, "volume_interstitial": 5000 * units.ml})
other.update({"plasma_flow": 100000 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
other.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
other.update({"num_cell": 1e13 * 0.5, "num_T": 1.1E+10, "num_NK": 3.1E+09})
other.update({"num_A": 10000, "num_B": 0})


############ model ############

def model(TCE, tumors, organs, connect_tumors = True):
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + tumors + organs]
  system = System(analytes, compartments)
  system.centrals = [plasma, lymph]
  system.tumors = tumors
  system.organs = organs
  
  for analyte in analytes:
    for central in centrals:
      system.set_volume(analyte, central["name"], central["volume"])
    for tumor in tumors:
      system.set_volume(analyte, tumor["name"], tumor["volume"] * tumor["volume_interstitial_proportion"])
    for organ in organs:
      system.set_volume(analyte, organ["name"], organ["volume_interstitial"])
  
  # whole-body clearance
  for compartment in compartments:
    for drug in drugs:
        system.add_flow(drug, compartment, None, system.get_volume(drug, compartment) * TCE["clearance"])
  
  # small forms plasma clearance
  for small in TCE["smalls"]:
    system.add_flow(small, "plasma", None, system.get_volume(drug, "plasma") * math.log(2)/(45 * units.MIN))
  
  for drug in drugs:
    # drug tumor flow
    for tumor in tumors:
      system.add_flow(drug, "plasma", tumor["name"], tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
      system.add_flow(drug, tumor["name"], "plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
    
    # drug organ flow
    for organ in organs:
      system.add_flow(drug, "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
      system.add_flow(drug, organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
      system.add_flow(drug, "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
  
  # exchange drugs between tumors if tumors are connected
  if connect_tumors:
    system.add_process(equilibrium([tumor["name"] for tumor in tumors], drugs))
  
  # mask cleavage
  system.add_process(TCE["cleavage_plasma"])
  system.add_process(TCE["cleavage_tumor"])
  
  # target binding
  for drug in drugs:
    off_C = TCE["off_C"]; on_C = {"n":TCE["off_C"] / TCE["affn_C"], "m":TCE["off_C"] / TCE["affm_C"]}[drug[0]]
    off_A = TCE["off_A"]; on_A = {"n":TCE["off_A"] / TCE["affn_A"], "m":TCE["off_A"] / TCE["affm_A"]}[drug[1]]
    off_B = TCE["off_B"]; on_B = TCE["off_B"] / TCE["aff_B"]
    avidity_target = TCE["avidity_target"]
    
    system.add_simple("plasma", ["C", f"{drug}"], [f"C-{drug}"], on_C, off_C)
    
    for organ in centrals + tumors + organs:
      system.add_simple(organ["name"], ["C", f"{drug}"], [f"C-{drug}"], on_C, off_C)
      
      system.add_simple(organ["name"], [f"{drug}", "A"], [f"{drug}-A"], on_A, off_A)
      system.add_simple(organ["name"], [f"{drug}", "B"], [f"{drug}-B"], on_B, off_B)
      system.add_simple(organ["name"], [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * avidity_target, off_B)
      system.add_simple(organ["name"], [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * avidity_target, off_A)
      
      for target in targets:
        system.add_simple(organ["name"], ["C", f"{drug}-{target}"], [f"C-{drug}-{target}"], on_C, off_C)
      
      for effector in effectors:
        system.add_simple(organ["name"], [f"{effector}-{drug}", "A"], [f"{effector}-{drug}-A"], on_A, off_A)
        system.add_simple(organ["name"], [f"{effector}-{drug}", "B"], [f"{effector}-{drug}-B"], on_B, off_B)
        system.add_simple(organ["name"], [f"{effector}-{drug}-A", "B"], [f"{effector}-{drug}-AB"], on_B * avidity_target, off_B)
        system.add_simple(organ["name"], [f"{effector}-{drug}-B", "A"], [f"{effector}-{drug}-AB"], on_A * avidity_target, off_A)
  
  # internalization
  system.add_process(TCE["internalization"])
  
  # initial concentrations
  for central in centrals:
    system.add_x("C", central["name"], 124000 * central["num_T"] / central["volume"] / units.avagadro)
    system.add_x("A", central["name"], central["conc_A"])
    system.add_x("B", central["name"], central["conc_B"])
  
  for tumor in tumors:
    system.add_x("C", tumor["name"], 124000 * tumor["density_T"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("A", tumor["name"], tumor["num_A"] * tumor["density_cell"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("B", tumor["name"], tumor["num_B"] * tumor["density_cell"] / tumor["volume_interstitial_proportion"] / units.avagadro)
  
  for organ in organs:
    system.add_x("C", organ["name"], 124000 * organ["num_T"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("A", organ["name"], organ["num_A"] * organ["num_cell"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("B", organ["name"], organ["num_B"] * organ["num_cell"] / organ["volume_interstitial"] / units.avagadro)
  
  return system



############# plot #############

def plot(system, name):
  #pickle.dump(system, open(f"{name}.pickle", "wb"))
  
  groups = [antigens_effector,
            antigens_target,
            drugs,
            dimers_effector,
            dimers_target,
            trimers]
  labels = ["effector", "target", "drug", "effector-drug", "drug-target", "effector-drug-target"]
  colors = ["tab:orange", "tab:blue", "black", "wheat", "skyblue", "tab:purple"]
  linestyles = ["solid", "solid", "solid", "solid", "solid", "solid", "solid"]
  system.plot(compartments = system.compartments, 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")
  
  groups = [[analyte for analyte in analytes if re.search("[mn][mn]", analyte)],
            [analyte for analyte in analytes if re.search("[n][mn]", analyte)],
            [analyte for analyte in analytes if re.search("[mn][n]", analyte)]]
  labels = ["(effector)-xx-(target)", "(effector)-nx-(target)", "(effector)-xn-(target)"]
  colors = ["black", "tab:orange", "tab:blue"]
  linestyles = ["solid", "dashed", "dashed"]
  system.plot(compartments = system.compartments, 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_drugs.png")
  
  groups = [["A"], ["B"],
            [analyte for analyte in analytes if re.match("[mn][mn]-A$", analyte)],
            [analyte for analyte in analytes if re.match("[mn][mn]-B$", analyte)],
            [analyte for analyte in analytes if re.match("[mn][mn]-AB$", analyte)],
            [analyte for analyte in analytes if re.search("-[mn][mn]-A", analyte)],
            [analyte for analyte in analytes if re.search("-[mn][mn]-B", analyte)],
            [analyte for analyte in analytes if re.search("-[mn][mn]-AB", analyte)]]
  labels = ["A", "B", "xxx-A", "xxx-B", "xxx-AB", "effector-xxx-A", "effector-xxx-B", "effector-xxx-AB"]
  colors = ["tab:blue", "tab:red", 
            "skyblue", "pink", "tab:purple", "skyblue", "pink", "tab:purple"]
  linestyles = ["solid", "solid", "dashed", "dashed", "dashed", "solid", "solid", "solid"]
  system.plot(compartments = system.compartments, 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_targets.png")
