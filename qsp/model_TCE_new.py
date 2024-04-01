from .qsp import *
import re
import itertools

### this model is mostly from ...

drugs = [f"{c}{a}{b}" for c in ("m", "n") for a in ("m", "n") for b in ("m", "n")]
targets = ["C", "A", "B", "AB"]
antigens = ["C", "A", "B"]
dimers = [f"{drug}-{target}" for drug in drugs for target in targets]
analytes = drugs + antigens + dimers



############ processes ############

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


def add_two_dicts(a, b):
  return dict(list(a.items()) + list(b.items()) + [(k, a[k] + b[k]) for k in set(b) & set(a)])

class transform:
  def __init__(self):
    self.system = None
    self.Qs = dict()
    
    self.analyteses_ = []
    self.analyteses_.append([analytes.index(f"{drug}") for drug in drugs])
    for target in targets:
      self.analyteses_.append([analytes.index(f"{drug}-{target}") for drug in drugs])
  
  def add(self, linker, reactant, products):
    self.system = None
    reactant_ = drugs.index(reactant)
    products_ = [drugs.index(product) for product in products] 
    
    for compartment, rate in linker:
      if compartment not in self.Qs:
        self.Qs[compartment] = np.zeros([len(drugs), len(drugs)])
      self.Qs[compartment][reactant_, reactant_] -= rate.number(1/units.h)
      self.Qs[compartment][reactant_, products_] += rate.number(1/units.h)
  
  def __add__(self, transform2): 
    buffer = transform()
    buffer.Qs = add_two_dicts(self.Qs, transform2.Qs)
    return buffer
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.Qs_ = {system.compartments.index(compartment):Q for compartment, Q in self.Qs.items()}
    
    t = t.number(units.h)
    for compartment_, Q in self.Qs_.items():
      for analytes_ in self.analyteses_:
        system.x[analytes_, compartment_] = system.x[analytes_, compartment_] @ expm(Q * t)


class internalization:
  def __init__(self, rates):
    self.system = None
    
    q = np.zeros(len(dimers))
    Q = np.zeros([len(dimers), len(analytes)])
    for target, products, rate in rates:
      idx_dimers = [dimers.index(f"{drug}-{target}") for drug in drugs if f"{drug}-{target}" in dimers]
      idx_products = [analytes.index(product) for product in products]
      q[idx_dimers] -= rate.number(1/units.h)
      for i in idx_dimers:
        np.add.at(Q, (i, idx_products), 1) # there may be repeated antigens
    
    self.q = q
    self.Q = Q
    self.idx_dimers = [analytes.index(dimer) for dimer in dimers]
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.compartments_ = [system.compartments.index(compartment) for compartment in system.compartments]
    
    t = t.number(units.h)
    for compartment_ in self.compartments_:
      delta_dimers = system.x[self.idx_dimers, compartment_] * (1 - np.exp(self.q * t))
      system.x[self.idx_dimers, compartment_] -= delta_dimers
      system.x[:, compartment_] += delta_dimers @ self.Q


############ drugs ############

GBR1302 = {}
GBR1302.update({"A": "HER2", "B": None})
GBR1302.update({"off_C": 10**-4 / units.s, "affn_C": 10 * units.nM, "affm_C": 1000 * units.nM})
GBR1302.update({"off_A": 10**-4 / units.s, "affn_A": 10 * units.nM, "affm_A": 1000 * units.nM})
GBR1302.update({"off_B": 10**-4 / units.s, "affn_B": math.inf * units.nM, "affm_B": math.inf * units.nM})
GBR1302.update({"avidity_effector": 1, "avidity_target": 1})
GBR1302.update({"clearance": math.log(2)/(70 * units.h)}); GBR1302["smalls"] = []
GBR1302["cleavage"] = None
GBR1302["internalization"] = internalization(rates_effector = [("C", ["C"], 0.1 / units.h)],
                                             rates_target = [("A", ["A"], 0.1 / units.h)])


############ model ############

def model(TCE, plasma, lymph, tumors, organs, connect_tumors = True):
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
  
  # target binding
  for drug in [f"{c}{a}{b}" for c in ("m", "n") for a in ("m", "n") for b in ("m", "n")]:
    off_C = TCE["off_C"]; on_C = {"n":TCE["off_C"] / TCE["affn_C"], "m":TCE["off_C"] / TCE["affm_C"]}[drug[0]]
    off_A = TCE["off_A"]; on_A = {"n":TCE["off_A"] / TCE["affn_A"], "m":TCE["off_A"] / TCE["affm_A"]}[drug[1]]
    off_B = TCE["off_B"]; on_B = {"n":TCE["off_B"] / TCE["affn_B"], "m":TCE["off_B"] / TCE["affm_B"]}[drug[2]]
    avidity_target = TCE["avidity_target"]
    
    for organ in centrals + tumors + organs:
      system.add_simple(organ["name"], ["C", f"{drug}"], [f"{drug}-C"], on_C, off_C)
      
      system.add_simple(organ["name"], [f"{drug}", "A"], [f"{drug}-A"], on_A, off_A)
      system.add_simple(organ["name"], [f"{drug}", "B"], [f"{drug}-B"], on_B, off_B)
      system.add_simple(organ["name"], [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * avidity_target, off_B)
      system.add_simple(organ["name"], [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * avidity_target, off_A)
  
  # mask cleavage
  if TCE["cleavage"] is not None:
    system.add_process(TCE["cleavage"])
  
  # internalization
  if TCE["internalization"] is not None:
    system.add_process(TCE["internalization"])
  
  # initial concentrations
  for central in centrals:
    system.add_x("C", central["name"], 124000 * central["num_T"] / central["volume"] / units.avagadro)
    if TCE['A'] is not None:
      system.add_x("A", central["name"], central[f"conc_{TCE['A']}"])
    if TCE['B'] is not None:
      system.add_x("B", central["name"], central[f"conc_{TCE['B']}"])
  
  for tumor in tumors:
    system.add_x("C", tumor["name"], 124000 * tumor["density_T"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    if TCE['A'] is not None:
      system.add_x("A", tumor["name"], tumor[f"num_{TCE['A']}"] * tumor["density_cell"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    if TCE['B'] is not None:
      system.add_x("B", tumor["name"], tumor[f"num_{TCE['B']}"] * tumor["density_cell"] / tumor["volume_interstitial_proportion"] / units.avagadro)
  
  for organ in organs:
    system.add_x("C", organ["name"], 124000 * organ["num_T"] / organ["volume_interstitial"] / units.avagadro)
    if TCE['A'] is not None:
      system.add_x("A", organ["name"], organ[f"num_{TCE['A']}"] * organ["num_cell"] / organ["volume_interstitial"] / units.avagadro)
    if TCE['B'] is not None:
      system.add_x("B", organ["name"], organ[f"num_{TCE['B']}"] * organ["num_cell"] / organ["volume_interstitial"] / units.avagadro)
  
  return system



############# plot #############

def plot(system, name):
  main_drugs = [f"{c}{a}{b}" for c in ("m", "n") for a in ("m", "n") for b in ("m", "n")]
  
  groups = [["P"],
            ["C"],
            ["A", "B"],
            main_drugs,
            ["p"],
            [f"p-P"],
            [f"{drug}-{target}" for drug in main_drugs for target in ["P", "PC"]],
            [f"{drug}-{target}" for drug in main_drugs for target in ["C", "PC"]],
            [f"{drug}-{target}" for drug in main_drugs for target in ["A", "B", "AB"]]]
  labels = ["PD1", "CD3", "target", "drug (mainbody)", "aPD1", "aPD1-PD1", "drug-PD1", "drug-CD3", "drug-target"]
  colors = ["tab:red", "tab:orange", "tab:blue", "black", "black", "pink", "pink", "wheat", "skyblue"]
  linestyles = ["solid", "solid", "solid", "solid", "dashed", "dashed", "solid", "solid", "solid"]
  system.plot(compartments = system.compartments, 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")
