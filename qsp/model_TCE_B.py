from .qsp import *
import re
import itertools

### this model is mostly from ...

cells = ["T", "C"]

solubles = ["A", "B"]
ligands = ["[T]CD3", "[C]A", "[C]B"]
drugs = [f"{c}{a}{b}" for c in ("m", "n") for a in ("m", "n") for b in ("m", "n")]
dimers = [f"{binding}-{drug}" for binding in ["A", "B", "AB", "[T]CD3", "[C]A", "[C]B", "[C]AB"] for drug in drugs]
analytes = solubles + ligands + drugs + dimers



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

linker = [("plasma", 0.07 / units.d), 
          ("lymph", 0.07 / units.d), 
          ("bone", 0.07 / units.d), 
          ("liver", 0.07 / units.d), 
          ("lung", 0.07 / units.d), 
          ("SI", 0.07 / units.d), 
          ("gallbladder", 0.07 / units.d)]

BD = {}
BD.update({"A": "CD19", "B": "BAFFR"})
BD.update({"off_C": 10**-4 / units.s, "affn_CD3": 10 * units.nM, "affm_CD3": 1000 * units.nM, "aff2d_CD3": None})
BD.update({"off_A": 10**-4 / units.s, "affn_A": 10 * units.nM, "affm_A": 1000 * units.nM, "aff2d_A": None})
BD.update({"off_B": 10**-4 / units.s, "affn_B": 10 * units.nM, "affm_B": 1000 * units.nM, "aff2d_B": None})
BD.update({"avidity": 20})
BD.update({"clearance": math.log(2)/(70 * units.h)})
BD["smalls"] = []
BD["internalization"] = internalization(rates = [("C", ["C"], 0.1 / units.h),
                                                 ("A", ["A"], 0.1 / units.h),
                                                 ("B", ["B"], 0.1 / units.h),
                                                 ("AB", ["A", "B"], 0.02 / units.h)])
BD["cleavage"] = transform()
for a, b in itertools.product(("m", "n"), ("m", "n")):
    BD["cleavage"].add(linker = linker, reactant = f"m{a}{b}", products = ["p", f"n{a}{b}"])
for c, b in itertools.product(("m", "n"), ("m", "n")):
    BD["cleavage"].add(linker = linker, reactant = f"{c}m{b}", products = [f"{c}n{b}"])
for c, a in itertools.product(("m", "n"), ("m", "n")):
    BD["cleavage"].add(linker = linker, reactant = f"{c}{a}m", products = [f"{c}{a}n"])


############ model ############

def model(TCE, plasma, lymph, tumors, organs, connect_tumors = True):
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + tumors + organs]
  system = System(compartments, analytes, cells)
  system.centrals = [plasma, lymph]
  system.tumors = tumors
  system.organs = organs
  
  for central in centrals:
    system.set_volume(central["name"], central["volume"])
  for tumor in tumors:
    system.set_volume(tumor["name"], tumor["volume"] * tumor["volume_interstitial_proportion"])
  for organ in organs:
    system.set_volume(organ["name"], organ["volume_interstitial"])
  
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
  for drug in drugs:
    off_CD3 = TCE["off_C"]; on_CD3 = {"n":TCE["off_CD3"] / TCE["affn_CD3"], "m":TCE["off_CD3"] / TCE["affm_CD3"]}[drug[0]]
    off_A = TCE["off_A"]; on_A = {"n":TCE["off_A"] / TCE["affn_A"], "m":TCE["off_A"] / TCE["affm_A"]}[drug[1]]
    off_B = TCE["off_B"]; on_B = {"n":TCE["off_B"] / TCE["affn_B"], "m":TCE["off_B"] / TCE["affm_B"]}[drug[2]]
    avidity_effector = TCE["avidity_effector"]
    avidity_target = TCE["avidity_target"]
    
    for organ in centrals + tumors + organs:
      system.add_simple(organ["name"], ["[T]CD3", f"{drug}"], [f"[T]CD3-{drug}"], on_C, off_C)
      
      system.add_simple(organ["name"], ["[C]A", f"{drug}"], [f"[C]A-{drug}"], on_A, off_A)
      system.add_simple(organ["name"], ["[C]B", f"{drug}", "B"], [f"[C]B-{drug}"], on_B, off_B)
      system.add_simple(organ["name"], ["[C]B", f"{drug}-A", "B"], [f"[C]AB-{drug}"], on_B * avidity_target, off_B)
      system.add_simple(organ["name"], ["[C]A", f"{drug}-B", "A"], [f"[C]AB-{drug}"], on_A * avidity_target, off_A)
  
  # mask cleavage
  if TCE["cleavage"] is not None:
    system.add_process(TCE["cleavage"])
  
  # internalization
  if TCE["internalization"] is not None:
    system.add_process(TCE["internalization"])
  
  # initial concentrations
  for central in centrals:
    system.add_c(central["name"], "T", central["num_T"] / central["volume"], ["CD3"], [124000])
    system.add_c(central["name"], "B", central["num_B"] / central["volume"], ["A", "B"], [20000, 10000])
  
  for tumor in tumors:
    system.add_c(tumor["name"], "T", tumor["density_T"] / tumor["volume_interstitial_proportion"], ["CD3"], [124000])
    system.add_c(tumor["name"], "B", tumor["density_B"] / tumor["volume_interstitial_proportion"], ["A", "B"], [20000, 10000])
  
  for organ in organs:
    system.add_c(central["name"], "T", organ["num_T"] / organ["volume_interstitial"], ["CD3"], [124000])
    system.add_c(central["name"], "B", organ["num_B"] / organ["volume_interstitial"], ["A", "B"], [20000, 10000])
  
  return system



############# plot #############

def plot(system, name):
  groups = [["C"],
            ["A", "B"],
            drugs,
            [f"{drug}-{target}" for drug in drugs for target in ["C"]],
            [f"{drug}-{target}" for drug in drugs for target in ["A", "B", "AB"]]]
  labels = ["CD3", "target", "drug", "drug-CD3", "drug-target"]
  colors = [ "tab:orange", "tab:blue", "black", "wheat", "skyblue"]
  linestyles = ["solid", "solid", "solid", "solid", "solid"]
  system.plot(compartments = system.compartments, 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")


############# demo ###############

