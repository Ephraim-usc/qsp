from .qsp import *
import re
import itertools

### this model is mostly from ...

cells = ["T", "B"]
areas = [200, 254]

antigens = ["[T]C", "[B]A", "H"] # H: hydroxyapatite, C: CD3, A: CD19
bindings = ["[T]C", "[B]A", "H"]
drugs = [f"{c}{a}{h}" for c in ("m", "n") for a in ("m", "n") for h in ("m", "n")]
dimers = [f"{binding}-{drug}" for binding in bindings for drug in drugs]
analytes = antigens + drugs + dimers



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
    for binding in bindings:
      self.analyteses_.append([analytes.index(f"{binding}-{drug}") for drug in drugs])
  
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
      self.Qs_ = {system.compartments.index(compartment):Q for compartment, Q in self.Qs.items() if compartment in system.compartments}
    
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
      idx_dimers = [dimers.index(f"{target}-{drug}") for drug in drugs if f"{target}-{drug}" in dimers]
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


contact_area_time = 10*units.s * math.pi*units.um**2
contact_freq = 4 * math.pi * (6.45 * units.um**2 / units.MIN) * 4*units.um # T cell diffusion rate according to https://pubmed.ncbi.nlm.nih.gov/29044117/ Figure 2D
contact_freqs = {"plasma": contact_freq * 10, "lymph": contact_freq * 0.1, "default": contact_freq}

class kill:
  def __init__(self, compartments, on2ds, contact_freq = contact_freqs,
               effector = "T", target = "B", contact_area_time = contact_area_time, synapse_efficiency = 0.01, damage = 0.5, regen = 0.1 / units.h):
    self.system = None
    self.compartments = compartments
    self.on2ds = on2ds # pandas data frame of unit um**2/s
    self.effector = effector
    self.target = target
    
    self.contact_freqs = [(contact_freqs[compartment].number(units.ml / units.h) if compartment in contact_freqs else contact_freqs["default"].number(units.ml / units.h)) for compartment in compartments]
    self.contact_area_time = contact_area_time.number(units.um**2 * units.s)
    self.synapse_efficiency = synapse_efficiency
    self.damage = damage
    self.regen = regen.number(1/units.h)
  
  def renormalize(self):
    deaths = self.hp <= 0
    alives = np.logical_not(deaths)
    for i in range(len(self.compartments_)):
      self.hp[deaths[:, i], i] = np.random.choice(self.hp[alives[:, i], i], deaths[:, i].sum())
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments if compartment in system.compartments]
      self.hp = np.ones([100000, len(self.compartments_)])

      ligands_effector = [ligand for ligand in self.on2ds.index.values if ligand in system.analytes]
      ligands_target = [ligand for ligand in self.on2ds.columns.values if ligand in system.analytes]
      
      self.effector_ = system.cells.index(self.effector)
      self.target_ = system.cells.index(self.target)
      self.ligands_effector_ = [system.analytes.index(ligand) for ligand in ligands_effector]
      self.ligands_target_ = [system.analytes.index(ligand) for ligand in ligands_target]
      self.on2ds_ = self.on2ds.loc[ligands_effector, ligands_target]
    
    t = t.number(units.h)
    contacts_expected = self.contact_freqs * system.c[self.effector_, self.compartments_] * t # average number of contacts with effector cells, for each target cell
    
    trimers = self.contact_area_time * np.array([system.y[self.ligands_effector_, compartment_] @ self.on2ds_ @ system.y[self.ligands_target_, compartment_] for compartment_ in self.compartments_])
    probs = 1 - (1 - self.synapse_efficiency)**trimers # probability that a contact would form a synapse
    
    contacts = np.stack([np.random.poisson(_, int(1e5)) for _ in contacts_expected], axis = 1)
    damages = np.random.binomial(contacts, probs) * self.damage
    self.hp = np.minimum(1.0, self.hp - damages + self.regen * t)
    
    deaths = (self.hp <= 0).mean(axis = 0)
    system.cell_death_(self.target_, self.compartments_, deaths)
    self.renormalize()


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
BD.update({"off_C": 10**-4 / units.s, "affn_C": 260 * units.nM, "affm_C": 26000 * units.nM})
BD.update({"off_A": 10**-4 / units.s, "affn_A": 1.49 * units.nM, "affm_A": 149 * units.nM})
BD.update({"off_H": 10**-4 / units.s, "affn_H": math.inf * units.nM, "affm_H": math.inf * units.nM})
BD.update({"on2dn_C": 1e-4 * units.um**2 / units.s, "on2dm_C": 1e-6 * units.um**2 / units.s})
BD.update({"on2dn_A": 1e-4 * units.um**2 / units.s, "on2dm_A": 1e-6 * units.um**2 / units.s})
BD.update({"clearance": math.log(2)/(120 * units.h)})
BD["smalls"] = []
BD["internalization"] = internalization(rates = [("[T]C", ["[T]C"], 0.1 / units.h),
                                                 ("[B]A", ["[B]A"], 0.1 / units.h)])
BD["cleavage"] = transform()
for a, h in itertools.product(("m", "n"), ("m", "n")):
    BD["cleavage"].add(linker = linker, reactant = f"m{a}{h}", products = [f"n{a}{h}"])
for c, h in itertools.product(("m", "n"), ("m", "n")):
    BD["cleavage"].add(linker = linker, reactant = f"{c}m{h}", products = [f"{c}n{h}"])
for c, a in itertools.product(("m", "n"), ("m", "n")):
    BD["cleavage"].add(linker = linker, reactant = f"{c}{a}m", products = [f"{c}{a}n"])


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
    off_H = TCE["off_H"]; on_H = {"n":TCE["off_H"] / TCE["affn_H"], "m":TCE["off_H"] / TCE["affm_H"]}[drug[2]]
    
    for organ in centrals + organs:
      system.add_simple(organ["name"], ["[T]C", f"{drug}"], [f"[T]C-{drug}"], on_C, off_C)
      system.add_simple(organ["name"], ["[B]A", f"{drug}"], [f"[B]A-{drug}"], on_A, off_A)
      system.add_simple(organ["name"], ["H", f"{drug}"], [f"H-{drug}"], on_H, off_H)
      
      system.add_simple(organ["name"], ["[T]C", f"H-{drug}"], [f"[T]C-{drug}", "H"], on_C, off_C)
      system.add_simple(organ["name"], ["[B]A", f"H-{drug}"], [f"[B]A-{drug}", "H"], on_A, off_A)
  
  ligands_effector = np.array(system.analytes)[np.array(system.ligands[0])]
  ligands_target = np.array(system.analytes)[np.array(system.ligands[1])]
  on2ds = pd.DataFrame(0, index = ligands_effector, columns = ligands_target) # in unit of um**2/s
  for drug in drugs:
    on2d_C = {"n":TCE["on2dn_C"], "m":TCE["on2dm_C"]}[drug[0]].number(units.um**2 / units.s)
    on2d_A = {"n":TCE["on2dn_A"], "m":TCE["on2dm_C"]}[drug[1]].number(units.um**2 / units.s)
    on2ds.loc[f"[T]C-{drug}", f"[B]A"] = on2d_A
    on2ds.loc[f"[T]C", f"[B]A-{drug}"] = on2d_A
  
  system.add_process(kill(compartments, on2ds))

  
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
              output = f"{name}_summary.png")


############# demo ###############
'''
from qsp import *
from qsp.human import *
from qsp.model_TCE_B import *

system = model(BD, plasma, lymph, [bone, lung, liver])
for _ in range(3):
  system.add_x("plasma", "nnn", 10 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "unmasked")
system.plot_cell(output = "unmasked.png")

system = model(BD, plasma, lymph, [bone, lung, liver])
for _ in range(3):
  system.add_x("plasma", "mmn", 100 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "masked")
system.plot_cell(output = "masked_cell.png")

TCE = BD.copy()
TCE.update({"off_H": 10**-4 / units.s, "affn_H": 1 * units.nM, "affm_H": 100 * units.nM})
system = model(TCE, plasma, lymph, [bone, lung, liver])
for _ in range(3):
  system.add_x("plasma", "mmn", 100 * units.nM)
  system.run(24 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "masked_HA")
system.plot_cell(output = "masked_HA.png")

'''
