from .qsp import *
import itertools
import sympy

############ constants ############

molecular_weight = 150000 * units.g/units.mol

def evalf_array(array, subs):
  return

def hill(x, EMAX, EC50, coef):
  return EMAX * (x**coef) / (x**coef + EC50**coef)

class process_equilibrium:
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

class process_transform:
  def __init__(self, cells, ligands):
    self.system = None
    self.ligands = ligands.copy()
    self.Qs = [ligand.Q for ligand in ligands]
    self.analytes_dict = {}
    self.analyteseses = [[[f"{ligand.name}:{state}" for state in ligand.states]] + \
                         [[f"{cell.name}:{binding}-{ligand.name}:{state}" for state in ligand.states] for cell in cells for binding in ligand.get_bindings(cell)] \
                         for ligand in ligands]
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.idxeseses = [[[system.analytes.index(analyte) for analyte in analytes] for analytes in analyteses] for analyteses in self.analyteseses]
    
    for idx_compartment in range(system.n_compartments):
      for idxeses, Q in zip(self.idxeseses, self.Qs):
        Q = np.array(Q.evalf(system.signals[idx_compartment]), dtype = np.float64)
        for idxes in idxeses:
          system.x[idxes, idx_compartment] = system.x[idxes, idx_compartment] @ expm(Q * t.number(units.h))







class Signal:
  def __init__(self):
    self.dict = {}
  
  def __getitem__(self, key):
    if key not in self.dict: 
      self.dict[key] = sympy.Symbol(key)
    return self.dict[key]

signal = Signal()


class Ligand:
  def __init__(self, name, n_sites, site_states, targets, clearance, smalls = None):
    self.name = name
    self.n_sites = n_sites
    self.site_states = site_states
    self.states = ["".join(tmp) for tmp in itertools.product(*site_states)]
    self.n_states = len(self.states)
    self.targets = targets
    self.clearance = clearance.number(1/units.h)
    
    self.V = np.ones([self.n_sites, self.n_sites]) # avidity
    self.Q = sympy.zeros(self.n_states, self.n_states) # transforming matrix, may contain signal variables
    self.smalls = smalls if smalls is not None else [] # states that are small molecules
  
  def set_avidity(self, site_A, site_B, avidity, symmetric = True):
    self.V[site_A, site_B] = avidity
    if symmetric:
      self.V[site_B, site_A] = avidity
  
  def add_transform(self, sites, state_from, state_to, rate):
    if type(sites) is not list:
      sites = [sites]
    
    site_states_from = self.site_states.copy()
    site_states_to = self.site_states.copy()
    for site, site_state_from, site_state_to in zip(sites, state_from, state_to):
      site_states_from[site] = [state_from]
      site_states_to[site] = [state_to]
    states_from = ["".join(tmp) for tmp in itertools.product(*site_states_from)]
    states_to = ["".join(tmp) for tmp in itertools.product(*site_states_to)]
    
    idxes_from = [self.states.index(tmp) for tmp in states_from]
    idxes_to = [self.states.index(tmp) for tmp in states_to]
    for idx_from, idx_to in zip(idxes_from, idxes_to):
      self.Q[idx_from, idx_from] -= rate.number(1/units.h)
      self.Q[idx_from, idx_to] += rate.number(1/units.h)
  
  def get_bindings(self, cell, state = None): # find all binding modes formed when a ligand binds to a cell
    targets_in_markers = [["_"] + [target for target in targets if target in cell.markers] for targets in self.targets]
    bindings = ["".join(complex) for complex in itertools.product(*targets_in_markers)][1:] # remove the all non-binding mode
    return bindings
  
  def print(self):
    Q = pd.DataFrame(np.array(self.Q), index = self.states, columns = self.states)
    print(Q)


X = Ligand(name = "X", 
           n_sites = 3, 
           site_states = [["n"], ["m", "n"], ["m", "n"]], 
           targets = [["P"], ["α", "R"], ["α", "R"]],
           clearance = math.log(2)/(70 * units.h))
X.add_transform(1, "m", "n", signal["enzyme"] * 0.01/units.h)
X.add_transform(2, "m", "n", signal["enzyme"] * 0.01/units.h)

IL2 = Ligand("IL2", 1, [["n"]], [["α", "R"]], math.log(2)/(70 * units.h), smalls = ["n"])


class Cell:
  def __init__(self, name, markers, initials, birth = None, death = None, prolif = None):
    self.name = name
    self.markers = markers
    self.initials = initials
    self.birth = birth.number(units.nM/units.h) if birth is not None else None
    self.death = death.number(1/units.h) if death is not None else None
    self.prolif = prolif.number(1/units.h) if prolif is not None else None

Treg = Cell("Treg", ["P", "α"], [30000, 300], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
nTh = Cell("nTh", [], [], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.002 / units.d)
aTh = Cell("aTh", ["P", "R"], [30000, 300], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Th = Cell("Th", ["P", "R"], [30000, 300], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
nTm = Cell("nTm", [], [], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.002 / units.d)
aTm = Cell("aTm", ["P", "α"], [30000, 1500], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Tm = Cell("Tm", ["P", "R"], [30000, 1500], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Teff = Cell("Teff", ["P", "α"], [30000, 1500], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Tex = Cell("Tex", ["P", "α"], [30000, 1500], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.1 / units.d)
NK = Cell("NK", ["α"], [30000, 3000], birth = signal["tumor"] * 1 * units.nM/units.d, death = 0.02 / units.d)





def PDSystem(organs, tumors, cells, ligands):

organs = [bone, lung]
tumors = [UT44]
cells = [Treg, nTh, aTh, Th, nTm, aTm, Tm, Teff, Tex, NK]
ligands = [X, IL2]

analytes_markers = [f"{cell.name}:{marker}" for cell in cells for marker in cell.markers]
analytes_ligands = [f"{ligand.name}:{state}" for ligand in ligands for state in ligand.states]
analytes_dimers = [f"{cell.name}:{binding}-{ligand.name}:{state}" for ligand in ligands for cell in cells for binding in ligand.get_bindings(cell) for state in ligand.states]
analytes = analytes_ligands + analytes_markers + analytes_dimers


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
for analyte in analytes_ligands:
  analyte = f"{ligand.name}:{state}"
  system.add_flow(analyte, compartment, None, system.get_volume(analyte, compartment) * ligand.clearance / units.h)

# small forms plasma clearance
for ligand in ligands:
  for state in ligand.smalls:
    analyte = f"{ligand.name}:{state}"
    system.add_flow(analyte, "plasma", None, system.get_volume(analyte, compartment) * math.log(2)/(45 * units.MIN))

for analyte in analytes_ligands:
  # drug tumor flow
  for tumor in tumors:
    system.add_flow(analyte, "plasma", tumor["name"], tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
    system.add_flow(analyte, tumor["name"], "plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
  
  # drug organ flow
  for organ in organs:
    system.add_flow(analyte, "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
    system.add_flow(analyte, organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
    system.add_flow(analyte, "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))

# exchange ligands between tumors if tumors are connected
if connect_tumors:
  system.add_process(equilibrium([tumor["name"] for tumor in tumors], analytes_ligands))

# mask cleavage
system.add_process(transform(compartments = [central["name"] for central in centrals] + [organ["name"] for organ in system.organs], rates = TCE["cleavage_plasma"]))
system.add_process(transform(compartments = [tumor["name"] for tumor in tumors], rates = TCE["cleavage_tumor"]))





















    

############ constants ############

molecular_weight = 150000 * units.g/units.mol

def hill(x, EMAX, EC50, coef):
  return EMAX * (x**coef) / (x**coef + EC50**coef)

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


proliferation_params = {}
proliferation_params["prolif_EMAX_8"] = (math.log(2) / (9.1 * units.h * 1.34)).number(1 / units.h)
proliferation_params["prolif_EMAX_4"] = (math.log(2) / (9.1 * units.h * 2.34)).number(1 / units.h)
proliferation_params["prolif_EMAX_nk"] = (math.log(2) / (9.1 * units.h * 1.22)).number(1 / units.h)
proliferation_params["prolif_EC50_8"] = 109
proliferation_params["prolif_EC50_4"] = 109
proliferation_params["prolif_EC50_nk"] = 109
proliferation_params["prolif_hill_8"] = 3.06
proliferation_params["prolif_hill_4"] = 1
proliferation_params["prolif_hill_nk"] = 1.3
proliferation_params["death_8"] = (math.log(2) / (70.95 * units.h)).number(1 / units.h)
proliferation_params["death_4"] = (math.log(2) / (138.55 * units.h)).number(1 / units.h)
proliferation_params["death_nk"] = (math.log(2) / (44.65 * units.h)).number(1 / units.h)
proliferation_params["births_8"] = {} # birth rates are initialized at runtime according to the initial number of immume cells
proliferation_params["births_4"] = {}
proliferation_params["births_nk"] = {}

migration_params = {}
migration_params["influx"] = (math.log(2) / (8.3 * units.h)).number(1 / units.h)
migration_params["marg_EMAX"] = None # natural marginalization rates are estimated according to equilibrium
migration_params["marg_EMAX"] = None
migration_params["marg_EMAX"] = None
migration_params["marg_EMAX_8"] = 20
migration_params["marg_EMAX_4"] = 8.24
migration_params["marg_EMAX_nk"] = 20
migration_params["marg_EC50_8"] = 350
migration_params["marg_EC50_4"] = 350
migration_params["marg_EC50_nk"] = 350
migration_params["marg_hill_8"] = 1
migration_params["marg_hill_4"] = 1
migration_params["marg_hill_nk"] = 1

class PD:
  def __init__(self, proliferation_params = proliferation_params, migration_params = migration_params):
    self.system = None
    self.params = {**proliferation_params, **migration_params}.copy()
    
    self.analytes = {}
    self.analytes["A"] = "A"; self.analytes["B"] = "B"
    self.analytes["8"] = "8"; self.analytes["C8"] = "C8"; self.analytes["R8"] = "R8"
    self.analytes["4"] = "4"; self.analytes["C4"] = "C4"; self.analytes["R4"] = "R4"
    self.analytes["nk"] = "nk"; self.analytes["Rnk"] = "Rnk"
    
    self.analytes["all_8"] = ["C8", "R8"] + [f"{effector}-{drug}" for effector in ["C8", "R8", "CR8"] for drug in drugs] + [f"{effector}-{drug}-{target}" for effector in ["C8", "R8", "CR8"] for drug in drugs for target in targets]
    self.analytes["all_4"] = ["C4", "R4"] + [f"{effector}-{drug}" for effector in ["C4", "R4", "CR4"] for drug in drugs] + [f"{effector}-{drug}-{target}" for effector in ["C4", "R4", "CR4"] for drug in drugs for target in targets]
    self.analytes["all_nk"] = ["Rnk"] + [f"Rnk-{drug}" for drug in drugs] + [f"Rnk-{drug}-{target}" for drug in drugs for target in targets]
    
    self.analytes["Rcomplexes_8"] = [f"{effector}-{drug}" for effector in ["R8", "CR8"] for drug in drugs] + [f"{effector}-{drug}-{target}" for effector in ["R8", "CR8"] for drug in drugs for target in targets]
    self.analytes["Rcomplexes_4"] = [f"{effector}-{drug}" for effector in ["R4", "CR4"] for drug in drugs] + [f"{effector}-{drug}-{target}" for effector in ["R4", "CR4"] for drug in drugs for target in targets]
    self.analytes["Rcomplexes_nk"] = [f"Rnk-{drug}" for drug in drugs] + [f"Rnk-{drug}-{target}" for drug in drugs for target in targets]
    
    self.analytes["trimers_8_A"] = [f"{effector}-{drug}-{target}" for effector in ["C8", "R8", "CR8"] for drug in drugs for target in ["A", "AB"]]
    self.analytes["trimers_8_B"] = [f"{effector}-{drug}-{target}" for effector in ["C8", "R8", "CR8"] for drug in drugs for target in ["B", "AB"]]
    self.analytes["trimers_4_A"] = [f"{effector}-{drug}-{target}" for effector in ["C4", "R4", "CR4"] for drug in drugs for target in ["A", "AB"]]
    self.analytes["trimers_4_B"] = [f"{effector}-{drug}-{target}" for effector in ["C4", "R4", "CR4"] for drug in drugs for target in ["B", "AB"]]
    self.analytes["trimers_nk_A"] = [f"Rnk-{drug}-{target}" for drug in drugs for target in ["A", "AB"]]
    self.analytes["trimers_nk_B"] = [f"Rnk-{drug}-{target}" for drug in drugs for target in ["B", "AB"]]
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.index_plasma = system.compartments.index("plasma")
      self.index_lymph = system.compartments.index("lymph")
      self.index_tumors = [system.compartments.index(tumor["name"]) for tumor in system.tumors]
      self.index_organs = [system.compartments.index(organ["name"]) for organ in system.organs]
      
      self.index = {}
      for key, value in self.analytes.items():
        if type(value) is list:
          self.index[key] = [system.analytes.index(analyte) for analyte in value]
        else:
          self.index[key] = system.analytes.index(value)
      
      for index_compartment in [self.index_lymph] + self.index_tumors + self.index_organs:
        self.params["births_8"][index_compartment] = system.x[self.index["8"], index_compartment] * self.params["death_8"]
        self.params["births_4"][index_compartment] = system.x[self.index["4"], index_compartment] * self.params["death_4"]
        self.params["births_nk"][index_compartment] = system.x[self.index["nk"], index_compartment] * self.params["death_nk"]
      
      self.params["marg_8"] = self.params["influx"] * system.V[self.index["8"], self.index_lymph] * system.x[self.index["8"], self.index_lymph] / (system.V[self.index["8"], self.index_plasma] * system.x[self.index["8"], self.index_plasma])
      self.params["marg_4"] = self.params["influx"] * system.V[self.index["4"], self.index_lymph] * system.x[self.index["4"], self.index_lymph] / (system.V[self.index["4"], self.index_plasma] * system.x[self.index["4"], self.index_plasma])
      self.params["marg_nk"] = self.params["influx"] * system.V[self.index["nk"], self.index_lymph] * system.x[self.index["nk"], self.index_lymph] / (system.V[self.index["nk"], self.index_plasma] * system.x[self.index["nk"], self.index_plasma])
    
    # body of the process
    index = self.index
    t = t.number(units.h)
    
    # proliferation
    for index_compartment in [self.index_lymph] + self.index_tumors + self.index_organs:
      prolif_8 = hill(system.x[index["Rcomplexes_8"], index_compartment].sum() / system.x[index["8"], index_compartment], self.params["prolif_EMAX_8"], self.params["prolif_EC50_8"], self.params["prolif_hill_8"])
      prolif_4 = hill(system.x[index["Rcomplexes_4"], index_compartment].sum() / system.x[index["4"], index_compartment], self.params["prolif_EMAX_4"], self.params["prolif_EC50_4"], self.params["prolif_hill_4"])
      prolif_nk = hill(system.x[index["Rcomplexes_nk"], index_compartment].sum() / system.x[index["nk"], index_compartment], self.params["prolif_EMAX_nk"], self.params["prolif_EC50_nk"], self.params["prolif_hill_nk"])
      
      # target antigens released from trimers when immune cells die
      system.x[index["A"], index_compartment] += system.x[index["trimers_8_A"], index_compartment].sum() * self.params["death_8"] * t
      system.x[index["B"], index_compartment] += system.x[index["trimers_8_B"], index_compartment].sum() * self.params["death_8"] * t
      system.x[index["A"], index_compartment] += system.x[index["trimers_4_A"], index_compartment].sum() * self.params["death_4"] * t
      system.x[index["B"], index_compartment] += system.x[index["trimers_4_B"], index_compartment].sum() * self.params["death_4"] * t
      system.x[index["A"], index_compartment] += system.x[index["trimers_nk_A"], index_compartment].sum() * self.params["death_nk"] * t
      system.x[index["B"], index_compartment] += system.x[index["trimers_nk_B"], index_compartment].sum() * self.params["death_nk"] * t
      
      system.x[index["all_8"], index_compartment] += system.x[index["all_8"], index_compartment] * (-self.params["death_8"]) * t
      system.x[index["all_4"], index_compartment] += system.x[index["all_4"], index_compartment] * (-self.params["death_4"]) * t
      system.x[index["all_nk"], index_compartment] += system.x[index["all_nk"], index_compartment] * (-self.params["death_nk"]) * t
      
      system.x[index["C8"], index_compartment] += (self.params["births_8"][index_compartment] + system.x[index["8"], index_compartment] * prolif_8) * num_C_8 * t
      system.x[index["R8"], index_compartment] += (self.params["births_8"][index_compartment] + system.x[index["8"], index_compartment] * prolif_8) * num_R_8 * t
      system.x[index["C4"], index_compartment] += (self.params["births_4"][index_compartment] + system.x[index["4"], index_compartment] * prolif_4) * num_C_4 * t
      system.x[index["R4"], index_compartment] += (self.params["births_4"][index_compartment] + system.x[index["4"], index_compartment] * prolif_4) * num_R_4 * t
      system.x[index["Rnk"], index_compartment] += (self.params["births_nk"][index_compartment] + system.x[index["nk"], index_compartment] * prolif_nk) * num_R_nk * t
      
      system.x[index["8"], index_compartment] += (self.params["births_8"][index_compartment] + system.x[index["8"], index_compartment] * (-self.params["death_8"] + prolif_8)) * t
      system.x[index["4"], index_compartment] += (self.params["births_4"][index_compartment] + system.x[index["4"], index_compartment] * (-self.params["death_4"] + prolif_4)) * t
      system.x[index["nk"], index_compartment] += (self.params["births_nk"][index_compartment] + system.x[index["nk"], index_compartment] * (-self.params["death_nk"] + prolif_nk)) * t
    
    # computing total amounts of migration
    marg_8 = hill(system.x[index["Rcomplexes_8"], self.index_plasma].sum() / system.x[index["8"], self.index_plasma], self.params["marg_EMAX_8"], self.params["marg_EC50_8"], self.params["marg_hill_8"])
    marg_4 = hill(system.x[index["Rcomplexes_4"], self.index_plasma].sum() / system.x[index["4"], self.index_plasma], self.params["marg_EMAX_4"], self.params["marg_EC50_4"], self.params["marg_hill_4"])
    marg_nk = hill(system.x[index["Rcomplexes_nk"], self.index_plasma].sum() / system.x[index["nk"], self.index_plasma], self.params["marg_EMAX_nk"], self.params["marg_EC50_nk"], self.params["marg_hill_nk"])
    
    migration_8 = system.V[index["8"], self.index_lymph] * system.x[index["8"], self.index_lymph] * self.params["influx"] * t - system.V[index["8"], self.index_plasma] * system.x[index["8"], self.index_plasma] * self.params["marg_8"] * (1 + marg_8) * t
    migration_4 = system.V[index["4"], self.index_lymph] * system.x[index["4"], self.index_lymph] * self.params["influx"] * t - system.V[index["4"], self.index_plasma] * system.x[index["4"], self.index_plasma] * self.params["marg_4"] * (1 + marg_4) * t
    migration_nk = system.V[index["nk"], self.index_lymph] * system.x[index["nk"], self.index_lymph] * self.params["influx"] * t - system.V[index["nk"], self.index_plasma] * system.x[index["nk"], self.index_plasma] * self.params["marg_nk"] * (1 + marg_nk) * t
    
    migration_all_8 = system.V[index["all_8"], self.index_lymph] * system.x[index["all_8"], self.index_lymph] * self.params["influx"] * t - system.V[index["all_8"], self.index_plasma] * system.x[index["all_8"], self.index_plasma] * self.params["marg_8"] * (1 + marg_8) * t
    migration_all_4 = system.V[index["all_4"], self.index_lymph] * system.x[index["all_4"], self.index_lymph] * self.params["influx"] * t - system.V[index["all_4"], self.index_plasma] * system.x[index["all_4"], self.index_plasma] * self.params["marg_4"] * (1 + marg_4) * t
    migration_all_nk = system.V[index["all_nk"], self.index_lymph] * system.x[index["all_nk"], self.index_lymph] * self.params["influx"] * t - system.V[index["all_nk"], self.index_plasma] * system.x[index["all_nk"], self.index_plasma] * self.params["marg_nk"] * (1 + marg_nk) * t
    
    # applying the migration
    system.x[index["8"], self.index_plasma] += migration_8 / system.V[index["8"], self.index_plasma]
    system.x[index["4"], self.index_plasma] += migration_4 / system.V[index["4"], self.index_plasma]
    system.x[index["nk"], self.index_plasma] += migration_nk / system.V[index["nk"], self.index_plasma]
    
    system.x[index["all_8"], self.index_plasma] += migration_all_8 / system.V[index["all_8"], self.index_plasma]
    system.x[index["all_4"], self.index_plasma] += migration_all_4 / system.V[index["all_4"], self.index_plasma]
    system.x[index["all_nk"], self.index_plasma] += migration_all_nk / system.V[index["all_nk"], self.index_plasma]

    system.x[index["8"], self.index_lymph] -= migration_8 / system.V[index["8"], self.index_lymph]
    system.x[index["4"], self.index_lymph] -= migration_4 / system.V[index["4"], self.index_lymph]
    system.x[index["nk"], self.index_lymph] -= migration_nk / system.V[index["nk"], self.index_lymph]
    
    system.x[index["all_8"], self.index_lymph] -= migration_all_8 / system.V[index["all_8"], self.index_lymph]
    system.x[index["all_4"], self.index_lymph] -= migration_all_4 / system.V[index["all_4"], self.index_lymph]
    system.x[index["all_nk"], self.index_lymph] -= migration_all_nk / system.V[index["all_nk"], self.index_lymph]



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
    
    self.analyteses = []
    self.analyteses.append(drugs)
    for effector in effectors:
      self.analyteses.append([f"{effector}-{drug}" for drug in drugs])
    for target in targets:
      self.analyteses.append([f"{effector}-{drug}" for drug in drugs])
    for effector in effectors:
      for target in targets:
        self.analyteses.append([f"{effector}-{drug}-{target}" for drug in drugs])
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.index_compartments = [system.compartments.index(compartment) for compartment in self.compartments]
      self.indexeses = [[system.analytes.index(analyte) for analyte in _analytes] for _analytes in self.analyteses]
    
    for index_compartment in self.index_compartments:
      for indexes in self.indexeses:
        system.x[indexes, index_compartment] = system.x[indexes, index_compartment] @ expm(self.Q * t.number(units.h))


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
      self.index_compartments = [system.compartments.index(compartment) for compartment in self.compartments]
      
      self.idx_dimers_effector = [system.analytes.index(analyte) for analyte in dimers_effector]
      self.idx_dimers_target = [system.analytes.index(analyte) for analyte in dimers_target]
      self.idx_antigens_effector = [system.analytes.index(analyte) for analyte in antigens_effector]
      self.idx_antigens_target = [system.analytes.index(analyte) for analyte in antigens_target]
    
    for index_compartment in self.index_compartments:
      delta_effector = system.x[self.idx_dimers_effector, index_compartment] * (1 - np.exp(self.q_effector * t.number(units.h)))
      system.x[self.idx_dimers_effector, index_compartment] -= delta_effector
      system.x[self.idx_antigens_effector, index_compartment] += delta_effector @ self.Q_effector
      
      delta_target = system.x[self.idx_dimers_target, index_compartment] * (1 - np.exp(self.q_target * t.number(units.h)))
      system.x[self.idx_dimers_target, index_compartment] -= delta_target
      system.x[self.idx_antigens_target, index_compartment] += delta_target @ self.Q_target



VIB = {}
VIB.update({"off_C": 10**-4 / units.s, "affn_C": 10 * units.nM, "affm_C": 200 * units.nM})
VIB.update({"off_R": 10**-4 / units.s, "affn_R": 1 * units.nM, "affm_R": 20 * units.nM})
VIB.update({"off_A": 10**-4 / units.s, "affn_A": 10 * units.nM, "affm_A": 200 * units.nM})
VIB.update({"off_B": 10**-4 / units.s, "aff_B": 10 * units.nM})
VIB.update({"avidity_effector": 19, "avidity_target": 19})
VIB.update({"clearance": math.log(2)/(70 * units.h)})
VIB["internalization_effector"] = [("C8", ["C8"], 0.1 / units.h), ("R8", ["R8"], 0.3 / units.h), ("CR8", ["C8", "R8"], 0.1 / units.h), ("C4", ["C4"], 0.1 / units.h), ("R4", ["R4"], 0.3 / units.h), ("CR4", ["C4", "R4"], 0.1 / units.h), ("Rnk", ["Rnk"], 0.3 / units.h)]
VIB["internalization_target"] = [("A", ["A"], 0.02 / units.h), ("B", ["B"], 0.02 / units.h), ("AB", ["A", "B"], 0.02 / units.h)]

VIBX = VIB.copy(); VIBX["smalls"] = []
VIBY = VIB.copy(); VIBY["smalls"] = ["mnn", "nnn"]

VIBX_I = VIBX.copy()
VIBX_I["cleavage_plasma"] = [("m..", "n..", 0.05 / units.d), (".m.", ".n.", 0.05 / units.d), ("..m", "..n", 0.05 / units.d)]
VIBX_I["cleavage_tumor"] = [("m..", "n..", 0.15 / units.d), (".m.", ".n.", 0.15 / units.d), ("..m", "..n", 0.15 / units.d)]

VIBX_II = VIBX.copy()
VIBX_II["cleavage_plasma"] = [("m..", "n..", 0.05 / units.d), (".mm", ".nn", 0.05 / units.d)]
VIBX_II["cleavage_tumor"] = [("m..", "n..", 0.15 / units.d), (".mm", ".nn", 0.15 / units.d)]

VIBY_I = VIBY.copy()
VIBY_I["cleavage_plasma"] = [("m..", "n..", 0.05 / units.d), (".m.", ".n.", 0.05 / units.d), ("..m", "..n", 0.05 / units.d)]
VIBY_I["cleavage_tumor"] = [("m..", "n..", 0.15 / units.d), (".m.", ".n.", 0.15 / units.d), ("..m", "..n", 0.15 / units.d)]

VIBY_II = VIBY.copy()
VIBY_II["cleavage_plasma"] = [("m..", "n..", 0.05 / units.d), (".mm", ".nn", 0.05 / units.d)]
VIBY_II["cleavage_tumor"] = [("m..", "n..", 0.15 / units.d), (".mm", ".nn", 0.15 / units.d)]



############ tumors ############

UT44 = {"name": "tumor"}
UT44.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
UT44.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
UT44.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
UT44.update({"diffusion": 10 * units.micrometer**2 / units.s})
UT44.update({"density_cell": 3e8 * 0.44 / units.ml, "density_8": 3e8 * 0.05 / units.ml, "density_4": 3e8 * 0.1 / units.ml, "density_nk": 3e8 * 0.02 / units.ml})
UT44.update({"num_A": 7e5, "num_B": 1.45e6})

FTC238 = {"name": "tumor"}
FTC238.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
FTC238.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
FTC238.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
FTC238.update({"diffusion": 10 * units.micrometer**2 / units.s})
FTC238.update({"density_cell": 3e8 * 0.44 / units.ml, "density_8": 3e8 * 0.05 / units.ml, "density_4": 3e8 * 0.1 / units.ml, "density_nk": 3e8 * 0.02 / units.ml})
FTC238.update({"num_A": 1e5, "num_B": 1e5})


############ organs ############

plasma = {"name": "plasma"}
plasma.update({"volume": 3126 * units.ml})
plasma.update({"num_8": 7.9E+09 * 0.33, "num_4": 7.9E+09 * 0.67, "num_nk": 1.6E+09})
plasma.update({"conc_A": 5 * units.ug/units.ml / (170 * units.kDa), "conc_B": 0 * units.nM})

lymph = {"name": "lymph"}
lymph.update({"volume": 274 * units.ml})
lymph.update({"num_8": 3.6E+11 * 0.33, "num_4": 3.6E+11 * 0.67, "num_nk": 6.7E+08})
lymph.update({"conc_A": 0 * units.nM, "conc_B": 0 * units.nM})

bone = {"name": "bone"}
bone.update({"volume_plasma": 224 * units.ml, "volume_interstitial": 1891 * units.ml})
bone.update({"plasma_flow": 2591 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
bone.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
bone.update({"num_cell": 4.77E+09 * 0.5, "num_8": 2.1E+10 * 0.33, "num_4": 2.1E+10 * 0.67, "num_nk": 3.3E+09})
bone.update({"num_A": 0, "num_B": 0})

lung = {"name": "lung"}
lung.update({"volume_plasma": 55 * units.ml, "volume_interstitial": 300 * units.ml})
lung.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
lung.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"num_cell": 2.36E+11 * 0.5, "num_8": 1.3E+10 * 0.33, "num_4": 1.3E+10 * 0.67, "num_nk": 7.2E+08})
lung.update({"num_A": 133439, "num_B": 0}) # num_B = 1019 from Liyuan

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"num_cell": 7.2e11 * 0.5, "num_8": 1.8E+10 * 0.33, "num_4": 1.8E+10 * 0.67, "num_nk": 8.1E+08})
SI.update({"num_A": 57075, "num_B": 39649})

other = {"name": "other"}
other.update({"volume_plasma": 1000 * units.ml, "volume_interstitial": 5000 * units.ml})
other.update({"plasma_flow": 100000 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
other.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
other.update({"num_cell": 1e13 * 0.5, "num_8": 1.1E+10 * 0.33, "num_4": 1.1E+10 * 0.67, "num_nk": 3.1E+09})
other.update({"num_A": 10000, "num_B": 0})


############ cells ############

Treg = {"name": "Treg"}
Treg["markers"] = ["P", "R"]
Treg["initials"] = {"P": 30000, "R": 300}
Treg["ligands"] = drugs
Treg["bindings"] = ["P", "R", "RR", "PR", "PRR"]
Treg["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
Treg["death"] = 0.01 / units.d
Treg["alpha"] = True

nTh = {"name": "nTh"}
nTh["markers"] = []
nTh["initials"] = {}
nTh["ligands"] = []
nTh["bindings"] = []
nTh["signals"] = {}
nTh["death"] = 0.002 / units.d
nTh["alpha"] = False

aTh = {"name": "aTh"}
aTh["markers"] = ["P", "R"]
aTh["initials"] = {"P": 30000, "R": 300}
aTh["ligands"] = drugs
aTh["bindings"] = ["P", "R", "RR", "PR", "PRR"]
aTh["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
aTh["death"] = 0.01 / units.d
aTh["alpha"] = False

Th = {"name": "Th"}
Th["markers"] = ["P", "R"]
Th["initials"] = {"P": 30000, "R": 300}
Th["ligands"] = drugs
Th["bindings"] = ["P", "R", "RR", "PR", "PRR"]
Th["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
Th["death"] = 0.01 / units.d
Th["alpha"] = False

nTm = {"name": "nTm"}
nTm["markers"] = []
nTm["initials"] = {}
nTm["ligands"] = []
nTm["bindings"] = []
nTm["signals"] = {}
nTm["death"] = 0.002 / units.d
nTm["alpha"] = False

aTm = {"name": "aTm"}
aTm["markers"] = ["P", "R"]
aTm["initials"] = {"P": 30000, "R": 1500}
aTm["ligands"] = drugs
aTm["bindings"] = ["P", "R", "RR", "PR", "PRR"]
aTm["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
aTm["death"] = 0.01 / units.d
aTm["alpha"] = True

Tm = {"name": "Tm"}
Tm["markers"] = ["P", "R"]
Tm["initials"] = {"P": 30000, "R": 1500}
Tm["ligands"] = drugs
Tm["bindings"] = ["P", "R", "RR", "PR", "PRR"]
Tm["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
Tm["death"] = 0.01 / units.d
Tm["alpha"] = False

Teff = {"name": "Teff"}
Teff["markers"] = ["P", "R"]
Teff["initials"] = {"P": 30000, "R": 1500}
Teff["ligands"] = drugs
Teff["bindings"] = ["P", "R", "RR", "PR", "PRR"]
Teff["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
Teff["death"] = 0.01 / units.d
Teff["alpha"] = True

Tex = {"name": "Tex"}
Tex["markers"] = ["P", "R"]
Tex["initials"] = {"P": 30000, "R": 1500}
Tex["ligands"] = drugs
Tex["bindings"] = ["P", "R", "RR", "PR", "PRR"]
Tex["signals"] = {"PD1":{"P":1, "PR":1, "PRR":1}, "IL2":{"R":1, "RR":2, "PR":1, "PRR":2}}
Tex["death"] = 0.1 / units.d
Tex["alpha"] = True

NK = {"name": "NK"}
NK["markers"] = ["R"]
NK["initials"] = {"R": 3000}
NK["ligands"] = drugs
NK["bindings"] = ["R", "RR"]
NK["signals"] = {"IL2":{"R":1, "RR":2}}
NK["death"] = 0.02 / units.d
NK["alpha"] = True

cells = [Treg, nTh, aTh, Th, nTm, aTm, Tm, Teff, Tex, NK]

############ model ############

def model(TCE, tumors, organs, connect_tumors = True):
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + tumors + organs]
  system = System(drugs, compartments, cells = cells)
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
  system.add_process(transform(compartments = [central["name"] for central in centrals] + [organ["name"] for organ in system.organs], rates = TCE["cleavage_plasma"]))
  system.add_process(transform(compartments = [tumor["name"] for tumor in tumors], rates = TCE["cleavage_tumor"]))
  
  # target binding
  for drug in drugs:
    off_C = TCE["off_C"]; on_C = {"n":TCE["off_C"] / TCE["affn_C"], "m":TCE["off_C"] / TCE["affm_C"]}[drug[0]]
    off_R = TCE["off_R"]; on_R = {"n":TCE["off_R"] / TCE["affn_R"], "m":TCE["off_R"] / TCE["affm_R"]}[drug[1]]
    off_A = TCE["off_A"]; on_A = {"n":TCE["off_A"] / TCE["affn_A"], "m":TCE["off_A"] / TCE["affm_A"]}[drug[2]]
    off_B = TCE["off_B"]; on_B = TCE["off_B"] / TCE["aff_B"]
    avidity_effector = TCE["avidity_effector"]
    avidity_target = TCE["avidity_target"]
    
    system.add_simple("plasma", ["C8", f"{drug}"], [f"C8-{drug}"], on_C, off_C)
    system.add_simple("plasma", ["C4", f"{drug}"], [f"C4-{drug}"], on_C, off_C)
    
    for organ in centrals + tumors + organs:
      system.add_simple(organ["name"], ["C8", f"{drug}"], [f"C8-{drug}"], on_C, off_C)
      system.add_simple(organ["name"], ["R8", f"{drug}"], [f"R8-{drug}"], on_R, off_R)
      system.add_simple(organ["name"], ["C8", f"R8-{drug}"], [f"CR8-{drug}"], on_C * avidity_effector, off_C)
      system.add_simple(organ["name"], ["R8", f"C8-{drug}"], [f"CR8-{drug}"], on_R * avidity_effector, off_R)
      system.add_simple(organ["name"], ["C4", f"{drug}"], [f"C4-{drug}"], on_C, off_C)
      system.add_simple(organ["name"], ["R4", f"{drug}"], [f"R4-{drug}"], on_R, off_R)
      system.add_simple(organ["name"], ["C4", f"R4-{drug}"], [f"CR4-{drug}"], on_C * avidity_effector, off_C)
      system.add_simple(organ["name"], ["R4", f"C4-{drug}"], [f"CR4-{drug}"], on_R * avidity_effector, off_R)
      system.add_simple(organ["name"], ["Rnk", f"{drug}"], [f"Rnk-{drug}"], on_R, off_R)
      
      system.add_simple(organ["name"], [f"{drug}", "A"], [f"{drug}-A"], on_A, off_A)
      system.add_simple(organ["name"], [f"{drug}", "B"], [f"{drug}-B"], on_B, off_B)
      system.add_simple(organ["name"], [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * avidity_target, off_B)
      system.add_simple(organ["name"], [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * avidity_target, off_A)
      
      for target in targets:
        system.add_simple(organ["name"], ["C8", f"{drug}-{target}"], [f"C8-{drug}-{target}"], on_C, off_C)
        system.add_simple(organ["name"], ["R8", f"{drug}-{target}"], [f"R8-{drug}-{target}"], on_R, off_R)
        system.add_simple(organ["name"], ["C8", f"R8-{drug}-{target}"], [f"CR8-{drug}-{target}"], on_C * avidity_effector, off_C)
        system.add_simple(organ["name"], ["R8", f"C8-{drug}-{target}"], [f"CR8-{drug}-{target}"], on_R * avidity_effector, off_R)
        system.add_simple(organ["name"], ["C4", f"{drug}-{target}"], [f"C4-{drug}-{target}"], on_C, off_C)
        system.add_simple(organ["name"], ["R4", f"{drug}-{target}"], [f"R4-{drug}-{target}"], on_R, off_R)
        system.add_simple(organ["name"], ["C4", f"R4-{drug}-{target}"], [f"CR4-{drug}-{target}"], on_C * avidity_effector, off_C)
        system.add_simple(organ["name"], ["R4", f"C4-{drug}-{target}"], [f"CR4-{drug}-{target}"], on_R * avidity_effector, off_R)
        system.add_simple(organ["name"], ["Rnk", f"{drug}-{target}"], [f"Rnk-{drug}-{target}"], on_R, off_R)
      
      for effector in effectors:
        system.add_simple(organ["name"], [f"{effector}-{drug}", "A"], [f"{effector}-{drug}-A"], on_A, off_A)
        system.add_simple(organ["name"], [f"{effector}-{drug}", "B"], [f"{effector}-{drug}-B"], on_B, off_B)
        system.add_simple(organ["name"], [f"{effector}-{drug}-A", "B"], [f"{effector}-{drug}-AB"], on_B * avidity_target, off_B)
        system.add_simple(organ["name"], [f"{effector}-{drug}-B", "A"], [f"{effector}-{drug}-AB"], on_A * avidity_target, off_A)
  
  # internalization
  system.add_process(internalization(compartments = system.compartments, rates_effector = TCE["internalization_effector"], rates_target = TCE["internalization_target"]))
  
  # initial concentrations
  for central in centrals:
    system.add_x("8", central["name"], central["num_8"] / central["volume"] / units.avagadro)
    system.add_x("4", central["name"], central["num_4"] / central["volume"] / units.avagadro)
    system.add_x("nk", central["name"], central["num_nk"] / central["volume"] / units.avagadro)
    
    system.add_x("C8", central["name"], num_C_8 * central["num_8"] / central["volume"] / units.avagadro)
    system.add_x("R8", central["name"], num_R_8 * central["num_8"] / central["volume"] / units.avagadro)
    system.add_x("C4", central["name"], num_C_4 * central["num_4"] / central["volume"] / units.avagadro)
    system.add_x("R4", central["name"], num_R_4 * central["num_4"] / central["volume"] / units.avagadro)
    system.add_x("Rnk", central["name"], num_R_nk * central["num_nk"] / central["volume"] / units.avagadro)
    system.add_x("A", central["name"], central["conc_A"])
    system.add_x("B", central["name"], central["conc_B"])
  
  for tumor in tumors:
    system.add_x("8", tumor["name"], tumor["density_8"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("4", tumor["name"], tumor["density_4"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("nk", tumor["name"], tumor["density_nk"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    
    system.add_x("C8", tumor["name"], num_C_8 * tumor["density_8"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("R8", tumor["name"], num_R_8 * tumor["density_8"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("C4", tumor["name"], num_C_4 * tumor["density_4"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("R4", tumor["name"], num_R_4 * tumor["density_4"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("Rnk", tumor["name"], num_R_nk * tumor["density_nk"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("A", tumor["name"], tumor["num_A"] * tumor["density_cell"] / tumor["volume_interstitial_proportion"] / units.avagadro)
    system.add_x("B", tumor["name"], tumor["num_B"] * tumor["density_cell"] / tumor["volume_interstitial_proportion"] / units.avagadro)
  
  for organ in organs:
    system.add_x("8", organ["name"], organ["num_8"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("4", organ["name"], organ["num_4"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("nk", organ["name"], organ["num_nk"] / organ["volume_interstitial"] / units.avagadro)
    
    system.add_x("C8", organ["name"], num_C_8 * organ["num_8"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("R8", organ["name"], num_R_8 * organ["num_8"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("C4", organ["name"], num_C_4 * organ["num_4"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("R4", organ["name"], num_R_4 * organ["num_4"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("Rnk", organ["name"], num_R_nk * organ["num_nk"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("A", organ["name"], organ["num_A"] * organ["num_cell"] / organ["volume_interstitial"] / units.avagadro)
    system.add_x("B", organ["name"], organ["num_B"] * organ["num_cell"] / organ["volume_interstitial"] / units.avagadro)
  
  # proliferation process
  system.add_process(PD())
  
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
  system.plot(groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")
  
  system.plot(groups = [["8"], ["4"], ["nk"]], linthresh = 1e-6, output = f"{name}_cells.png")
  
  groups = [[analyte for analyte in analytes if re.search("[mn][mn][mn]", analyte)],
            [analyte for analyte in analytes if re.search("[n][mn][mn]", analyte)],
            [analyte for analyte in analytes if re.search("[mn][n][mn]", analyte)],
            [analyte for analyte in analytes if re.search("[mn][mn][n]", analyte)]]
  labels = ["(effector)-xxx-(target)", "(effector)-nxx-(target)", "(effector)-xnx-(target)", "(effector)-xxn-(target)",]
  colors = ["black", "tab:orange", "tab:red", "tab:blue"]
  linestyles = ["solid", "dashed", "dashed", "dashed"]
  system.plot(groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_drugs.png")

  '''
  groups = [["C"], ["R"], ["Rnk"],
            [analyte for analyte in analytes if re.search("C-[mn][mn][mn]$", analyte)],
            [analyte for analyte in analytes if re.search("R-[mn][mn][mn]$", analyte)],
            [analyte for analyte in analytes if re.search("CR-[mn][mn][mn]$", analyte)],
            [analyte for analyte in analytes if re.search("Rnk-[mn][mn][mn]$", analyte)],
            [analyte for analyte in analytes if re.search("C-[mn][mn][mn]-", analyte)],
            [analyte for analyte in analytes if re.search("R-[mn][mn][mn]-", analyte)],
            [analyte for analyte in analytes if re.search("CR-[mn][mn][mn]-", analyte)],
            [analyte for analyte in analytes if re.search("Rnk-[mn][mn][mn]-", analyte)]]
  labels = ["C", "R", "Rnk", "C-xxx", "R-xxx", "CR-xxx", "Rnk-xxx", "C-xxx-target", "R-xxx-target", "CR-xxx-target", "Rnk-xxx-target"]
  colors = ["tab:orange", "tab:blue", "tab:red", 
            "wheat", "skyblue", "tab:green", "pink", "wheat", "skyblue", "tab:green", "pink"]
  linestyles = ["solid", "solid", "solid", "dashed", "dashed", "dashed", "dashed", "solid", "solid", "solid", "solid"]
  system.plot(groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_effectors.png")
  '''
  
  groups = [["A"], ["B"],
            [analyte for analyte in analytes if re.match("[mn][mn][mn]-A$", analyte)],
            [analyte for analyte in analytes if re.match("[mn][mn][mn]-B$", analyte)],
            [analyte for analyte in analytes if re.match("[mn][mn][mn]-AB$", analyte)],
            [analyte for analyte in analytes if re.search("-[mn][mn][mn]-A", analyte)],
            [analyte for analyte in analytes if re.search("-[mn][mn][mn]-B", analyte)],
            [analyte for analyte in analytes if re.search("-[mn][mn][mn]-AB", analyte)]]
  labels = ["A", "B", "xxx-A", "xxx-B", "xxx-AB", "effector-xxx-A", "effector-xxx-B", "effector-xxx-AB"]
  colors = ["tab:blue", "tab:red", 
            "skyblue", "pink", "tab:purple", "skyblue", "pink", "tab:purple"]
  linestyles = ["solid", "solid", "dashed", "dashed", "dashed", "solid", "solid", "solid"]
  system.plot(groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_targets.png")
