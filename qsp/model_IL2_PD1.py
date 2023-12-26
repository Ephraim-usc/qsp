from .qsp import *
import itertools
import sympy

############ constants ############

molecular_weight = 150000 * units.g/units.mol

def evalf_array(expr, subs_array, n):
  array_of_subs = [{key:value[i] for key, value in subs_array.items()} for i in range(n)]
  return np.array([expr.evalf(subs = subs) for subs in array_of_subs])

def hill(x, EMAX, EC50, coef):
  return EMAX * (x**coef) / (x**coef + EC50**coef)

def add_cell(system, cell, compartment, value):
  value = value.number(units.nM)
  idx_cell = system.cells.index(cell)
  idx_compartment = system.compartments.index(compartment)
  idx_markers = [system.analytes.index(f"{cell.name}:{marker}") for marker in cell.markers]
  system.c[idx_cell, idx_compartment] += value
  system.x[idx_markers, idx_compartment] += np.array(cell.initials) * value
  

class process_compute_cellular_signals:
  def __init__(self, cells, ligands):
    self.system = None
    self.cells = cells
    self.ligands = ligands
    
    self.signal_dimers_analytes = {}
    for cell in cells:
      self.signal_dimers_analytes[cell.name] = {}
      for marker in cell.markers:
        self.signal_dimers_analytes[cell.name][marker] = cell.get_cellular_signal_dimers(marker, ligands)
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      self.signal_dimers_idxes = {}
      for cell in self.cells:
        self.signal_dimers_idxes[cell.name] = {}
        for marker in cell.markers:
          self.signal_dimers_idxes[cell.name][marker] = [system.analytes.index(analyte) for analyte in self.signal_dimers_analytes[cell.name][marker]]
    
    signals = {}
    for cell in self.cells:
      idx_cell = system.cells.index(cell)
      signals[cell.name] = {}
      for marker in cell.markers:
        idx_signal_dimers = self.signal_dimers_idxes[cell.name][marker]
        sum_dimers = system.x[idx_signal_dimers, :].sum(axis = 0)
        sum_cells = system.c[idx_cell, :]
        signals[cell.name][SIGNALS_CEL[marker]] = np.divide(sum_dimers, sum_cells, out = np.zeros_like(sum_dimers), where = sum_cells!=0)
    system.signals_cellular = signals


class process_cell_dynamics:
  def __init__(self, cells, ligands):
    self.system = None
    self.cells = cells
    
    self.all_analytes = {}
    self.marker_analytes = {}
    for cell in cells:
      self.all_analytes[cell.name] = cell.get_all_analytes(ligands)
      self.marker_analytes[cell.name] = [f"{cell.name}:{marker}" for marker in cell.markers]
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      self.all_analytes_idxes = {}
      self.marker_analytes_idxes = {}
      for cell in self.cells:
        self.all_analytes_idxes[cell.name] = [system.analytes.index(analyte) for analyte in self.all_analytes[cell.name]]
        self.marker_analytes_idxes[cell.name] = [system.analytes.index(analyte) for analyte in self.marker_analytes[cell.name]]
    
    # computing pair-wise cell ratios
    for cell_1 in self.cells:
      idx_cell_1 = system.cells.index(cell_1)
      for cell_2 in self.cells:
        idx_cell_2 = system.cells.index(cell_2)
        sum_1 = system.c[idx_cell_1, :]
        sum_2 = system.c[idx_cell_2, :]
        system.signals_env[SIGNALS_ENV[f"{cell_1.name}_per_{cell_2.name}"]] = np.divide(sum_1, sum_2, out = np.zeros_like(sum_1), where = sum_2!=0)
    
    t = t.number(units.h)
    for cell in self.cells:
      idx_cell = system.cells.index(cell)
      idx_all_analytes = self.all_analytes_idxes[cell.name]
      idx_marker_analytes = self.marker_analytes_idxes[cell.name]
      
      signals = {**system.signals_cellular[cell.name], **system.signals_env}
      deaths = evalf_array(cell.death, signals, system.n_compartments).astype(float) if isinstance(cell.death, sympy.Expr) else cell.death
      diffs = evalf_array(cell.diff, signals, system.n_compartments).astype(float) if isinstance(cell.diff, sympy.Expr) else cell.diff
      prolifs = evalf_array(cell.prolif, signals, system.n_compartments).astype(float) if isinstance(cell.prolif, sympy.Expr) else cell.prolif
      births = evalf_array(cell.birth, signals, system.n_compartments).astype(float) if isinstance(cell.birth, sympy.Expr) else cell.birth
      minus = 1 - np.exp(- (deaths + diffs) * t)
      plus = births * t + system.c[idx_cell, :] * (np.exp(prolifs * t) - 1)
      
      system.x[idx_all_analytes, :] -= system.x[idx_all_analytes, :] * minus
      system.x[idx_marker_analytes, :] += np.outer(cell.initials, plus)
      system.c[idx_cell, :] -= system.c[idx_cell, :] * minus
      system.c[idx_cell, :] += plus
      
      if cell.diff_cell is not None:
        idx_cell_dest = system.cells.index(cell.diff_cell)
        idx_marker_analytes_dest = self.marker_analytes_idxes[cell.diff_cell.name]
        
        plus_dest = minus * diffs / (diffs + deaths) * cell.diff_copy
        system.x[idx_marker_analytes_dest, :] += np.outer(cell.diff_cell.initials, plus_dest)
        system.c[idx_cell_dest, :] += plus_dest
      


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
        Q = np.array(Q.evalf(subs = system.signals[idx_compartment]), dtype = np.float64)
        for idxes in idxeses:
          system.x[idxes, idx_compartment] = system.x[idxes, idx_compartment] @ expm(Q * t.number(units.h))







class Signal:
  def __init__(self, category):
    self.category = category
    self.dict = {}
  
  def __getitem__(self, key):
    if key not in self.dict:
      self.dict[key] = sympy.Symbol(f"{self.category}_{key}")
    return self.dict[key]

SIGNALS_ENV = Signal("env")
SIGNALS_CEL = Signal("cel")



class Ligand:
  def __init__(self, name, affs_offs, masks, clearance, smalls = None): # currently only consider mask/naked states
    self.name = name
    
    self.n_sites = len(affs_offs)
    self.affs_offs = affs_offs
    self.targets = [list(tmp.keys()) for tmp in affs_offs]
    
    self.masks = masks
    self.site_states = [["n"] if mask is None else ["m", "n"] for mask in masks]
    self.states = ["".join(tmp) for tmp in itertools.product(*self.site_states)]
    self.n_states = len(self.states)
    
    self.clearance = clearance.number(1/units.h)
    
    self.V = np.ones([self.n_sites, self.n_sites]) # avidity
    self.Q = sympy.zeros(self.n_states, self.n_states) # transforming matrix, may contain signal variables
    self.smalls = smalls if smalls is not None else [] # states that are small molecules
  
  def set_avidity(self, site_A, site_B, avidity, symmetric = True):
    self.V[site_A, site_B] = avidity
    if symmetric:
      self.V[site_B, site_A] = avidity
  
  def add_cleavage(self, sites, rate):
    if type(sites) is not list:
      sites = [sites]
    
    site_states_from = self.site_states.copy()
    site_states_to = self.site_states.copy()
    for site in sites:
      site_states_from[site] = ["m"]
      site_states_to[site] = ["n"]
    states_from = ["".join(tmp) for tmp in itertools.product(*site_states_from)]
    states_to = ["".join(tmp) for tmp in itertools.product(*site_states_to)]
    
    idxes_from = [self.states.index(tmp) for tmp in states_from]
    idxes_to = [self.states.index(tmp) for tmp in states_to]
    for idx_from, idx_to in zip(idxes_from, idxes_to):
      self.Q[idx_from, idx_from] -= rate.number(1/units.h)
      self.Q[idx_from, idx_to] += rate.number(1/units.h)
  
  def get_bindings(self, cell): # find all binding modes formed when a ligand binds to a cell
    site_bindings = [["_"] + [target for target in targets if target in cell.markers] for targets in self.targets]
    bindings = ["".join(complex) for complex in itertools.product(*site_bindings)][1:] # remove the all non-binding mode
    return bindings
  
  def get_reactions(self, cell):
    site_bindings = [["_"] + [target for target in targets if target in cell.markers] for targets in self.targets]
    
    reactions = []
    for state in self.states:
      for site in range(self.n_sites):
        site_bindings_from = site_bindings.copy()
        site_bindings_from[site] = ["_"]
        bindings_from = ["".join(complex) for complex in itertools.product(*site_bindings_from)]
        for marker in site_bindings[site][1:]:
          site_bindings_to = site_bindings.copy()
          site_bindings_to[site] = [marker]
          bindings_to = ["".join(complex) for complex in itertools.product(*site_bindings_to)]
          for binding_from, binding_to in zip(bindings_from, bindings_to):
            if binding_from == len(binding_from) * "_":
              analyte_from = f"{self.name}:{state}"
            else:
              analyte_from = f"{cell.name}:{binding_from}-{self.name}:{state}"
            analyte_to = f"{cell.name}:{binding_to}-{self.name}:{state}"
            
            aff, off = self.affs_offs[site][marker]
            if state[site] == "m":
              aff *= self.masks[site]
            aff /= max([(self.V[i, site] if binding_from[i] != "_" else 1) for i in range(self.n_sites)]) # take the maximum avidity from currently bound sites
            
            reactions.append(([analyte_from, f"{cell.name}:{marker}"], [analyte_to], aff, off))
    return reactions
  
  def print(self):
    Q = pd.DataFrame(np.array(self.Q), index = self.states, columns = self.states)
    print(Q)

X = Ligand(name = "X", 
           affs_offs = [{"P": (1 * units.nM, 1e-4 / units.s)}, \
                      {"α": (0.01 * units.nM, 2e-4 / units.s), "R": (1 * units.nM, 1e-4 / units.s)}, \
                      {"α": (0.01 * units.nM, 2e-4 / units.s), "R": (1 * units.nM, 1e-4 / units.s)}],
           masks = [None, 20, 20],
           clearance = math.log(2)/(70 * units.h))
X.add_cleavage(1, signals["enzyme"] * 0.01/units.h)
X.add_cleavage(2, signals["enzyme"] * 0.01/units.h)
X.set_avidity(0, 1, 20)
X.set_avidity(0, 2, 20)
X.set_avidity(1, 2, 20)

IL2 = Ligand(name = "IL2", 
             affs_offs = [{"α": (0.01 * units.nM, 2e-4 / units.s), "R": (1 * units.nM, 1e-4 / units.s)}],
             masks = [None],
             clearance = math.log(2)/(70 * units.h), 
             smalls = ["n"])


class Cell:
  def __init__(self, name, markers, initials, birth = None, death = None, prolif = None, diff = None, diff_copy = 1, diff_cell = None):
    self.name = name
    self.markers = markers
    self.initials = initials
    self.birth = birth.number(units.nM/units.h) if birth is not None else 0.0
    self.death = death.number(1/units.h) if death is not None else 0.0
    self.prolif = prolif.number(1/units.h) if prolif is not None else 0.0
    self.diff = diff.number(1/units.h) if diff is not None else 0.0
    self.diff_copy = diff_copy
    self.diff_cell = diff_cell
  
  def get_all_analytes(self, ligands):
    buffer = [f"{self.name}:{marker}" for marker in self.markers]
    buffer += self.get_all_dimers(ligands)
    return buffer
  
  def get_all_dimers(self, ligands):
    buffer = []
    for ligand in ligands:
      buffer += [f"{self.name}:{binding}-{ligand.name}:{state}" for binding in ligand.get_bindings(self) for state in ligand.states]
    return buffer
  
  def get_cellular_signal_dimers(self, marker, ligands): # dimers with multiple markers bound are repeated
    buffer = []
    for ligand in ligands:
      bindings = ligand.get_bindings(self); counts = [binding.count(marker) for binding in bindings]
      bindings_with_repeats = [binding for binding, count in zip(bindings, counts) for i in range(count)]
      buffer += [f"{self.name}:{binding}-{ligand.name}:{state}" for binding in bindings_with_repeats for state in ligand.states]
    return buffer

tumor_cell_total_density = 3e8 / units.ml / units.avagadro
TREG_RATIO = 0.1
Treg = Cell("Treg", ["P", "α"], [30000, 300], 
            death = SIGNALS_ENV["tumor"] * 0.01 / units.d, 
            birth = SIGNALS_ENV["tumor"] * tumor_cell_total_density * 0.1 * TREG_RATIO * (0.01 / units.d))
Th = Cell("Th", ["P", "R"], [30000, 300], 
          death = SIGNALS_ENV["tumor"] * 0.01 / units.d, 
          birth = SIGNALS_ENV["tumor"] * tumor_cell_total_density * 0.1 * (1-TREG_RATIO) * (0.01 / units.d))
Tm = Cell("Tm", ["P", "R"], [30000, 1500], 
          birth = SIGNALS_ENV["tumor"] * 1 * units.nM/units.d, 
          death = 0.01 / units.d)
Teff = Cell("Teff", ["P", "α"], [30000, 1500],
            birth = SIGNALS_ENV["tumor"] * 0.01 / units.d * tumor_cell_total_density * 0.05,
            death = SIGNALS_ENV["tumor"] * 0.01 / units.d,
            prolif = SIGNALS_ENV["tumor"] * 0.01 / units.d * (SIGNALS_CEL["α"] - SIGNALS_CEL["P"]),
            diff = SIGNALS_ENV["tumor"] * 0.1 / units.d * (SIGNALS_ENV["Treg_per_Teff"] + SIGNALS_CEL["P"] + SIGNALS_CEL["α"]),
            diff_cell = Tex)
Tex = Cell("Tex", ["P", "α"], [30000, 1500],
           death = 0.1 / units.d)
NK = Cell("NK", ["α"], [3000],
          birth = SIGNALS_ENV["tumor"] * 0.02 / units.d * tumor_cell_total_density * 0.02,
          death = SIGNALS_ENV["tumor"] * 0.02 / units.d)


'''
Treg = Cell("Treg", ["P", "α"], [30000, 300], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
nTh = Cell("nTh", [], [], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.002 / units.d)
aTh = Cell("aTh", ["P", "R"], [30000, 300], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Th = Cell("Th", ["P", "R"], [30000, 300], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
nTm = Cell("nTm", [], [], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.002 / units.d)
aTm = Cell("aTm", ["P", "α"], [30000, 1500], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Tm = Cell("Tm", ["P", "R"], [30000, 1500], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.01 / units.d)
Teff = Cell("Teff", ["P", "α"], [30000, 1500], birth = signals["tumor"] * 1 * units.nM/units.d, prolif = 0.01 * signals["α"] / units.d, death = 0.01 / units.d)
Tex = Cell("Tex", ["P", "α"], [30000, 1500], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.1 / units.d)
NK = Cell("NK", ["α"], [30000, 3000], birth = signals["tumor"] * 1 * units.nM/units.d, death = 0.02 / units.d)
'''




def PDSystem(organs, tumors, cells, ligands):

organs = [bone, lung]
tumors = [UT44]
cells = [Treg, Th, Teff, Tex, NK]
ligands = [X, IL2]

analytes_markers = [f"{cell.name}:{marker}" for cell in cells for marker in cell.markers]
analytes_ligands = [f"{ligand.name}:{state}" for ligand in ligands for state in ligand.states]
analytes_dimers = [f"{cell.name}:{binding}-{ligand.name}:{state}" for ligand in ligands for cell in cells for binding in ligand.get_bindings(cell) for state in ligand.states]
analytes = analytes_ligands + analytes_markers + analytes_dimers

centrals = [plasma, lymph]
compartments = [organ["name"] for organ in centrals + tumors + organs]
system = System(analytes, compartments)

system.cells = cells
system.centrals = [plasma, lymph]
system.tumors = tumors
system.organs = organs
system.signals_env = {}
system.signals_env[SIGNALS_ENV["tumor"]] = [1 if organ["name"] in [tumor["name"] for tumor in tumors] else 0 for organ in centrals + tumors + organs]
system.signals_env[SIGNALS_ENV["enzyme"]] = [organ["enzyme"] for organ in centrals + tumors + organs]
system.c = np.zeros([len(cells), system.n_compartments], dtype = float) # in units.nM


for analyte in analytes:
  for central in centrals:
    system.set_volume(analyte, central["name"], central["volume"])
  for tumor in tumors:
    system.set_volume(analyte, tumor["name"], tumor["volume"] * tumor["volume_interstitial_proportion"])
  for organ in organs:
    system.set_volume(analyte, organ["name"], organ["volume_interstitial"])

# whole-body clearance
for compartment in compartments:
  for ligand in ligands:
    for state in ligand.states:
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
system.add_process(process_transform(cells, ligands))

# binding reactions
for compartment in compartments:
  for cell in cells:
    for ligand in ligands:
      for reactants, products, aff, off in ligand.get_reactions(cell):
        system.add_simple(compartment, reactants, products, off/aff, off)


# add cells
for central in centrals:
  add_cell(system, Treg, central["name"], central["num_4"] * TREG_RATIO / central["volume"] / units.avagadro)
  add_cell(system, Th, central["name"], central["num_4"] * (1-TREG_RATIO) / central["volume"] / units.avagadro)
  add_cell(system, Teff, central["name"], central["num_8"] / central["volume"] / units.avagadro)
  add_cell(system, NK, central["name"], central["num_nk"] / central["volume"] / units.avagadro)

for organ in organs:
  add_cell(system, Treg, organ["name"], organ["num_4"] * TREG_RATIO / central["volume"] / units.avagadro)
  add_cell(system, Th, organ["name"], organ["num_4"] * (1-TREG_RATIO) / central["volume"] / units.avagadro)
  add_cell(system, Teff, organ["name"], organ["num_8"] / central["volume"] / units.avagadro)
  add_cell(system, NK, organ["name"], organ["num_nk"] / central["volume"] / units.avagadro)

for tumor in tumors:
  add_cell(system, Treg, tumor["name"], tumor["density_4"] * TREG_RATIO / units.avagadro)
  add_cell(system, Th, tumor["name"], tumor["density_4"] * (1-TREG_RATIO) / units.avagadro)
  add_cell(system, Teff, tumor["name"], tumor["density_8"] / units.avagadro)
  add_cell(system, NK, tumor["name"], tumor["density_nk"] / units.avagadro)






    

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




############ tumors ############

UT44 = {"name": "tumor"}
UT44.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
UT44.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
UT44.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
UT44.update({"diffusion": 10 * units.micrometer**2 / units.s})
UT44.update({"density_cell": 3e8 * 0.44 / units.ml, "density_8": 3e8 * 0.05 / units.ml, "density_4": 3e8 * 0.1 / units.ml, "density_nk": 3e8 * 0.02 / units.ml})
UT44.update({"num_A": 7e5, "num_B": 1.45e6})
UT44.update({"enzyme": 15})

FTC238 = {"name": "tumor"}
FTC238.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
FTC238.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
FTC238.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
FTC238.update({"diffusion": 10 * units.micrometer**2 / units.s})
FTC238.update({"density_cell": 3e8 * 0.44 / units.ml, "density_8": 3e8 * 0.05 / units.ml, "density_4": 3e8 * 0.1 / units.ml, "density_nk": 3e8 * 0.02 / units.ml})
FTC238.update({"num_A": 1e5, "num_B": 1e5})
FTC238.update({"enzyme": 15})

############ organs ############

plasma = {"name": "plasma"}
plasma.update({"volume": 3126 * units.ml})
plasma.update({"num_8": 7.9E+09 * 0.33, "num_4": 7.9E+09 * 0.67, "num_nk": 1.6E+09})
plasma.update({"conc_A": 5 * units.ug/units.ml / (170 * units.kDa), "conc_B": 0 * units.nM})
plasma.update({"enzyme": 1})

lymph = {"name": "lymph"}
lymph.update({"volume": 274 * units.ml})
lymph.update({"num_8": 3.6E+11 * 0.33, "num_4": 3.6E+11 * 0.67, "num_nk": 6.7E+08})
lymph.update({"conc_A": 0 * units.nM, "conc_B": 0 * units.nM})
lymph.update({"enzyme": 1})

bone = {"name": "bone"}
bone.update({"volume_plasma": 224 * units.ml, "volume_interstitial": 1891 * units.ml})
bone.update({"plasma_flow": 2591 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
bone.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
bone.update({"num_cell": 4.77E+09 * 0.5, "num_8": 2.1E+10 * 0.33, "num_4": 2.1E+10 * 0.67, "num_nk": 3.3E+09})
bone.update({"num_A": 0, "num_B": 0})
bone.update({"enzyme": 1})

lung = {"name": "lung"}
lung.update({"volume_plasma": 55 * units.ml, "volume_interstitial": 300 * units.ml})
lung.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
lung.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"num_cell": 2.36E+11 * 0.5, "num_8": 1.3E+10 * 0.33, "num_4": 1.3E+10 * 0.67, "num_nk": 7.2E+08})
lung.update({"num_A": 133439, "num_B": 0}) # num_B = 1019 from Liyuan
lung.update({"enzyme": 1})

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"num_cell": 7.2e11 * 0.5, "num_8": 1.8E+10 * 0.33, "num_4": 1.8E+10 * 0.67, "num_nk": 8.1E+08})
SI.update({"num_A": 57075, "num_B": 39649})
SI.update({"enzyme": 1})

other = {"name": "other"}
other.update({"volume_plasma": 1000 * units.ml, "volume_interstitial": 5000 * units.ml})
other.update({"plasma_flow": 100000 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
other.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
other.update({"num_cell": 1e13 * 0.5, "num_8": 1.1E+10 * 0.33, "num_4": 1.1E+10 * 0.67, "num_nk": 3.1E+09})
other.update({"num_A": 10000, "num_B": 0})
other.update({"enzyme": 10})


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
