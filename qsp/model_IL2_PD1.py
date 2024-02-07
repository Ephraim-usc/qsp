from .qsp import *
import itertools
import sympy

############ tumors ############

UT44 = {"name": "tumor"}
UT44.update({"volume": 170 * units.ul, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
UT44.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
UT44.update({"capillary_radius": 10 * units.um, "capillary_permeability": 3e-7 * units.cm/units.s})
UT44.update({"diffusion": 10 * units.um**2 / units.s})
UT44.update({"density_cell": 3e8 * 0.44 / units.ml, "density_8": 3e8 * 0.05 / units.ml, "density_4": 3e8 * 0.1 / units.ml, "density_nk": 3e8 * 0.02 / units.ml})
UT44.update({"num_A": 7e5, "num_B": 1.45e6})
UT44.update({"enzyme": 15})

FTC238 = {"name": "tumor"}
FTC238.update({"volume": 170 * units.ul, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
FTC238.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
FTC238.update({"capillary_radius": 10 * units.um, "capillary_permeability": 3e-7 * units.cm/units.s})
FTC238.update({"diffusion": 10 * units.um**2 / units.s})
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
bone.update({"vascular_reflection": 0.85, "lymphatic_reflection": 0.2})
bone.update({"num_cell": 4.77E+09 * 0.5, "num_8": 2.1E+10 * 0.33, "num_4": 2.1E+10 * 0.67, "num_nk": 3.3E+09})
bone.update({"num_A": 0, "num_B": 0})
bone.update({"enzyme": 1})

lung = {"name": "lung"}
lung.update({"volume_plasma": 55 * units.ml, "volume_interstitial": 300 * units.ml})
lung.update({"vascular_reflection": 0.95, "lymphatic_reflection": 0.2})
lung.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"num_cell": 2.36E+11 * 0.5, "num_8": 1.3E+10 * 0.33, "num_4": 1.3E+10 * 0.67, "num_nk": 7.2E+08})
lung.update({"num_A": 133439, "num_B": 0}) # num_B = 1019 from Liyuan
lung.update({"enzyme": 1})

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"vascular_reflection": 0.9, "lymphatic_reflection": 0.2})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"num_cell": 7.2e11 * 0.5, "num_8": 1.8E+10 * 0.33, "num_4": 1.8E+10 * 0.67, "num_nk": 8.1E+08})
SI.update({"num_A": 57075, "num_B": 39649})
SI.update({"enzyme": 1})

other = {"name": "other"}
other.update({"volume_plasma": 1000 * units.ml, "volume_interstitial": 5000 * units.ml})
other.update({"plasma_flow": 100000 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
other.update({"vascular_reflection": 0.95, "lymphatic_reflection": 0.2})
other.update({"num_cell": 1e13 * 0.5, "num_8": 1.1E+10 * 0.33, "num_4": 1.1E+10 * 0.67, "num_nk": 3.1E+09})
other.update({"num_A": 10000, "num_B": 0})
other.update({"enzyme": 10})


############ constants ############

molecular_weight = 150000 * units.g/units.mol

def evalf_array(expr, subs_array, n):
  array_of_subs = [{key:value[i] for key, value in subs_array.items()} for i in range(n)]
  return np.array([expr.evalf(subs = subs) for subs in array_of_subs])

def np_safe_divide(a, b):
  return np.divide(a, b, out = np.zeros_like(a), where = b!=0)

def hill(x, EC50, EMAX = 1, coef = 1):
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
        signals[cell.name][SIGNALS_CEL[marker]] = np_safe_divide(sum_dimers, sum_cells)
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
        system.signals_env[SIGNALS_ENV[f"{cell_1.name}_per_{cell_2.name}"]] = np_safe_divide(system.c[idx_cell_1, :], system.c[idx_cell_2, :])
    
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

      deaths = np.maximum(0.0, deaths)
      diffs = np.maximum(0.0, diffs)
      prolifs = np.maximum(0.0, prolifs)
      births = np.maximum(0.0, births)
      
      minus = 1 - np.exp(- (deaths + diffs) * t)
      plus = births * t + system.c[idx_cell, :] * (np.exp(prolifs * t) - 1)
      
      system.x[idx_all_analytes, :] -= system.x[idx_all_analytes, :] * minus
      system.x[idx_marker_analytes, :] += np.outer(cell.initials, plus)
      system.c[idx_cell, :] -= system.c[idx_cell, :] * minus
      system.c[idx_cell, :] += plus
      
      if cell.diff_cell is not None:
        idx_cell_dest = system.cells.index(cell.diff_cell)
        idx_marker_analytes_dest = self.marker_analytes_idxes[cell.diff_cell.name]
        
        plus_dest = system.c[idx_cell, :] * minus * np_safe_divide(diffs, diffs + deaths) * cell.diff_copy
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
    
    t = t.number(units.h)
    for idxeses, Q in zip(self.idxeseses, self.Qs): # for each ligand
      Qs = evalf_array(Q, system.signals_env, system.n_compartments).astype(float) # Q matrices for each compartment
      for idxes in idxeses: # for each class of analytes
        for idx_compartment in range(system.n_compartments):
          system.x[idxes, idx_compartment] = system.x[idxes, idx_compartment] @ expm(Qs[idx_compartment] * t)




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
  
  def get_binding_reactions(self, cell):
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
X.add_cleavage(1, SIGNALS_ENV["enzyme"] * 0.01/units.d)
X.add_cleavage(2, SIGNALS_ENV["enzyme"] * 0.01/units.d)
X.set_avidity(0, 1, 20)
X.set_avidity(0, 2, 20)
X.set_avidity(1, 2, 20)

IL2 = Ligand(name = "IL2", 
             affs_offs = [{"α": (0.01 * units.nM, 2e-4 / units.s), "R": (1 * units.nM, 1e-4 / units.s)}],
             masks = [None],
             clearance = math.log(2)/(70 * units.h), 
             smalls = ["n"])


class Cell:
  def __init__(self, name, markers, initials, ints, int_mode = "geomean", birth = None, death = None, prolif = None, diff = None, diff_copy = 1, diff_cell = None):
    self.name = name
    self.markers = markers
    self.initials = initials
    self.ints = [internalization.number(1/units.h) for internalization in ints]
    self.int_mode = int_mode
    
    self.birth = birth.number(units.nM/units.h) if birth is not None else 0.0
    self.death = death.number(1/units.h) if death is not None else 0.0
    self.prolif = prolif.number(1/units.h) if prolif is not None else 0.0
    self.diff = diff.number(1/units.h) if diff is not None else 0.0
    self.diff_copy = diff_copy
    self.diff_cell = diff_cell
  
  def get_bindings(self, ligand): # find all binding modes formed when a ligand binds to a cell
    site_bindings = [["_"] + [target for target in targets if target in self.markers] for targets in ligand.targets]
    bindings = ["".join(complex) for complex in itertools.product(*site_bindings)][1:] # remove the all non-binding mode
    return bindings
  
  def get_all_analytes(self, ligands):
    buffer = [f"{self.name}:{marker}" for marker in self.markers]
    for ligand in ligands:
      buffer += [f"{self.name}:{binding}-{ligand.name}:{state}" for binding in self.get_bindings(ligand) for state in ligand.states]
    return buffer
  
  def get_cellular_signal_dimers(self, marker, ligands): # dimers with multiple markers bound are repeated
    buffer = []
    for ligand in ligands:
      bindings = self.get_bindings(ligand); counts = [binding.count(marker) for binding in bindings]
      bindings_with_repeats = [binding for binding, count in zip(bindings, counts) for i in range(count)]
      buffer += [f"{self.name}:{binding}-{ligand.name}:{state}" for binding in bindings_with_repeats for state in ligand.states]
    return buffer
  
  def get_internalization_reactions(self, ligand):
    ints = {marker:internalization for marker, internalization in zip(self.markers, self.ints)}
    bindings = self.get_bindings(ligand)
    productses = [[f"{self.name}:{marker}" for marker in binding if marker != "_"] for binding in bindings]
    rateses = [[ints[marker] for marker in binding if marker != "_"] for binding in bindings]
    if self.int_mode == "geomean":
      rates = [np.array(tmp).prod()**(1.0/len(tmp)) for tmp in rateses]
    
    reactions = []
    for binding, products, rate in zip(bindings, productses, rates):
      for state in ligand.states:
        reactions.append(([f"{self.name}:{binding}-{ligand.name}:{state}"], products, rate/units.h))
    return reactions
    
    
# currently assuming T cells not proliferating in lymph, and equivalently shrink the dosage by ~10 times.
# PD1/PDL1的T细胞失活效果的EC50是每细胞1000个，Hill coefficient是2
# PD1/PDL1的抑制T细胞杀伤效果的EC50是每细胞250个，Hill coefficient是2
# PD1/PDL1的亲和力是7.2uM=7200nM，而肿瘤细胞表达160万PDL1（对应900nM），因此因此PD1占据率大约为10%
tumor_cell_total_density = 3e8 / units.ml / units.avagadro
TREG_RATIO = 0.1
Treg = Cell("Treg", ["P", "α"], [30000, 300], [0.05/units.h, 2.0/units.h],
            birth = SIGNALS_ENV["tumor"] * 0.01 / units.d * tumor_cell_total_density * 0.1 * TREG_RATIO,
            death = SIGNALS_ENV["tumor"] * 0.01 / units.d,
            prolif = SIGNALS_ENV["tumor"] * 0.5 / units.d * (0.05 - hill(30000 - SIGNALS_CEL["P"], 10000, EMAX = 0.1) + hill(SIGNALS_CEL["α"], 100, coef = 1, EMAX = 0.5)))
Th = Cell("Th", ["P", "R"], [30000, 300], [0.05/units.h, 2.0/units.h],
          birth = SIGNALS_ENV["tumor"] * 0.01 / units.d * tumor_cell_total_density * 0.1 * (1-TREG_RATIO),
          death = SIGNALS_ENV["tumor"] * 0.01 / units.d,
          prolif = SIGNALS_ENV["tumor"] * 0.5 / units.d * (0.05 - hill(30000 - SIGNALS_CEL["P"], 10000, EMAX = 0.1) + hill(SIGNALS_CEL["R"], 100, coef = 1, EMAX = 0.5)),
          diff = SIGNALS_ENV["tumor"] * 0.1 / units.d * hill(SIGNALS_CEL["R"], EC50 = 100, coef = 1.0),
          diff_cell = Treg)
Tm = Cell("Tm", ["P", "R"], [30000, 1500], [0.05/units.h, 2.0/units.h])
Tex = Cell("Tex", ["P", "α"], [60000, 1500], [0.05/units.h, 2.0/units.h],
           death = 0.1 / units.d)
Teff = Cell("Teff", ["P", "α"], [60000, 1500], [0.05/units.h, 2.0/units.h],
            birth = SIGNALS_ENV["tumor"] * 0.01 / units.d * tumor_cell_total_density * 0.05,
            death = SIGNALS_ENV["tumor"] * 0.01 / units.d,
            prolif = SIGNALS_ENV["tumor"] * 1.386 / units.d * (0.05 - hill(60000 - SIGNALS_CEL["P"], 10000, EMAX = 0.1) + hill(SIGNALS_CEL["α"], 100, coef = 3.1, EMAX = 0.5)),
            diff = SIGNALS_ENV["tumor"] * 0.5 / units.d * (hill(SIGNALS_ENV["Treg_per_Teff"], 1, EMAX = 0.05) + hill(60000 - SIGNALS_CEL["P"], 10000, EMAX = 0.1) + hill(SIGNALS_CEL["α"], 100, coef = 3.1, EMAX = 0.5)),
            diff_cell = Tex)
NK = Cell("NK", ["α"], [3000], [2.0/units.h],
          birth = SIGNALS_ENV["tumor"] * 0.01 / units.d * tumor_cell_total_density * 0.02,
          prolif = SIGNALS_ENV["tumor"] * 1.512 / units.d * hill(SIGNALS_CEL["α"], 100, coef = 1.3, EMAX = 0.5),
          death = SIGNALS_ENV["tumor"] * 0.01 / units.d)


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
  analytes_markers = [f"{cell.name}:{marker}" for cell in cells for marker in cell.markers]
  analytes_ligands = [f"{ligand.name}:{state}" for ligand in ligands for state in ligand.states]
  analytes_dimers = [f"{cell.name}:{binding}-{ligand.name}:{state}" for ligand in ligands for cell in cells for binding in ligand.get_bindings(cell) for state in ligand.states]
  analytes = analytes_ligands + analytes_markers + analytes_dimers
  
  centrals = [plasma, lymph]
  compartments = [organ["name"] for organ in centrals + tumors + organs]
  system = System(analytes, compartments, cells = cells)
  
  system.cells = cells
  system.centrals = [plasma, lymph]
  system.tumors = tumors
  system.organs = organs
  system.signals_env = {}
  system.signals_env[SIGNALS_ENV["tumor"]] = [1 if organ["name"] in [tumor["name"] for tumor in tumors] else 0 for organ in centrals + tumors + organs]
  system.signals_env[SIGNALS_ENV["enzyme"]] = [organ["enzyme"] for organ in centrals + tumors + organs]
  
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
  if len(tumors) > 1 and connect_tumors:
    system.add_process(equilibrium([tumor["name"] for tumor in tumors], analytes_ligands))
  
  # mask cleavage
  system.add_process(process_transform(cells, ligands))
  
  # binding reactions
  for compartment in compartments:
    for cell in cells:
      for ligand in ligands:
        for reactants, products, aff, off in ligand.get_binding_reactions(cell):
          system.add_simple(compartment, reactants, products, off/aff, off)
  
  # internalization
  for compartment in compartments:
    for cell in cells:
      for ligand in ligands:
        for reactants, products, rate in cell.get_internalization_reactions(ligand):
          system.add_simple(compartment, reactants, products, rate)
  
  # add cells
  for central in centrals:
    add_cell(system, Treg, central["name"], central["num_4"] * TREG_RATIO / central["volume"] / units.avagadro)
    add_cell(system, Th, central["name"], central["num_4"] * (1-TREG_RATIO) / central["volume"] / units.avagadro)
    add_cell(system, Tm, central["name"], central["num_8"] / central["volume"] / units.avagadro)
    add_cell(system, NK, central["name"], central["num_nk"] / central["volume"] / units.avagadro)
  
  for organ in organs:
    add_cell(system, Treg, organ["name"], organ["num_4"] * TREG_RATIO / central["volume"] / units.avagadro)
    add_cell(system, Th, organ["name"], organ["num_4"] * (1-TREG_RATIO) / central["volume"] / units.avagadro)
    add_cell(system, Tm, organ["name"], organ["num_8"] / central["volume"] / units.avagadro)
    add_cell(system, NK, organ["name"], organ["num_nk"] / central["volume"] / units.avagadro)
  
  for tumor in tumors:
    add_cell(system, Treg, tumor["name"], tumor["density_4"] * TREG_RATIO / units.avagadro)
    add_cell(system, Th, tumor["name"], tumor["density_4"] * (1-TREG_RATIO) / units.avagadro)
    add_cell(system, Teff, tumor["name"], tumor["density_8"] / units.avagadro)
    add_cell(system, NK, tumor["name"], tumor["density_nk"] / units.avagadro)
  
  # process of cell dynamics
  system.add_process(process_compute_cellular_signals(cells, ligands))
  system.add_process(process_cell_dynamics(cells, ligands))

  return system
