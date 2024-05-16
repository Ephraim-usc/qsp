from .qsp import *
import itertools

### this model is mostly from ...

cells = ["T"]
areas = [200]

antigens = ["[T]P"]
bindings = ["[T]P"]
drugs = ["n"]
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


