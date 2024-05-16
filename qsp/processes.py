class Equilibrium:
  def __init__(self):
    self.system = None
    self.data = dict()
  
  def add(analyte, compartments):
    if analyte not in self.data:
      self.data[analyte] = []
    for compartment in compartments:
      if compartment not in self.data[analyte]:
        self.data[analyte].append(compartment)
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      self.data_ = dict()
      for analyte, compartments in self.data.items():
        if analyte in system.analytes:
          analyte_ = system.analytes.index(analyte)
        else:
          continue
        self.data_[analyte_] = []
        for compartment in compartments:
          if compartment in system.compartments:
            compartment_ = system.compartments.index(compartment)
          else:
            continue
          self.data_[analyte_].append(compartment_)
    
    for analyte_, compartments_ in self.data_.items():
      x = system.x[analyte_, compartments_]
      volumes = system.V[analyte_, compartments_]
      system.x[analyte_, compartments_] = np.average(x, weights = volumes)


class Transform:
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
