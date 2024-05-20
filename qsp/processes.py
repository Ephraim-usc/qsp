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
    self.data = dict()
  
  def add(self, compartment, reactant, products, rate):
    pass


def add_cleavage(system, linker, drug_source, drug_dests, bindings):
  for compartment, rate in linker:
    if compartment not in system.compartments:
      continue
    system.add_transform(compartment, drug_source, drug_dests, rate)
    for binding in bindings:
      analyte_source = f"{binding}-{drug_source}"
      analyte_dests = [f"{binding}-{drug_dest}" for drug_dest in drug_dests]
      system.add_transform(compartment, analyte_source, analyte_dests, rate)

def add_internalization(system, binding_source, drug_sources, analyte_dests, rate):
  for compartment in system.compartments:
    for drug_source in drug_sources:
      analyte_source = f"{binding_source}-{drug_source}"
      system.add_transform(compartment, analyte_source, analyte_dests, rate)
