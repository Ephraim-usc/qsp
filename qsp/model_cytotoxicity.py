from .qsp import *

class Drug:
  def __init__(self, num_binding_sites):
    self.num_binding_sites = num_binding_sites
  
  def add_affinity(self, ligand, aff, off):
    pass
  
  def add_cis(self, ligand1, ligand2, on2D):
    pass
  
  def add_trans(self, ligand1, ligand2, on2D):
    pass

X = Drug(3)
X.add_affinity(0, "CD3", 1*units.nM, 1e-4/units.s)
X.add_affinity(1, "EGFR", 10*units.nM, 1e-4/units.s)
X.add_affinity(2, "CA9", 20*units.nM, 1e-4/units.s)
X.add_cis(1, 2, )





def Cytotoxicity(system, cell_effector, cell_target, params):
  pass





class Cell:
  def __init__(self, radius = None, area = None, markers = None, copies = None):
    self.radius = radius
    if area is None:
      self.area = math.pi * radius**2
    
    self.markers = markers
    self.copies = copies


Tcell = Cell(radius = 5*units.um**2, markers = ["CD3"], copies = [50000])
HT29 = Cell(radius = 5*units.um**2, markers = ["EGFR", "CA9"], copies = [80000, 80000])


system = System(compartment = "container", cells = [Tcell, HT29], drugs = [X], solutes = [])
system.add_process()


