from .qsp import *

class Drug:
  def __init__(self, num_binding_sites):
    self.num_binding_sites = num_binding_sites
  
  # affinity定义的是自由药物（即状态为空）时结合各ligand的能力，以及脱落至自由状态的速率
  def add_affinity(self, binding_site, ligand, aff, off):
    pass
  
  # cis定义的是已经在细胞表面时，结合同一细胞各ligand的能力，以及解离（之后仍在细胞上）的速率。这是avidity、cross-arm binding等现象的广义描述
  # 二维浓度的单位是 1/um^2，因此二维on rate的单位是 um^2/s
  def add_cis(self, binding_site, ligand, on2D, off, complex = None):
    pass
  
  # trans定义的是已经在细胞表面时，如果此细胞与另一细胞接触，结合对方细胞各ligand的能力，以及解离的速率。这是形成trimer等过程的广义描述。
  def add_trans(self, binding_site, ligand, on2D, off, complex = None):
    pass

  def finalize()


X = Drug(3)
X.add_affinity(0, "CD3", 1*units.nM, 1e-4/units.s)
X.add_affinity(1, "EGFR", 10*units.nM, 1e-4/units.s)
X.add_affinity(2, "CA9", 20*units.nM, 1e-4/units.s)
X.add_cis(1, "EGFR", 1 * units.mm**2/units.s, state = ",,CA9")
X.add_cis(2, "CA9", 1 * units.mm**2/units.s, state = ",EGFR,")
X.add_trans(0, "CD3", 1 * units.mm**2/units.s)
X.add_trans(1, "EGFR", 1 * units.mm**2/units.s)
X.add_trans(2, "CA9", 1 * units.mm**2/units.s)



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


