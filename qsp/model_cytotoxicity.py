from .qsp import *

class Drug:
  def __init__(self, num_binding_sites):
    self.num_binding_sites = num_binding_sites
  
  def add_affinity(self, ligand, aff, aff2D, off):
    pass
  
  def add_cis(self, ligand1, ligand2, on2D):
    pass
  
  def add_trans(self, ligand1, ligand2, on2D):
    pass

drug = Drug(3)
drug.add_affinity(0, "CD3", 1*units.nM, 1e-4/units.s)
drug.add_affinity(1, "EGFR", 10*units.nM, 1e-4/units.s)
drug.add_affinity(2, "CA9", 20*units.nM, 1e-4/units.s)
drug.add_cis(1, 2, )


nM = mole/l = 1e-15 mole/um**3
mole/um**2 = 1e15 nM * um

(1/um**2) / (1/um**2) / (1/um**2) / s = um**2/s
