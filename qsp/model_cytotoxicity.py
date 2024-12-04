from .qsp import *

class Drug:
  def __init__(self, num_binding_sites):
    self.num_binding_sites = num_binding_sites
  
  def add_affinity(self, ligand, aff, off):
    buffer = transform([], [], [])
    buffer.Qs =  add_two_dicts(self.Qs, transform2.Qs)
    return buffer

drug = Drug(3)
drug.add_affinity("CD3", 1*units.nM, 1e-4/units.s)
drug.add_affinity("EGFR", 10*units.nM, 1e-4/units.s)
drug.add_affinity("CA9", 20*units.nM, 1e-4/units.s)
drug.add_avidity("CA9", 20*units.nM, 1e-4/units.s)
