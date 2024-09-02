from .qsp import *

def flowcytom(x, num_A, num_B, aff_A, off_A, int_A, aff_B, off_B, int_B, CAB, t, cell_density = 1e6 / units.ml):
  compartments = ["container"]
  system = System(compartments, analytes)
  system.set_volume("container", 100*units.ul)
  
  system.add_simple("container", ["A", "drug"], ["A-drug"], forward = off_A/aff_A, backward = off_A)
  system.add_simple("container", ["B", "drug"], ["B-drug"], forward = off_B/aff_B, backward = off_B)
  system.add_simple("container", ["A", "B-drug"], ["AB-drug"], forward = off_A/aff_A * CAB, backward = off_A)
  system.add_simple("container", ["B", "A-drug"], ["AB-drug"], forward = off_B/aff_B * CAB, backward = off_B)
  
  system.add_x("container", "A", num_A * cell_density / units.avagadro)
  system.add_x("container", "B", num_B * cell_density / units.avagadro)
  system.add_x("container", "drug", x)
  
  system.add_flow("A-drug", "container", None, int_A * 100*units.ul)
  system.add_flow("B-drug", "container", None, int_A * 100*units.ul)
  system.add_flow("AB-drug", "container", None, min(int_A, int_B) * 100*units.ul)
  
  system.run(t)
  
  nums = [(x*units.nM / cell_density * units.avagadro).number(1) for x in system.x[:,0]]
  return {analyte:num for analyte, num in zip(system.analytes, nums)}
