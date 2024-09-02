from .qsp import *

def flowcytom(x, num_A, num_B, aff_A, off_A, aff_B, off_B, CAB, t, cell_density = 1e6 / units.ml):
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
  system.run(t)
  return {analyte:x for analyte, x in zip(system.analytes, system.x[:,0])}
