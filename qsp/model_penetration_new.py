from .qsp import *
from qsp.processes import *
from qsp.human import *
from qsp.tumors import *


def model(num_antigen, aff, off, int_rate, half_life,
          tumor_capillary_radius = 10 * units.um, tumor_capillary_permeability = 3e-7 * units.cm/units.s, tumor_diffusion = 10 * units.um**2 / units.s, tumor_cell_density = 3e8 * 0.44 / units.ml,
          tumor_capillary = 100*units.cm, tumor_layer_depth = 10*units.um, tumor_num_layers = 100):
  analytes = ["antigen", "drug", "antigen-drug"]
  centrals = [plasma, lymph]
  organs = [bone, lung, liver, SI, other]
  tumors = [FTC238.copy() for _ in range(10)]
  for i in range(10):
    tumors[i]["name"] = f"tumor_{i}"
  
  compartments = [organ["name"] for organ in centrals + organs + tumors]
  system = System(compartments, analytes)
  
  # define volumes
  for central in centrals:
    system.set_volume(central["name"], central["volume"])
  for organ in organs:
    system.set_volume(organ["name"], organ["volume_interstitial"])
  for i in range(tumor_num_layers):
    radius = tumor_capillary_radius + i * tumor_layer_depth
    area = tumor_capillary * (2 * math.pi * radius)
    system.set_volume(f"tumor_{i}", area * tumor_layer_depth)
  
  # organ distribution
  for organ in organs:
    system.add_flow("drug", "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
    system.add_flow("drug", organ["name"], "lymph", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
    system.add_flow("drug", "lymph", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
  
  # tumor distribution
  system.add_flow("drug", "plasma", "tumor_0", tumor_capillary * (2 * math.pi * tumor_capillary_radius) * tumor_capillary_permeability)
  system.add_flow("drug", "tumor_0", "plasma", tumor_capillary * (2 * math.pi * tumor_capillary_radius) * tumor_capillary_permeability)
  
  for i in range(tumor_num_layers - 1):
    radius = tumor_capillary_radius + (i + 1) * tumor_layer_depth
    area = tumor_capillary * (2 * math.pi * radius)
    system.add_flow("drug", f"tumor_{i}", f"tumor_{i+1}", area/tumor_layer_depth * tumor_diffusion)
    system.add_flow("drug", f"tumor_{i+1}", f"tumor_{i}", area/tumor_layer_depth * tumor_diffusion)
  
  # bulk clearance
  for compartment in compartments:
    system.add_flow(drug, compartment, None, system.get_volume(compartment) * math.log(2)/half_life)
  
  # binding and internalizing kinetics
  for compartment in compartments:
    system.add_simple(compartment, ["antigen", "drug"], ["antigen-drug"], aff*off, off)
    system.add_transform(compartment, "antigen-drug", ["antigen"], rate = int_rate)
  
  # initial concentrations
  for i in range(tumor_num_layers):
    system.set_x("antigen", f"tumor_{i}", num_antigen * tumor_cell_density / units.avagadro)


############# demo ###############
'''
from qsp import *
from qsp.processes import *
from qsp.human import *
from qsp.tumors import *
from qsp.model_PD1 import *

system = model(num_antigen = 10000, aff = 1*units.nM, off = 1e-4*units.s, int_rate = 0.2/units.h, half_life = 80*units.h)
system.add_x("plasma", "n", 100 * units.nM)
system.run(168 * units.h, t_step = 1/6 * units.h)

coverages = np.array([x[-1, -1] / (x[-1, -1] + x[0, -1]) for t, x, c in system.history])
coverage = coverages.mean()

'''
