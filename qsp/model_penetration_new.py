from .qsp import *
from qsp.human import *


def model(num_antigen, aff, off, int_rate, half_life,
          tumor_capillary_radius = 10 * units.um, tumor_capillary_permeability = 3e-7 * units.cm/units.s, tumor_diffusion = 10 * units.um**2 / units.s, tumor_cell_density = 3e8 * 0.44 / units.ml,
          tumor_capillary = 100*units.cm, tumor_layer_depth = 10*units.um, tumor_num_layers = 100):
  analytes = ["antigen", "drug", "antigen-drug"]
  centrals = [plasma, lymph]
  organs = [bone, lung, liver, SI, other]
  
  compartments = [organ["name"] for organ in centrals + organs] + [f"tumor_{i}" for i in range(tumor_num_layers)]
  system = System(compartments, analytes)
  system.params = locals().copy()
  
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
    system.add_flow("drug", compartment, None, system.get_volume(compartment) * math.log(2)/half_life)
  
  # binding and internalizing kinetics
  for compartment in compartments:
    system.add_simple(compartment, ["antigen", "drug"], ["antigen-drug"], off/aff, off)
    system.add_transform(compartment, "antigen-drug", ["antigen"], rate = int_rate)
  
  # initial concentrations
  for i in range(tumor_num_layers):
    system.set_x(f"tumor_{i}", "antigen", num_antigen * tumor_cell_density / units.avagadro)
  
  return system


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_penetration(system, compartments, labels, colors, linestyles = None, title = None, output = "tmp.png"):
  if compartments is None:
    compartments = system.compartments
  compartments = [system.compartments.index(compartment) for compartment in compartments]
  
  if linestyles is None:
    linestyles = ["solid"] * len(compartments)

  SF = units.nM * units.avagadro / system.params["cell_density"]
  Xmax = max([t for t, x, c in system.history])
  Ymax = max([x[2, compartment]*SF for t, x, c in system.history for compartment in compartments])
  
  fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (4, 4))
  for compartment, label, color, linestyle in zip(compartments, labels, colors, linestyles):
    X = [t for t, x, c in system.history]
    Y = [x[2, compartment]*SF for t, x, c in system.history]
    RATIO = max(Y) / Ymax
    ax.plot(X, Y, label = f"{label}, avg={RATIO * 100:.2}%", color = color)
  
  if Xmax > 100:
    ax.set_xticks([0, 24, 48, 72, 96, 120, 144, 168])
    ax.set_xticklabels(["0", "1", "2", "3", "4", "5", "6", "7"])
  else:
    ax.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
  
  ax.set_xlim(0, Xmax)
  if Xmax > 100:
    ax.set_xlabel("time (d)")
  else:
    ax.set_xlabel("time (h)")
  ax.set_ylim(0, Ymax)
  ax.grid(axis = "y", color = "grey", linewidth = 1)
  if title:
    ax.set_title(title)
  ax.legend(loc = "upper right", prop={'size': 6})
  
  if output is None:
    fig.show()
  else:
    fig.savefig(output, dpi = 300)
    plt.close(fig)


############# demo ###############
'''
from qsp import *
from qsp.model_penetration_new import *

system = model(num_antigen = 10000, aff = 1*units.nM, off = 1e-4/units.s, int_rate = 0.2/units.h, half_life = 80*units.h)
system.add_x("plasma", "drug", 10 * units.nM)
system.run(168 * units.h, t_step = 1/6 * units.h, verbose = True)

compartments = ["plasma", "tumor_0", "tumor_5", "tumor_10", "tumor_15", "tumor_20"]
plot_penetration(system, 
                 compartments = ["plasma", "tumor_0", "tumor_5", "tumor_10", "tumor_15", "tumor_20"], 
                 labels = ["plasma", "0um", "50um", "100um", "150um", "200um"], 
                 colors = ["black", "red", "orange", "gold", "green", "blue"],
                 output = "tmp.png")


system.plot(compartments = ["tumor_0", "tumor_10", "tumor_20", "tumor_30", "tumor_40"], 
              #groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"penetration.png")

'''
