from .qsp import *

### this model is mostly from ...

targets = ["A", "B", "AB"]
drugs = ["m", "n"]
dimers = [f"{drug}-{target}" for drug in drugs for target in targets] + [f"C-{drug}" for drug in drugs]
trimers = [f"C-{drug}-{target}" for drug in drugs for target in targets]
analytes = ["C", "A", "B"] + drugs + dimers + trimers



############ constants ############

molecular_weight = 150000 * units.g/units.mol

class intratumoral_equilibrium:
  def __init__(self):
    self.system = None
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      self.compartments_ = [system.compartments.index(tumor["name"]) for tumor in system.tumors]
      self.analytes_ = [system.analytes.index(f"{drug}") for drug in drugs]
    
    for analyte_ in self.analytes_:
      x = system.x[analyte_, self.compartments_]
      volumes = system.V[analyte_, self.compartments_]
      system.x[analyte_, self.compartments_] = np.average(x, weights = volumes)


############ host ############

mouse = {}
mouse.update({"volume_plasma": 1.26 * units.ml})

human = {}
human.update({"volume_plasma": 2877 * units.ml})


############ drug ############

class cleavage:
  def __init__(self, compartments, rate):
    self.system = None
    self.compartments = compartments
    
    Q = np.zeros([2, 2])
    rate = rate.number(1/units.h)
    Q[1, 0] += rate; Q[0, 0] -= rate
    self.Q = Q
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      if callable(self.compartments):
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments(system) if compartment in system.compartments]
      else:
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments if compartment in system.compartments]
      
      self.analyteses_ = []
      self.analyteses_.append([system.analytes.index(f"{drug}") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"{drug}-A") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"{drug}-B") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"{drug}-AB") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-A") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-B") for drug in drugs])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-AB") for drug in drugs])
    
    for compartment in self.compartments_:
      for analytes_ in self.analyteses_:
        x = system.x[analytes_, compartment]
        system.x[analytes_, compartment] = expm(self.Q * t.number(units.h)) @ x

'''
VIB4 = {}
VIB4.update({"off_C": 8.09e-3 / units.s, "affn_C": 6.56e-8 * units.molar, "affm_C": 3.059e-6 * units.molar})
VIB4.update({"off_A": 3e-3 / units.s, "affn_A": 2e-8 * units.molar, "affm_A": 3.34e-7 * units.molar})
VIB4.update({"off_B": 4.138e-4 / units.s, "aff_B": 1.7e-9 * units.molar})
VIB4.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIB4.update({"avidity": 69})
VIB4.update({"clearance": math.log(2)/(70 * units.h)}); VIB4["smalls"] = ["a"]
VIB4.update({"internalization_Tcell": 0.1 / units.h, "internalization_tumor": 0.02 / units.h, "internalization_organ": 0.02 / units.h})
VIB4["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = 0.0527 / units.d, 
                                  rate_A = 0.0527 / units.d)
VIB4["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = 0.1783 / units.d, 
                                 rate_A = 0.1783 / units.d)

VIBX = {}
VIBX.update({"off_C": 10**-3.67417 / units.s, "affn_C": 10**-9.75451 * units.molar, "affm_C": 10**-7.75451 * units.molar})
VIBX.update({"off_A": 10**-5 / units.s, "affn_A": 10**-8 * units.molar, "affm_A": 10**-6 * units.molar})
VIBX.update({"off_B": 10**-5 / units.s, "aff_B": 10**-8.71011 * units.molar})
VIBX.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIBX.update({"avidity": 19})
VIBX.update({"clearance": math.log(2)/(70 * units.h)}); VIBX["smalls"] = ["a"]
VIBX.update({"internalization_Tcell": 0.1 / units.h, "internalization_tumor": 0.02 / units.h, "internalization_organ": 0.02 / units.h})
VIBX["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = 0.05 / units.d, 
                                  rate_A = 0.01546 / units.d)
VIBX["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = 0.15 / units.d, 
                                 rate_A = 0.04638 / units.d)

VIBY = {}
VIBY.update({"off_C": 10**-3.59947 / units.s, "affn_C": 10**-9.73151 * units.molar, "affm_C": 10**-7.73151 * units.molar})
VIBY.update({"off_A": 10**-5 / units.s, "affn_A": 10**-8 * units.molar, "affm_A": 10**-6 * units.molar})
VIBY.update({"off_B": 10**-4.82697 / units.s, "aff_B": 10**-8.82898 * units.molar})
VIBY.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIBY.update({"avidity": 19})
VIBY.update({"clearance": math.log(2)/(70 * units.h)}); VIBY["smalls"] = ["mn", "nn", "a"]
VIBY.update({"internalization_Tcell": 0.1 / units.h, "internalization_tumor": 0.02 / units.h, "internalization_organ": 0.02 / units.h})
VIBY["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = 0.05 / units.d, 
                                  rate_A = 0.05 / units.d)
VIBY["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = 0.15 / units.d, 
                                 rate_A = 0.15 / units.d)
'''

VIBZ = {}
VIBZ.update({"off_C": 10**-3.59947 / units.s, "affn_C": 10**-9.73151 * units.molar, "affm_C": 10**-7.73151 * units.molar})
VIBZ.update({"off_A": 10**-5 / units.s, "affn_A": 10**-8 * units.molar, "affm_A": 10**-6 * units.molar})
VIBZ.update({"off_B": 10**-4.82697 / units.s, "aff_B": 10**-8.82898 * units.molar})
VIBZ.update({"avidity": 19})
VIBZ.update({"clearance": math.log(2)/(70 * units.h)}); VIBZ["smalls"] = ["n"]
VIBZ.update({"internalization_Tcell": 0.1 / units.h, "internalization_tumor": 0.02 / units.h, "internalization_organ": 0.02 / units.h})
VIBZ["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], rate = 0.05 / units.d)
VIBZ["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], rate = 0.15 / units.d)


############ tumors ############

UT44 = {"name": "tumor"}
UT44.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
UT44.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
UT44.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
UT44.update({"diffusion": 10 * units.micrometer**2 / units.s})
UT44.update({"cell_density": 3e8 / units.ml, "T_cell_density": 3e7 / units.ml})
UT44.update({"num_A": 7e5, "num_B": 1.45e6})

FTC238 = {"name": "tumor"}
FTC238.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
FTC238.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
FTC238.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
FTC238.update({"diffusion": 10 * units.micrometer**2 / units.s})
FTC238.update({"cell_density": 3e8 / units.ml, "T_cell_density": 3e7 / units.ml})
FTC238.update({"num_A": 1e5, "num_B": 1e5})


############ organs ############

other = {"name": "other"}
other.update({"volume_plasma": 500 * units.ml, "volume_interstitial": 3000 * units.ml})
other.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
other.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
other.update({"cell_density": 1e8 / units.ml, "T_cell_density": 3e6 / units.ml})
other.update({"num_A": 100000, "num_B": 0})

lung = {"name": "lung"}
lung.update({"volume_plasma": 55 * units.ml, "volume_interstitial": 300 * units.ml})
lung.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
lung.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"cell_density": 1e8 / units.ml, "T_cell_density": 3e7 / units.ml})
lung.update({"num_A": 133439, "num_B": 0}) # num_B = 1019 from Liyuan

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"cell_density": 1e8 / units.ml, "T_cell_density": 3e6 / units.ml})
SI.update({"num_A": 57075, "num_B": 39649})


############ model ############

def model(host, TCE, tumors, organs, connect_tumors = True):
  compartments = ["plasma"] + [organ["name"] for organ in tumors + organs]
  system = System(analytes, compartments)
  system.tumors = tumors
  system.organs = organs
  
  for analyte in analytes:
    system.set_volume(analyte, "plasma", host["volume_plasma"])
    for tumor in tumors:
      system.set_volume(analyte, tumor["name"], tumor["volume"] * tumor["volume_interstitial_proportion"])
    for organ in organs:
      system.set_volume(analyte, organ["name"], organ["volume_interstitial"])
  
  # whole-body clearance
  for compartment in compartments:
    for drug in drugs:
        system.add_flow(drug, compartment, None, system.get_volume(drug, compartment) * TCE["clearance"])
  
  # small forms plasma clearance
  for small in TCE["smalls"]:
    system.add_flow(small, "plasma", None, system.get_volume(drug, "plasma") * math.log(2)/(45 * units.MIN))
  
  for drug in drugs:
    # tumor flow
    for tumor in tumors:
      system.add_flow(drug, "plasma", tumor["name"], tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
      system.add_flow(drug, tumor["name"], "plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
    
    # organ flow
    for organ in organs:
      system.add_flow(drug, "plasma", organ["name"], organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["vascular_reflection"]))
      system.add_flow(drug, organ["name"], "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - organ["lymphatic_reflection"]))
  
  # exchange drugs between tumors if tumors are connected
  if connect_tumors:
    system.add_process(intratumoral_equilibrium())
  
  # mask cleavage
  system.add_process(TCE["cleavage_plasma"])
  system.add_process(TCE["cleavage_tumor"])
  
  # target binding
  for drug in drugs:
    off_C = TCE["off_C"]; on_C = {"n":TCE["off_C"] / TCE["affn_C"], "m":TCE["off_C"] / TCE["affm_C"]}[drug]
    off_A = TCE["off_A"]; on_A = {"n":TCE["off_A"] / TCE["affn_A"], "m":TCE["off_A"] / TCE["affm_A"]}[drug]
    off_B = TCE["off_B"]; on_B = TCE["off_B"] / TCE["aff_B"]
    avidity = TCE["avidity"]
    
    system.add_simple("plasma", ["C", f"{drug}"], [f"C-{drug}"], on_C, off_C)
    
    for organ in tumors + organs:
      system.add_simple(organ["name"], [f"{drug}", "A"], [f"{drug}-A"], on_A, off_A)
      system.add_simple(organ["name"], [f"{drug}", "B"], [f"{drug}-B"], on_B, off_B)
      system.add_simple(organ["name"], [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * avidity, off_B)
      system.add_simple(organ["name"], [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * avidity, off_A)
      
      system.add_simple(organ["name"], ["C", f"{drug}"], [f"C-{drug}"], on_C, off_C)
      system.add_simple(organ["name"], ["C", f"{drug}-A"], [f"C-{drug}-A"], on_C, off_C)
      system.add_simple(organ["name"], ["C", f"{drug}-B"], [f"C-{drug}-B"], on_C, off_C)
      system.add_simple(organ["name"], ["C", f"{drug}-AB"], [f"C-{drug}-AB"], on_C, off_C)
      
      system.add_simple(organ["name"], [f"C-{drug}", "A"], [f"C-{drug}-A"], on_A, off_A)
      system.add_simple(organ["name"], [f"C-{drug}", "B"], [f"C-{drug}-B"], on_B, off_B)
      system.add_simple(organ["name"], [f"C-{drug}-A", "B"], [f"C-{drug}-AB"], on_B * avidity, off_B)
      system.add_simple(organ["name"], [f"C-{drug}-B", "A"], [f"C-{drug}-AB"], on_A * avidity, off_A)

  # internalization
  for drug in drugs:
    for compartment in compartments:
      system.add_simple(compartment, [f"C-{drug}"], ["C"], TCE["internalization_Tcell"])
    for tumor in tumors:
      system.add_simple(tumor["name"], [f"{drug}-A"], ["A"], TCE["internalization_tumor"])
      system.add_simple(tumor["name"], [f"{drug}-B"], ["B"], TCE["internalization_tumor"])
      system.add_simple(tumor["name"], [f"{drug}-AB"], ["A", "B"], TCE["internalization_tumor"])
    for organ in organs:
      system.add_simple(organ["name"], [f"{drug}-A"], ["A"], TCE["internalization_organ"])
      system.add_simple(organ["name"], [f"{drug}-B"], ["B"], TCE["internalization_organ"])
      system.add_simple(organ["name"], [f"{drug}-AB"], ["A", "B"], TCE["internalization_organ"])
  
  
  # initial concentrations
  T_cell_density_blood = 1e6 / units.ml
  system.add_x("C", "plasma", 124000 * T_cell_density_blood / units.avagadro)
  
  for tumor in tumors:
    system.add_x("C", tumor["name"], 124000 * tumor["T_cell_density"] / units.avagadro)
    system.add_x("A", tumor["name"], tumor["num_A"] * tumor["cell_density"] / units.avagadro)
    system.add_x("B", tumor["name"], tumor["num_B"] * tumor["cell_density"] / units.avagadro)
  
  for organ in organs:
    system.add_x("C", organ["name"], 124000 * organ["T_cell_density"] / units.avagadro)
    system.add_x("A", organ["name"], organ["num_A"] * organ["cell_density"] / units.avagadro)
    system.add_x("B", organ["name"], organ["num_B"] * organ["cell_density"] / units.avagadro)
  
  return system


############ plot ############

def plot(system, name):
  #pickle.dump(system, open(f"{name}.pickle", "wb"))
  
  targets = ["A", "B", "AB"]
  groups = [["C"],
            ["A", "B"],
            drugs,
            [f"C-{drug}" for drug in drugs],
            [f"{drug}-{target}" for drug in drugs for target in targets],
            [f"C-{drug}-{target}" for drug in drugs for target in targets]]
  labels = ["C", "target", "cap", "drug", "C-drug", "drug-target", "C-drug-target"]
  colors = ["tab:orange", "tab:blue", "black", "wheat", "skyblue", "tab:purple"]
  linestyles = ["solid", "solid", "solid", "solid", "solid", "solid"]
  system.plot(compartments = ["plasma"] + [tumor["name"] for tumor in system.tumors] + [organ["name"] for organ in system.organs], 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")
  
  groups = [drugs + [f"{drug}-{target}" for drug in drugs for target in targets] + [f"C-{drug}" for drug in drugs] + [f"C-{drug}-{target}" for drug in drugs for target in targets],
            ["n"] + [f"n-{target}" for target in targets] + [f"C-n"] + [f"C-n-{target}" for target in targets]]
  labels = ["(C-)x(-target)", "(C-)n(-target)"]
  colors = ["tab:green", "wheat"]
  linestyles = ["solid", "dashed"]
  system.plot(compartments = ["plasma"] + [tumor["name"] for tumor in system.tumors] + [organ["name"] for organ in system.organs], 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_drugs.png")
