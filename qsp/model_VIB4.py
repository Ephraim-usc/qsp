from .qsp import *

### this model is mostly from ...

targets = ["A", "B", "AB"]
drugs = [f"{c}{a}" for c in ["m", "n"] for a in ["m", "n", "a"]]
dimers = [f"{drug}-{target}" for drug in drugs for target in targets] + [f"C-{drug}" for drug in drugs]
trimers = [f"C-{drug}-{target}" for drug in drugs for target in targets]
analytes = ["C", "A", "B", "a"] + drugs + dimers + trimers



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
  def __init__(self, compartments, rate_C, rate_A):
    self.system = None
    self.compartments = compartments
    
    Q = np.zeros([7, 7])
    rate_C = rate_C.number(1/units.h)
    rate_A = rate_A.number(1/units.h)
    Q[[3,4,5], [0,1,2]] += rate_C; Q[[0,1,2], [0,1,2]] -= rate_C
    Q[[2,5], [0,3]] += rate_A; Q[[0,3], [0,3]] -= rate_A
    Q[[6], [0,3]] += rate_A
    self.Q = Q
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      if callable(self.compartments):
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments(system) if compartment in system.compartments]
      else:
        self.compartments_ = [system.compartments.index(compartment) for compartment in self.compartments if compartment in system.compartments]
      
      self.analyteses_ = []
      self.analyteses_.append([system.analytes.index(f"{drug}") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"{drug}-A") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"{drug}-B") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"{drug}-AB") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"C-{drug}") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-A") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-B") for drug in drugs] + [system.analytes.index("a")])
      self.analyteses_.append([system.analytes.index(f"C-{drug}-AB") for drug in drugs] + [system.analytes.index("a")])
    
    for compartment in self.compartments_:
      for analytes_ in self.analyteses_:
        x = system.x[analytes_, compartment]
        system.x[analytes_, compartment] = expm(self.Q * t.number(units.h)) @ x



'''
VIB1 = {}
VIB1.update({"off_C": 1.3e-2 / units.s, "affn_C": 9.034e-8 * units.molar, "affm_C": 8.545e-7 * units.molar})
VIB1.update({"off_A": 1e-2 / units.s, "affn_A": 2.14e-7 * units.molar, "affm_A": 9.163e-7 * units.molar})
VIB1.update({"off_B": 4.6e-4 / units.s, "aff_B": 1.1e-9 * units.molar})
VIB1.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIB1.update({"avidity": 69})
VIB1.update({"clearance": math.log(2)/(40 * units.h)})
VIB1.update({"internalization_Tcell": 0.5 / units.h, "internalization_tumor": 0.1 / units.h, "internalization_organ": 0.05 / units.h})
VIB1["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = 0.0527 / units.d, 
                                  rate_A = 0.0527 / units.d)
VIB1["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = 0.1783 / units.d, 
                                 rate_A = 0.1783 / units.d)

VIB5 = {}
VIB5.update({"off_C": 1e-2 / units.s, "affn_C": 1.26e-8 * units.molar, "affm_C": 4.16e-7 * units.molar})
VIB5.update({"off_A": 3e-3 / units.s, "affn_A": 2e-8 * units.molar, "affm_A": 3.34e-7 * units.molar})
VIB5.update({"off_B": 4e-4 / units.s, "aff_B": 1.7e-9 * units.molar})
VIB5.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIB5.update({"avidity": 69})
VIB5.update({"clearance": math.log(2)/(40 * units.h)})
VIB5.update({"internalization_Tcell": 0.5 / units.h, "internalization_tumor": 0.1 / units.h, "internalization_organ": 0.05 / units.h})
VIB5["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = 0.0527 / units.d, 
                                  rate_A = 0.0527 / units.d)
VIB5["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = 0.1783 / units.d, 
                                 rate_A = 0.1783 / units.d)
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
VIBX.update({"off_C": 3.0720e-4 / units.s, "affn_C": 2.4589e-10 * units.molar, "affm_C": 2.4589e-8 * units.molar})
VIBX.update({"off_A": 1.2059e-5 / units.s, "affn_A": 1.0232e-9 * units.molar, "affm_A": 1.0232e-7 * units.molar})
VIBX.update({"off_B": 1.2344e-4 / units.s, "aff_B": 1e-8 * units.molar})
VIBX.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIBX.update({"avidity": 19})
VIBX.update({"clearance": math.log(2)/(70 * units.h)}); VIBX["smalls"] = ["a"]
VIBX.update({"internalization_Tcell": 0.1 / units.h, "internalization_tumor": 0.02 / units.h, "internalization_organ": 0.02 / units.h})
VIBX["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = 0.05 / units.d, 
                                  rate_A = 0.01687 / units.d)
VIBX["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = 0.15 / units.d, 
                                 rate_A = 0.05061 / units.d)

VIBY = {}
VIBY.update({"off_C": 7.0734e-4 / units.s, "affn_C": 3.2902e-10 * units.molar, "affm_C": 3.2902e-8 * units.molar})
VIBY.update({"off_A": 1.5704e-5 / units.s, "affn_A": 1e-8 * units.molar, "affm_A": 1e-6 * units.molar})
VIBY.update({"off_B": 3.5120e-5 / units.s, "aff_B": 1e-8 * units.molar})
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

JANUX = {}
JANUX.update({"off_C": 1.36e-2 / units.s, "affn_C": 2e-9 * units.molar, "affm_C": 1.55e-6 * units.molar}) #
JANUX.update({"off_A": 1.37e-2 / units.s, "affn_A": 3e-9 * units.molar, "affm_A": 5.51e-7 * units.molar}) #
JANUX.update({"off_B": 4.138e-4 / units.s, "aff_B": math.inf * units.molar})
JANUX.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
JANUX.update({"avidity": 69})
JANUX.update({"clearance": math.log(2)/(100 * units.h)}); JANUX["smalls"] = ["nm", "nn", "na", "a"]
JANUX.update({"internalization_Tcell": 0.1 / units.h, "internalization_tumor": 0.02 / units.h, "internalization_organ": 0.02 / units.h})
JANUX["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ['name'] for organ in system.organs], 
                                  rate_C = 0.0527 / units.d, 
                                  rate_A = 0.0527 / units.d)
JANUX["cleavage_tumor"] = cleavage(lambda system: [tumor['name'] for tumor in system.tumors], 
                                 rate_C = 0.1783 / units.d, 
                                 rate_A = 0.1783 / units.d)


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
    for drug in drugs + ["a"]:
        system.add_flow(drug, compartment, None, system.get_volume(drug, compartment) * TCE["clearance"])
  
  # small forms plasma clearance
  for small in TCE["smalls"]:
    system.add_flow(small, "plasma", None, system.get_volume(drug, "plasma") * math.log(2)/(45 * units.MIN))
  
  for drug in drugs + ["a"]:
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
  
  for compartment in compartments:
    system.add_simple(compartment, ["mn", "a"], ["ma"], TCE["off_a"] / TCE["aff_a"], TCE["off_a"])
    system.add_simple(compartment, ["nn", "a"], ["na"], TCE["off_a"] / TCE["aff_a"], TCE["off_a"])
    system.add_simple(compartment, ["C-mn", "a"], ["ma"], TCE["off_a"] / TCE["aff_a"], TCE["off_a"])
    system.add_simple(compartment, ["C-nn", "a"], ["na"], TCE["off_a"] / TCE["aff_a"], TCE["off_a"])
  
  # target binding
  for drug in drugs:
    off_C = TCE["off_C"]; on_C = TCE["off_C"] / TCE["affn_C"] if drug[0] == "n" else TCE["off_C"] / TCE["affm_C"]
    off_A = TCE["off_A"]; on_A = TCE["off_A"] / TCE["affn_A"] if drug[1] == "n" else TCE["off_A"] / TCE["affm_A"]
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
