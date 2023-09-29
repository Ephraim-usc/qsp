from .qsp import *

### this model is mostly from ...

antigens = ["C", "A", "B", "a"]
drugs = [f"{c}{a}" for c in ["m", "n"] for a in ["m", "n", "a"]]
C_conjugates = [f"{drug}-{Ag}" for drug in drugs for Ag in ["A", "B", "AB"]]
T_conjugates = [f"C-{drug}" for drug in drugs]
trimers = [f"C-{drug}-{Ag}" for drug in drugs for Ag in ["A", "B", "AB"]]
analytes = antigens + drugs + C_conjugates + T_conjugates + trimers



############ constants ############

molecular_weight = 150000 * units.g/units.mol

class intratumoral_equilibrium:
  def __init__(self):
    self.system = None
  
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      
      self.compartments_ = [system.compartments.index(f"{tumor['name']}_interstitial") for tumor in system.tumors]
      self.analytes_ = [system.analytes.index(f"{drug}") for drug in drugs]
    
    for analyte_ in self.analytes_:
      x = system.x[analyte_, self.compartments_]
      volumes = system.V[analyte_, self.compartments_]
      system.x[analyte_, self.compartments_] = np.average(x, weights = volumes)


############ host ############

mouse = {}
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})

human = {}
human.update({"volume_central": 2877 * units.ml, "volume_peripheral": 2857 * units.ml})
human.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})


############ drug ############

class cleavage:
  def __init__(self, compartments, rate_C, rate_A):
    self.system = None
    self.compartments = compartments
    
    Q = np.zeros([7, 7])
    rate_C = rate_C.number(1/units.h)
    rate_A = rate_A.number(1/units.h)
    rate_B = rate_B.number(1/units.h)
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


VIB4 = {}
VIB4.update({"off_C": 8.09e-3 / units.s, "affn_C": 6.56e-8 * units.molar, "affm_C": 3.059e-6 * units.molar})
VIB4.update({"off_A": 3e-3 / units.s, "affn_A": 2e-8 * units.molar, "affm_A": 3.34e-7 * units.molar})
VIB4.update({"off_B": 4.138e-4 / units.s, "aff_B": 1.7e-9 * units.molar})
VIB4.update({"off_a": 8.09e-3 / units.s, "aff_a": 1e-9 * units.molar})
VIB4.update({"avidity": 20})
VIB4.update({"clearance": math.log(2)/(40 * units.h)})
VIB4["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [f"{organ['name']}_interstitial" for organ in system.organs], 
                                  rate_C = 0.0527 / units.d, 
                                  rate_A = 0.0527 / units.d)
VIB4["cleavage_tumor"] = cleavage(lambda system: [f"{tumor['name']}_interstitial" for tumor in system.tumors], 
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
other.update({"cell_density": 1e8 / units.ml, "T_cell_density": 3e6 / units.ml})
other.update({"num_A": 100000, "num_B": 0})

lung = {"name": "lung"}
lung.update({"volume_plasma": 55 * units.ml, "volume_interstitial": 300 * units.ml})
lung.update({"plasma_flow": 181913 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"cell_density": 1e8 / units.ml, "T_cell_density": 3e7 / units.ml})
lung.update({"num_A": 133439, "num_B": 0}) # num_B = 1019 from Liyuan

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"cell_density": 1e8 / units.ml, "T_cell_density": 3e6 / units.ml})
SI.update({"num_A": 57075, "num_B": 39649})


############ model ############

def model(host, TCE, tumors, organs, connect_tumors = False):
  compartments = ["plasma"] + [f"{organ['name']}_{tissue}" for organ in tumors + organs for tissue in ["plasma", "interstitial"]]
  system = System(analytes, compartments)
  system.tumors = tumors
  system.organs = organs
  
  for analyte in analytes:
    system.set_volume(analyte, "plasma", host["volume_central"])
    for tumor in tumors:
      system.set_volume(analyte, f"{tumor['name']}_plasma", tumor["volume"] * tumor["volume_plasma_proportion"])
      system.set_volume(analyte, f"{tumor['name']}_interstitial", tumor["volume"] * tumor["volume_interstitial_proportion"])
    for organ in organs:
      system.set_volume(analyte, f"{organ['name']}_plasma", organ["volume_plasma"])
      system.set_volume(analyte, f"{organ['name']}_interstitial", organ["volume_interstitial"])
  
  
  for drug in drugs:
    # whole-body clearance
    for compartment in compartments:
      system.add_flow(drug, compartment, None, system.get_volume(drug, compartment) * TCE["clearance"])
    
    # tumor flow
    for tumor in tumors:
      system.add_flow(drug, "plasma", f"{tumor['name']}_plasma", tumor["volume"] * tumor["plasma_flow_density"])
      system.add_flow(drug, f"{tumor['name']}_plasma", "plasma", tumor["volume"] * tumor["plasma_flow_density"])
      system.add_flow(drug, f"{tumor['name']}_plasma", f"{tumor['name']}_interstitial", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
      system.add_flow(drug, f"{tumor['name']}_interstitial", f"{tumor['name']}_plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
    
    # organ flow
    for organ in organs:
      system.add_flow(drug, "plasma", f"{organ['name']}_plasma", organ["plasma_flow"])
      system.add_flow(drug, f"{organ['name']}_plasma", "plasma", organ["plasma_flow"])
      system.add_flow(drug, f"{organ['name']}_plasma", f"{organ['name']}_interstitial", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - host["vascular_reflection"]))
      system.add_flow(drug, f"{organ['name']}_interstitial", "plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - host["lymphatic_reflection"]))
  
  
  # exchange drugs between tumors if tumors are connected
  if connect_tumors:
    system.add_process(intratumoral_equilibrium())
  
  
  # mask cleavage
  system.add_process(TCE["cleavage_plasma"])
  system.add_process(TCE["cleavage_tumor"])
  
  for compartment in compartments:
    system.add_simple(compartment, ["mn", "a"], ["ma"], TCE["off_a"] / TCE["affn_a"], TCE["off_a"])
    system.add_simple(compartment, ["nn", "a"], ["na"], TCE["off_a"] / TCE["affn_a"], TCE["off_a"])
  
  
  # target binding
  for drug in drugs:
    off_C = TCE["off_C"]; on_C = TCE["off_C"] / TCE["affn_C"] if drug[0] == "n" else TCE["off_C"] / TCE["affm_C"]
    off_A = TCE["off_A"]; on_A = TCE["off_A"] / TCE["affn_A"] if drug[1] == "n" else TCE["off_A"] / TCE["affm_A"]
    off_B = TCE["off_B"]; on_B = TCE["off_B"] / TCE["aff_B"]
    avidity = TCE["avidity"]
    
    system.add_simple("plasma", ["C", f"{drug}"], [f"C-{drug}"], on_C, off_C)
    
    for organ in tumors + organs:
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}", "A"], [f"{drug}-A"], on_A, off_A)
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}", "B"], [f"{drug}-B"], on_B, off_B)
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * avidity, off_B)
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * avidity, off_A)
      
      system.add_simple(f"{organ['name']}_interstitial", ["C", f"{drug}"], [f"C-{drug}"], on_C, off_C)
      system.add_simple(f"{organ['name']}_interstitial", ["C", f"{drug}-A"], [f"C-{drug}-A"], on_C, off_C)
      system.add_simple(f"{organ['name']}_interstitial", ["C", f"{drug}-B"], [f"C-{drug}-B"], on_C, off_C)
      system.add_simple(f"{organ['name']}_interstitial", ["C", f"{drug}-AB"], [f"C-{drug}-AB"], on_C, off_C)
      
      system.add_simple(f"{organ['name']}_interstitial", [f"C-{drug}", "A"], [f"C-{drug}-A"], on_A, off_A)
      system.add_simple(f"{organ['name']}_interstitial", [f"C-{drug}", "B"], [f"C-{drug}-B"], on_B, off_B)
      system.add_simple(f"{organ['name']}_interstitial", [f"C-{drug}-A", "B"], [f"C-{drug}-AB"], on_B * avidity, off_B)
      system.add_simple(f"{organ['name']}_interstitial", [f"C-{drug}-B", "A"], [f"C-{drug}-AB"], on_A * avidity, off_A)
  
  
  # initial concentrations
  T_cell_density_blood = 1e6 / units.ml
  system.add_x("C", "plasma", 124000 * T_cell_density_blood / units.avagadro)
  
  for tumor in tumors:
    system.add_x("C", f"{tumor['name']}_interstitial", 124000 * tumor["T_cell_density"] / units.avagadro)
    system.add_x("A", f"{tumor['name']}_interstitial", tumor["num_A"] * tumor["cell_density"] / units.avagadro)
    system.add_x("B", f"{tumor['name']}_interstitial", tumor["num_B"] * tumor["cell_density"] / units.avagadro)
  
  for organ in organs:
    system.add_x("C", f"{organ['name']}_interstitial", 124000 * organ["T_cell_density"] / units.avagadro)
    system.add_x("A", f"{organ['name']}_interstitial", organ["num_A"] * organ["cell_density"] / units.avagadro)
    system.add_x("B", f"{organ['name']}_interstitial", organ["num_B"] * organ["cell_density"] / units.avagadro)
  
  return system
