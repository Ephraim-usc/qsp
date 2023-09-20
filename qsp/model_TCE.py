from .qsp import *

### this model is mostly from ...

CD3s = ["CD3eff", "CD3reg"]
targets = ["A", "B"]
drugs = ["mm", "mn", "nm", "nn"]
C_conjugates = [f"{drug}-{Ag}" for drug in drugs for Ag in ["A", "B", "AB"]]
T_conjugates = [f"{CD3}-{drug}" for drug in drugs for CD3 in ["CD3eff", "CD3reg"]]
trimers = [f"{CD3}-{drug}-{Ag}" for drug in drugs for Ag in ["A", "B", "AB"] for CD3 in ["CD3eff", "CD3reg"]]
analytes = CD3s + targets + drugs + C_conjugates + T_conjugates + trimers


############ constants ############

molecular_weight = 150000 * units.g/units.mol

def Hill(EMAX, EC50, coefficient, x):
  if coefficient == 1:
    return EMAX/(1 + (x/EC50))
  else:
    return EMAX/(1 + (x/EC50) ** coefficient)


############ host ############

class nonlinear_clearance:
  def __init__(self, EMAX, EC50):
    self.system = None
    self.EMAX = EMAX.number(units.nM*units.ml/ units.h)
    self.EC50 = EC50.number(units.nM)
  def __call__(self, system, t):
    if self.system is not system:
      self.system = system
      self.analytes = [system.analytes.index(analyte) for analyte in ["mm", "mn", "nm", "nn"]]
      self.compartment = system.compartments.index("central")
    x = system.x[self.analytes, self.compartment]
    s = x.sum()
    if s == 0:
      return
    rate = - Hill(self.EMAX, self.EC50, 1, s) * (x/s)
    system.x[self.analytes, self.compartment] += rate * t.number(units.h)

mouse = {}
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"distribution": 4.82 * units.ml/units.d})
mouse.update({"clearance": 0.334 * units.ml/units.d})
mouse.update({"nonlinear_clearance": nonlinear_clearance(0.518 * units.microgram/units.d / molecular_weight, 0.366 * units.microgram/units.ml / molecular_weight)})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})

human = {}
human.update({"volume_central": 2877 * units.ml, "volume_peripheral": 2857 * units.ml})
human.update({"distribution": 384 * units.ml/units.d})
human.update({"clearance": 167 * units.ml/units.d})
human.update({"nonlinear_clearance": nonlinear_clearance(114 * units.microgram/units.d / molecular_weight, 0.078 * units.microgram/units.ml / molecular_weight)})
human.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})


############ drug ############

R72 = {}
R72.update({"off_CD3": 0.106 / units.h, "off_A": 0.106 / units.h, "off_B": 0.106 / units.h})
R72.update({"on_CD3": R72["off_CD3"] / 90 * units.nM, 
            "on_A": R72["off_A"] / 203 * units.nM, 
            "on_B": R72["off_B"] / 1.07 * units.nM})
R72.update({"avidity": 1000})
R72.update({"breath_CD3": 90/846, "breath_A": 203/916, "breath_B": 1})
R72.update({"cleavage_plasma_CD3": 0.0527 / units.d, "cleavage_plasma_A": 0.0527 / units.d, "cleavage_plasma_B": 0.0527 / units.d})
R72.update({"cleavage_tumor_CD3": 0.1783 / units.d, "cleavage_tumor_A": 0.1783 / units.d, "cleavage_tumor_B": 0.1783 / units.d})

R77 = {}
R77.update({"off_CD3": 0.106 / units.h, "off_A": 0.106 / units.h, "off_B": 0.106 / units.h})
R77.update({"on_CD3": R72["off_CD3"] / 25 * units.nM, 
            "on_A": R72["off_A"] / 11 * units.nM, 
            "on_B": R72["off_B"] / 189 * units.nM})
R77.update({"avidity": 1000})
R77.update({"breath_CD3": 25/527, "breath_A": 11/243, "breath_B": 1})
R77.update({"cleavage_plasma_CD3": 0.0527 / units.d, "cleavage_plasma_A": 0.0527 / units.d, "cleavage_plasma_B": 0.0527 / units.d})
R77.update({"cleavage_tumor_CD3": 0.1783 / units.d, "cleavage_tumor_A": 0.1783 / units.d, "cleavage_tumor_B": 0.1783 / units.d})


############ tumors ############

MC38 = {"name": "tumor"}
MC38.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
MC38.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
MC38.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
MC38.update({"diffusion": 10 * units.micrometer**2 / units.s})
MC38.update({"cell_density": 3e8 / units.ml})
MC38.update({"num_A": 3000, "num_B": 1000000})


############ organs ############

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_endosomal": 1.93 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"cell_density": 1e8 / units.ml})
SI.update({"num_A": 1000, "num_B": 100000})


############ model ############

def model(host, TCE, tumor, organs):
  compartments = ["central", "peripheral", "tumor_plasma", "tumor_interstitial"] + [f"{organ['name']}_{tissue}" for organ in organs for tissue in ["plasma", "interstitial"]]
  system = System(analytes, compartments)
  
  for analyte in analytes:
    system.set_volume(analyte, "central", host["volume_central"])
    system.set_volume(analyte, "peripheral", host["volume_peripheral"])
    system.set_volume(analyte, "tumor_plasma", tumor["volume"] * tumor["volume_plasma_proportion"])
    system.set_volume(analyte, "tumor_interstitial", tumor["volume"] * tumor["volume_interstitial_proportion"])
    for organ in organs:
      system.set_volume(analyte, f"{organ['name']}_plasma", organ["volume_plasma"])
      system.set_volume(analyte, f"{organ['name']}_interstitial", organ["volume_interstitial"])
  
  # distribution and clearance
  for drug in drugs:
    system.add_flow(drug, "central", "peripheral", host["distribution"])
    system.add_flow(drug, "peripheral", "central", host["distribution"])
    system.add_flow(drug, "central", None, host["clearance"])
  system.add_process(host["nonlinear_clearance"])
  
  # organ flow
  for drug in drugs:
    for organ in organs:
      system.add_flow(drug, "central", f"{organ['name']}_plasma", organ["plasma_flow"])
      system.add_flow(drug, f"{organ['name']}_plasma", "central", organ["plasma_flow"])
      system.add_flow(drug, "central", f"{organ['name']}_plasma", organ["plasma_flow"] * organ["lymphatic_flow_ratio"])
      system.add_flow(drug, f"{organ['name']}_plasma", f"{organ['name']}_interstitial", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - host["vascular_reflection"]))
      system.add_flow(drug, f"{organ['name']}_interstitial", "central", organ["plasma_flow"] * organ["lymphatic_flow_ratio"] * (1 - host["lymphatic_reflection"]))
  
  # tumor flow
  for drug in drugs:
    system.add_flow(drug, "central", "tumor_plasma", tumor["volume"] * tumor["plasma_flow_density"])
    system.add_flow(drug, "tumor_plasma", "central", tumor["volume"] * tumor["plasma_flow_density"])
    system.add_flow(drug, "tumor_plasma", "tumor_interstitial", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
    system.add_flow(drug, "tumor_interstitial", "tumor_plasma", tumor["volume"] * tumor["volume_plasma_proportion"] * (2 / tumor["capillary_radius"]) * tumor["capillary_permeability"])
  
  
  # mask cleavage
  system.add_simple("central", ["mm"], ["nm"], TCE["cleavage_plasma_CD3"])
  system.add_simple("central", ["mm"], ["mn"], TCE["cleavage_plasma_A"])
  system.add_simple("central", ["mn"], ["nn"], TCE["cleavage_plasma_CD3"])
  system.add_simple("central", ["nm"], ["nn"], TCE["cleavage_plasma_A"])
  
  system.add_simple("tumor_interstitial", ["mm"], ["nm"], TCE["cleavage_tumor_CD3"])
  system.add_simple("tumor_interstitial", ["mm"], ["mn"], TCE["cleavage_tumor_A"])
  system.add_simple("tumor_interstitial", ["mn"], ["nn"], TCE["cleavage_tumor_CD3"])
  system.add_simple("tumor_interstitial", ["nm"], ["nn"], TCE["cleavage_tumor_A"])
  
  
  # target binding
  for drug in drugs:
    on_CD3, on_A, on_B = TCE["on_CD3"], TCE["on_A"], TCE["on_B]
    if drug[0] == "m":
      on_CD3 *= TCE["breath_CD3"]
    if drug[1] == "m":
      on_A *= TCE["breath_A"]
    
    for CD3 in ["CD3eff", "CD3reg"]:
      system.add_simple("central", [f"{CD3}", f"{drug}"], [f"{CD3}-{drug}"], on_CD3, TCE["off_CD3"])
    
    for organ in [tumor] + organs:
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}", "A"], [f"{drug}-A"], on_A, TCE["off_A"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}", "B"], [f"{drug}-B"], on_B, TCE["off_B"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * TCE["avidity"], TCE["off_B"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * TCE["avidity"], TCE["off_A"])
      
      for CD3 in ["CD3eff", "CD3reg"]:
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}", f"{drug}"], [f"{CD3}-{drug}"], on_CD3, TCE["off_CD3"])
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}", f"{drug}-A"], [f"{CD3}-{drug}-A"], on_CD3, TCE["off_CD3"])
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}", f"{drug}-B"], [f"{CD3}-{drug}-B"], on_CD3, TCE["off_CD3"])
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}", f"{drug}-AB"], [f"{CD3}-{drug}-AB"], on_CD3, TCE["off_CD3"])
        
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}-{drug}", "A"], [f"{CD3}-{drug}-A"], on_A, TCE["off_A"])
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}-{drug}", "B"], [f"{CD3}-{drug}-B"], on_B, TCE["off_B"])
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}-{drug}-A", "B"], [f"{CD3}-{drug}-AB"], on_B * TCE["avidity"], TCE["off_B"])
        system.add_simple(f"{organ['name']}_interstitial", [f"{CD3}-{drug}-B", "A"], [f"{CD3}-{drug}-AB"], on_A * TCE["avidity"], TCE["off_A"])
  
  
  # initial concentrations
  Treg_density_blood = 138070411 / (6 * units.l)
  Treg_density_lymph = 107544935 / (0.021 * units.l)
  Treg_density_peripheral = 2.317e11 / (62.24 * units.l)
  Treg_density_tumor = 9303338 / (1 * units.l)
  
  system.add_x("CD3eff", "central", 124000 * (Treg_density_blood*0.6) / units.avagadro)
  system.add_x("CD3reg", "central", 124000 * Treg_density_blood / units.avagadro)
  
  system.add_x("CD3eff", "peripheral", 124000 * (Treg_density_peripheral*0.6) / units.avagadro)
  system.add_x("CD3reg", "peripheral", 124000 * Treg_density_peripheral / units.avagadro)
  
  system.add_x("CD3eff", "tumor_interstitial", 124000 * (4.3 * 400/units.microliter) / units.avagadro)
  system.add_x("CD3reg", "tumor_interstitial", 124000 * (4.3 * 600/units.microliter) / units.avagadro)
  system.add_x("A", "tumor_interstitial", tumor["num_A"] * tumor["cell_density"] / units.avagadro)
  system.add_x("B", "tumor_interstitial", tumor["num_B"] * tumor["cell_density"] / units.avagadro)
  for organ in organs:
    system.add_x("CD3eff", f"{organ['name']}_interstitial", 124000 * (400/units.microliter) / units.avagadro)
    system.add_x("CD3reg", f"{organ['name']}_interstitial", 124000 * (600/units.microliter) / units.avagadro)
    system.add_x("A", f"{organ['name']}_interstitial", organ["num_A"] * organ["cell_density"] / units.avagadro)
    system.add_x("B", f"{organ['name']}_interstitial", organ["num_B"] * organ["cell_density"] / units.avagadro)
  
  return system
