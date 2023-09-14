from .qsp import *

### this model is mostly from ...

CD3s = ["CD3eff", "CD3reg"]
targets = ["A", "B"]
drugs = ["bimasked", "leftmasked", "rightmasked", "unmasked"]
C_conjugates = [f"{drug}-{Ag}" for drug in drugs for Ag in ["A", "B", "AA", "AB", "BB"]]
T_conjugates = [f"{CD3}-{drug}" for drug in drugs for CD3 in ["CD3eff", "CD3reg"]]
trimers = [f"{CD3}-{drug}-{Ag}" for drug in drugs for Ag in ["A", "B", "AA", "AB", "BB"] for CD3 in ["CD3eff", "CD3reg"]]
analytes = cells + targets + drugs + C_conjugates + T_conjugates + trimers


############ constants ############

molecular_weight = 150000 * units.g/units.mol
tumor_cell_density = 3e8 / units.ml

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
      self.analytes = [system.analytes.index(analyte) for analyte in ["bimasked", "leftmasked", "rightmasked", "unmasked"]]
      self.compartment = system.compartments.index("central")
    x = system.x[self.analytes, self.compartment]
    s = x.sum()
    if s == 0:
      return
    rate = - Hill(EMAX, EC50, 1, s) * (x/s)
    system.x[self.analytes, self.compartment] += rate * t.number(units.h)

mouse = {}
mouse.update({"volume_central": 1.26 * units.ml, "volume_peripheral": 0.819 * units.ml})
mouse.update({"distribution": 4.82 * units.ml/units.d})
mouse.update({"clearance": 0.334 * units.ml/units.d})
mouse.update({"nonlinear_clearance": nonlinear_clearance(0.518 * units.microgram/units.d / molecular_weight, 0.366 * units.microgram/units.ml / molecular_weight)})
mouse.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
mouse.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
mouse.update({"FcRn": 49.8 * units.micromolar, "FcRn-on": 0.0806 * 1/units.nM/units.h, "FcRn-off": 6.55 / units.h})

human = {}
human.update({"volume_central": 2877 * units.ml, "volume_peripheral": 2857 * units.ml})
human.update({"distribution": 384 * units.ml/units.d})
human.update({"clearance": 167 * units.ml/units.d})
human.update({"nonlinear_clearance": nonlinear_clearance(114 * units.microgram/units.d / molecular_weight, 0.078 * units.microgram/units.ml / molecular_weight)})
human.update({"vascular_reflection": 0.842, "lymphatic_reflection": 0.2})
human.update({"endosomal_pinocytosis": 0.0366 / units.h, "endosomal_degradation": 42.9 / units.h, "vascular_recycle": 0.715})
human.update({"FcRn": 49.8 * units.micromolar, "FcRn-on": 0.792 * 1/units.nM/units.h, "FcRn-off": 23.9 / units.h})


############ antigens ############

V1 = {}
V1.update({"on_CD3": 0.34 * 1/units.nM/units.d, "on_A": 0.34 * 1/units.nM/units.d, "on_B": 0.34 * 1/units.nM/units.d})
V1.update({"avidity": 1000})
V1.update({"off_CD3": 0.106 / units.h, "off_A": 0.106 / units.h, "off_B": 0.106 / units.h})
V1.update({"breath_CD3_l": 0.01, "breath_A_l": 1, "breath_B_l": 1})
V1.update({"breath_CD3_r": 1, "breath_A_r": 0.01, "breath_B_r": 0.01})
V1.update({"cleavage_plasma_l": 0.0527 / units.d, "cleavage_plasma_r": 0.0527 / units.d})
V1.update({"cleavage_tumor_l": 0.1783 / units.d, "cleavage_tumor_r": 0.1783 / units.d})

############ antigens ############

EGFR = {"num_tumor": None, "num_SI": None}
CAIX = {"num_tumor": None, "num_SI": None}


############ tumors ############

MC38 = {}
MC38.update({"volume": 170 * units.microliter, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
MC38.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
MC38.update({"capillary_radius": 10 * units.micrometer, "capillary_permeability": 3e-7 * units.cm/units.s})
MC38.update({"diffusion": 10 * units.micrometer**2 / units.s})


############ organs ############

SI = {"name": "SI"}
SI.update({"volume_plasma": 6.15 * units.ml, "volume_endosomal": 1.93 * units.ml, "volume_interstitial": 67.1 * units.ml})
SI.update({"plasma_flow": 12368 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})




############ model ############

def model(host, TCE, A, B, tumor, organs):
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
  system.add_simple("central", ["bimasked"], ["rightmasked"], TCE["cleavage_plasma_l"])
  system.add_simple("central", ["bimasked"], ["leftmasked"], TCE["cleavage_plasma_r"])
  system.add_simple("central", ["leftmasked"], ["unmasked"], TCE["cleavage_plasma_l"])
  system.add_simple("central", ["rightmasked"], ["unmasked"], TCE["cleavage_plasma_r"])
  
  system.add_simple("tumor_interstitial", ["bimasked"], ["rightmasked"], TCE["cleavage_tumor_l"])
  system.add_simple("tumor_interstitial", ["bimasked"], ["leftmasked"], TCE["cleavage_tumor_r"])
  system.add_simple("tumor_interstitial", ["leftmasked"], ["unmasked"], TCE["cleavage_tumor_l"])
  system.add_simple("tumor_interstitial", ["rightmasked"], ["unmasked"], TCE["cleavage_tumor_r"])
  
  
  # target binding
  for drug in ["bimasked", "leftmasked", "rightmasked", "unmasked"]:
    if drug == "bimasked":
      on_CD3 = TCE["on_A"] * TCE["breath_CD3_l"] * TCE["breath_CD3_r"]
      on_A = TCE["on_A"] * TCE["breath_A_l"] * TCE["breath_A_r"]
      on_B = TCE["on_B"] * TCE["breath_B_l"] * TCE["breath_B_r"]
    if drug == "leftmasked":
      on_CD3 = TCE["on_A"] * TCE["breath_CD3_l"]
      on_A = TCE["on_A"] * TCE["breath_A_l"]
      on_B = TCE["on_B"] * TCE["breath_B_l"]
    if drug == "rightmasked":
      on_CD3 = TCE["on_A"] * TCE["breath_CD3_r"]
      on_A = TCE["on_A"] * TCE["breath_A_r"]
      on_B = TCE["on_B"] * TCE["breath_B_r"]
    if drug == "unmasked":
      on_CD3 = TCE["on_A"]
      on_A = TCE["on_A"]
      on_B = TCE["on_B"]
    
    # in organs
    for organ in organs:
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}", "A"], [f"{drug}-A"], on_A, TCE["off_A"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}", "B"], [f"{drug}-B"], on_B, TCE["off_B"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-A", "A"], [f"{drug}-AA"], on_A * TCE["avidity"], TCE["off_A"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * TCE["avidity"], TCE["off_B"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * TCE["avidity"], TCE["off_A"])
      system.add_simple(f"{organ['name']}_interstitial", [f"{drug}-B", "B"], [f"{drug}-BB"], on_B * TCE["avidity"], TCE["off_B"])
    
    # in tumor
    system.add_simple("tumor_interstitial", [f"{drug}", "A"], [f"{drug}-A"], on_A, TCE["off_A"])
    system.add_simple("tumor_interstitial", [f"{drug}", "B"], [f"{drug}-B"], on_B, TCE["off_B"])
    system.add_simple("tumor_interstitial", [f"{drug}-A", "A"], [f"{drug}-AA"], on_A * TCE["avidity"], TCE["off_A"])
    system.add_simple("tumor_interstitial", [f"{drug}-A", "B"], [f"{drug}-AB"], on_B * TCE["avidity"], TCE["off_B"])
    system.add_simple("tumor_interstitial", [f"{drug}-B", "A"], [f"{drug}-AB"], on_A * TCE["avidity"], TCE["off_A"])
    system.add_simple("tumor_interstitial", [f"{drug}-B", "B"], [f"{drug}-BB"], on_B * TCE["avidity"], TCE["off_B"])
    
    for CD3 in ["CD3eff", "CD3reg"]:
      system.add_simple("tumor_interstitial", [f"{CD3}", f"{drug}"], [f"{CD3}-{drug}"], on_CD3, TCE["off_CD3"])
      system.add_simple("tumor_interstitial", [f"{CD3}", f"{drug}-A"], [f"{CD3}-{drug}-A"], on_CD3, TCE["off_CD3"])
      system.add_simple("tumor_interstitial", [f"{CD3}", f"{drug}-B"], [f"{CD3}-{drug}-B"], on_CD3, TCE["off_CD3"])
      system.add_simple("tumor_interstitial", [f"{CD3}", f"{drug}-AA"], [f"{CD3}-{drug}-AA"], on_CD3, TCE["off_CD3"])
      system.add_simple("tumor_interstitial", [f"{CD3}", f"{drug}-AB"], [f"{CD3}-{drug}-AB"], on_CD3, TCE["off_CD3"])
      system.add_simple("tumor_interstitial", [f"{CD3}", f"{drug}-BB"], [f"{CD3}-{drug}-BB"], on_CD3, TCE["off_CD3"])
      
      system.add_simple("tumor_interstitial", [f"{CD3}-{drug}", "A"], [f"{CD3}-{drug}-A"], on_A, TCE["off_A"])
      system.add_simple("tumor_interstitial", [f"{CD3}-{drug}", "B"], [f"{CD3}-{drug}-B"], on_B, TCE["off_B"])
      system.add_simple("tumor_interstitial", [f"{CD3}-{drug}-A", "A"], [f"{CD3}-{drug}-AA"], on_A * TCE["avidity"], TCE["off_A"])
      system.add_simple("tumor_interstitial", [f"{CD3}-{drug}-A", "B"], [f"{CD3}-{drug}-AB"], on_B * TCE["avidity"], TCE["off_B"])
      system.add_simple("tumor_interstitial", [f"{CD3}-{drug}-B", "A"], [f"{CD3}-{drug}-AB"], on_A * TCE["avidity"], TCE["off_A"])
      system.add_simple("tumor_interstitial", [f"{CD3}-{drug}-B", "B"], [f"{CD3}-{drug}-BB"], on_B * TCE["avidity"], TCE["off_B"])
  
  return system
