from qsp import *
from qsp.model_VIB4 import *
import pickle

double = FTC238.copy()
double["name"] = "double"

single = FTC238.copy()
single["num_B"] = 0
single["name"] = "single"

def plot(system, name):
  #pickle.dump(system, open(f"{name}.pickle", "wb"))
  
  targets = ["A", "B", "AB"]
  groups = [["C"],
            ["A", "B"],
            ["a"],
            drugs,
            [f"C-{drug}" for drug in drugs],
            [f"{drug}-{target}" for drug in drugs for target in targets],
            [f"C-{drug}-{target}" for drug in drugs for target in targets]]
  labels = ["C", "target", "cap", "drug", "C-drug", "drug-target", "C-drug-target"]
  colors = ["tab:orange", "tab:blue", "black", "black", "wheat", "skyblue", "tab:purple"]
  linestyles = ["solid", "solid", "dotted", "solid", "solid", "solid", "solid"]
  system.plot(compartments = ["plasma"] + [f"{tumor['name']}_interstitial" for tumor in system.tumors] + [f"{organ['name']}_interstitial" for organ in system.organs], 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")
  
  groups = [drugs + [f"{drug}-{target}" for drug in drugs for target in targets] + [f"C-{drug}" for drug in drugs] + [f"C-{drug}-{target}" for drug in drugs for target in targets],
            [f"n{a}" for a in ["m","n"]] + [f"n{a}-{target}" for a in ["m","n"] for target in targets] + [f"C-n{a}" for a in ["m","n"]] + [f"C-n{a}-{target}" for a in ["m","n"] for target in targets],
            [f"{c}n" for c in ["m","n"]] + [f"{c}n-{target}" for c in ["m","n"] for target in targets] + [f"C-{c}n" for c in ["m","n"]] + [f"C-{c}n-{target}" for c in ["m","n"] for target in targets],
            ["a"]]
  labels = ["(C-)xxx(-target)", "(C-)nx(-target)", "(C-)xn(-target)", "cap"]
  colors = ["tab:green", "wheat", "skyblue", "black"]
  linestyles = ["solid", "dashed", "dashed", "dotted"]
  system.plot(compartments = ["plasma"] + [f"{tumor['name']}_interstitial" for tumor in system.tumors] + [f"{organ['name']}_interstitial" for organ in system.organs], 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_drugs.png")
  
  groups = [["A"], ["B"], drugs,
            [f"{drug}-A" for drug in drugs] + [f"C-{drug}-A" for drug in drugs],
            [f"{drug}-B" for drug in drugs] + [f"C-{drug}-B" for drug in drugs],
            [f"{drug}-AB" for drug in drugs] + [f"C-{drug}-AB" for drug in drugs]]
  labels = ["A", "B", "(C-)drug", "(C-)drug-A", "(C-)drug-B", "(C-)drug-AB"]
  colors = ["tab:blue", "tab:red", "black", "skyblue", "salmon", "tab:purple"]
  system.plot(compartments = ["plasma"] + [f"{tumor['name']}_interstitial" for tumor in system.tumors] + [f"{organ['name']}_interstitial" for organ in system.organs],
              groups = groups, labels = labels, colors = colors,
              output = f"{name}_targets.png")


params = [[0.05, 0.15, 0.05, 0.15, x] for x in [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]] + \
         [[x, 0.15, 0.05, 0.15, 1e-9] for x in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]] + \
         [[0.05, x, 0.05, 0.15, 1e-9] for x in [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]] + \
         [[0.05, 0.15, x, 0.15, 1e-9] for x in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]] + \
         [[0.05, 0.15, 0.05, x, 1e-9] for x in [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]] + \
         [[0.05*x, 0.15*x, 0.05, 0.15, 1e-9] for x in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2]] + \
         [[0.05, 0.15, 0.05*x, 0.15*x, 1e-9] for x in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2]]

results = pd.DataFrame(columns = ["rate_C_plasma", "rate_C_tumor", "rate_A_plasma", "rate_A_tumor", "aff_a", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for rate_C_plasma, rate_C_tumor, rate_A_plasma, rate_A_tumor, aff_a in params:
  TCE = VIB4.copy()
  TCE["aff_a"] = aff_a * units.molar
  TCE["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [f"{organ['name']}_interstitial" for organ in system.organs], 
                                  rate_C = rate_C_plasma / units.d, 
                                  rate_A = rate_A_plasma / units.d)
  TCE["cleavage_tumor"] = cleavage(lambda system: [f"{tumor['name']}_interstitial" for tumor in system.tumors], 
                                 rate_C = rate_C_tumor / units.d, 
                                 rate_A = rate_A_tumor / units.d)
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([rate_C_plasma, rate_C_tumor, rate_A_plasma, rate_A_tumor, aff_a, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])
  print(results)

results.to_csv("results.csv")



results = pd.DataFrame(columns = ["aff_a", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for aff_a in [1e-12, 1e-11]:
  TCE = VIB4.copy()
  TCE["aff_a"] = aff_a * units.molar
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([aff_a, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])

print(results)

results = pd.DataFrame(columns = ["rate_C_plasma", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for rate_C_plasma in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
  TCE = VIB4.copy()
  TCE["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [f"{organ['name']}_interstitial" for organ in system.organs], 
                                  rate_C = rate_C_plasma / units.d, 
                                  rate_A = 0.0527 / units.d)
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([rate_C_plasma, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])

print(results)

results = pd.DataFrame(columns = ["rate_C_tumor", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for rate_C_tumor in [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]:
  TCE = VIB4.copy()
  TCE["cleavage_tumor"] = cleavage(lambda system: [f"{tumor['name']}_interstitial" for tumor in system.tumors], 
                                 rate_C = rate_C_tumor / units.d, 
                                 rate_A = 0.1783 / units.d)
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([rate_C_tumor, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])

print(results)

results = pd.DataFrame(columns = ["rate_A_plasma", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for rate_A_plasma in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
  TCE = VIB4.copy()
  TCE["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [f"{organ['name']}_interstitial" for organ in system.organs], 
                                  rate_C = 0.0527 / units.d, 
                                  rate_A = rate_A_plasma / units.d)
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([rate_A_plasma, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])

print(results)

results = pd.DataFrame(columns = ["rate_A_tumor", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for rate_A_tumor in [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]:
  TCE = VIB4.copy()
  TCE["cleavage_tumor"] = cleavage(lambda system: [f"{tumor['name']}_interstitial" for tumor in system.tumors], 
                                 rate_C = 0.1783 / units.d, 
                                 rate_A = rate_A_tumor / units.d)
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([rate_A_tumor, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])

print(results)
