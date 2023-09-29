from qsp import *
from qsp.model_VIB4 import *

double = FTC238.copy()
double["name"] = "double"

single = FTC238.copy()
single["num_B"] = 0
single["name"] = "single"

def plot(system, name):
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
            [f"n{a}{b}" for a in ["m","n"] for b in ["m","n"]] + [f"n{a}{b}-{target}" for a in ["m","n"] for b in ["m","n"] for target in targets] + [f"C-n{a}{b}" for a in ["m","n"] for b in ["m","n"]] + [f"C-n{a}{b}-{target}" for a in ["m","n"] for b in ["m","n"] for target in targets],
            [f"{c}n{b}" for c in ["m","n"] for b in ["m","n"]] + [f"{c}n{b}-{target}" for c in ["m","n"] for b in ["m","n"] for target in targets] + [f"C-{c}n{b}" for c in ["m","n"] for b in ["m","n"]] + [f"C-{c}n{b}-{target}" for c in ["m","n"] for b in ["m","n"] for target in targets],
            ["a"]]
  labels = ["(C-)xxx(-target)", "(C-)nx(-target)", "(C-)xn(-target)", "cap"]
  colors = ["tab:green", "wheat", "skyblue", "salmon", "black"]
  linestyles = ["solid", "dashed", "dashed", "dashed", "dotted"]
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


results = pd.DataFrame(columns = ["aff_a", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for aff_a in [1e-10, 1e-9, 1e-8, 1e-7]:
  TCE = VIB4.copy()
  TCE["aff_a"] = aff_a * units.molar
  
  system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  plot(system, f"{aff_a}")
  
  summary = system.summary(trimers)
  results.loc[results.shape[0]] = np.array([aff_a, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])
