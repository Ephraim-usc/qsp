from qsp import *
from qsp.model_VIB4 import *
from qsp.optimizer import *

tumor_AB = FTC238.copy()
tumor_AB["name"] = "tumor_AB"

tumor_A = FTC238.copy()
tumor_A["num_B"] = 0
tumor_A["name"] = "tumor_A"

tumor_B = FTC238.copy()
tumor_B["num_A"] = 0
tumor_B["name"] = "tumor_B"


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
  system.plot(compartments = ["plasma"] + [tumor["name"] for tumor in system.tumors] + [organ["name"] for organ in system.organs], 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_summary.png")
  
  groups = [drugs + [f"{drug}-{target}" for drug in drugs for target in targets] + [f"C-{drug}" for drug in drugs] + [f"C-{drug}-{target}" for drug in drugs for target in targets],
            [f"n{a}" for a in ["m","n"]] + [f"n{a}-{target}" for a in ["m","n"] for target in targets] + [f"C-n{a}" for a in ["m","n"]] + [f"C-n{a}-{target}" for a in ["m","n"] for target in targets],
            [f"{c}n" for c in ["m","n"]] + [f"{c}n-{target}" for c in ["m","n"] for target in targets] + [f"C-{c}n" for c in ["m","n"]] + [f"C-{c}n-{target}" for c in ["m","n"] for target in targets],
            [f"{c}a" for c in ["m","n"]] + [f"{c}a-{target}" for c in ["m","n"] for target in targets] + [f"C-{c}a" for c in ["m","n"]] + [f"C-{c}a-{target}" for c in ["m","n"] for target in targets]]
  labels = ["(C-)xxx(-target)", "(C-)nx(-target)", "(C-)xn(-target)", "(C-)xa(-target)"]
  colors = ["tab:green", "wheat", "skyblue", "salmon"]
  linestyles = ["solid", "dashed", "dashed", "dashed"]
  system.plot(compartments = ["plasma"] + [tumor["name"] for tumor in system.tumors] + [organ["name"] for organ in system.organs], 
              groups = groups, labels = labels, colors = colors, linestyles = linestyles,
              output = f"{name}_drugs.png")


results = pd.DataFrame(columns = ["affn_C", "affm_C", "tumor_AB", "tumor_A", "tumor_B", "lung", "SI"])
for affn_C in 10**np.arange(-10, -6, 0.4):
  for affm_C in 10**np.arange(-10, -6, 0.4):
    TCE = VIBY.copy()
    TCE["affn_C"] = affn_C * units.molar; TCE["affm_C"] = affm_C * units.molar
    system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
    system.add_x("mm", "plasma", 1 * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
    summary = system.summary(trimers)["average"]
    results.loc[results.shape[0]] = np.array([affn_C, affm_C, summary.loc["tumor_AB"], summary.loc["tumor_A"], summary.loc["tumor_B"], summary.loc["lung"], summary.loc["SI"]])
    results.to_csv("aff_C_map.csv")


results = pd.DataFrame(columns = ["affn_A", "affn_B", "tumor_AB", "tumor_A", "tumor_B", "lung", "SI"])
for affn_A in 10**np.arange(-10, -6, 0.4):
  for affn_B in 10**np.arange(-10, -6, 0.4):
    TCE = VIBY.copy()
    TCE["affn_A"] = affn_A * units.molar; TCE["affm_A"] = affn_A*100 * units.molar
    TCE["affn_B"] = affn_B * units.molar; TCE["affm_B"] = affn_B*100 * units.molar
    system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
    system.add_x("mm", "plasma", 1 * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
    summary = system.summary(trimers)["average"]
    results.loc[results.shape[0]] = np.array([affn_A, affn_B, summary.loc["tumor_AB"], summary.loc["tumor_A"], summary.loc["tumor_B"], summary.loc["lung"], summary.loc["SI"]])
    results.to_csv("aff_A_B_map.csv")


results = pd.DataFrame(columns = ["off_a", "aff_a", "tumor_AB", "tumor_A", "tumor_B", "lung", "SI"])
for off_a in 10**np.arange(-5, -2, 0.3):
  for aff_a in 10**np.arange(-10, -6, 0.4):
    TCE = VIBY.copy()
    TCE["off_a"] = off_a / units.s; TCE["aff_a"] = aff_a * units.molar
    system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
    system.add_x("mm", "plasma", 1 * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
    summary = system.summary(trimers)["average"]
    results.loc[results.shape[0]] = np.array([off_a, aff_a, summary.loc["tumor_AB"], summary.loc["tumor_A"], summary.loc["tumor_B"], summary.loc["lung"], summary.loc["SI"]])
    results.to_csv("aff_a_map.csv")


results = pd.DataFrame(columns = ["avidity", "tumor_AB", "tumor_A", "tumor_B", "lung", "SI"])
for avidity in [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]:
    TCE = VIBY.copy()
    TCE["avidity"] = avidity
    system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
    system.add_x("mm", "plasma", 1 * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
    summary = system.summary(trimers)["average"]
    results.loc[results.shape[0]] = np.array([avidity, summary.loc["tumor_AB"], summary.loc["tumor_A"], summary.loc["tumor_B"], summary.loc["lung"], summary.loc["SI"]])
    results.to_csv("avidity.csv")


results = pd.DataFrame(columns = ["num_A", "num_B", "tumor_AB", "tumor_A", "tumor_B", "lung", "SI"])
for num_A in 10**np.arange(3, 7, 0.4):
  for num_B in 10**np.arange(3, 7, 0.4):
    TCE = VIBY.copy()
    tumor_AB_ = tumor_AB.copy(); tumor_AB_["num_A"] = num_A; tumor_AB_["num_B"] = num_B
    tumor_A_ = tumor_A.copy(); tumor_A_["num_A"] = num_A
    tumor_B_ = tumor_B.copy(); tumor_B_["num_B"] = num_B
    system = model(human, TCE, [tumor_AB_, tumor_A_, tumor_B_], [other, lung, SI], connect_tumors = True)
    system.add_x("mm", "plasma", 1 * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
    summary = system.summary(trimers)["average"]
    results.loc[results.shape[0]] = np.array([num_A, num_B, summary.loc["tumor_AB"], summary.loc["tumor_A"], summary.loc["tumor_B"], summary.loc["lung"], summary.loc["SI"]])
    results.to_csv("num_A_B.csv")


results = pd.DataFrame(columns = ["volume_AB", "tumor_AB", "tumor_A", "tumor_B", "lung", "SI"])
for volume_AB in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    TCE = VIBY.copy()
    tumor_AB_ = tumor_AB.copy(); tumor_AB_["volume"] *= volume_AB
    system = model(human, TCE, [tumor_AB_, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
    system.add_x("mm", "plasma", 1 * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
    summary = system.summary(trimers)["average"]
    results.loc[results.shape[0]] = np.array([volume_AB, summary.loc["tumor_AB"], summary.loc["tumor_A"], summary.loc["tumor_B"], summary.loc["lung"], summary.loc["SI"]])
    results.to_csv("volume_AB.csv")











def ratio(affn_C, affn_A, aff_B, off_C, off_A, off_B, halflife, cleavage_C_plasma, cleavage_A_plasma):
  TCE = BEST.copy()
  TCE["affn_C"] = (10**affn_C) * units.molar; TCE["affm_C"] = (10**affn_C) * 100 * units.molar; TCE["off_C"] = (10**off_C) / units.s
  TCE["affn_A"] = (10**affn_A) * units.molar; TCE["affm_A"] = (10**affn_A) * 100 * units.molar; TCE["off_A"] = (10**off_A) / units.s
  TCE["aff_B"] = (10**aff_B) * units.molar; TCE["off_B"] = (10**off_B) / units.s
  TCE.update({"clearance": math.log(2)/(halflife * units.h)})
  TCE["cleavage_plasma"] = cleavage(lambda system: ["plasma"] + [organ["name"] for organ in system.organs], 
                                  rate_C = cleavage_C_plasma / units.d, 
                                  rate_A = cleavage_A_plasma / units.d)
  TCE["cleavage_tumor"] = cleavage(lambda system: [tumor["name"] for tumor in system.tumors], 
                                 rate_C = cleavage_C_plasma*3 / units.d, 
                                 rate_A = cleavage_A_plasma*3 / units.d)
  
  system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 1 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  
  plot(system, f"{affn_C}_{affn_A}_{aff_B}_{off_C}_{off_A}_{off_B}_{halflife}_{cleavage_C_plasma}_{cleavage_A_plasma}")
  summary = system.summary(trimers)["average"]
  return summary["tumor_A"] / summary["lung"]


names = ["affn_C", "affn_A", "aff_B", "off_C", "off_A", "off_B", "halflife", "cleavage_C_plasma", "cleavage_A_plasma"]
values = [-8, -9, -8, -3, -3, -3, 40, 0.05, 0.05]
limits = [(-10, -7), (-10, -7), (-10, -7), (-5, -2), (-5, -2), (-5, -2), (10, 200), (0.01, 0.1), (0.01, 0.1)]
search = Search(names, values, limits)
maximize(ratio, search)






X = VIB4.copy(); Y = VIB4.copy()
X.update({"off_A": 3e-4 / units.s, "affn_A": 5e-10 * units.molar, "affm_A": 1e-8 * units.molar})
Y.update({"off_A": 7e-3 / units.s, "affn_A": 2.5e-8 * units.molar, "affm_A": 5e-7 * units.molar})

for interrnalization_tumor in [0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2]:
  TCE = VIB4.copy()
  TCE["internalization_tumor"] = interrnalization_tumor / units.h
  system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 10 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  plot(system, f"VIB4_internalization={interrnalization_tumor}")

for interrnalization_tumor in [0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2]:
  TCE = X.copy()
  TCE["internalization_tumor"] = interrnalization_tumor / units.h
  system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 10 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  plot(system, f"X_internalization={interrnalization_tumor}")

for interrnalization_tumor in [0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2]:
  TCE = Y.copy()
  TCE["internalization_tumor"] = interrnalization_tumor / units.h
  system = model(human, TCE, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
  system.add_x("mm", "plasma", 10 * units.nM)
  system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
  plot(system, f"Y_internalization={interrnalization_tumor}")













params = [[0.05, 0.15, 0.05, 0.15, x] for x in [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]] + \
         [[x, 0.15, 0.05, 0.15, 1e-9] for x in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]] + \
         [[0.05, x, 0.05, 0.15, 1e-9] for x in [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]] + \
         [[0.05, 0.15, x, 0.15, 1e-9] for x in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]] + \
         [[0.05, 0.15, 0.05, x, 1e-9] for x in [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]] + \
         [[0.05*x, 0.15*x, 0.05, 0.15, 1e-9] for x in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2]] + \
         [[0.05, 0.15, 0.05*x, 0.15*x, 1e-9] for x in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2]]

params = [[0.05, 0.15, 0.05, 0.15, 1e-9, x] for x in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]]

results = pd.DataFrame(columns = ["rate_C_plasma", "rate_C_tumor", "rate_A_plasma", "rate_A_tumor", "aff_a", "half_life", "avg_trimer_double", "avg_trimer_single", "avg_trimer_lung"])
for rate_C_plasma, rate_C_tumor, rate_A_plasma, rate_A_tumor, aff_a, half_life in params:
  TCE = VIB4.copy()
  TCE["aff_a"] = aff_a * units.molar
  TCE["clearance"] = math.log(2)/(40 * units.h)
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
  results.loc[results.shape[0]] = np.array([rate_C_plasma, rate_C_tumor, rate_A_plasma, rate_A_tumor, aff_a, half_life, summary.loc["double_interstitial", "average"], summary.loc["single_interstitial", "average"], summary.loc["lung_interstitial", "average"]])
  print(results)

results.to_csv("results.csv")




system = model(human, VIB4, [tumor_AB, tumor_A, tumor_B], [other, lung, SI], connect_tumors = True)
system.add_x("mm", "plasma", 1 * units.nM)
system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "VIB4_1nM")

system = model(human, VIB1, [double, single], [other, lung, SI], connect_tumors = True)
system.add_x("mm", "plasma", 1 * units.nM)
system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "VIB1_1nM")

system = model(human, VIB5, [double, single], [other, lung, SI], connect_tumors = True)
system.add_x("mm", "plasma", 1 * units.nM)
system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "VIB5_1nM")

TCE = JANUX.copy()
TCE["clearance"] = {"mm": math.log(2)/(100 * units.h),
                    "mn": math.log(2)/(100 * units.h),
                    "nm": math.log(2)/(0.25 * units.h),
                    "nn": math.log(2)/(0.25 * units.h),
                    "a": math.log(2)/(0.25 * units.h)}
system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
system.add_x("mm", "plasma", 0.01 * units.nM)
system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "JANUX_100h_0.25h_0.01nM")

TCE = JANUX.copy()
TCE["clearance"] = {"mm": math.log(2)/(100 * units.h),
                    "mn": math.log(2)/(100 * units.h),
                    "nm": math.log(2)/(100 * units.h),
                    "nn": math.log(2)/(100 * units.h),
                    "a": math.log(2)/(0.25 * units.h)}
system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
system.add_x("mm", "plasma", 0.01 * units.nM)
system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "JANUX_100h_100h_0.01nM")

TCE = JANUX.copy()
TCE["clearance"] = {"mm": math.log(2)/(40 * units.h),
                    "mn": math.log(2)/(40 * units.h),
                    "nm": math.log(2)/(40 * units.h),
                    "nn": math.log(2)/(40 * units.h),
                    "a": math.log(2)/(0.25 * units.h)}
system = model(human, TCE, [double, single], [other, lung, SI], connect_tumors = True)
system.add_x("mm", "plasma", 0.01 * units.nM)
system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)
plot(system, "JANUX_40h_40h_0.01nM")
