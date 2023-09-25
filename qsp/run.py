from qsp import *
from qsp.model_TCE import *

pd.option_context(
                  'display.max_rows', None,
                  'display.max_columns', None,
                  'display.precision', 6,
                 )

for dose in [0.1, 1, 10, 100]:
  for construct, construct_name in zip([R72, R77], ["R72", "R77"]):
    system = model(human, construct, FTC238, [other, lung, SI])
    system.add_x("mmm", "plasma", dose * units.nM)
    system.run(300 * units.h, t_step = 1/60 * units.h, t_record = 1 * units.h)

    groups = [["C"],
              ["A", "B"],
              drugs,
              [f"C-{drug}" for drug in drugs],
              [f"{drug}-{target}" for drug in drugs for target in ["A", "B", "AB"]],
              [f"C-{drug}-{target}" for drug in drugs for target in ["A", "B", "AB"]]]
    labels = ["C", "target", "drug", "C-drug", "drug-target", "C-drug-target"]
    colors = ["tab:orange", "tab:blue", "black", "wheat", "skyblue", "tab:purple"]
    system.plot(compartments = ["plasma", "tumor_interstitial", "other_interstitial", "lung_interstitial", "SI_interstitial"],
                groups = groups, labels = labels, colors = colors,
                output = f"{dose}_{construct_name}_summary.png")
    
    '''
    groups = [[f"{drug}"] +
              [f"{CD3}-{drug}" for CD3 in ["CD3eff", "CD3reg"]] +
              [f"{drug}-{target}" for target in ["A", "B", "AB"]] +
              [f"{CD3}-{drug}-{target}" for CD3 in ["CD3eff", "CD3reg"] for target in ["A", "B", "AB"]]
              for drug in drugs]
    labels = [f"(CD3-){drug}(-target)" for drug in drugs]
    colors = ["tab:green"] * 4
    linestyles = ["solid", (0, (1, 1)), "dashed", "dotted"]
    system.plot(compartments = ["plasma", "tumor_interstitial", "other_interstitial", "lung_interstitial", "SI_interstitial"],
                groups = groups, labels = labels, colors = colors, linestyles = linestyles,
                output = f"{dose}_{construct_name}_drugs.png")
    '''

    groups = [["A"], ["B"], drugs,
              [f"{drug}-A" for drug in drugs] + [f"C-{drug}-A" for drug in drugs],
              [f"{drug}-B" for drug in drugs] + [f"C-{drug}-B" for drug in drugs],
              [f"{drug}-AB" for drug in drugs] + [f"C-{drug}-AB" for drug in drugs]]
    labels = ["A", "B", "(C-)drug", "(C-)drug-A", "(C-)drug-B", "(C-)drug-AB"]
    colors = ["tab:blue", "tab:red", "black", "skyblue", "salmon", "tab:purple"]
    system.plot(compartments = ["tumor_interstitial", "other_interstitial", "lung_interstitial", "SI_interstitial"],
                groups = groups, labels = labels, colors = colors,
                output = f"{dose}_{construct_name}_targets.png")
