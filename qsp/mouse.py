plasma = {"name": "plasma"}
plasma.update({"volume": 0.944 * units.ml})
plasma.update({"num_T": 7.9E+09 * (0.944/3126), "num_NK": 1.6E+09 * (0.944/3126)})
plasma.update({"conc_A": 5 * units.ug/units.l / (170 * units.kDa), "conc_B": 0 * units.nM})

lymph = {"name": "lymph"}
lymph.update({"volume": 0.113 * units.ml})
lymph.update({"num_T": 3.6E+11 * (0.113/274), "num_NK": 6.7E+08 * (0.113/274)})
lymph.update({"conc_A": 0 * units.nM, "conc_B": 0 * units.nM})

bone = {"name": "bone"}
bone.update({"volume_plasma": 0.0621 * units.ml, "volume_interstitial": 0.525 * units.ml})
bone.update({"plasma_flow": 15.2 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
bone.update({"vascular_reflection": 0.85, "lymphatic_reflection": 0.2})
bone.update({"num_cell": 4.77E+09 * 0.5 * (0.525/1891), "num_T": 2.1E+10 * (0.525/1891), "num_NK": 3.3E+09 * (0.525/1891)})
bone.update({"num_A": 0, "num_B": 0})

lung = {"name": "lung"}
lung.update({"volume_plasma": 0.0295 * units.ml, "volume_interstitial": 0.0384 * units.ml})
lung.update({"vascular_reflection": 0.95, "lymphatic_reflection": 0.2})
lung.update({"plasma_flow": 373 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
lung.update({"num_cell": 2.36E+11 * 0.5 * (0.0384/300), "num_T": 1.3E+10 * (0.0384/300), "num_NK": 7.2E+08 * (0.0384/300)})
lung.update({"num_A": 1000, "num_B": 0}) # num_B = 1019 from Liyuan

liver = {"name": "liver"}
liver.update({"volume_plasma": 0.164 * units.ml, "volume_interstitial": 0.385 * units.ml})
liver.update({"vascular_reflection": 0.85, "lymphatic_reflection": 0.2})
liver.update({"plasma_flow": 10.3 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
liver.update({"num_cell": 3.90E+11 * 0.5 * (0.0385/429), "num_T": 7.9E+09 * (0.0385/429), "num_NK": 5.5E+09 * (0.0385/429)})
liver.update({"num_A": 10000, "num_B": 0})

SI = {"name": "SI"}
SI.update({"volume_plasma": 0.0116 * units.ml, "volume_interstitial": 0.127 * units.ml})
SI.update({"vascular_reflection": 0.9, "lymphatic_reflection": 0.2})
SI.update({"plasma_flow": 58.1 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
SI.update({"num_cell": 7.2e11 * 0.5 * (0.127/67.1), "num_T": 1.8E+10 * (0.127/67.1), "num_NK": 8.1E+08 * (0.127/67.1)})
SI.update({"num_A": 1000, "num_B": 10000})

kidney = {"name": "kidney"}
kidney.update({"volume_plasma": 0.0289 * units.ml, "volume_interstitial": 0.0788 * units.ml})
kidney.update({"vascular_reflection": 0.9, "lymphatic_reflection": 0.2})
kidney.update({"plasma_flow": 68.5 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
kidney.update({"num_cell": 1.06E+11 * 0.5 * (0.0788/49.8), "num_T": 496504047 * (0.0788/49.8), "num_NK": 1023295100 * (0.0788/49.8)}) # T anc NK numbers are calculated as density * volume
kidney.update({"num_A": 10000, "num_B": 0})

gallbladder = {"name": "gallbladder"}
gallbladder.update({"volume_plasma": 0.0289/4 * units.ml, "volume_interstitial": 0.0788/4 * units.ml})
gallbladder.update({"vascular_reflection": 0.9, "lymphatic_reflection": 0.2})
gallbladder.update({"plasma_flow": 68.5/4 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
gallbladder.update({"num_cell": 1.06E+11 * 0.5/4 * (0.0788/49.8), "num_T": 496504047/4 * (0.0788/49.8), "num_NK": 1023295100/4 * (0.0788/49.8)}) # T anc NK numbers are calculated as density * volume
gallbladder.update({"num_A": 10000, "num_B": 100000})

other = {"name": "other"}
other.update({"volume_plasma": 1 * units.ml, "volume_interstitial": 5 * units.ml})
other.update({"plasma_flow": 100 * units.ml/units.h, "lymphatic_flow_ratio": 0.002})
other.update({"vascular_reflection": 0.95, "lymphatic_reflection": 0.2})
other.update({"num_cell": 1e13 * 0.5 * (1/1000), "num_T": 1.1E+10 * (1/1000), "num_NK": 3.1E+09 * (1/1000)})
other.update({"num_A": 10000, "num_B": 0})
