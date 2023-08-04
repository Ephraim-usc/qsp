organs = ["heart", "lung", "muscle", "skin", "adipose", "bone", "brain", "kidney", "liver", "SI", "LI", "pancreas", "thymus", "spleen", "other"]
tissues = ["plasma", "BC", "interstitial", "endosomal", "membrane", "cellular"]
centrals = ["plasma", "BC", "lymph"]
compartments = centrals + [f"{organ}_{tissue}" for organ in organs for tissue in tissues]

volumes_mouse = np.array([0.944, 0.773, 0.113,
                          0.00585, 0.00479, 0.0217, 0.000760, 0.0217, 0.119,
                          0.0295, 0.0241, 0.0384, 0.00102, 0.0384, 0.111,
                          0.249, 0.204, 1.47, 0.0566, 1.47, 9.34,
                          0.188, 0.154, 1.66, 0.0251, 1.66, 3.00,
                          0.0218, 0.0178, 0.337, 0.00991, 0.337, 1.60,
                          0.0621, 0.0508, 0.525, 0.0141, 0.525, 2.17,
                          0.0107, 0.00873, 0.0873, 0.00243, 0.0873, 0.376,
                          0.0289, 0.0236, 0.0788, 0.00263, 0.0788, 0.391,
                          0.164, 0.134, 0.385, 0.00963, 0.385, 1.23,
                          0.0116, 0.00950, 0.127, 0.00364, 0.127, 0.577,
                          0.0050, 0.00409, 0.0545, 0.00157, 0.0545, 0.248,
                          0.00534, 0.00437, 0.0169, 0.000485, 0.0169, 0.0699,
                          0.0005, 0.000405, 0.00153, 0.00005, 0.00153, 0.00653, 
                          0.0154, 0.0126, 0.0254, 0.000635, 0.0254, 0.0730,
                          0.0195, 0.0160, 0.0797, 0.00233, 0.0797, 0.348,
                         ]) * units.ml

volumes_human = np.array([3126, 2558, 274,
                          13.1, 10.8, 48.8, 1.71, 48.8, 267,
                          55.0, 45.0, 300, 5.00, 300, 595,
                          662, 541, 3910, 150, 3910, 24815,
                          127, 104, 1125, 17.0, 1125, 2035,
                          148, 121, 2289, 67.3, 2289, 10840,
                          224, 183, 1891, 50.8, 1891, 7817,
                          31.9, 26.1, 261, 7.25, 261, 1124,
                          18.2, 14.9, 49.8, 1.66, 49.8, 247,
                          183, 149, 429, 10.7, 429, 1371,
                          6.15, 5.03, 67.1, 1.93, 67.1, 305,
                          8.74, 7.15, 95.3, 2.74, 95.3, 434,
                          5.70, 4.66, 18.0, 0.518, 18.0, 74.7,
                          0.353, 0.288, 1.09, 0.0321, 1.09, 4.65,
                          26.8, 21.9, 44.3, 1.11, 44.3, 127,
                          204, 167, 831, 24.3, 831, 3626
                         ]) * units.ml

plasma_flows_mouse = np.array([36.5, 373, 86.1, 27.8, 13.4, 15.2, 11.8, 68.5, 10.3, 58.1, 17.3, 6.24, 1.19, 8.18, 10.9]) * units.ml/units.h
plasma_flows_human = np.array([7752, 181913, 33469, 11626, 11233, 2591, 21453, 36402, 13210, 12368, 12867, 3056, 353, 6343, 5521]) * units.ml/units.h

BC_flows_mouse = np.array([29.9, 305, 70.5, 22.8, 11.0, 12.4, 9.64, 56.1, 8.40, 47.5, 14.1, 5.10, 0.97, 6.70, 8.91]) * units.ml/units.h
BC_flows_human = np.array([6342, 148838, 27383, 9512, 9191, 2120, 17553, 29784, 10808, 10120, 10527, 2500, 289, 5189, 4517]) * units.ml/units.h

vascular_reflections = np.array([0.95, 0.95, 0.95, 0.95, 0.95, 0.85, 0.99, 0.9, 0.85, 0.9, 0.95, 0.9, 0.9, 0.85, 0.95])
lymphatic_reflection = 0.2

mouse = {}
mouse.update({f"volume_{compartment}":x for compartment, x in zip(compartments, volumes_mouse)})
mouse.update({f"plasma_flow_{organ}":x for organ, x in zip(organs, plasma_flows_mouse)})
mouse.update({f"BC_flow_{organ}":x for organ, x in zip(organs, BC_flows_mouse)})
mouse.update({f"vascular_reflection_{organ}":x for organ, x in zip(organs, vascular_reflectionse)})
mouse.update({"lymphatic_reflection":lymphatic_reflection})

human = {}
human.update({f"volume_{compartment}":x for compartment, x in zip(compartments, volumes_human)})
human.update({f"plasma_flow_{organ}":x for organ, x in zip(organs, plasma_flows_human)})
human.update({f"BC_flow_{organ}":x for organ, x in zip(organs, BC_flows_human)})
human.update({f"vascular_reflection_{organ}":x for organ, x in zip(organs, vascular_reflectionse)})
human.update({"lymphatic_reflection":lymphatic_reflection})



