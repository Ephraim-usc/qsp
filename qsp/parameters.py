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

plasma_flows = {"heart": 36.5 * units.ml/units.h, 
                "lung": 373 * units.ml/units.h,
                "muscle": 86.1 * units.ml/units.h,
                "skin": 27.8 * units.ml/units.h,
                "adipose": 13.4 * units.ml/units.h,
                "bone": 15.2 * units.ml/units.h,
                "brain": 11.8 * units.ml/units.h,
                "kidney": 68.5 * units.ml/units.h,
                "liver": 10.3 * units.ml/units.h,
                "SI": 58.1 * units.ml/units.h,
                "LI": 17.3 * units.ml/units.h,
                "pancreas": 6.24 * units.ml/units.h,
                "thymus": 1.19 * units.ml/units.h,
                "spleen": 8.18 * units.ml/units.h,
                "other": 10.9 * units.ml/units.h}

vascular_reflection_coefficients = {"heart": 0.95,
                                   "lung": 0.95,
                                   "muscle": 0.95,
                                   "skin": 0.95,
                                   "adipose": 0.95,
                                   "bone": 0.85,
                                   "brain": 0.99,
                                   "kidney": 0.9,
                                   "liver": 0.85,
                                   "SI": 0.9,
                                   "LI": 0.95,
                                   "pancreas": 0.9,
                                   "thymus": 0.9,
                                   "spleen": 0.85,
                                   "other": 0.95}

lymphatic_reflection_coefficient = 0.2
