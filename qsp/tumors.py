UT44 = {"name": "tumor"}
UT44.update({"volume": 170 * units.ul, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
UT44.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
UT44.update({"capillary_radius": 10 * units.um, "capillary_permeability": 3e-7 * units.cm/units.s})
UT44.update({"diffusion": 10 * units.um**2 / units.s})
UT44.update({"density_cell": 3e8 * 0.44 / units.ml, "density_T": 3e8 * 0.15 / units.ml, "density_NK": 3e8 * 0.02 / units.ml})
UT44.update({"num_A": 7e5, "num_B": 1.45e6})

FTC238 = {"name": "tumor"}
FTC238.update({"volume": 170 * units.ul, "volume_plasma_proportion": 0.07, "volume_interstitial_proportion": 0.55})
FTC238.update({"plasma_flow_density": 12.7 / units.h, "lymphatic_flow_ratio": 0.002})
FTC238.update({"capillary_radius": 10 * units.um, "capillary_permeability": 3e-7 * units.cm/units.s})
FTC238.update({"diffusion": 10 * units.um**2 / units.s})
FTC238.update({"density_cell": 3e8 * 0.44 / units.ml, "density_T": 3e8 * 0.15 / units.ml, "density_NK": 3e8 * 0.02 / units.ml})
FTC238.update({"num_A": 1e5, "num_B": 1e5})
