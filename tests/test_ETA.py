from mordred import Calculator, ExtendedTopochemicalAtom
from rdkit import Chem
from numpy.testing import assert_almost_equal

references = {
    # Roy et. al. Quantitative Structure-Activity Relationships in Drug Design, Predictive Toxicology, and Risk Assessment. p. 65
    'O=Cc1c(Cl)cccc1': {
        "alpha": 4.547,
        "alpha'": 0.505,
        "shape_p": 0.230,
        "shape_y": 0.220,
        "shape_x": 0.000,
        "dAlpha_A": 0.005,
        "dAlpha_B": 0.000,
        "eta": 3.536,
        "eta_R": 9.998,
        "eta_F'": 0.718,
        "eta_L": 1.543,
        "eta_RL": 4.343,
        "eta_FL": 2.800,
        "eta_BR'": 0.018,
        "epsilon_1": 0.661,
        "epsilon_2": 0.861,
        "epsilon_3": 0.433,
        "epsilon_4": 0.553,
        "epsilon_5": 0.861,
        "dEpsilon_A": 0.228,
        "dEpsilon_B": 0.108,
        "dEpsilon_C": -0.119,
        "dEpsilon_D": 0.000,
        "beta_s'": 0.556,
        "beta_ns'": 0.889,
        "beta'": 1.444,
        "dBeta'": 0.333,
        "beta_ns_d'": 0.056,
        "psi_1": 0.586,
        "dPsi_A": 0.128,
        "dPsi_B": 0.000,
    },
    'n1cc(O)ccc1': {
        "alpha": 3.233,
        "alpha'": 0.462,
        "shape_p": 0.103,
        "shape_y": 0.155,
        "shape_x": 0.000,
        "dAlpha_A": 0.000,
        "dAlpha_B": 0.038,
        "eta": 2.048,
        "eta_R": 6.626,  # original:6.267. use sum of [eta_R]_i
        "eta_F'": 0.654,
        "eta_L": 1.073,
        "eta_RL": 3.394,
        "eta_FL": 2.321,
        "eta_BR'": 0.015,
        "epsilon_1": 0.631,
        "epsilon_2": 0.867,
        "epsilon_3": 0.433,
        "epsilon_4": 0.548,
        "epsilon_5": 0.796,
        "dEpsilon_A": 0.197,
        "dEpsilon_B": 0.083,
        "dEpsilon_C": -0.115,
        "dEpsilon_D": 0.071,
        "beta_s'": 0.607,
        "beta_ns'": 0.929,
        "beta'": 1.536,
        "dBeta'": 0.321,
        "beta_ns_d'": 0.071,
        "psi_1": 0.533,
        "dPsi_A": 0.181,
        "dPsi_B": 0.000,
    },
# Roy, et. al. Understanding the Basics of QSAR for Applications in Pharmaceutical Sciences and Risk Assessment. p.145
    'O=Cc1ccccc1': {
        "alpha'": 0.479,
        "shape_p": 0.087,
        "shape_y": 0.130,
        "shape_x": 0.000,
        "dAlpha_A": 0.000,
        "dAlpha_B": 0.021,
        "eta": 2.576,
        "eta_R": 8.023,
        "eta_F'": 0.681,
        "eta_L": 1.303,
        "eta_BR'": 0.009,
        "epsilon_1": 0.583,
        "epsilon_2": 0.796,
        "epsilon_3": 0.433,
        "epsilon_4": 0.498,
        "epsilon_5": 0.796,
        "dEpsilon_A": 0.150,
        "dEpsilon_B": 0.085,
        "dEpsilon_C": -0.065,
        "dEpsilon_D": 0.000,
        "beta_s'": 0.531,
        "beta_ns'": 0.938,
        "beta'": 1.469,
        "dBeta'": 0.406,
        "beta_ns_d'": 0.000,
        "psi_1": 0.602,
        "dPsi_A": 0.112,
        "dPsi_B": 0.000,
    },
}

def test_ETA():
    calc = Calculator(ExtendedTopochemicalAtom)

    for smi, desireds in references.items():
        mol = Chem.MolFromSmiles(smi)
        actuals = {str(d): v for d, v in calc(mol)}

        for name, desired in desireds.items():
            yield assert_almost_equal, actuals['ETA_' + name], desired, 2, '{} of {}'.format(name, smi)
