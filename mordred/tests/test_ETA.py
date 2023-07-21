from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred import Calculator, ExtendedTopochemicalAtom

references = {
    # Roy et. al. Quantitative Structure-Activity Relationships in Drug Design,
    #     Predictive Toxicology, and Risk Assessment. p. 65
    "O=Cc1c(Cl)cccc1": {
        "ETA_alpha": 4.547,
        "AETA_alpha": 0.505,  # noqa: S001
        "ETA_shape_p": 0.230,
        "ETA_shape_y": 0.220,
        "ETA_shape_x": 0.000,  # noqa: S001
        "ETA_dAlpha_A": 0.005,  # noqa: S001
        "ETA_dAlpha_B": 0.000,
        "ETA_eta": 3.536,
        "ETA_eta_R": 9.998,
        "AETA_eta_F": 0.718,  # noqa: S001
        "ETA_eta_L": 1.543,
        "ETA_eta_RL": 4.343,
        "ETA_eta_FL": 2.800,  # noqa: S001
        "AETA_eta_BR": 0.018,  # noqa: S001
        "ETA_epsilon_1": 0.661,
        "ETA_epsilon_2": 0.861,
        "ETA_epsilon_3": 0.433,
        "ETA_epsilon_4": 0.553,
        "ETA_epsilon_5": 0.861,
        "ETA_dEpsilon_A": 0.228,  # noqa: S001
        "ETA_dEpsilon_B": 0.108,
        "ETA_dEpsilon_C": -0.119,
        "ETA_dEpsilon_D": 0.000,
        "AETA_beta_s": 0.556,  # noqa: S001
        "AETA_beta_ns": 0.889,  # noqa: S001
        "AETA_beta": 1.444,  # noqa: S001
        "AETA_dBeta": 0.333,
        "AETA_beta_ns_d": 0.056,  # noqa: S001
        "ETA_psi_1": 0.586,
        "ETA_dPsi_A": 0.128,  # noqa: S001
        "ETA_dPsi_B": 0.000,
    },
    "n1cc(O)ccc1": {
        "ETA_alpha": 3.233,
        "AETA_alpha": 0.462,  # noqa: S001
        "ETA_shape_p": 0.103,
        "ETA_shape_y": 0.155,
        "ETA_shape_x": 0.000,  # noqa: S001
        "ETA_dAlpha_A": 0.000,  # noqa: S001
        "ETA_dAlpha_B": 0.038,
        "ETA_eta": 2.048,
        "ETA_eta_R": 6.626,  # original:6.267. use sum of [eta_R]_i
        "AETA_eta_F": 0.654,  # noqa: S001
        "ETA_eta_L": 1.073,
        "ETA_eta_RL": 3.394,
        "ETA_eta_FL": 2.321,  # noqa: S001
        "AETA_eta_BR": 0.015,  # noqa: S001
        "ETA_epsilon_1": 0.631,
        "ETA_epsilon_2": 0.867,
        "ETA_epsilon_3": 0.433,
        "ETA_epsilon_4": 0.548,
        "ETA_epsilon_5": 0.796,
        "ETA_dEpsilon_A": 0.197,  # noqa: S001
        "ETA_dEpsilon_B": 0.083,
        "ETA_dEpsilon_C": -0.115,
        "ETA_dEpsilon_D": 0.071,
        "AETA_beta_s": 0.607,  # noqa: S001
        "AETA_beta_ns": 0.929,  # noqa: S001
        "AETA_beta": 1.536,  # noqa: S001
        "AETA_dBeta": 0.321,
        "AETA_beta_ns_d": 0.071,  # noqa: S001
        "ETA_psi_1": 0.533,
        "ETA_dPsi_A": 0.181,  # noqa: S001
        "ETA_dPsi_B": 0.000,
    },
    # Roy, et. al. Understanding the Basics of QSAR for Applications
    #     in Pharmaceutical Sciences and Risk Assessment. p.145
    "O=Cc1ccccc1": {  # noqa: S001
        "AETA_alpha": 0.479,
        "ETA_shape_p": 0.087,
        "ETA_shape_y": 0.130,
        "ETA_shape_x": 0.000,  # noqa: S001
        "ETA_dAlpha_A": 0.000,  # noqa: S001
        "ETA_dAlpha_B": 0.021,
        "ETA_eta": 2.576,
        "ETA_eta_R": 8.023,
        "AETA_eta_F": 0.681,  # noqa: S001
        "ETA_eta_L": 1.303,
        "AETA_eta_BR": 0.009,  # noqa: S001
        "ETA_epsilon_1": 0.583,
        "ETA_epsilon_2": 0.796,
        "ETA_epsilon_3": 0.433,
        "ETA_epsilon_4": 0.498,
        "ETA_epsilon_5": 0.796,
        "ETA_dEpsilon_A": 0.150,  # noqa: S001
        "ETA_dEpsilon_B": 0.085,
        "ETA_dEpsilon_C": -0.065,
        "ETA_dEpsilon_D": 0.000,
        "AETA_beta_s": 0.531,  # noqa: S001
        "AETA_beta_ns": 0.938,  # noqa: S001
        "AETA_beta": 1.469,  # noqa: S001
        "AETA_dBeta": 0.406,
        "AETA_beta_ns_d": 0.000,  # noqa: S001
        "ETA_psi_1": 0.602,
        "ETA_dPsi_A": 0.112,  # noqa: S001
        "ETA_dPsi_B": 0.000,
    },
}


def test_ETA():
    calc = Calculator(ExtendedTopochemicalAtom)

    for smi, desireds in references.items():
        mol = Chem.MolFromSmiles(smi)
        actuals = {str(d): v for d, v in zip(calc.descriptors, calc(mol))}

        for name, desired in desireds.items():
            assert (
                assert_almost_equal(
                    actuals[name], desired, 2, "{} of {}".format(name, smi)
                )
                is None
            )
