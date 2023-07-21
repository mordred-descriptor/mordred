from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred import Calculator, DistanceMatrix, AdjacencyMatrix, is_missing

# Balaban, A. T.; Ciubotariu, D.; Medeleanu, M.
# Topological indices and real number vertex invariants based on graph eigenvalues or eigenvectors.
# J. Chem. Inf. Comput. Sci. 1991, 31 (4), 517-523 DOI: 10.1021/ci00004a014.

descs = "          VE1_A     VE3_A       VE1_D     VE3_D       VR1_A       VR1_D".strip().split()
data = """
CC                 1.41421   -1.26286    1.41421   -1.26286     1.41421     1.41421
CCC                1.70711   -0.66917    1.71563   -0.66419   4:3.36358   4:3.72243
CCCC               1.9465    -0.25026    1.97417   -0.23614     5.89199   3:6.52546
CC(C)C             1.93185   -0.25781    1.97226   -0.23711   4:5.5836    3:6.90091
CCCCC              2.1547     0.0745     2.20361    0.09695   4:8.98667   3:9.73947
CC(C)CC            2.13099    0.06344    2.20196    0.0962    4:8.62989  3:10.15828
CC(C)(C)C          2.12132    0.05889    2.20397    0.09711   4:8.00002  4:10.74143
CCCCCC             2.3419     0.34014    2.4118     0.36955  4:12.62883  3:13.31649
CCC(C)CC           2.3008     0.32243    2.40851    0.36818  4:12.26119    13.88003
CC(C)CCC           2.31281    0.32764    2.41168    0.3695   4:12.39112  3:13.67976
CC(C)C(C)C         2.3094     0.32616    2.41209    0.36967  3:11.52993    14.14868
CC(C)(C)CC         2.2855     0.31576    2.4111     0.36926  4:11.63736  4:14.40727
CCCCCCC            2.51367    0.56507    2.60364    0.60024  3:16.80002  4:17.22301
CCC(CC)CC          2.44949    0.5392     2.59754    0.59789    16.3923   3:18.05924
CCC(C)CCC          2.45839    0.54283    2.60089    0.59918  4:16.73877  3:17.78548
CC(C)CCCC          2.48138    0.55214    2.60499    0.60075  3:16.81935  4:17.51358
CC(C)C(C)CC        2.45983    0.54342    2.60267    0.59986  4:15.71166  4:18.19395
CCC(C)(C)CC        2.42986    0.53116    2.6005     0.59903  3:15.79398  4:18.5519
CC(C)CC(C)C        2.5        0.55962    2.60668    0.6014     15.31371  3:17.86571
CC(C)(C)CCC        2.42522    0.52925    2.60382    0.60031  3:16.77771  4:18.20072
CC(C)(C)C(C)C      2.45369    0.54092    2.60656    0.60136  4:14.73031  i:15.84788
CCCCCCCC           2.67347    0.76023    2.78244    0.80018  3:21.48306  4:21.43353
CCC(C)CCCC         2.60405    0.73392    2.78102    0.79967  3:22.12919  3:21.93652
CCC(CC)CCC         2.58816    0.7278     2.77621    0.79794  3:21.55154  4:22.33865
CCC(C)C(C)CC       2.59508    0.73047    2.77988    0.79926  3:20.49061  3:22.59672
CCC(C)(CC)CC     4:2.55992    0.71683    2.77683    0.79817  3:20.39206  3:23.11878
CCCC(C)CCC         2.59808    0.73163    2.77872    0.79885  3:21.90424  3:22.06447
CC(C)CCCCC         2.63927    0.74736    2.78488    0.80106  3:21.88178  3:21.65563
CC(C)C(CC)CC       2.59417    0.73012    2.77888    0.79891  4:20.30079    22.71023
CC(C)C(C)CCC       2.59179    0.72921    2.78154    0.79986  3:21.29856    22.38291
CC(C)CC(C)CC       2.63932    0.74738    2.78375    0.80066  4:19.93766  1:22.24115
CCC(C)(C)CCC       2.55245    0.71391  1:2.77026  2:0.7958   3:21.69885  0:22.65171
CC(C)CCC(C)C       2.68328    0.7639   3:2.78755  3:0.80202  4:19.35733  2:21.91039
CC(C)C(C)C(C)C     2.61804    0.73928    2.78509    0.80114  3:19.13133  3:22.76268
CC(C)C(C)(C)CC     2.58138    0.72518    2.7826     0.80024    19.39011  4:23.197
CC(C)(C)CCCC       2.54201    0.70981    2.78428  0:0.88085  3:23.97131    22.25066
CC(C)(C)C(C)CC     2.58279    0.72573    2.7838     0.80067    19.83745  3:23.03447
CC(C)(C)CC(C)C     2.61183    0.73691    2.78783    0.80212  4:19.55118  3:22.60561
CC(C)(C)C(C)(C)C   2.6026     0.73337    2.7892     0.80261  4:17.88167  3:23.56702
""".strip().split(
    "\n"
)


def parse_reference(a):
    if a[:2] == "i:":
        return 0, None
    elif a[1] == ":":
        return int(a[0]), float(a[2:])

    return 5, float(a)


def test_VEA():
    calc = Calculator([AdjacencyMatrix, DistanceMatrix])

    for line in data:
        line = line.strip().split()

        smi = line[0]
        mol = Chem.MolFromSmiles(smi)

        desireds = dict(zip(descs, map(parse_reference, line[1:])))
        actuals = {str(k): v for k, v in zip(calc.descriptors, calc(mol))}

        for desc in descs:
            actual = actuals[desc]
            decimal, desired = desireds[desc]
            if desired is None:
                continue

            assert not is_missing(actual), actual

            assert (
                assert_almost_equal(
                    actual,
                    desired,
                    decimal,
                    "{} of {}".format(desc, smi),
                )
            ) is None
