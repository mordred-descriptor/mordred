from __future__ import print_function

from numpy.testing import assert_almost_equal
from rdkit import Chem

from mordred.SA import SA

_sa_testing_data = {'Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]': '3.166',
                    'Cn1cc(NC=O)cc1C(=O)Nc1cc(C(=O)Nc2cc(C(=O)NCCC(N)=[NH2+])n(C)c2)n(C)c1': '3.328',
                    'OC(c1ccncc1)c1ccc(OCC[NH+]2CCCC2)cc1': '3.822', 'CC(C(=O)[O-])c1ccc(-c2ccccc2)cc1': '2.462',
                    'C[NH+](C)CC(O)Cn1c2ccc(Br)cc2c2cc(Br)ccc21': '3.577', 'NC(=[NH2+])NCC1COc2ccccc2O1': '3.290',
                    'CCC(C)(C)[NH2+]CC(O)COc1ccccc1C#N': '3.698', 'CC12CCC3C(CCC4CC(=O)CCC43C)C1CCC2=O': '3.912',
                    'COc1ccc(OC(=O)N(CC(=O)[O-])Cc2ccc(OCCc3nc(-c4ccccc4)oc3C)cc2)cc1': '2.644',
                    'COc1ccccc1OC(=O)c1ccccc1': '1.342', 'CC(C)CC[NH2+]CC1COc2ccccc2O1': '3.701',
                    'CN1CCN(C(=O)OC2c3nccnc3C(=O)N2c2ccc(Cl)cn2)CC1': '3.196',
                    'CCC1(c2ccccc2)C(=O)N(COC)C(=O)N(COC)C1=O': '2.759', 'Nc1ccc(S(=O)(=O)Nc2ccccc2)cc1': '1.529',
                    'O=C([O-])CCCNC(=O)NC1CCCCC1': '2.493', 'CCC(C)C(C(=O)OC1CC[N+](C)(C)CC1)c1ccccc1': '3.399',
                    'CCC(C)SSc1ncc[nH]1': '3.983', 'CC[N+](C)(CC)CCOC(=O)C(O)(c1cccs1)C1CCCC1': '3.471',
                    'CC12CCC3C(CCC4CC(=O)CCC43C)C1CCC2O': '3.994', 'CC12CCC3C4CCC(=O)C=C4CCC3C1CCC2O': '4.056',
                    'OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O': '4.282',
                    'CC(C)CC(CC[NH+](C(C)C)C(C)C)(C(N)=O)c1ccccn1': '4.092',
                    'C=CC1(C)CC(=O)C2(O)C(C)(O1)C(OC(C)=O)C(OC(=O)CC[NH+](C)C)C1C(C)(C)CCC(O)C12C': '5.519',
                    'C=CC[NH+]1CCCC1CNC(=O)c1cc(S(N)(=O)=O)cc(OC)c1OC': '4.286', 'CC(=O)OC1C[NH+]2CCC1CC2': '5.711',
                    'CC1(O)CCC2C3CCC4=CC(=O)CCC4(C)C3CCC21C': '4.022',
                    'CC(=O)OC1(C(C)=O)CCC2C3C=C(Cl)C4=CC(=O)C5CC5C4(C)C3CCC21C': '4.827',
                    'C#CC1(O)CCC2C3CCc4cc(OC)ccc4C3CCC21C': '3.810',
                    'C=CC1(C)CC(OC(=O)CSCC[NH+](CC)CC)C2(C)C3C(=O)CCC3(CCC2C)C(C)C1O': '6.200',
                    'O=C([O-])C(=O)Nc1nc(-c2ccc3c(c2)OCCO3)cs1': '2.594',
                    'CC[NH+]1CCCC1CNC(=O)C(O)(c1ccccc1)c1ccccc1': '3.950',
                    'CC(C)(OCc1nn(Cc2ccccc2)c2ccccc12)C(=O)[O-]': '2.573',
                    'Cc1nnc(C(C)C)n1C1CC2CCC(C1)[NH+]2CCC(NC(=O)C1CCC(F)(F)CC1)c1ccccc1': '5.316',
                    'Nc1ncnc2c1ncn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])C(O)C1O': '5.290',
                    'O=C([O-])CNC(=O)c1ccccc1': '2.097', 'c1ccc(OCc2ccc(CCCN3CCOCC3)cc2)cc1': '1.702',
                    'CC=CC1=C(C(=O)[O-])N2C(=O)C(NC(=O)C(N)c3ccc(O)cc3)C2SC1': '4.042',
                    'C[NH+]1CCCC1COc1cccnc1': '4.510', 'O=C([O-])C(O)CC(O)C(O)CO': '4.398',
                    'C[NH+]1CCCC1Cc1c[nH]c2ccc(CCS(=O)(=O)c3ccccc3)cc12': '3.921',
                    'C(=Cc1ccccc1)C[NH+]1CCN(C(c2ccccc2)c2ccccc2)CC1': '2.973', 'CC(c1ccccc1)N(C)C=O': '2.562',
                    'CC(=O)C1CCC2C3CCC4CC(C)(O)CCC4(C)C3CCC12C': '4.279', 'COc1cc(O)c(C(=O)c2ccccc2)c(O)c1': '1.868',
                    'COc1c2c(cc3c1C(O)N(C)CC3)OCO2': '3.183', 'CCC(C(=O)[O-])c1ccc(CC(C)C)cc1': '2.827',
                    'O=C([O-])C1[NH+]=C(c2ccccc2)c2cc(Cl)ccc2NC1(O)O': '4.011', 'OCC(O)COc1ccc(Cl)cc1': '2.102',
                    'NC(=O)NC(=O)C(Cl)c1ccccc1': '2.455', 'OC(c1ccccc1)(c1ccccc1)C1C[NH+]2CCC1CC2': '4.530',
                    'C[NH2+]CC(C)c1ccccc1': '3.471', 'Clc1cccc(Cl)c1N=C1NCCO1': '3.267',
                    '[NH3+]C(Cc1ccccc1)C(=O)CCl': '3.251', 'CC(C)Cn1cnc2c1c1ccccc1nc2N': '2.230',
                    'CC(O)CN(C)c1ccc(NN)nn1': '3.193', 'CC1(O)CCC2C3CCC4=CC(=O)CCC4=C3C=CC21C': '4.461',
                    'CCC(C(=O)[O-])c1ccc(-c2ccccc2)cc1': '2.505',
                    'CC(=O)OCC1OC(n2ncc(=O)[nH]c2=O)C(OC(C)=O)C1OC(C)=O': '3.832',
                    'Cn1cc(C(=O)c2cccc3ccccc32)cc1C(=O)[O-]': '2.456',
                    'Cc1cccc(-c2nc3ccccc3c(Nc3ccc4[nH]ncc4c3)n2)n1': '2.358', 'O=C([O-])C1CC2CCCCC2[NH2+]1': '5.422',
                    'CC(=O)OCC(=O)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C': '4.187', 'O=C(OCc1ccccc1)C(O)c1ccccc1': '2.038',
                    'CC(=O)OCC(=O)C1(O)CCC2C3CCC4=CC(=O)C=CC4(C)C3C(O)CC21C': '4.394',
                    'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1': '2.144',
                    'COCc1cccc(CC(O)C=CC2C(O)CC(=O)C2CCSCCCC(=O)OC)c1': '3.934',
                    'CCC(=O)N(c1ccccc1)C1CC[NH+](C(C)Cc2ccccc2)CC1': '3.582',
                    'CCOC(=O)Nc1ccc2c(c1)N(C(=O)CCN1CCOCC1)c1ccccc1S2': '2.446',
                    'O=C([O-])Cc1cc(=O)[nH]c(=O)[nH]1': '3.258', 'NC(=O)C([NH3+])Cc1c[nH]c2ccccc12': '3.224',
                    'NC(=O)C([NH3+])Cc1ccc(O)cc1': '3.280', 'O=C(c1cc2ccccc2o1)N1CCN(Cc2ccccc2)CC1': '1.799',
                    'O=C(CO)C(O)C(O)CO': '3.473', 'CC(Cc1ccccc1)NC(=O)C([NH3+])CCCC[NH3+]': '3.967',
                    'C[NH+]1CCC(c2c(O)cc(=O)c3c(O)cc(-c4ccccc4Cl)oc2-3)C(O)C1': '4.616',
                    'CN(C)c1ccc(O)c2c1CC1CC3C([NH+](C)C)C(=O)C(C(N)=O)=C(O)C3(O)C(=O)C1=C2O': '4.713',
                    'Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)CO)c2cc1C': '3.791',
                    'C[NH+]1C2CCC1CC(OC(=O)c1c[nH]c3ccccc13)C2': '4.892', 'Cc1ccccc1NC(=O)C(C)[NH+]1CCCC1': '3.809',
                    'O=S(=O)([O-])CCN1CCOCC1': '2.776', 'C[NH+]1CCN(CC(=O)N2c3ccccc3C(=O)Nc3cccnc32)CC1': '3.379',
                    'CCCCCC=CCC=CCCCCCCCC(=O)[O-]': '2.805', 'CC(CC([NH3+])C(=O)[O-])C(=O)[O-]': '5.690',
                    'CC1c2cccc(O)c2C(=O)C2=C(O)C3(O)C(O)=C(C(N)=O)C(=O)C([NH+](C)C)C3C(O)C21': '5.069',
                    'Cc1cc2nc3nc([O-])[nH]c(=O)c3nc2cc1C': '3.079', 'CC1=CC(C)C2(CO)COC(c3ccc(O)cc3)C1C2C': '4.749',
                    'CC[NH+]1CCC(=C2c3ccccc3CCc3ccccc32)C1C': '3.925'}


def test_SA():
    for smiles in _sa_testing_data.keys():
        mol = Chem.MolFromSmiles(smiles)
        yield assert_almost_equal(float(_sa_testing_data[smiles]), SA()(mol), 3)
