# ------------------------------------------------
# built-ins
import unittest

# local
from tests.utils.base import TestPostgapBase

# ------------------------------------------------

# TODO: Expected (target, disease) pairs from experts
EXPECTED_ASSOCIATIONS = [
    ('ENSG00001', 'EFO_00001'),
    ('ENSG00002', 'EFO_00001')
]

# TODO: Expected (snp, gene) pairs from experts
# extracted by qaqc.GoldStandardExamples.py
EXPECTED_ASSIGNMENTS = {
    'rs10872142': {'FRK'},
    'rs10889356': {'DOCK7, ANGPTL3'},
    'rs10941679': {'MRPS30', 'FGF10'},
    'rs111961716': {'ABHD8', 'ANKLE1'},
    'rs12740374': {'SORT1'},
    'rs1421085': {'ARID5B', 'IRX3', 'IRX5'},
    'rs1427407': {'BCL11A'},
    'rs2277862': {'CPNE1'},
    'rs356168': {'SNCA'},
    'rs378854': {'MYC', 'PVT1'},
    'rs4442975': {'IGFBP5'},
    'rs4500960': {'TBR1'},
    'rs554219': {'CCND1'},
    'rs56069439': {'ABHD8', 'ANKLE1'},
    'rs7463708': {'PCAT1'},
    'rs75915166': {'CCND1'},
    'rs7606173': {'BCL11A'},
    'rs78540526': {'CCND1'},
    'rs9926296': {'EDN1'}
}


class TestPostgapGoldStandards(TestPostgapBase):
    def test_each_expected_association_is_present(self):
        self.skipTest('CHECK ALL GOLD STANDARD EXPECTED ASSOCIATIONS')


if __name__ == '__main__':
    unittest.main()
