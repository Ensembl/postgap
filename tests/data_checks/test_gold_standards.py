# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

# TODO: Expected (target, disease) pairs from experts
EXPECTED_ASSOCIATIONS = [
    ('ENSG00001', 'EFO_00001'),
    ('ENSG00002', 'EFO_00001')
]

class TestPostgapGoldStandards(TestPostgapBase):

    def test_each_expected_association_is_present(self):
        self.skipTest('CHECK ALL GOLD STANDARD EXPECTED ASSOCIATIONS')


if __name__ == '__main__':
    unittest.main()
