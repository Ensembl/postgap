# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

class TestPValueFiltering(TestPostgapBase):

    def test_one_pvalue_per_gwas_pmid_AND_efo_id(self):
        self.skipTest('ONE PVALUE PER GWAS PMID AND EFO ID')

    def test_largest_pvalue_of_set_per_gwas_pmid_AND_efo_id(self):
        self.skipTest('LARGEST PVALUE OF SET PER GWAS PMID AND EFO ID')


if __name__ == '__main__':
    unittest.main()
