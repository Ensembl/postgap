# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

class TestGWASCatalogCoverage(TestPostgapBase):

    def test_each_gwas_efo_covered(self):
        self.skipTest('EACH GWAS EFO ID COVERED IN POSTGAP OUTPUT')

    def test_each_gwas_snp_covered(self):
        self.skipTest('EACH GWAS SNP ID COVERED IN POSTGAP OUTPUT')


if __name__ == '__main__':
    unittest.main()
