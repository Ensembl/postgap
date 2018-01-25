# ------------------------------------------------
# built-ins
import unittest

# local
from .base import TestPostgapBase
# ------------------------------------------------

class TestPostgapPerLdSnpANDGwasSnp(TestPostgapBase):

    def setUp(self):
        self.per_ld_snp_and_gwas_snp = self.pg.groupby(['ld_snp_rsID', 'gwas_snp'])

    def test_each_ld_snp_rsID_and_gwas_snp_pair_has_unique_r2(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_ld_snp_and_gwas_snp.r2
        )


if __name__ == '__main__':
    unittest.main()
