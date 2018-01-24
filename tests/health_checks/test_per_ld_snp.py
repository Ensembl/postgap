# ------------------------------------------------
# built-ins
import unittest

# local
from .base import TestPostgapBase
# ------------------------------------------------

class TestPostgapPerLdSnp(TestPostgapBase):

    def setUp(self):
        self.per_ld_snp = self.pg.groupby('ld_snp_rsID')

    def test_each_ld_snp_rsID_has_unique_chrom(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_ld_snp.chrom
        )

    def test_each_ld_snp_rsID_has_unique_pos(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_ld_snp.pos
        )

    def test_each_ld_snp_rsID_has_unique_GRCh38_chrom(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_ld_snp.GRCh38_chrom
        )

    def test_each_ld_snp_rsID_has_unique_GRCh38_pos(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_ld_snp.GRCh38_pos
        )


if __name__ == '__main__':
    unittest.main()
