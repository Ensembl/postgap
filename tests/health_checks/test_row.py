# ------------------------------------------------
# built-ins
import unittest

# pipped

# local
from .base import TestPostgapBase
# ------------------------------------------------

class TestPostgapRow(TestPostgapBase):

    def test_canary(self):
        self.assertTrue(True)

    def test_score_range(self):
        self.assert_series_in_range(self.pg.score, 0, 6)

    def test_r2_range(self):
        self.assert_series_in_range(self.pg.r2, 0.7, 1.0)

    def test_gwas_pval_range(self):
        self.assert_series_in_range(self.pg.gwas_pvalue, 0.0, 1.0)


if __name__ == '__main__':
    unittest.main()
