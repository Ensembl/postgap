# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

class TestPostgapPerGeneANDLdSnp(TestPostgapBase):

    def setUp(self):
        self.per_gene_and_ld_snp = self.pg.groupby(['gene_id', 'ld_snp_rsID'])

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_VEP(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.VEP
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_GTEx(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.GTEx
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_PCHiC(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.PCHiC
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_DHS(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.DHS
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_Fantom5(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.Fantom5
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_Nearest(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.Nearest
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_Regulome(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.Regulome
        )

    def test_each_gene_id_and_ld_snp_rsID_pair_has_unique_score(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene_and_ld_snp.score
        )


if __name__ == '__main__':
    unittest.main()
