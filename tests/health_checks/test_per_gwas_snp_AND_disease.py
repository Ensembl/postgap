# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

class TestPostgapPerGwasSnpANDDisease(TestPostgapBase):

    def setUp(self):
        self.per_gwas_snp_and_disease = self.pg.groupby(['gwas_snp', 'disease_efo_id'])

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_pvalue(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_pvalue')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_pvalue_description(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_pvalue_description')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_odds_ratio(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_odds_ratio')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_beta(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_beta')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_size(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_size')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_pmid(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_pmid')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_study(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_study')

    def test_each_gwas_snp_and_disease_efo_id_pair_has_unique_gwas_reported_trait(self):
        self.skipTest('CHECK FOR UNIQUENESS OF gwas_reported_trait')


if __name__ == '__main__':
    unittest.main()
