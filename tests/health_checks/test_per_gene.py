# ------------------------------------------------
# built-ins
import unittest

# local
from .base import TestPostgapBase
# ------------------------------------------------

class TestPostgapPerGene(TestPostgapBase):

    def setUp(self):
        self.per_gene = self.pg.groupby('gene_id')

    def test_each_gene_id_has_unique_gene_symbol(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene.gene_symbol
        )

    def test_each_gene_id_has_unique_gene_chrom(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene.gene_chrom
        )

    def test_each_gene_id_has_unique_gene_tss(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene.gene_tss
        )


if __name__ == '__main__':
    unittest.main()
