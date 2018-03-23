# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
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

    def test_each_gene_id_has_unique_GRCh38_gene_chrom(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene.GRCh38_gene_chrom
        )

    def test_each_gene_id_has_unique_GRCh38_gene_pos(self):
        self.assert_groupby_series_is_unique_per_group(
            self.per_gene.GRCh38_gene_pos
        )


if __name__ == '__main__':
    unittest.main()
