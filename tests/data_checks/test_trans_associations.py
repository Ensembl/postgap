# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

class TestTransAssociations(TestPostgapBase):
    '''
    Consider first association per group only, as chromosomes
    being consistent within groups is tested elsewhere.
    '''

    def setUp(self):
        self.per_gene_and_ld_snp = self.pg.groupby(['gene_id', 'ld_snp_rsID'])

    def assert_no_trans_associations(self, chrom_field, gene_chrom_field):
        firsts = self.per_gene_and_ld_snp.head(1)
        chroms_match = (firsts[gene_chrom_field] == firsts[chrom_field])
        all_chroms_match = chroms_match.all()
        first_exception = None
        if (not all_chroms_match):
            fields = ['ld_snp_rsID', 'gene_id', 'gene_symbol', chrom_field, gene_chrom_field]
            first_exception = firsts[~chroms_match][fields].head(1).to_string(index=False)
        self.assertTrue(all_chroms_match, first_exception)

    def test_trans_associations_filtered(self):
        self.assert_no_trans_associations('chrom', 'gene_chrom')

    def test_trans_associations_filtered_GRCh38(self):
        self.assert_no_trans_associations('GRCh38_chrom', 'GRCh38_gene_chrom')


if __name__ == '__main__':
    unittest.main()
