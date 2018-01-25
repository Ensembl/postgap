# ------------------------------------------------
# built-ins
import unittest

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

# GRCh37 MHC region chr6:28,477,797-33,448,354
# (see https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
MHC_START_GRCH37 = 28477797
MHC_END_GRCH37 = 33448354

# GRCh38 MHC region chr6:28,510,120-33,480,577
# (see https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
MHC_START_GRCH38 = 28510120
MHC_END_GRCH38 = 33480577

MHC_CHROM = 6

class TestPostgapPerGene(TestPostgapBase):

    def test_canary(self):
        self.assertTrue(True)

    def test_mhc_region_filtered_gene_tss(self):
        self.assert_series_not_in_range(self.pg[self.pg.gene_chrom == MHC_CHROM].gene_tss,
                                        MHC_START_GRCH37, MHC_END_GRCH37)

    def test_mhc_region_filtered_GRCh38_gene_pos(self):
        self.assert_series_not_in_range(self.pg[self.pg.GRCh38_gene_chrom == MHC_CHROM].GRCh38_gene_pos,
                                        MHC_START_GRCH38, MHC_END_GRCH38)

    def test_mhc_region_filtered_pos(self):
        self.assert_series_not_in_range(self.pg[self.pg.chrom == MHC_CHROM].pos,
                                        MHC_START_GRCH37, MHC_END_GRCH37)

    def test_mhc_region_filtered_GRCh38_pos(self):
        self.assert_series_not_in_range(self.pg[self.pg.GRCh38_chrom == MHC_CHROM].GRCh38_pos,
                                        MHC_START_GRCH38, MHC_END_GRCH38)


if __name__ == '__main__':
    unittest.main()
