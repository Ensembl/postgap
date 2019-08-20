# ------------------------------------------------
# built-ins
import unittest

# pipped

# local
from utils.base import TestPostgapBase
# ------------------------------------------------

class TestPostgapRow(TestPostgapBase):

    def test_canary(self):
        self.assertTrue(True)

    # ld_snp_rsID
    def test_col_format_ld_snp_rsID(self):
        self.assert_series_valid_snp_id(self.pg.ld_snp_rsID)

    # chrom
    def test_col_format_chrom(self):
        self.assert_series_valid_chrom(self.pg.chrom)

    # pos
    def test_col_format_pos(self):
        self.assert_series_valid_genomic_coord(self.pg.pos)

    # GRCh38_chrom
    def test_col_format_GRCh38_chrom(self):
        self.assert_series_valid_chrom(self.pg.GRCh38_chrom)

    # GRCh38_pos
    def test_col_format_GRCh38_pos(self):
        self.assert_series_valid_genomic_coord(self.pg.GRCh38_pos)

    # afr_maf
    def test_col_format_afr_maf(self):
        self.assert_series_in_range(self.pg.afr, 0.0, 1.0, allow_na=True)

    # amr_maf
    def test_col_format_amr_maf(self):
        self.assert_series_in_range(self.pg.amr, 0.0, 1.0, allow_na=True)

    # eas_maf
    def test_col_format_eas_maf(self):
        self.assert_series_in_range(self.pg.eas, 0.0, 1.0, allow_na=True)

    # eur_maf
    def test_col_format_eur_maf(self):
        self.assert_series_in_range(self.pg.eur, 0.0, 1.0, allow_na=True)

    # sas_maf
    def test_col_format_sas_maf(self):
        self.assert_series_in_range(self.pg.sas, 0.0, 1.0, allow_na=True)

    # gene_id
    def test_col_format_gene_id(self):
        self.assert_series_valid_gene_id(self.pg.gene_id)

    # gene_chrom
    def test_col_format_gene_chrom(self):
        self.assert_series_valid_chrom(self.pg.gene_chrom)

    # gene_tss
    def test_col_format_gene_tss(self):
        self.assert_series_valid_genomic_coord(self.pg.gene_tss)

    # GRCh38_gene_chrom
    def test_col_format_GRCh38_gene_chrom(self):
        self.assert_series_valid_chrom(self.pg.GRCh38_gene_chrom)

    # GRCh38_gene_pos
    def test_col_format_GRCh38_gene_pos(self):
        self.assert_series_valid_genomic_coord(self.pg.GRCh38_gene_pos)

    # disease_name

    # disease_efo_id

    # r2
    def test_col_format_r2(self):
        self.assert_series_in_range(self.pg.r2, 0.7, 1.0)

    # cluster_id

    # gwas_source
    def test_col_format_gwas_source(self):
        self.assert_series_valid_gwas_source(self.pg.gwas_source)

    # gwas_snp
    def test_col_format_gwas_snp(self):
        self.assert_series_valid_snp_id(self.pg.gwas_snp)

    # gwas_pvalue
    def test_col_format_gwas_pval(self):
        self.assert_series_in_range(self.pg.gwas_pvalue, 0.0, 1.0)

    # GTEx
    def test_col_format_GTEx(self):
        self.assert_series_in_range(self.pg.GTEx, 0.0, 1.0)

    # VEP
    # Fantom5
    def test_col_format_Fantom5(self):
        self.assert_series_in_range(self.pg.Fantom5, 0.0, 1.0)

    # DHS
    def test_col_format_DHS(self):
        self.assert_series_in_range(self.pg.DHS, 0.0, 1.0)

    # PCHiC
    def test_col_format_PCHiC(self):
        self.assert_series_in_range(self.pg.PCHiC, 0.0, 1.0)


if __name__ == '__main__':
    unittest.main()
