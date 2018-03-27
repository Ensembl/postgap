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
        self.assert_series_in_range(self.pg.afr_maf, 0.0, 1.0, allow_na=True)

    # amr_maf
    def test_col_format_amr_maf(self):
        self.assert_series_in_range(self.pg.amr_maf, 0.0, 1.0, allow_na=True)

    # eas_maf
    def test_col_format_eas_maf(self):
        self.assert_series_in_range(self.pg.eas_maf, 0.0, 1.0, allow_na=True)

    # eur_maf
    def test_col_format_eur_maf(self):
        self.assert_series_in_range(self.pg.eur_maf, 0.0, 1.0, allow_na=True)

    # sas_maf
    def test_col_format_sas_maf(self):
        self.assert_series_in_range(self.pg.sas_maf, 0.0, 1.0, allow_na=True)

    # gene_symbol
    def test_col_format_gene_symbol(self):
        self.skipTest('DEFINE FORMAT FOR gene_symbol')

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
    def test_col_format_disease_efo_id(self):
        self.assert_series_valid_efo_id(self.pg.disease_efo_id)

    # score
    def test_col_format_score(self):
        self.skipTest('DEFINE FORMAT FOR score')

    # rank
    def test_col_format_rank(self):
        self.skipTest('DEFINE FORMAT FOR rank')

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

    # gwas_pvalue_description
    def test_col_format_gwas_pvalue_description(self):
        self.skipTest('DEFINE FORMAT FOR gwas_pvalue_description')

    # gwas_odds_ratio
    def test_col_format_gwas_odds_ratio(self):
        self.skipTest('DEFINE FORMAT FOR gwas_odds_ratio')

    # gwas_beta
    def test_col_format_gwas_beta(self):
        self.skipTest('DEFINE FORMAT FOR gwas_beta')

    # gwas_size
    def test_col_format_gwas_size(self):
        self.skipTest('DEFINE FORMAT FOR gwas_size')

    # gwas_pmid
    def test_col_format_gwas_pmid(self):
        self.skipTest('DEFINE FORMAT FOR gwas_pmid')
        
    # gwas_study
    def test_col_format_gwas_study(self):
        self.skipTest('DEFINE FORMAT FOR gwas_study')

    # gwas_reported_trait
    def test_col_format_gwas_reported_trait(self):
        self.skipTest('DEFINE FORMAT FOR gwas_reported_trait')

    # ls_snp_is_gwas_snp
    def test_col_format_ls_snp_is_gwas_snp(self):
        self.skipTest('DEFINE FORMAT FOR ls_snp_is_gwas_snp')

    # vep_terms
    def test_col_format_vep_terms(self):
        self.skipTest('DEFINE FORMAT FOR vep_terms')

    # vep_sum
    def test_col_format_vep_sum(self):
        self.skipTest('DEFINE FORMAT FOR vep_sum')

    # vep_mean
    def test_col_format_vep_mean(self):
        self.skipTest('DEFINE FORMAT FOR vep_mean')

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

    # Nearest
    def test_col_format_Nearest(self):
        self.skipTest('DEFINE FORMAT FOR Nearest')

    # Regulome
    def test_col_format_Regulome(self):
        self.skipTest('DEFINE FORMAT FOR Regulome')

    # VEP_reg
    def test_col_format_VEP_reg(self):
        self.skipTest('DEFINE FORMAT FOR VEP_reg')

    # vep_max_score
    def test_col_format_VEP_max_score(self):
        self.skipTest('DEFINE FORMAT FOR VEP_max_score')

    # TODO: The following fields may be generated by Open Targets;
    #       probably should NOT test in Ensembl/postgap

    # # fg_score
    # def test_col_format_fg_score(self):
    #     self.skipTest('DEFINE FORMAT FOR fg_score')

    # # snp_gene_dist
    # def test_col_format_snp_gene_dist(self):
    #     self.skipTest('DEFINE FORMAT FOR snp_gene_dist')


if __name__ == '__main__':
    unittest.main()
