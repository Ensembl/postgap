# ------------------------------------------------
# built-ins
import unittest
# ------------------------------------------------

VALID_CHROMOSOMES = [*[str(chr) for chr in range(23)], 'X', 'Y']
VALID_GWAS_SOURCES = ['GWAS Catalog']
VALID_SNP_ID_REGEX = '^rs\d+$'


class TestPostgapBase(unittest.TestCase):
    """
    Base class for postgap tests. Provides common utility methods.
    """

    def __init__(self, test_name, postgap=None):
        super(TestPostgapBase, self).__init__(test_name)
        self.pg = postgap

    def assert_series_in_range(self, series, low, high):
        """
        Check if all values in a `pandas.Series` are in the range [low, high].
        """
        between = series.between(low, high)
        all_between = between.all()
        self.assertTrue(all_between,
                        series[~between].head(1).to_string(index=False))

    def assert_valid_snp_id(self, series):
        """
        Check if all values in a `pandas.Series` are valid SNP ids.
        """
        valid_snp_id = series.str.match(VALID_SNP_ID_REGEX)
        all_valid_snp_id = valid_snp_id.all()
        self.assertTrue(all_valid_snp_id,
                        series[~valid_snp_id].head(1).to_string(index=False))

    def assert_series_valid_chrom(self, series):
        """
        Check if all values in a `pandas.Series` are valid chromosomes.
        """
        series_str = series.apply(str)
        chrs = series_str.unique()
        all_valid = all(chr in VALID_CHROMOSOMES for chr in chrs)
        invalid_freqs = series_str[~series_str.isin(VALID_CHROMOSOMES)].value_counts()
        self.assertTrue(all_valid, invalid_freqs)

    def assert_series_valid_gwas_source(self, series):
        """
        Check if all values in a `pandas.Series` are valid GWAS sources.
        """
        sources = series.unique()
        all_valid = all(s in VALID_GWAS_SOURCES for s in sources)
        invalid_freqs = series[~series.isin(VALID_GWAS_SOURCES)].value_counts()
        self.assertTrue(all_valid, invalid_freqs)
        